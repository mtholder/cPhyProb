#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

from cPhyProb import get_logger, ComputationalResourceFlags
_LOG = get_logger(__name__)

from cPhyProb.ccore.dsct_model import cstate_set_lookup_ctor, cphyprob_init, cphyprob_free
from cPhyProb.asrv import RateHetManager

class LeafDataCodingEnum:
    STATE_CODE_ARRAY = 0
    PARTIAL_LIKELIHOOD_ARRAY = 1




def _partial_ambig_check_or_convert(do_convert, row, num_states):
    """If `do_convert` is True, then the function will return a new list in which
    any ambiguity code is changed to the `num_states` state code.

    If `do_convert` is False, then the function will return None if the `row` contains
    a partial ambiguity code.  Otherwise `row` will be returned
    """
    if do_convert:
        cr = []
        for cell in row:
            if cell > num_states or cell < 0:
                cr.append(num_states)
            else:
                cr.append(cell)
        return cr
    for cell in row:
        if cell > num_states or cell < 0:
            return None
    return row

class PhyloSubsetLikelihoodCalculator(object):
    def __init__(self, char_matrix=None, model=None, discrete_char_type=None, **kwargs):
        """
            `char_matrix` should be a list.
        It must be either a list of strings or a list of lists (if it is strings,
        then it must a `discrete_char_type` argument should be passed in to
        allow for translation).

            `discrete_char_type` or `model` can supply the number of states.

            `model` will be the source of num_rate_categories if it is None, then
        a `num_rate_categories` argument can be used.

        other kwargs and their defaults:
            `convert_partial_ambig_to_missing` = False
            `num_partial_storage_per_internal_edge` = 1
            `num_extra_partial_storage` = 0
            `max_num_internal_edges` = None
            `num_eigen_storage ` =  1
            `num_prob_mat_per_internal_edge` = 1
            `num_prob_mat_per_terminal_edge` = 1
            `num_extra_prob_mat_storage` = 0
            `rescale_every` = 10
            `rooted` = False
            `num_rate_categories` = 1 only used if a model is not passed in

        """
        self.initialized = False
        kwarg_defaults = { "convert_partial_ambig_to_missing": False,
                          "num_partial_storage_per_internal_edge": 1,
                          "num_extra_partial_storage": 0,
                          "prob_mat_per_edge": 1,
                          "max_num_internal_edges": None,
                          "num_eigen_storage": 1,
                          "num_prob_mat_per_internal_edge": 1,
                          "num_prob_mat_per_terminal_edge": 1,
                          "num_extra_prob_mat_storage": 0,
                          "rescale_every": 10,
                          "rooted": False,
                          }
        for k in kwargs.keys():
            if not k in kwarg_defaults:
                raise ValueError("Argument '%s' is not supported" % k)
        kwarg_defaults.update(kwargs)
        for k, v in kwarg_defaults.iteritems():
            setattr(self, k, v)
        if char_matrix is None:
            mat = []
        else:
            mat = char_matrix
        self.num_leaves = len(mat)
        self.num_patterns = 0
        dt = discrete_char_type
        if dt is not None:
            state_code_lookup = dt.state_code_lookup
            self.num_states = dt.num_states
        else:
            num_states = kwargs.get("num_states")
            if model and (not num_states):
                num_states = model.num_states
            state_code_lookup = cstate_set_lookup_ctor(num_states, expansions)
            self.num_states = num_states
        if model:
            if model.num_states != self.num_states:
                raise ValueError("The number of states in the model must agree with the character type")
            self.num_rate_categories = model.num_rate_categories
            self.num_prob_models = model.num_prob_models
        else:
            self.num_rate_categories = kwargs.get("num_rate_categories", 1)
        self.rate_het = kwargs.get("rate_het")
        if not self.rate_het:
            self.rate_het = RateHetManager(rates=[1.0])

        enc_data = []
        if self.num_leaves > 0:
            first_el = mat[0]
            if isinstance(first_el, list):
                enc_data = mat
                for row in mat:
                    if not isinstance(row, list):
                        raise ValueError("All rows of the matrix must be of the same type")
            else:
                if not dt:
                    raise ValueError("discrete_char_type must be passed as an argument if the matrix rows are not lists of coded states")
                enc_data = [dt.to_indices(i) for i in mat]
            self.num_patterns = len(first_el)

        # The format of the leaf data is not known at this point.
        # Depending on whether or not leaf data has partial ambiguities that must
        #   be maintained, we must decide on whether the data can go in state
        #   code array or a partial.
        # Here we compile lists of the raw data for each storage type and a
        #   mapping (leaf_index_2_destination) indicating the location of the data)
        # leaf_data_matrix has two lists (one for each format type).  Each list
        #   contains elements that are tuples of (leaf index, list of data as integers)
        self.leaf_data_matrix = [[], []]
        self.leaf_index_2_destination = []
        leaf_index = 0
        _LOG.debug("enc_data = %s " % str(enc_data))
        while enc_data:
            row = enc_data.pop(0)
            cr = _partial_ambig_check_or_convert(self.convert_partial_ambig_to_missing, row, self.num_states)
            if cr is None:
                dest_index = LeafDataCodingEnum.PARTIAL_LIKELIHOOD_ARRAY
                cr = row
            else:
                dest_index = LeafDataCodingEnum.STATE_CODE_ARRAY

            dest_list = self.leaf_data_matrix[dest_index]
            dest_list_pos = len(dest_list)
            self.leaf_index_2_destination.append((dest_index, dest_list_pos))
            dest_list.append((leaf_index, cr))

            leaf_index += 1

        self.num_leaves = leaf_index

        # here we figure out how much memory we'll need to allocate
        self.num_state_code_arrays = len(self.leaf_data_matrix[LeafDataCodingEnum.STATE_CODE_ARRAY])
        if self.rooted:
            self.num_internal_edges = self.num_leaves - 2
        else:
            self.num_internal_edges = self.num_leaves - 3

        if self.max_num_internal_edges is None:
            self.max_num_internal_edges = self.num_internal_edges
        elif self.num_internal_edges > self.max_num_internal_edges:
            self.num_internal_edges = self.max_num_internal_edges

        self.num_edges = self.num_internal_edges + self.num_leaves
        self.num_terminal_edges = self.num_leaves

        self.num_leaves_as_partials = len(self.leaf_data_matrix[LeafDataCodingEnum.PARTIAL_LIKELIHOOD_ARRAY])
        self.num_partial_structs = (self.num_internal_edges + self.num_extra_partial_storage)*self.num_partial_storage_per_internal_edge
        self.num_partial_structs += self.num_leaves_as_partials

        self.num_prob_mat_arrays = self.num_extra_partial_storage
        self.num_prob_mat_arrays += self.num_prob_mat_per_internal_edge*self.num_internal_edges
        self.num_prob_mat_arrays += self.num_prob_mat_per_terminal_edge*self.num_terminal_edges
        self.total_num_prob_mats = self.num_rate_categories*self.num_prob_mat_arrays

        # We assume that we'll need a tree-wide rescaling array as well as the per-node rescalers
        # We add 2 instead of 1 so that the division rounds up
        # We multiply by num_partial_storage_per_internal_edge because this gives
        #   us a hint of how many sets of partials will be in use.
        self.num_rescalings_multipliers = (2 + int(self.num_edges/float(self.rescale_every)))*self.num_partial_storage_per_internal_edge
        #ldobjs = [cleaf_data_ctor(i, state_set_lookup, n_categ) for i in enc_data]

        self.pattern_weights = kwargs.get("pattern_weights")
        if not self.pattern_weights:
            self.pattern_weights = [1.0]*self.num_patterns

        self.impl_handle = None

        self.resource_identifier = 0 # this is used to tell beagle what computational resoure we want to run on.
        self.resource_flags = ComputationalResourceFlags.COMPLEX & ComputationalResourceFlags.LOG_SCALER

    def activate(self):
        if not self.initialized:
            self.allocate_storage()

    def allocate_storage(self):
        if not self.impl_handle is None:
            self.free_storage()
        arg_list = [self.num_leaves,
                    self.num_patterns,
                    self.num_states,
                    self.num_rate_categories,
                    self.num_state_code_arrays,
                    self.num_partial_structs,
                    self.num_prob_models, # number of q matrices
                    self.total_num_prob_mats,
                    self.num_eigen_storage,
                    self.num_rescalings_multipliers,
                    self.resource_identifier,
                    self.resource_flags]
        _LOG.debug("Calling cphyprob_init(%s)" % ', '.join([str(i) for i in arg_list]))
        self.impl_handle = cphyprob_init(*arg_list)
        self.initialized = True
    def free_storage(self):
        if self.initialized:
            self.initialized = False
            cphyprob_free(self.impl_handle)
    def __del__(self):
        self.free_storage()

class PhyloLikelihoodCalculator(object):
    def __init__(self, subsets=None, **kwargs):
        if subsets is None:
            self.model_subsets = []
        else:
            self.model_subsets = list(subsets)
        if kwargs:
            self.add_matrix(**kwargs)

    def add_matrix(self, **kwargs):
        pass

    def activate(self):
        for mod in self.model_subsets:
            mod.activate()

    def allocate_storage(self):
        for mod in self.model_subsets:
            mod.allocate_storage()

    def free_storage(self):
        for mod in self.model_subsets:
            mod.free_storage()

    def __del__(self):
        self.free_storage()


'''
"""Classes for encapsulating datatypes, and factory for encapsulating data
structures for characters at the leaves of a tree and the conditional
likelihood array for a tree.
"""
from cPhyProb.ccore.dsct_model import cleaf_data_ctor, cla_ctor, full_la_ctor
from cPhyProb.prob_mat import ProbMatrixArray
from cPhyProb.ccore.dsct_model import calc_ln_L_across_term, \
    calc_ln_L_across_internal, calc_anc_from_two_tips, \
    calc_anc_from_one_tip_one_intern, calc_anc_from_two_internals, \
    add_leaf_to_cla, add_internal_to_cla, full_la_set_freqs
from cPhyProb.ccore.dsct_model import get_cla_vals, get_full_la_vals

def print_clas(cla_struct):
    clas, rescalers = get_cla_vals(cla_struct)
    print "clas=", clas
    print "rescalers=", rescalers

def print_full_la(full_la):
    lnls, freqs, clas, rescalers = get_full_la_vals(full_la)
    print "site lnL=", lnls
    print "freqs   =", freqs
    print "clas    =", clas
    print "rescale =", rescalers

def calc_cla(par, leaves, internals, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold=10):
    par_struct = par.calc_struct
    n_leaves = len(leaves)
    n_internals = len(internals)
    if n_leaves > 0:
        l0, p0 = leaves[0].refresh_structs()
        if n_leaves > 1:
            l1, p1 = leaves[1].refresh_structs()
            calc_anc_from_two_tips(par_struct, l0, p0.raw_prob_mat, l1, p1.raw_prob_mat, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
            #print_clas(par_struct)
            for leaf in leaves[2:]:
                l1, p1 = leaf.refresh_structs()
                add_leaf_to_cla(par_struct, l1, p1.raw_prob_mat, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
                #print_clas(par_struct)
            for internal in internals:
                l1, p1 = internal.refresh_structs()
                add_internal_to_cla(par_struct, l1, p1.raw_prob_mat, beg_subset, end_subset, beg_cat, end_cat, 0, 1, rescale_threshold)
                #print_clas(par_struct)
            return
    if n_leaves + n_internals == 0:
        raise ValueError("Must have at least 2 neighbors to call calc_cla")
    if n_leaves == 0:
        l0, p0 = internals[0].refresh_structs()
        f = calc_anc_from_two_internals
        curr_internal = 1
    else:
        f = calc_anc_from_one_tip_one_intern
        curr_internal = 0
    internal = internals[curr_internal]
    l1, p1 = internal.refresh_structs()
    f(par_struct, l0, p0.raw_prob_mat, l1, p1.raw_prob_mat, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
    #print_clas(par_struct)
    for internal in internals[1+curr_internal:]:
        l1, p1 = internal.refresh_structs()
        add_internal_to_cla(par_struct, l1, p1.raw_prob_mat, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
        #print_clas(par_struct)

def calc_lnL(full_la, focal_node, leaves, internals, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold=10):
    n_leaves = len(leaves)
    n_internals = len(internals)
    focal_node_struct = focal_node.calc_struct
    if n_leaves + n_internals < 3:
        raise ValueError("Must have at least 3 neighbors to call calc_lnL")
    if n_leaves > 0:
        other, other_p = leaves[0].refresh_structs()
        calc_cla(focal_node, leaves[1:], internals, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
        lnL = calc_ln_L_across_term(full_la, other, other_p.raw_prob_mat, focal_node_struct, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
    else:
        other, other_p = internals[0].refresh_structs()
        calc_cla(focal_node, leaves, internals[1:], beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
        lnL = calc_ln_L_across_internal(full_la, other, other_p.raw_prob_mat, focal_node_struct, beg_subset, end_subset, beg_cat, end_cat, rescale_threshold)
    #print_full_la(full_la)
    return lnL

_INVALID_CALC_TEMPS = 0
_VALID_CALC_TEMPS = 1
_VALID_IF_REJECTED = 2
_RECALC_NEEDED = 3

class NodeCalcStruct(object):
    pass

class EdgeCalcStruct(object):
    def __init__(self, calc_struct, pmat_arr, model=None, alt_pmat_arr=None,  calc_key=None):
        """Creates a wrapper around the c structs used to calculate and store the
        conditional likelihoods for one edge.

         * `calc_struct` is the struct that holds the results of calculations using
            Felsenstein's pruning algorithm.  For terminal edges this will be a
            leaf_data struct. For internal edges this will be a cla (conditional
            likelihood array).  In either case, this opaque structure is used in
            `calc_cla` and `calc_lnL` functions
        """
        self.calc_struct = calc_struct
        self.pmat_arr = pmat_arr
        self.alt_pmat_arr = alt_pmat_arr
        self.model = model
        self.brlen = None
        self.alt_brlen = None
        self.pmat_calc = False
        self.alt_pmat_calc = False
        self.calc_key = calc_key

    def _check_brlen(self, b):
        if b < 0.0:
            raise ValueError("Branch length cannot be negative")

    def revert_brlen(self):
        if self.alt_pmat_arr:
            self.alt_pmat_arr, self.pmat_arr = self.pmat_arr, self.alt_pmat_arr
            self.alt_pmat_calc, self.pmat_calc = self.alt_pmat_calc, self.pmat_calc
        else:
            self.pmat_calc = False
        self.alt_brlen, self.brlen = self.brlen, self.alt_brlen

    def set_brlen(self, b, calc_key=None):
        self._check_brlen(b)
        self.pmat_calc = False
        self.brlen = b

    def propose_model_change(self):
        b = self.brlen
        self.revert_brlen()
        self.set_brlen(b)

    def propose_brlen(self, b):
        self._check_brlen(b)
        self.revert_brlen()
        self.set_brlen(b)

    def refresh_structs(self, force_recalc=False):
        """Returns a (leaf_data_type, pmat_array_type) for the leaf with the
        current branch length.
        """
        if force_recalc or not self.pmat_calc:
            if self.brlen is None:
                raise ValueError("branch length must be set before "\
                                 "refresh_leaf_data can be called")
            self.pmat_arr.calc_probs(self.brlen)
            self.pmat_calc = True
        return self.calc_struct, self.pmat_arr

class PMatStorage:
    """An class that acts like an enumeration of the ways for allocating
    ProbMatrixArray
    """
    EXTERNAL_PMATS = 0 # at least two pmats will need to created for the likelihood calculations
    ONE_PER_EDGE = 1 # every EdgeCalcStruct gets its own PMatArray - enables more
                     #   caching of temporaries, but requires more memory
    STORE_BACKUPS = 2 # allocated more pmats than branches, to avoid recalculation
                     # of the same values (if a set of parameters is rejected).



def get_tree_decorators(disc_char_mat, model, n_pmats=1, **kwargs):
    """Takes:
     * `disc_char_mat` -- the character data as a container of iterables,
     * `model` -- a DiscreteCharType object for the data.
     * `n_categ` -- the number of categories in the model (>1 for a mixture
        model).
     * `kwargs` (see below)
    Returns a tuple of the form:
        ([leaf EdgeCalcStruct], [internal EdgeCalcStruct], `full_la_obj`)
    if n_pmats > 0, or of the form:
        ([cleaf_data_obj], [cla_obj], `full_la_obj`)
    if n_pmats == 0

    The order of the leaf EdgeCalcStruct will correspond to the order of the
        rows of `disc_char_mat`
    The internal EdgeCalcStruct will have the necessary memory allocated, but
        will be unintialized (and assignable to any internal).
    The `full_la_obj` is needed for calculating the lnLikelihood from the
        "effective root" of the tree.

    `n_pmats` is the number of ProbMatArray objects to be assigned to
        the EdgeCalcStruct objects.
        If `n_pmats` == 0 then cleaf_data_obj and cla_obj will be returned "raw"
        rather than wrapped in EdgeCalcStruct objects.


    The following `kwargs` can modify the behavior:
     * `encode_chars` -- defaults to `True` if `False` then `disc_char_mat` is
        assumed to have already been translated to the appropriate integer
        codes (and each row of `disc_char_mat` must be a list).
     * `pattern_weights` -- a list of real weights for each pattern. If
        omitted, then each pattern is assumed to have a weight of 1
     * `n_cla_per_edge` -- defaults to 1, can be set higher for caching clas in
        multiple "directions" on the tree
     * `extra_clas` -- if sent augments the count of clas to be allocated
     * `rooted` -- create enough cla's for a rooted tree.
     * `n_clas` -- if sent overrides the number of cla objects calculated from
        the number of leaves (causes `rooted`, `n_cla_per_edge`, `extra_clas`
        to be ignored)
     """
    disc_char_type = model.char_type
    n_tips = len(disc_char_mat)
    n_clas = kwargs.get("n_clas")
    if n_clas is None:
        n_clas = n_tips - 2
        if kwargs.get("rooted"):
            n_clas += 1
        n_clas += kwargs.get("extra_clas", 0)
        if n_clas < 1:
            n_clas = 1
    if kwargs.get("encode_chars", True):
        enc_data = [disc_char_type.to_indices(row) for row in disc_char_mat]
    else:
        enc_data = disc_char_mat
    if not enc_data:
        return [], [], None
    n_sites = len(enc_data[0])
    for i in xrange(1, n_tips):
        if n_sites != len(enc_data[i]):
            m = "The number of characters differs between row 0 (%d) and row %d (%d." % (i, n_sites, len(enc_data[i]))
            raise ValueError(m)
    state_set_lookup = disc_char_type.get_cstate_set_lookup()
    num_states = disc_char_type.get_n_states()
    n_categ = model.n_prob_matrices
    ldobjs = [cleaf_data_ctor(i, state_set_lookup, n_categ) for i in enc_data]
    clas =  [cla_ctor(n_sites, num_states, n_categ) for i in range(n_clas)]
    full_la =  full_la_ctor(n_sites, num_states, n_categ)
    if n_pmats > 0:
        if n_pmats == 1:
            ldobjs = [EdgeCalcStruct(i, ProbMatrixArray(model=model), model=model) for i in ldobjs]
            clas = [EdgeCalcStruct(i, ProbMatrixArray(model=model), model=model) for i in clas]
        elif n_pmats == 2:
            ldobjs = [EdgeCalcStruct(i, ProbMatrixArray(model=model), model=model, alt_pmat_arr=ProbMatrixArray(model=model)) for i in ldobjs]
            clas = [EdgeCalcStruct(i, ProbMatrixArray(model=model), model=model, alt_pmat_arr=ProbMatrixArray(model=model)) for i in clas]
        else:
            raise ValueError("Currently a maximum 2 ProbMatArrays per edge can be requested")
    mp = list(model.get_model_probs())
    ml = model.model_list
    sf = [list(i.state_freqs) for i in ml]
    assert len(sf) == len(mp)
    full_la_set_freqs(full_la, mp, sf)
    return ldobjs, clas, full_la

'''

################################################################################
# cPhyProb is a package implementing some probability calculations used in
#   calculating likelihoods on phylogenies.
#
# Copyright (C) 2005-2007  Mark Holder mtholder@gmail.com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU  General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
# You should have received a copy of the GNU  General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
################################################################################
