#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Classes for models of the evolution of discrete characters.
"""
import copy
from itertools import izip
from cPhyProb import get_logger
from cPhyProb import Parameter, MutableFloatParameter, FloatParameter
_LOG = get_logger(__name__)
from cPhyProb.discrete_char_type import DNAType, CreateDNAType


PurePythonImpl = True
class PurePythonModelImpl(object):
    def __init__(self, num_states):
        self.num_state  = num_states
        self.q_mat = None
    def set_q_mat(self, q_mat):
        self.q_mat = q_mat


def _create_model_impl(num_states):
    if PurePythonImpl:
        return PurePythonModelImpl(num_states)
    from cPhyProb.ccore.dsct_model import cdsctm_ctor, cdsctm_set_q_mat
    return cdsctm_ctor(num_states)


if PurePythonImpl:
    def cdsctm_set_q_mat(model, q_mat):
        model.set_q_mat(q_mat)
    
class SiteModel(object):
    def __init__(self, sub_models, asrv=None, mixture_proportions=None):
        if isinstance(sub_models, list):
            self.sub_models = sub_models
        else:
            self.sub_models = [sub_models]
        self.asrv = None
        if len(self.sub_models) == 1:
            if mixture_proportions is not None:
                raise ValueError("mixture_proportions cannot be used when only one sub_model is specified")
        else:
            if not self.sub_models:
                raise ValueError("At least one `sub_model` must be specified")
            if asrv is not None:
                raise ValueError("Among site rate variation cannot be specified if a mixture of model is used (to accommodate both rate heterogeneity and mixtures, creat SiteModel objects with rate heterogeneity and mix add them as sub_models to a containing SiteModel)")
            if mixture_proportions is None:
                raise ValueError("mixture_proportions must be specified when more than one sub_model is used")
            
        # all sub_models must have the same number of states
        self._num_states = self.sub_models[0].num_states
        for sm in self.sub_models[1:]:
            assert(sm.num_states == self._num_states)
        
    def get_num_eigen_solutions(self):
        return sum([i.get_num_eigen_solutions() for i in self.sub_models])

    def get_num_states(self):
        return self._num_states
    
    def get_num_categories(self):
        t = sum([i.num_categories for i in self.sub_models])
        if self.asrv:
            return self.asrv.num_categories * t
        return t

    def get_num_prob_models(self):
        return sum([i.get_num_prob_models() for i in self.sub_models])

    get_num_rate_categories = get_num_categories
    
    num_rate_categories = property(get_num_rate_categories)
    
    num_states = property(get_num_states)
    num_prob_models = property(get_num_prob_models)
    
class RevDiscreteModel(object):
    """Calculates a transition probability matrix for a discrete state 
    character when given a branch length.
    """
    # pylint: disable-msg=R0903
    def __init__(self, r_mat=None, r_upper=None, equil_freq=None, char_type=None, params=None):
        """Initializes a model. 
        
        `r_mat` must be a square matrix (list of lists of floats).  Diagonal
            elements will be ignored. The columns of this matrix will be 
            multiplied by their corresponding equilibirium frequencies.
        `r_upper` can be used instead of `r_mat.` It should correspond to the
            upper triangle of the R-matrix.
        `equil_freq` must be None a list of the equilibrium frequencies.
        """
        self.char_type = char_type
        if r_upper:
            if r_mat:
                raise ValueError("r_mat or r_upper cannot both be specified")
            r_mat = _r_upper_to_r_mat(r_upper) 
        elif not r_mat:
            raise ValueError("Either r_mat or r_upper must be given")
        self._num_states = len(r_mat)
        if self._num_states < 2:
            raise ValueError("The number of states must be > 1.")
        if self.char_type:
            if self.char_type.num_states != self._num_states:
                raise ValueError("The number of states in the char_type must agree with other args")
        if r_mat:
            self._verify_r_mat(r_mat)
        else:
            raise ValueError("Either r_mat or r_upper must be specified")
        if equil_freq:
            self._verify_eq_freq_len(equil_freq)
        self._dsctm = _create_model_impl(self._num_states)

        priv_mat = []
        priv_q_mat = []
        for row in r_mat:
            priv_row = list(row)
            for n, cell in enumerate(row):
                if not isinstance(cell, Parameter):
                    priv_row[n] = FloatParameter(cell)
            priv_mat.append(tuple(priv_row))
            priv_q_mat.append([float(i) for i in priv_row])
        self._r_mat = tuple(priv_mat)
        self._q_mat = priv_q_mat
        
        equal_freq = 1.0/self._num_states
        self._state_freqs = [equal_freq] * self._num_states
        if equil_freq:
            self.set_state_freqs(equil_freq)
        else:
            self._recalc_q_mat()
        self._ti_probs = None
        if params:
            self._parameters = tuple(params)
        else:    
            self._parameters = ()
        for i in self._parameters:
            i.add_listener(self)
    def get_num_eigen_solutions(self):
        return 1
    def get_num_states(self):
        return self._num_states
    def get_num_categories(self):
        return 1
    def get_num_prob_models(self):
        return 1

    num_categories = property(get_num_categories)
    num_states = property(get_num_states)
    num_prob_models = property(get_num_prob_models)

    def get_all_submodels(self):
        return [self]
        
    def _get_c_model_list(self):
        """Returns a list to references to the (opaque) C-struct that this 
        object wraps.
        """
        return [self._dsctm]

    def _create_prob_matrix_obj(self):
        """Factory for ProbMatrix or ProbMatrixArray object with `self` as its 
        model."""
        return ProbMatrix(model=self)

    def get_probs(self, br_len):
        """Returns the transition probability matrix as a list of list of floats
        for `br_len.`
        """
        if self._ti_probs is None:
            self._ti_probs = self._create_prob_matrix_obj()
        return self._ti_probs.get_probs(br_len)

    def get_q_mat(self):
        "Returns a deep copy of the values in the Q-matrix."
        q = []
        for row in self._q_mat:
            q.append(list(row))
        return q

    def get_r_mat(self):
        "Accessor for R-matrix."
        return self._r_mat

    def get_r_upper(self):
        "Accessor for R-matrix."
        return _r_mat_to_r_upper(self._r_mat)
    
    def get_state_freqs(self):
        "Accessor for state frequencies."
        return self._state_freqs

    def set_r_upper(self, r_upper):
        "Sets the R-Matrix based on `r_upper` -- the upper triangle."
        r_mat = _r_upper_to_r_mat(r_upper)
        self.set_r_mat(r_mat)

    def set_r_mat(self, r_mat):
        """`r_mat` must be a square matrix (list of lists of floats).  Diagonal
        elements will be ignored. The columns of this matrix will be multiplied
        by their corresponding equilibirium frequencies.
        """
        self._verify_r_mat(r_mat)
        for row_n, row in enumerate(r_mat):
            internal_row = self._r_mat[row_n]
            for col_n, cell in enumerate(row):
                internal_cell = internal_row[col_n]
                print "repr(internal_cell) =", repr(internal_cell)
                if float(cell) != internal_cell.value:
                    internal_cell.value = cell
        self._recalc_q_mat()
        
    def set_state_freqs(self, equil_freq):
        """Replaces the equilibrium state frequencies.
        
        `equil_freq` must be a list of non-negative floats with length equal to 
            (or one less than) the number of states. If less, the last element
            is obtained by subtraction. The frequencies must sum to 1.0 (within
            rounding error).
        """
        self._verify_eq_freq_len(equil_freq)
        tmp = copy.deepcopy(self._state_freqs)
        try:
            sum_f = 0.0
            for i, val in enumerate(equil_freq):
                if val < 0.0:
                    if val < -1.0e-06: 
                        raise ValueError("State frequencies cannot be negative")
                    val = 0.0
                prev_sum_f = sum_f
                sum_f += val
                if sum_f > 1.0:
                    if (i + 1) != self._num_states or sum_f > 1.00001:
                        raise ValueError("Sum of state frequencies cannot "\
                                         "exceed 1.0")
                    val = 1.0 - prev_sum_f
                self._state_freqs[i] = float(val)
            if len(equil_freq) == (self._num_states + 1):
                self._state_freqs[-1] = 1.0 - sum_f
        except:
            self._state_freqs = tmp
            raise
        self._recalc_q_mat()
    
    def _recalc_q_mat(self):
        """Uses self._state_freqs and self._r_mat to refresh self._q_mat.
        
        As required by the current version of cdsctm_set_q_mat, the Q-matrix
            is rescaled such that each row sums to 0.0 and the weighted 
            average of the diagonal elements is -1.0 (with the weight being the
            associated equilibrium frequency)."""
        _LOG.debug("in self._recalc_q_mat")
        w_mat_sum = 0.0
        for i, rows in enumerate(izip(self._r_mat, self._q_mat)):
            r_row, q_row = rows
            row_sum = 0.0
            for j, to_state in enumerate(izip(r_row, self._state_freqs)):
                r_val, stf = [float(p) for p in to_state]
                if i != j:
                    v = r_val * stf
                    assert(r_val >= 0.0)
                    assert(stf >= 0.0)
                    q_row[j] = v
                    row_sum += v
            q_row[i] = -row_sum
            w_mat_sum += row_sum*self._state_freqs[i]
        for q_row in self._q_mat:
            for i in xrange(len(q_row)):
                q_row[i] /= w_mat_sum
        _LOG.debug("calling cdsctm_set_q_mat...")
        cdsctm_set_q_mat(self._dsctm, self._q_mat)
        _LOG.debug("back from cdsctm_set_q_mat...")
    
    def _verify_eq_freq_len(self, equil_freq):
        "Raises a ValueError if `equil_freq` is not the correct length"
        nef = len(equil_freq)
        if nef != self._num_states:
            if nef != (1 + self._num_states):
                raise ValueError("The number of states in the r_mat and "\
                             "equil_freq must agree.")

    def _verify_r_mat(self, r_mat):
        """Raises a ValueError if `r_mat` is not the correct shape or has 
        negative values (other than on the diagonal).
        """
        for i, row in enumerate(r_mat):
            if len(row) != self.num_states:
                raise ValueError("The R-matrix must be square.")
            for j, val in enumerate(row):
                if i != j:
                    if val < 0.0:
                        raise ValueError("Off-diagonal elements of the "\
                                         "R-matrix cannot be negative.")
                    if abs(val - float(r_mat[j][i])) > 1.0e-6:
                        raise ValueError("R-matrix must be symmetric")
    def get_num_states(self):
        "Returns the number of states"
        return self._num_states

    def get_num_prob_matrices(self):
        return  1

    def get_adjusted_brlens_list(self, b):
        return [b]

    def get_model_list(self):
        return [self]

    def get_model_probs(self):
        return (1.0,)

    def get_parameters(self):
        return self._parameters

    model_probs = property(get_model_probs)
    num_prob_matrices = property(get_num_prob_matrices)
    q_mat = property(get_q_mat)
    r_mat = property(get_r_mat, set_r_mat)
    r_upper = property(get_r_upper, set_r_upper)
    state_freqs = property(get_state_freqs, set_state_freqs)
    model_list = property(get_model_list)


if False:
    def Kimura2ParameterModel(kappa):
        dna = CreateDNAType()
        if kappa is None:
            kappa = MutableFloatParameter(1.0)
        elif not isinstance(kappa, Parameter):
            kappa = MutableFloatParameter(kappa)
        return RevDiscreteModel(r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna, params=[kappa])

class Kimura2ParameterModel(RevDiscreteModel):
    """A RevDiscreteModel instance with the Kimura (1980) parameters."""
    def __init__(self, kappa=None):
        dna = DNAType()
        if kappa is None:
            kappa = 1.0
        if not isinstance(kappa, Parameter):
            kappa = FloatParameter(kappa)
        self._kappa = kappa
        _LOG.debug("Created dna type in Kimura2ParameterModel")
        RevDiscreteModel.__init__(self, r_upper=[[1.0, kappa, 1.0], [1.0, kappa], [1.0]], char_type=dna, params=[kappa])
        _LOG.debug("called RevDiscreteModel.__init__")
        self.set_kappa(kappa)   
        _LOG.debug("called Kimura2ParameterModel.set_kappa called")
    def get_kappa(self):
        return self._kappa
    def set_kappa(self, k):
        f = float(k)
        if f < 0.0:
            raise ValueError("kappa must be greater than 0.0")
        self._kappa.set_value = k
        r_upper = [[1.0, k, 1.0], [1.0, k], [1.0]]
        self.r_upper = r_upper
    kappa = property(get_kappa, set_kappa)


class JukesCantorModel(RevDiscreteModel):
    """A RevDiscreteModel instance with no free parameters"""
    def __init__(self, kappa=None):
        dna = DNAType()
        _LOG.debug("Created dna type in JukesCantorModel")
        RevDiscreteModel.__init__(self, r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna, params=[])
        _LOG.debug("called RevDiscreteModel.__init__")




def _r_upper_to_r_mat(r_upper):
    """Convert the upper triangle of a symmetric matrix to the full matrix with 
    0.0 on the diagonal.
    """
    len_upper = len(r_upper)
    r_mat = []
    blank = [0.0]
    for i in xrange(len_upper):
        next_row = blank + r_upper[i]
        blank.append(0.0)
        for j in xrange(i):
            next_row[j] = r_upper[j][i-j-1]
        r_mat.append(next_row)
    for j in xrange(len_upper):
        blank[j] = r_upper[j][len_upper-j-1]
    r_mat.append(blank)
    return r_mat

def _r_mat_to_r_upper(r_mat):
    """Returns the upper triangle (list of lists) of a square matrix `r_mat`"""
    r_upper = []
    for i in xrange(len(r_mat)):
        r_upper.append(r_mat[i][i+1:])
    return r_upper



oldcode='''
from cPhyProb.prob_mat import ProbMatrixArray, ProbMatrix
from cPhyProb.ccore.dsct_model import cdsctm_ctor, cdsctm_set_q_mat


class SiteModel(object):
    def __init__(self, submodels, model_probs, char_type, ):
        self._submodels = submodels
        self._model_probs = model_probs
        self.char_type = char_type

    def get_model_probs(self):
        return self._model_probs

    def set_model_probs(self, m):
        self._model_probs = m

    def get_model_list(self):
        s = []
        for i in self._submodels:
            s.extend(i.get_model_list())
        return s
        
    def get_num_prob_matrices(self):
        return sum([i.get_num_prob_matrices() in self._submodels])

    def get_adjusted_brlens_list(self, b):
        s = []
        for i in self._submodels:
            s.extend(i.get_adjusted_brlens_list(b))
        return s

    def get_all_submodels(self):
        s = []
        for i in self._submodels:
            s.extend(i.get_all_submodels())
        return s

    def _get_c_model_list(self):
        s = []
        for i in self._submodels:
            s.extend(i._get_c_model_list())
        return s

    def _create_prob_matrix_obj(self):
        return ProbMatrixArray(model=self)
    model_probs = property(get_model_probs, set_model_probs)
    model_list = property(get_model_list)

def JukesCantorModel():
    """Returns a RevDiscreteModel instance with the Jukes-Cantor (1969)
    parameters.
    """
    dna = DNAType()
    return RevDiscreteModel(r_upper = [[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna)





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
'''
