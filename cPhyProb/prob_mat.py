#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Classes for models of the evolution of discrete characters.
"""
import copy
from itertools import izip

from cPhyProb.ccore.dsct_model import cpmat_array_ctor, calc_pmat_array, cpmat_array_get_mat_list


def slice_to_iterable(s):
    """Returns an xrange object that will iterate through the indices specified
    by the slice.
    """
    start = s.start or 0
    step = s.step
    if not step:
        return xrange(start, s.stop)
    return xrange(start, s.stop, step)

def normalize_getitem_key(key, n):
    """Returns a list of non-negative indices for a `key` on a sequence-like
    container with a size of `n`
    
    `key` can be an int, slice of ints, or an iterable of ints
    Exceptions are raised for incorrect types or out-of-range indices.
    Negative indices are converted to the non-negative equivalent.
    """
    if isinstance(key, slice):
        inds = slice_to_iterable(slice)
    elif isinstance(key, int):
        i = int(key)
        if i != key:
            raise TypeError("int or slice expected")
        inds = [i]
    else:
        inds = key # any iterable of ints 
    mat_inds = []
    for i in inds:
        if i < 0:
            i += n
        if i >= n  or i < 0:
            raise IndexError("index out of range")
        mat_inds.append(i)
    return mat_inds

class ProbMatrixArray(object):
    def __init__(self, **kwargs):
        """Creates a new array of transition probability matrices.
        
        `kwargs` must contain:
            *`model` can be passed in (which makes the model arg
                optional in calls to `calc_probs()`, or
        OR
            * `n_states` the number of states in each probability matrix.
          AND:
            * `n_mat` is the number of matrices, 
        """
        self.model = kwargs.get("model")
        self._n_mat = kwargs.get("n_mat")
        self._n_states = kwargs.get("n_states")
        if self.model is not None:
            n_mat = self.model.n_prob_matrices
            if self._n_mat is None:
                self._n_mat = n_mat
            elif self._n_mat != n_mat:
                raise ValueError("model.n_prob_matrices and n_mat should not differ.")
            mns = self.model.n_states
            if self._n_states is None:
                self._n_states = mns
            elif mns != self._n_states:
                raise ValueError("model.n_states and n_states should not differ.")
            mod_l = self.model.get_model_list()
        else:
            if self._n_mat is None:
                raise ValueError("model or n_mat must be specified.")
            if self._n_states is None:
                raise ValueError("A model or n_states must be specified.")
            mod_l = [None] * self._n_mat
        n_mat = self._n_mat
        self._mat_arr = cpmat_array_ctor(n_mat, self._n_states)
        self.calculated = [False] * n_mat
        self._last_out = None
        as_list = []
        for i in xrange(n_mat):
            p = ProbMatrix(pmat_array=self, pos=i, model=mod_l[i])
            as_list.append(p)
        self._prob_matrix_list = tuple(as_list)
        self.pos_list = range(n_mat)

    def get_n_mat(self):
        "returns the number of probability matrices in the ProbMatrixArray"
        return self._n_mat

    def get_n_states(self):
        """returns the number of states in each probability matrix in the 
        ProbMatrixArray"""
        return self._n_states

    def get_prob_matrix_tuple(self):
        """Returns the tuple of ProbMatrix that correspond to each matrix in the
        array.
        """
        return self._prob_matrix_list

    def __getitem__(self, key):
        """Returns the matrices with the indices indicated by the element(s) in
        `key`.
        """
        mat_inds = normalize_getitem_key(key, self.n_mat)
            
        return [ProbMatrix(pmat_array=self, pos=i) for i in mat_inds]

    def calc_probs(self, br_len, model=None, mat_n=None):
        """Calculates the transitition probability matrices for the branch
        length `br_len`
        
        If `self.model` is not None, then `model` must be sent (`model` will be 
            used instead of `self.model` if both are not None).
        `mat_n` can be an index of slice only some of the matrices should be 
            recalculated.  if `mat_n` is None, then all matrices are 
            recalculated.
        """
        if br_len < 0.0:
            raise ValueError("`br_len` cannot be negative")
        n_mat = self.n_mat
        if mat_n is not None:
            mat_inds = normalize_getitem_key(mat_n, n_mat)
            mat_inds.sort()
            if mat_inds == self.pos_list:
                mat_inds = []
        else:
            mat_inds = []

        m = model is None and self.model or model
        if m is None:
            mat_list = self.matrices
            model_list = []
            br_lens = []
            for i, mat in enumerate(mat_list):
                if mat.model is None:
                    if not mat_inds or i in mat_inds:
                       raise ValueError("`model` must be specified")
                model_list.extend(mat.model._get_c_model_list())
                br_lens.extend(m.get_adjusted_brlens_list(br_len))
        else:
            model_list = m._get_c_model_list()
            br_lens = m.get_adjusted_brlens_list(br_len)
        self._last_out = None
        if len(mat_inds) == 0:
            calc_pmat_array(self._mat_arr, model_list, br_lens)
            for i in range(n_mat):
                self.calculated[i] = True
        else:
            for i in mat_inds:
                calc_pmat(self._mat_arr, i, model_list[i], br_lens[i])
                self.calculated[i] = True
    def get_raw_prob_mat(self):
        return self._mat_arr
    def _get_last_output(self):
        """Returns the transition probabilities from the last calculation as
        a list.
        """
        if self._last_out is None:
            self._last_out = cpmat_array_get_mat_list(self._mat_arr)
        return self._last_out
    
    def probs_as_list(self, key=None):
        """ Returns the specified probability matrices as a 3D matrix (list of
        list of list of floats).
        
        `key` can be ommitted for all matrices or can be a valid index or slice.
        """
        s = key is None and slice(self.n_mat) or key
        inds = normalize_getitem_key(s, self.n_mat)
        all_mats = self._get_last_output()
        to_ret = []
        for i in inds:
            if not self.calculated[i]:
                raise ValueError("Transition probability matrix must be "\
                                 "calculated before it is requested.")
            to_ret.append(all_mats[i])
        return to_ret

    def get_pmat_array_and_pos(self):
        """Returns this ProbMatrixArray object the list of indices for the 
        element ProbMatrix objects.
        """
        return self.pmat_arr, self.pos_list

    def get_probs(self, br_len, model=None, mat_n=None):
        """Convenience function equivalent to calling `calc_probs` then
        `probs_as_list`.
        """
        self.calc_probs(br_len, model, mat_n)
        return self.probs_as_list(mat_n)
    raw_prob_mat = property(get_raw_prob_mat)
    n_states = property(get_n_states)
    n_mat = property(get_n_mat)
    matrices = property(get_prob_matrix_tuple)

class ProbMatrix(object):
    def __init__(self, **kwargs):
        """Creates a new Transition probability object (and an ProbMatrixArray
        if needed.
        
        `kwargs` must contain:
         * `pmat_array`: a reference to a ProbMatrixArray object, and
         * `pos`: (the position of this ProbMat in the array).
         Or:
         * `n_states`: the # of character states.  In this case a 
            ProbMatrixArray object of length 1 will be created.
        `kwargs` can also contain a `model` arg (which makes the model arg
            optional in calls to `calc_probs()`
        """
        self.model = kwargs.get("model")
        self.pmat_arr = kwargs.get("pmat_array")
        if self.pmat_arr is not None:
            self.pos = kwargs.get("pos")
            if self.pos is None:
                raise ValueError("`pos` expected")
            self._n_states = self.pmat_arr.n_states
            if self.model and self.model.n_states != self._n_states:
                raise ValueError("pmat_array.n_states and "\
                                 "model.n_states must agree.")
        else:
            self._n_states = kwargs.get("n_states")
            if self._n_states is None:
                if self.model is not None:
                    self._n_states = self.model.n_states
                else:
                    raise ValueError("`n_states` or `pmat_array` expected")
            elif self.model is not None:
                if self.model.n_state != self._n_states:
                    raise ValueError("n_states and model.n_states must agree.")
            self.pmat_arr = ProbMatrixArray(n_mat=1,
                                            n_states=self._n_states,
                                            model=self.model)
            self.pos = 0
        self.pos_list =  [self.pos]

    def get_pmat_array_and_pos(self):
        """Returns the ProbMatrixArray object that the matrix belongs to a the 
        position in that array as an list with an index in it. (other versions
        of this function may contain lists with several indices
        """
        return self.pmat_arr, self.pos_list

    def probs_as_list(self):
        """Returns the the probability matrix as a list of lists of floats.
        """
        return self.pmat_arr.probs_as_list(self.pos)[0]

    def __getitem__(self, key):
        "Returns the specified row(s) in the probability matrix"
        l = self.probs_as_list()
        return l[key]

    def calc_probs(self, br_len, model=None):
        m = model is None and self.model or None
        self.pmat_arr.calc_probs(br_len, m, self.pos)
    
    def get_probs(self, br_len, model=None):
        """Convenience function equivalent to calling `calc_probs` then
        `probs_as_list`.
        """
        self.calc_probs(br_len, model)
        return self.probs_as_list()

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
