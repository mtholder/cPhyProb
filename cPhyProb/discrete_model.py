#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Classes for models of the evolution of discrete characters.
"""
import copy
from itertools import izip
from cPhyProb.prob_mat import ProbMatrixArray, ProbMatrix
from cPhyProb.discrete_char_type import DNAType
from cPhyProb.ccore.dsct_model import cdsctm_ctor, cdsctm_get_dimen, \
    cdsctm_set_q_mat, casrvo_ctor, casrvo_set_shape,\
    casrvo_get_rates, casrvo_get_shape, cpmat_array_ctor,\
    cpmat_array_get_mat_list, calc_pmat_array, calc_pmat

_MIN_GAMMA_SHAPE = 1.0e-07
_LOW_RATE_MIN_GAMMA_SHAPE = 1.0e-200

class RevDiscreteModel(object):
    """Calculates a transition probability matrix for a discrete state 
    character when given a branch length.
    """
    # pylint: disable-msg=R0903
    def __init__(self, r_mat=None, r_upper=None, equil_freq=None, char_type=None):
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
        self._n_states = len(r_mat)
        if self._n_states < 2:
            raise ValueError("The number of states must be > 1.")
        if self.char_type:
            if self.char_type.n_states != self._n_states:
                raise ValueError("The number of states in the char_type must agree with other args")
        if r_mat:
            self._verify_r_mat(r_mat)
        else:
            raise ValueError("Either r_mat or r_upper must be specified")
        if equil_freq:
            self._verify_eq_freq_len(equil_freq)
        self._dsctm = cdsctm_ctor(self._n_states)
        self._r_mat = copy.deepcopy(r_mat)
        equal_freq = 1.0/self._n_states
        self._state_freqs = [equal_freq] * self._n_states
        self._q_mat = copy.deepcopy(self._r_mat)
        if equil_freq:
            self.set_state_freqs(equil_freq)
        else:
            self._recalc_q_mat()
        self._ti_probs = None

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
        "Accessor for Q-matrix."
        return self._q_mat

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
        self._r_mat = copy.deepcopy(r_mat)
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
                    if (i + 1) != self._n_states or sum_f > 1.00001:
                        raise ValueError("Sum of state frequencies cannot "\
                                         "exceed 1.0")
                    val = 1.0 - prev_sum_f
                self._state_freqs[i] = float(val)
            if len(equil_freq) == (self._n_states + 1):
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
        w_mat_sum = 0.0
        for i, rows in enumerate(izip(self._r_mat, self._q_mat)):
            r_row, q_row = rows
            row_sum = 0.0
            for j, to_state in enumerate(izip(r_row, self._state_freqs)):
                r_val, stf = to_state
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
        cdsctm_set_q_mat(self._dsctm, self._q_mat)
    
    def _verify_eq_freq_len(self, equil_freq):
        "Raises a ValueError if `equil_freq` is not the correct length"
        nef = len(equil_freq)
        if nef != self._n_states:
            if nef != (1 + self._n_states):
                raise ValueError("The number of states in the r_mat and "\
                             "equil_freq must agree.")

    def _verify_r_mat(self, r_mat):
        """Raises a ValueError if `r_mat` is not the correct shape or has 
        negative values (other than on the diagonal).
        """
        for i, row in enumerate(r_mat):
            if len(row) != self.n_states:
                raise ValueError("The R-matrix must be square.")
            for j, val in enumerate(row):
                if i != j:
                    if val < 0.0:
                        raise ValueError("Off-diagonal elements of the "\
                                         "R-matrix cannot be negative.")
                    if abs(val - float(r_mat[j][i])) > 1.0e-6:
                        raise ValueError("R-matrix must be symmetric")
    def get_n_states(self):
        "Returns the number of states"
        return self._n_states

    def get_n_prob_matrices(self):
        return  1

    def get_adjusted_brlens_list(self, b):
        return [b]

    def get_model_list(self):
        return [self]

    def get_model_probs(self):
        return (1.0,)

    model_probs = property(get_model_probs)
    n_prob_matrices = property(get_n_prob_matrices)
    n_states = property(get_n_states)
    q_mat = property(get_q_mat)
    r_mat = property(get_r_mat, set_r_mat)
    r_upper = property(get_r_upper, set_r_upper)
    state_freqs = property(get_state_freqs, set_state_freqs)
    model_list = property(get_model_list)

class RevDiscreteMixtureModel(RevDiscreteModel):
    def __init__(self, submodels, model_probs, char_type):
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
        
    def get_n_prob_matrices(self):
        return sum([i.get_n_prob_matrices() in self._submodels])

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

def JukesCantor():
    """Returns a RevDiscreteModel instance with the Jukes-Cantor (1969)
    parameters.
    """
    dna = DNAType()
    return RevDiscreteModel(r_upper = [[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna)

class Kimura2Parameter(RevDiscreteModel):
    """A RevDiscreteModel instance with the Kimura (1980) parameters."""
    def __init__(self, kappa):
        dna = DNAType()
        RevDiscreteModel.__init__(self, r_upper = [[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna)
        self.set_kappa(kappa)
    def get_kappa(self):
        return self._kappa
    def set_kappa(self, k):
        f = float(k)
        if f < 0.0:
            raise ValueError("kappa must be greater than 0.0")
        self._kappa = k
        r_upper = [[1.0, k, 1.0], [1.0, k], [1.0]]
        self.r_upper = r_upper
    kappa = property(get_kappa, set_kappa)

class RateHetType:
    "Enumeration of the kinds of among-site rate heterogeneity."
    # pylint: disable-msg=R0903,W0232
    GAMMA_EQ_CAT_MEDIAN = 0 # keep this 0 (or coordinate with the ccore.ASRVObj
    GAMMA_EQ_CAT_MEAN = 1
    MAX_RATE_HET_TYPE = 1

class RateHetManager(object):
    "Calculates the rate multipliers for different rate categories."
    def __init__(self, 
                 shape=0.5, 
                 n_cat=4, 
                 distrib_code=RateHetType.GAMMA_EQ_CAT_MEAN):
        """Currently just supports the discrete approximation of 
        gamma-distributed rates.
        """
        if distrib_code > RateHetType.MAX_RATE_HET_TYPE:
            raise ValueError("Illegal value for distrib_code")
        if n_cat < 1 or int(n_cat) != n_cat:
            raise ValueError("n_cat must be a positive integer")
        if shape <= 0.0:
            raise ValueError("shape must be greater than 0.0")
        self._distrib_code = distrib_code
        self._n_cat = n_cat
        self._shape = shape
        s = max(_MIN_GAMMA_SHAPE, shape)
        self._asrv = casrvo_ctor(s, n_cat, distrib_code)

    def get_n_cat(self):
        "Returns the number of rate categories."
        return self._n_cat

    def get_rates(self):
        "Returns a list with a rate for each category."
        if self._shape > _MIN_GAMMA_SHAPE:
            return casrvo_get_rates(self._asrv)
        rate_approx = [_LOW_RATE_MIN_GAMMA_SHAPE ] * self._n_cat
        rate_approx[self._n_cat - 1] = float(self._n_cat)
        return rate_approx

    def get_shape(self):
        "Returns the shape parameter for a Gamma distribution over rates."
        return self._shape
        #return casrvo_get_shape(self._asrv)

    def set_shape(self, shape):
        "Sets the shape parameter for a Gamma distribution over rates."
        self._shape = shape
        if self._shape > _MIN_GAMMA_SHAPE:
            casrvo_set_shape(self._asrv, shape)

    n_cat = property(get_n_cat)
    rates = property(get_rates)
    shape = property(get_shape, set_shape)


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
