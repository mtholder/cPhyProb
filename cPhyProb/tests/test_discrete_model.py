#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests of RevDiscreteModel"
import unittest
from cPhyProb.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.discrete_model import RevDiscreteModel, _r_upper_to_r_mat, \
    _r_mat_to_r_upper

class TRDSModelTest(unittest.TestCase):
    def test_bad(self):
        self.assertRaises(ValueError, RevDiscreteModel, [])
        self.assertRaises(ValueError, RevDiscreteModel, [[0,1],[.5,0]])
    def test_init(self):
        a = RevDiscreteModel([["ignored", 2.],[2.,""]])
        assert_mat_eq(self, a.q_mat, [[-1.0, 1.0],[1.0,-1.0]])

    def test_r_upper_to_r_mat(self):
        cases = [([[1.0]],
                  [[0.0, 1.0],
                   [1.0,0.0]]),
                 ([[1.0,1.0],
                   [1.0]],
                  [[0.0, 1.0, 1.0],
                   [1.0, 0.0, 1.0],
                   [1.0, 1.0, 0.0]]),
                 ([[1.0,2.0],
                   [3.0]],
                  [[0.0, 1.0, 2.0],
                   [1.0, 0.0, 3.0],
                   [2.0, 3.0, 0.0]]),
                ]
        for t_case in cases:
            rup, rmat = t_case
            r_mat_gen = _r_upper_to_r_mat(rup)
            r_up_gen = _r_mat_to_r_upper(rmat)
            assert_mat_eq(self, r_mat_gen, rmat)
            assert_mat_eq(self, r_up_gen, rup)
    def test_set_q_mat(self):
        a = RevDiscreteModel(r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0],])
        oth = 1.0/3.0
        assert_list_eq(self, a.state_freqs, [0.25]*4)
        assert_mat_eq(self, a.q_mat, [[-1.0, oth, oth, oth],
                                       [oth, -1.0, oth, oth],
                                       [oth, oth, -1.0, oth],
                                       [oth, oth, oth, -1.0]])
        a.r_upper = [[1.0, 2.0, 1.0], [1.0, 2.0], [1.0],] # kimura kappa=2
        assert_list_eq(self, a.state_freqs, [0.25]*4)
        assert_mat_eq(self, a.q_mat, [[-1.0, .25, .5, .25],
                                      [.25, -1.0, .25, .5],
                                      [.5, .25, -1.0, .25],
                                      [.25, .5, .25, -1.0]])
        a.state_freqs = [0.2, 0.3, 0.15, 0.35]
        assert_list_eq(self, a.state_freqs, [0.2, 0.3, 0.15, 0.35])
        exp = [[-0.9547738, 0.301507, 0.3015075, 0.351758], 
               [0.20100502, -1.05527638, 0.15075376, 0.7035175], 
               [0.40201005, 0.301507537, -1.05527638, 0.35175879], 
               [0.2010050, 0.603015075, 0.150753768, -0.954773869]]
        assert_mat_eq(self, a.q_mat, exp)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(TRDSModelTest)

# pylint: disable-msg=C0103
def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
    """
    return additional_tests()

if __name__ == "__main__":
    unittest.main()

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
