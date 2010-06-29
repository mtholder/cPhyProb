#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests of among-site rate heterogeneity code"
import unittest
from cPhyProb.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.discrete_model import RateHetManager, RateHetType
class ASRVTest(unittest.TestCase):

    def test_bad(self):
        bad_rht = RateHetType.MAX_RATE_HET_TYPE + 1
        self.assertRaises(ValueError, RateHetManager, distrib_code=bad_rht)
        self.assertRaises(ValueError, RateHetManager, shape=0.0)
        self.assertRaises(ValueError, RateHetManager, shape=-10.0)
        self.assertRaises(ValueError, RateHetManager, n_cat=0)
        self.assertRaises(ValueError, RateHetManager, n_cat=1.1)
        RateHetManager(n_cat=1)
        self.assertRaises(ValueError, RateHetManager, n_cat="a")

    def test_mean(self):
        rhtype = RateHetType.GAMMA_EQ_CAT_MEAN
        rhm = RateHetManager(shape=0.5, n_cat=4, distrib_code=rhtype)
        expected = [0.03338775, 0.25191592, 0.82026848, 2.89442785]
        assert_list_eq(self, rhm.rates, expected)
        self.assertEqual(rhm.n_cat, 4)
        self.assertTrue(abs(rhm.shape-0.5) < 1.0e-5)

        rhm.set_shape(4.2)
        expected = [0.46502269, 0.78185321, 1.08432009, 1.66880400]
        assert_list_eq(self, rhm.rates, expected)
        self.assertEqual(rhm.n_cat, 4)
        self.assertTrue(abs(rhm.shape-4.2) < 1.0e-5)

        rhm = RateHetManager(shape=0.53, n_cat=10, distrib_code=rhtype)
        expected = [0.00680095, 0.04407109, 0.11615143, 0.22625412, 0.38203227,
                    0.59769745, 0.90010946, 1.34639481, 2.09367223, 4.28681618]
        assert_list_eq(self, rhm.rates, expected)
        self.assertEqual(rhm.n_cat, 10)
        self.assertTrue(abs(rhm.shape-0.53) < 1.0e-5)
        
    def test_median(self):
        rhtype = RateHetType.GAMMA_EQ_CAT_MEDIAN
        rhm = RateHetManager(shape=0.5, n_cat=4, distrib_code=rhtype)
        expected = [0.02907775, 0.28071453, 0.92477307, 2.76543465]
        assert_list_eq(self, rhm.rates, expected)

        rhm.set_shape(4.2)
        expected = [0.49655820, 0.79929783, 1.10263121, 1.60151276]
        assert_list_eq(self, rhm.rates, expected)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(ASRVTest)

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
