#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests datatyp internals"
import unittest
from cPhyProb.tests.util import *
x='''
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.phy_calc import get_tree_decorators, calc_cla, calc_lnL
from cPhyProb.discrete_model import RevDiscreteModel, _r_upper_to_r_mat, \
    _r_mat_to_r_upper, JukesCantor
from cPhyProb.prob_mat import ProbMatrixArray
from cPhyProb.discrete_char_type import DNAType

class LikelihoodTest(unittest.TestCase):
    def test_jc_three_tax(self):
        mod = JukesCantor()
        seq_data = ["ACGT", "ACGC", "ACGT"]
        leaf_data, clas, full_la = get_tree_decorators(seq_data, mod)
        leaf0, leaf1, leaf2 = leaf_data
        internal_d = clas[0]
        leaf0.set_brlen(0.1)
        leaf1.set_brlen(0.07)
        leaf2.set_brlen(0.05)
        lnL = calc_lnL(full_la, internal_d, leaves=leaf_data, internals=(), beg_subset=0, end_subset=1, beg_cat=0, end_cat=1)
        self.assertTrue(abs(10.12296 + lnL) < 1e-05)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(LikelihoodTest)

# pylint: disable-msg=C0103
def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
    """
    return additional_tests()

if __name__ == "__main__":
    unittest.main()
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
