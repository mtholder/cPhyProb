#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests datatype internals"
import unittest
from cPhyProb.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.discrete_char_type import DiscreteCharType, DNAType
from cPhyProb.discrete_model import JukesCantor
from cPhyProb.phy_calc import get_tree_decorators
class TreeDecoratorTest(unittest.TestCase):
    def test_bad_dna(self):
        jc = JukesCantor()
        inds = ["AC", "A"]
        self.assertRaises(ValueError, get_tree_decorators, inds, jc, 1)

    def test_dna(self):
        jc = JukesCantor()
        inds = ["AC", "AG"]
        ldo, cla, f = get_tree_decorators(inds, jc, 1)
        self.assertEquals(len(ldo), 2)
        self.assertEquals(len(cla), 1)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(TreeDecoratorTest)

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
