#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests datatyp internals"
import unittest
from cPhyProb.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.discrete_char_type import DiscreteCharType, DNAType
class DiscreteCharTypeTest(unittest.TestCase):
    def test_DNA(self):
        d = DNAType()

    def test_bad(self):
        self.assertRaises(ValueError, DiscreteCharType, "")
        self.assertRaises(ValueError, DiscreteCharType, "ACA")
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("A","A"),))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W","AT"),("W","AT")))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W","AT"),("Y","CK")))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W","AT"),("Y","CT")), missing="W")
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W","AT"),("Y","CT")), missing="?", aliases=[("A","A")])
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W","AT"),("Y","CT")), missing="?", aliases=[("a","a")], ignore_case=False)

    def test_dna(self):
        d = DNAType()
        inds = d.to_indices("ACNGT-WAYKCSXBDVMRH")
        self.assertEquals(inds, [0, 1, 4, 2, 3, 4, 8, 0, 6, 10, 1, 9, 4, 14, 13, 11, 7, 5, 12])
        self.assertEquals(d.num_states, 4)
        self.assertEquals(d.states, ('A', 'C', 'G', 'T'))
        self.assertEquals(d.all_symbols, ('A', 'C', 'G', 'T', 'N', 'R', 'Y', 'M', 'W', 'S', 'K', 'V', 'H', 'D', 'B'))
        self.assertEquals(d.ignore_case, True)
        self.assertEquals(d.symbol_to_ind, {'A': 0, 'C': 1, 'B': 14, 'D': 13, 'G': 2, 'H': 12, 'K': 10, 'M': 7, 'N': 4, 'S': 9, 'R': 5, 'T': 3, 'W': 8, 'V': 11, 'Y': 6, 'X': 4, '-': 4, '?': 4})
        self.assertEquals(d.state_code_lookup, ((0,), (1,), (2,), (3,), (0, 1, 2, 3), (0, 2), (1, 3), (0, 1), (0, 3), (1, 2), (2, 3), (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))
        self.assertTrue(d.state_code_lookup is not None)
        self.assertEquals(d.partial_ambiguity_indices, tuple(range(5,15)))

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(DiscreteCharTypeTest)

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
