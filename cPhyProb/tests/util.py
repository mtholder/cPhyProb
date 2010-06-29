#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
from itertools import izip
import unittest
def assert_mat_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two matrices."
    for ret_row, exp_row in izip(returned, expected):
        assert_list_eq(self, ret_row, exp_row)

def assert_list_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two lists."
    for ret_val, exp_val in izip(returned, expected):
        self.assertAlmostEqual(ret_val, exp_val, places=5)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestSuite([])

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
