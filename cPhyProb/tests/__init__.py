#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Test cases and suites for cProbPhy
"""

__all__ = [
    "util",
    "test_asrv",
    "test_discrete_char_type",
    #"test_discrete_model",
    #"test_likelihood",
    #"test_ti_prob",
    #"test_tree_decorators",
    ]

import unittest
# pylint: disable-msg=C0111,W0401,W0611

def even_more_tests(all_suites):
    "finds test from module introspection"
    return
    for i in __all__:
        module = __import__("cPhyProb.tests.%s" % i)
        print module
        tests_mod = getattr(module, "tests")
        sub_test_mod = getattr(tests_mod, i)
        suite = sub_test_mod.additional_tests()
        all_suites.append(suite)

def additional_tests():
    """Creates a unittest.TestSuite from all of the modules in `cPhyProb.tests`
    
    \todo uncommenting even_more_tests line results in test from "setup.py test"
        being run 3 times each.  I don't know why. (even with it commented out
        they are being run twice
    """
    all_suites = []
    even_more_tests(all_suites)
    return unittest.TestSuite(all_suites)

def test_all():
    "Runs all of the unittests in `cPhyProb.tests`"
    runner = unittest.TextTestRunner()
    runner.run(additional_tests())
    
if __name__ == "__main__":
    test_all()

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
