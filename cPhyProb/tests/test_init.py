#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests of among-site rate heterogeneity code"
import unittest
from cPhyProb import get_logger, ComputationalResourceFlags
_LOG = get_logger(__name__)
from cPhyProb.ccore.dsct_model import get_num_computational_resources, get_resource_info

class IntegrationTest(unittest.TestCase):

    def test_new_api(self):
        r = get_num_computational_resources()
        _LOG.debug("%d resources" % r)
        self.assertTrue(r > 0)
        for i in range(r):
            t = get_resource_info(i)
            _LOG.debug(ComputationalResourceFlags.resource_descrip_to_str(t))
        self.assertRaises(IndexError, get_resource_info, r)
        self.assertRaises(IndexError, get_resource_info, -r-1)
        self.assertRaises(TypeError, get_resource_info, "0")
            

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(IntegrationTest)

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
