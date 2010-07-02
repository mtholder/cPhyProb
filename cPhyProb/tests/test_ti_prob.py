#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
'''
"unit tests of among-site rate heterogeneity code"
import unittest
from cPhyProb.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from cPhyProb.prob_mat import ProbMatrixArray, ProbMatrix
from cPhyProb.discrete_model import RevDiscreteModel, JukesCantor, Kimura2Parameter
class ATest(unittest.TestCase):

    def test_jc(self):
        jc = JukesCantor()
        m = jc.get_probs(0.01)
        ti = tv = 0.003311999
        s = 1.0 - ti - 2.0*tv
        assert_mat_eq(self, m, [[s, tv, ti, tv],
                                [tv, s, tv, ti],
                                [ti, tv, s, tv],
                                [tv, ti, tv, s],])
        
    def test_kimura(self):
        k2p = Kimura2Parameter(1.0)
        m = k2p.get_probs(0.01)
        ti = tv = 0.003311999
        s = 1.0 - ti - 2.0*tv
        assert_mat_eq(self, m, [[s, tv, ti, tv],
                                [tv, s, tv, ti],
                                [ti, tv, s, tv],
                                [tv, ti, tv, s],])
        k2p.kappa = 10.0
        m = k2p.get_probs(0.01)
        ti = 0.008252
        tv = 0.000832
        s = 1.0 - ti - 2.0*tv
        assert_mat_eq(self, m, [[s, tv, ti, tv],
                                [tv, s, tv, ti],
                                [ti, tv, s, tv],
                                [tv, ti, tv, s],])
    def test_gtr(self):
        f = [0.23, 0.31, 0.21]
        s = sum(f)
        f.append(1.0-s)
        gtr = RevDiscreteModel(r_upper=[[1.0, 3.0, 1.2], [2.0, 4.1],[0.8]], 
                               equil_freq=f)
        m = gtr.get_probs(0.01)
        assert_mat_eq(self, m, [[0.99204347826086958, 0.0019969565217391304, 0.0040317391304347822, 0.001926086956521739], 
                                [0.0014816129032258064, 0.98928709677419346, 0.0026880645161290323, 0.0065419354838709672], 
                                [0.0044157142857142858, 0.0039680952380952384, 0.9903238095238095, 0.0012942857142857144], 
                                [0.0017719999999999999, 0.0081119999999999994, 0.0010872, 0.98902800000000002]])
    def test_pmat_alloc(self):
        jc = JukesCantor()
        pma = ProbMatrixArray(model=jc)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(ATest)

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
