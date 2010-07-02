#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests of among-site rate heterogeneity code"
import unittest
from cPhyProb.tests.util import *
from cPhyProb.asrv import RateHetManager, RateHetType, GammaRateHetManager
from cPhyProb.discrete_model import Kimura2Parameter, SiteModel
from cPhyProb.phy_calc import PhyloSubsetLikelihoodCalculator, PhyloLikelihoodCalculator
from cPhyProb.discrete_char_type import DNAType

class IntegrationTest(unittest.TestCase):

    def test_new_api(self):
        char_mat_inp = ["ACGT", "AGGT"]
        kimura_model = Kimura2Parameter()
        site_model = SiteModel(kimura_model)
        subsetCalc = PhyloSubsetLikelihoodCalculator(char_matrix=char_mat_inp,
                            discrete_char_type=DNAType(),
                            convert_partial_ambig_to_missing=False,
                            num_partial_storage_per_internal_edge=2,
                            num_extra_partial_storage=2,
                            max_num_internal_edges=3,
                            model=site_model,
                            num_eigen_storage = 2*site_model.get_num_eigen_solutions(),
                            num_prob_mat_per_internal_edge=1,
                            num_prob_mat_per_terminal_edge=1,
                            num_extra_prob_mat_storage=2,
                            rescale_every=10
                            )
        likelihoodCalc = PhyloLikelihoodCalculator(subsets=[subsetCalc])
        likelihoodCalc.activate()
    def x():

        calcContext.activate()
        root = Node(children=[Node(index=n, char_vector=i, edge_length = 0.1) for n, i in enumerate(calcContext.char_matrix)])
        calcContext.decorate_nodes(iter(root), attribute="likelihood_blob", attr_as_list=False)
        kappa = Kimura.parameters()[0]
        kappa.value = 2.0
        like1 = root.likelihood_blob.get_likelihood()
       
        kappa.value = 3.0
        like2 = root.likelihood_blob.get_likelihood()
        cache_receipt2 = calcContext.create_cache(model=True, partials=True, prob_mat=True)
        kappa.value = 4.0
        like3 = root.likelihood_blob.get_likelihood()
        calcContext.restore_snapshot2(cache_receipt2)
        like2again = root.likelihood_blob.get_likelihood()
        self.assertEquals(like2, like2again)
        self.assertEquals(kappa.value, 3.0)
        calcContext.flush_caches(cache_receipt2)

        root.child[0].likelihood_blob.edge_length = 0.2
        calcContext.set_edge_length(root.child[0], 0.2)
        
        x = """def invalidate(curr_nd, toward_root_nd):
	neighbors = curr_nd.adjacent_nodes()
	if i in neighbors:
		if i is not toward_root_nd:
			p = i.partials(curr_nd):
			if p:
				invalidate(i, curr_nd)


	for i in nd.partials(toward_root_nd):
		i.is_dirty = True"""
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
