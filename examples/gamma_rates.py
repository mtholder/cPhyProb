#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

import sys
from optparse import OptionParser
from cPhyProb.discrete_model import RateHetManager


parser = OptionParser()
parser.add_option("-s", "--shape", dest="shape", default=0.5,
                  type="float",
                  help="Shape of the Gamma distribution (mean is 1.0).")
parser.add_option("-n", "--ncat", dest="n_cat", default=4,
                  type="int",
                  help="The number of categories")
(options, args) = parser.parse_args()

fmt_s = "%2d   %-8g"
if options.shape <= 0.0:
    sys.exit("Gamma shape parameter must be greater than 0")
if options.n_cat < 2:
    print(fmt_s % (1, 1.0))
    sys.exit(0)

rhm = RateHetManager(shape=options.shape, n_cat=options.n_cat)
rates = rhm.rates
for i in range(options.n_cat):
    print(fmt_s % (1 + i, rates[i]))

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
