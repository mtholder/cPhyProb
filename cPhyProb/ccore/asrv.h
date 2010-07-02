/**
 * Among-site rate variation calculation in C.  Does not use BEAGLE.
 */
#if ! defined(ASRV_H)
#define ASRV_H
#ifdef __cplusplus
extern "C" 
{
#endif




#include "cphyprob_defs.h"


/**
 * stores an array of doubles and an array of frequencies representing the
 * a probability for each element in the array.
 * An example usage is for storing an array of rates for n_categories and the
 * probability that a site would belong to each of the rate categories.
 */
typedef struct {
	PyObject_HEAD
	unsigned n;
	int style; /* enum facet for to indicate type -- this is a hack
				  0 ="use median for rate categ"
				  1 ="use mean for rate categ"
			   */
	double param; /*shape param of the gamma distribution -- this is a hack*/
	double * val;
	double * freq;
} ASRVObj;



ASRVObj* asrv_obj_new(unsigned dim, int style, double param);
void asrv_obj_dtor(ASRVObj* asrh);



void internal_asrv_set_shape(ASRVObj *asrh, double val);






#ifdef __cplusplus
}
/* end of extern c bit */
#endif






#endif
/*
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
*/


