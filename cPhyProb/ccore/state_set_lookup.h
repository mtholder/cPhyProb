/**
 * Datastructure needed for encoding of characters into indices in C. 
 * Does not use BEAGLE.
 */
#if ! defined(STATE_SET_LOOKUP_H)
#define STATE_SET_LOOKUP_H
#ifdef __cplusplus
extern "C" 
{
#endif




#include "cphyprob_defs.h"

typedef struct {
	PyObject_HEAD
	int ** state_lookup; /* maps an state_set_index to the array indicating the # of states and then listing them */
	unsigned nStates; /*the number of "fundamental" states.*/
	unsigned n_state_sets; /*the number of combination of states*/
} StateSetLookupStruct;



StateSetLookupStruct* state_set_lookup_new(unsigned nStates, unsigned n_state_sets, int ** p);
void state_set_lookup_dtor(StateSetLookupStruct* dsct_model);



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

#include "asrv.h"
