/**
 *
 * Public interface of cPhyprob
 */
#if ! defined(C_PHY_PROB_H)
#define C_PHY_PROB_H

#ifdef __cplusplus
extern "C"
{
#endif


#include "cphyprob_defs.h"
#include "cphyprob_inlines.c"
#include "asrv.h"
#include "state_set_lookup.h"
#include "dsct_model.h"


struct DSCTModelObj;

/* Each LikeCalculatorInstance contains the internal storage needed to
    calculate likelihoods for one subset in a (possibly) partitioned matrix.

    Multiple Models may be associate with the instance, but all model share the
    same number of states.

*/
typedef struct LikeCalculatorInstance {
	unsigned int numLeaves;
	unsigned long numPatterns;
	unsigned int num_states;
	unsigned int num_rate_categories;
	unsigned int num_state_code_arrays;
	unsigned int numPartialStructs;
	unsigned int num_prob_models;
	unsigned int numProbMats;
	unsigned int numEigenStorage;
	unsigned int num_rescalings_multipliers;
	int resource_arg;
	long resource_flag;

#	if defined(USE_BEAGLE_LIB)
		int beagleInstanceIndex; /* handle of the beagle instance */
#	else
		double *** probMatArray; /* numProbMats * num_states * num_states */
		double *** probMatAliasArray; /* a vector of length numProbMats which aliases elements of probMatArray */
#	endif
	DSCTModelObj ** probModelArray;
	EigenSolutionStruct ** eigenSolutionStructs;

} LikeCalculatorInstance;






#if 0
// The CLA layout is [sites] x [rates] x [states]
//	------------------------------------------
//	4 bases, 3 sites, 5 rate categories
//
//	          11111111112222222222333333333344444444445555555555
//	012345678901234567890123456789012345678901234567890123456789
//	------------------------------------------------------------
//	012301230123012301230123012301230123012301230123012301230123 <== base b
//	000011112222333344440000111122223333444400001111222233334444 <== rate r
//	000000000000000000001111111111111111111122222222222222222222 <== site s
//	indexing:             ^ this is s*4*5 + r*4 + b = 22

/**
 * Python bindings
 */


typedef struct {
	unsigned total_n_sites;
	unsigned n_categ_or_subs; /* if n_categ_arr is 0L, then this is the number of categories in the mixture. If n_categ_arr is not NULL, then this is the number of subsets (the length of the n_categ_arr array)*/
	unsigned * n_categ_arr; /*array of the number of mixture components in each of the subsets (0L if there is no partitioning).*/
	unsigned * subset_component_offset; /* array of the sum of the number of mixture components that precede each subset (0L if there is no partitioning).*/
	unsigned * subset_char_offsets; /* The number of characters that precede each subset */
	unsigned * subset_cla_offsets; /* The number of sites*categories that precede each subset */
	unsigned max_n_categ;
} CategSubsetObj;


typedef struct {
	PyObject_HEAD
	int * ssind; /* index of the state set for each character pattern. */
	double *** pmat_columns; /* [subset_ind][to_state_set_ind][(categ_ind*nStates) + from_state_ind] */
	int ** state_lookup; /* maps an state_set_index to the array indicating the # of states and then listing them */
	unsigned nStates;
	unsigned n_state_sets;
	/* # categories does not need to be stored at the leaves
		unsigned n_categ;
	*/
	StateSetLookupStruct * state_set_lookup; /*if NULL, then state_lookup is owned by this struct, if not NULL then it will alias that field. */
	void * pmat_obj_ptr; /* used to check validity of data (pointer to the pmat_obj used in calculations)*/
	int calcTimeStamp; /* used to check validity of data ("timestamp" of the pmat_obj used in calculations)*/
} LeafDataObj;

#if defined(BUILDING_FOR_PYTHON)
	staticforward PyTypeObject leaf_data_type;
#endif

typedef double cla_float_t; /*Determines the precision of the likelihood */

typedef struct {
	PyObject_HEAD
	cla_float_t * cla; /* the conditional likelihood array */
#	if SEPARATE_RESCALER_ARRAY
		rescale_history_t * n_rescalings; /* NULL if no rescaling has been done*/
		rescale_history_t * rescale_swap_space; /* holder for n_rescalings array if no rescaling has been done */
#	endif
	unsigned n_edges_since_rescaling_check;
	/*The following are mainly used for error checking */
	unsigned len_per_categ; /* # of elements in in cla for every categ and site (either nStates or nStates + 1); */
	unsigned nStates;
	CategSubsetObj cso;
} CLAObj;

#if defined(BUILDING_FOR_PYTHON)
	staticforward PyTypeObject cla_type;
#endif

#if defined(PRINTING_LOTS) && PRINTING_LOTS
	void print_cla(CLAObj * par);
#endif

typedef struct {
	PyObject_HEAD
	/* the conditional likelihood array for the tree
		\TODO: this does not need to be stored, but it makes the root likelihood function easier to write
		 	(because we can delegate a call to a function that we have already written).
	*/
	CLAObj * full_cla;
	/* An array of length 2 * n_sites (plus scratch at the end).
		element 2i will contain the ln-Likelihood for pattern i and (2i+1) will be
		the pattern weight for that pattern.  Weights will usually be the count of
		the number of times that the site was seen, but they may be floats.
	*/
	cla_float_t * pat_lnL_wts;
	/* matrix with a row for each category.
	Each row consists of the state frequencies (conditional on the site
	belonging to that categ) and the category frequency.
	So, the row[nStates] element is the frequency of that categ
	*/
	double ** state_categ_freqs;
	/* len (2*nStates * n_categ) *ALIAS* to the end of the pat_lnL_wts array*/
	double * scratch;
	cla_float_t lnLikelihood; /*weighted sum over all patterns*/
	/*The following are mainly used for error checking */
	unsigned nStates;
	CategSubsetObj cso;
} FullLAObj;

#if defined(BUILDING_FOR_PYTHON)
	staticforward PyTypeObject full_la_type;
#endif


typedef struct {
	PyObject_HEAD
	StateSetLookupStruct * sharedStateSetLookupStruct;
	LeafDataObj ** leafData;
	CLAObj ** clas;
	FullLAObj * treeLike;
	ASRVObj * asrv;
	DSCTModelObj ** model;
	PMatArrayObj ** pmats;
	unsigned nLeafData;
	unsigned nPartials;
	unsigned nPMat;
	unsigned nModels;
} LikeStructsBundle;



void cla_dtor(CLAObj* cla_obj);
void full_la_dtor(FullLAObj* cla_obj);
void leaf_data_dtor(LeafDataObj* leaf_data);
void freeLikeStructFields(LikeStructsBundle *);
void zeroLikeStructFields(LikeStructsBundle *);


CLAObj* cla_new(unsigned n_sites, unsigned nStates, unsigned n_categ);
CLAObj* cla_ss_partitioned_new(unsigned n_sites, unsigned nStates);
DSCTModelObj* dsctModelNew(unsigned dim);
FullLAObj* full_la_new(unsigned n_sites, unsigned nStates, unsigned n_categ);
FullLAObj* full_la_ss_partitioned_new(unsigned n_sites, unsigned nStates);
LeafDataObj* leaf_data_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl);
LeafDataObj* leaf_data_ss_partitioned_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl);




int do_add_internal_to_cla(const CLAObj* fd, const  PMatArrayObj *fmat, CLAObj* par,  const CalcContext context);
int do_add_leaf_to_cla(LeafDataObj* fd, const PMatArrayObj *fmat, CLAObj* par, const CalcContext context);
int do_internal_cla(const CLAObj * fd, const PMatArrayObj *fmat, const CLAObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context);
int do_one_tip_cla(const CLAObj * fd, const PMatArrayObj *fmat, LeafDataObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context);


int do_two_tip_cla(LeafDataObj* fd, const PMatArrayObj *fmat, LeafDataObj* sd, const  PMatArrayObj *smat, CLAObj* par, const CalcContext);
double internal_edge_ln_likelihood(const CLAObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj *, const CalcContext context);
double terminal_edge_ln_likelihood(LeafDataObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj *, const CalcContext context);
double ln_likelihood(FullLAObj *full_la);
double partitioned_ln_likelihood(FullLAObj *full_la, double * subsetLnL);


#endif


long create_likelihood_calc_instance(unsigned int numLeaves,
									 unsigned long numPatterns,
									 unsigned int num_states,
									 unsigned int num_rate_categories,
									 unsigned int num_state_code_arrays,
									 unsigned int numPartialStructs,
									 unsigned int num_prob_models,
									 unsigned int numProbMats,
									 unsigned int numEigenStorage,
									 unsigned int num_rescalings_multipliers,
									 int resource_arg,
									 long resource_flag);

void free_likelihood_calc_instance(long handle);

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
