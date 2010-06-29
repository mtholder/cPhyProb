/**
 *
 * Extension adapted from Alex Martelli's example code in the Python Cookbook 17.1
 *
 * ASRV code from the PAML package (ref below):
 *
 *
 * Yang, Z.  2007.  PAML 4: a program package for phylogenetic analysis by maximum
 *	likelihood.  Molecular Biology and Evolution 24: 1586-1591
 *	(http://abacus.gene.ucl.ac.uk/software/paml.html)
 *
 *
 *  Eigensystem code from MrBayes 3.2
 *  by Fredrik Ronquist, John P. Huelsenbeck, and Paul van der Mark
 *  updated from the CVS Repository: mrbayes from the root:
 *   :pserver:anonymous@mrbayes.cvs.sourceforge.net:/cvsroot/mrbayes
 *  on 2007-Oct-05
 *
 *  Copyright 2002-2007
 *  	John P. Huelsenbeck
 *  	Fredrik Ronquist
 *  	Paul van der Mark
 *
 *	Some of that code was written (or translated from FORTRAN) by David Swofford
 *
 * Other code Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
 * mtholder@gmail.com
 * (see bottom of file)
 *
 *
 */
#if ! defined(DSCT_MODEL_H)
#define DSCT_MODEL_H

#ifdef __cplusplus
extern "C" 
{
#endif



#if defined(NO_INLINE) && NO_INLINE
#	define INLINE
#else
#	define INLINE inline
#endif



#if defined(BUILDING_FOR_PYTHON)

#	include "Python.h"

#else

#	define PyObject_HEAD char pyobjectHead;
#	define PyExc_ValueError 1	
#	define PyExc_RuntimeError 2
#	define PyExc_IndexError 3
#	define PyExc_TypeError 4
#	define PyObject_New(a,b) ((a*) malloc(sizeof(a)));
#	define PyObject_Del(a) (free((a)));
#	define Py_DECREF(a)
#	define Py_INCREF(a)

	void PyErr_NoMemory();
	void PyErr_SetString(int, const char *);
	
	
#endif 


#if defined(PRINTING_LOTS) && PRINTING_LOTS

#	include <stdio.h>
#	define PRINTF(v) (printf(v))
#	define PRINTF1(v,a) (printf(v,a))
#	define PRINTF2(v,a,aa) (printf(v,a,aa))
#	define PRINTF3(v,a,aa,aaa) (printf(v,a,aa,aaa))

#else

#	define PRINTF(v)
#	define PRINTF1(v,a)
#	define PRINTF2(v,a,aa)
#	define PRINTF3(v,a,aa,aaa)

#endif



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
	PyObject_HEAD
	int ** state_lookup; /* maps an state_set_index to the array indicating the # of states and then listing them */
	unsigned n_states; /*the number of "fundamental" states.*/
	unsigned n_state_sets; /*the number of combination of states*/
} StateSetLookupStruct;

#if defined(BUILDING_FOR_PYTHON)
	staticforward PyTypeObject state_set_lookup_type;
#endif 

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
	double *** pmat_columns; /* [subset_ind][to_state_set_ind][(categ_ind*n_states) + from_state_ind] */
	int ** state_lookup; /* maps an state_set_index to the array indicating the # of states and then listing them */
	unsigned n_states;
	unsigned n_state_sets;
	/* # categories does not need to be stored at the leaves
		unsigned n_categ; 
	*/
	StateSetLookupStruct * sslookup; /*if NULL, then state_lookup is owned by this struct, if not NULL then it will alias that field. */
	void * pmat_obj_ptr; /* used to check validity of data (pointer to the pmat_obj used in calculations)*/
	int calc_time_stamp; /* used to check validity of data ("timestamp" of the pmat_obj used in calculations)*/
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
	unsigned len_per_categ; /* # of elements in in cla for every categ and site (either n_states or n_states + 1); */
	unsigned n_states;
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
	So, the row[n_states] element is the frequency of that categ
	*/
	double ** state_categ_freqs;
	/* len (2*n_states * n_categ) *ALIAS* to the end of the pat_lnL_wts array*/
	double * scratch; 
	cla_float_t lnLikelihood; /*weighted sum over all patterns*/
	/*The following are mainly used for error checking */
	unsigned n_states;
	CategSubsetObj cso;
} FullLAObj;

#if defined(BUILDING_FOR_PYTHON)
	staticforward PyTypeObject full_la_type;
#endif

/* Define the type that corresponds to a model (holds the Q-Matrix and other
	temporary fields used to speed up calculation).
*/
typedef struct {
	PyObject_HEAD
	unsigned dim;
	double **q_mat;
	int eigen_calc_dirty; /*0 if the eigen calculations are up-to-date*/
	double * eigen_values;
	double * im_eigen_values;
	double ** eigen_vectors;
	double ** inv_eigen_vectors;
	double * cijk; /*len dim^3 */
	double ** work_mat; /*dim by dim*/
	double * d_work;  /*len dim*/
	int * i_work; /*len dim*/
	double *** mat_del_handle;
  	double * arr_del_handle;
} DSCTModelObj;
#if defined(BUILDING_FOR_PYTHON)
	/* Forward declare the python type object*/
	staticforward PyTypeObject dsct_model_type;
#endif
/* Define the type that corresponds to a model (holds the Q-Matrix and other
	temporary fields used to speed up calculation).
*/
typedef struct {
	PyObject_HEAD
	unsigned n_mat;
	unsigned n_states;
	double *** p_mat; /* [categ_ind][from_state][to_state] */
	double * first_el; /*alias to the first element*/
	DSCTModelObj ** model_aliases; /*aliases to the model used for this pmat (may not be stable over repeated invocations -- these are used as scratch) */
	double * brlen_aliases; /*copies of the branch length used for this pmat (may not be stable over repeated invocations -- these are used as scratch) */
	int calc_time_stamp;
} PMatArrayObj;
#if defined(BUILDING_FOR_PYTHON)
	/* Forward declare the python type object*/
	staticforward PyTypeObject pmat_array_type;
#endif
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
#if defined(BUILDING_FOR_PYTHON)
	/* Forward declare the python type object*/
	staticforward PyTypeObject asrv_type;
#endif


double **allocateDblMatrix(unsigned n_rows, unsigned n_cols);
double ***allocateDbl3DMatrix(unsigned n_mat, unsigned n_rows, unsigned n_cols);
void freeDblMatrix(double **);
void freeDbl3DMatrix(double ***);
void cpmat_array_dtor(PMatArrayObj* pmat_array_obj);
void asrv_obj_dtor(ASRVObj* asrh);
StateSetLookupStruct* sslookup_new(unsigned n_states, unsigned n_state_sets, int ** p);
void sslookup_dtor(StateSetLookupStruct* dsct_model);
LeafDataObj* leaf_data_ss_partitioned_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl);
LeafDataObj* leaf_data_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl);
void leaf_data_dtor(LeafDataObj* leaf_data);
CLAObj* cla_new(unsigned n_sites, unsigned n_states, unsigned n_categ);
CLAObj* cla_ss_partitioned_new(unsigned n_sites, unsigned n_states);
void cla_dtor(CLAObj* cla_obj);
FullLAObj* full_la_new(unsigned n_sites, unsigned n_states, unsigned n_categ);
FullLAObj* full_la_ss_partitioned_new(unsigned n_sites, unsigned n_states);

void full_la_dtor(FullLAObj* cla_obj);
DSCTModelObj* dsct_model_new(unsigned dim);
void cdsctm_dtor(DSCTModelObj* dsct_model);
PMatArrayObj* cpmat_array_new(unsigned n_mat, unsigned n_states);
ASRVObj* asrv_obj_new(unsigned dim, int style, double param);
void internal_asrv_set_shape(ASRVObj *asrh, double val);

typedef struct {
	unsigned categ_beg;
	unsigned categ_end;
	unsigned subset_beg;
	unsigned subset_end;
	unsigned rescale_threshold;
} CalcContext;
int configure_context(int ibeg_subset, int iend_subset, unsigned n_subsets, int ibeg_cat, int iend_cat, int irescale_thresh, unsigned n_cat, CalcContext * context);




int do_add_internal_to_cla(const CLAObj* fd, const  PMatArrayObj *fmat, CLAObj* par,  const CalcContext context);
int do_add_leaf_to_cla(LeafDataObj* fd, const PMatArrayObj *fmat, CLAObj* par, const CalcContext context);
int do_internal_cla(const CLAObj * fd, const PMatArrayObj *fmat, const CLAObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context);
int do_one_tip_cla(const CLAObj * fd, const PMatArrayObj *fmat, LeafDataObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context);
int do_asrv_pmat_array_calc(PMatArrayObj *p, DSCTModelObj * mod, ASRVObj *asrv, double edgeLen);
int do_pmat_array_calc(PMatArrayObj *p, unsigned start_ind, unsigned end_ind);
int do_pmat_calc(double **p, DSCTModelObj *mod, double branch_len);
int do_two_tip_cla(LeafDataObj* fd, const PMatArrayObj *fmat, LeafDataObj* sd, const  PMatArrayObj *smat, CLAObj* par, const CalcContext);
double internal_edge_ln_likelihood(const CLAObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj *, const CalcContext context);
double terminal_edge_ln_likelihood(LeafDataObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj *, const CalcContext context);
double ln_likelihood(FullLAObj *full_la);
double partitioned_ln_likelihood(FullLAObj *full_la, double * subsetLnL);
int copy_cla(CLAObj * dest, const CLAObj * source);
int next_calc_stamp();

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
