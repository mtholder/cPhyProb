#if 0






typedef struct {
	unsigned categ_beg;
	unsigned categ_end;
	unsigned subset_beg;
	unsigned subset_end;
	unsigned rescale_threshold;
} CalcContext;
int configure_context(int ibeg_subset, int iend_subset, unsigned n_subsets, int ibeg_cat, int iend_cat, int irescale_thresh, unsigned n_cat, CalcContext * context);




int copy_cla(CLAObj * dest, const CLAObj * source);
int next_calc_stamp(void);

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


void freeLikeStructFields(LikeStructsBundle *); 
void zeroLikeStructFields(LikeStructsBundle *); 

#endif



#if 0
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#if defined(STANDALONE_C_PROG) || (defined(PRINTING_LOTS) && PRINTING_LOTS)
#	include <stdio.h>
#endif



double **allocateDblMatrix(unsigned n_rows, unsigned n_cols);
double ***allocateDbl3DMatrix(unsigned nMatrices, unsigned n_rows, unsigned n_cols);
void freeDblMatrix(double **);
void freeDbl3DMatrix(double ***);










#if defined(ALL_IN_ONE) || !defined(BUILDING_FOR_PYTHON)







int new_cso_ss_partitioned(CategSubsetObj * cso, unsigned n_sites, unsigned nStates);

void cso_dtor(CategSubsetObj * cso);

#if !defined(BUILDING_FOR_PYTHON)
	void PyErr_NoMemory() {}
	void PyErr_SetString(int v, const char *c)
	{
		CPHYPROB_DEBUG_PRINTF2("Error: code %d\nmessage: %s\n", v, c);
	}
#endif



/**/
#define K_LN_IGNORABLE_ADDITION 41

/* check that a maximum of  1 form of rescaling to avoid underflow has been
	specified
This section set SEPARATE_RESCALER_ARRAY to the appropriat value:
	0 if the ln(rescale factor) is intercalated after the last state, or
	1 if an array that stores ln(rescale factor) for each site/categ is stored
		in the cla
*/
#if defined(CONSTANT_STEP_RESCALING) && CONSTANT_STEP_RESCALING
	/* In this form of rescaling we calculate a constant number to bump all
	likelihoods up by (won't take the highest likeilhood up to 1.0), but allows
	us to store the rescalings as integer count.
	*/
#	if defined(INT_INCR_RESCALING) && INT_INCR_RESCALING
#		error "INT_INCR_RESCALING and CONSTANT_STEP_RESCALING cannot both be specified"
#	elif defined(TO_EXACTLY_ONE_RESCALING) && TO_EXACTLY_ONE_RESCALING
#		error "TO_EXACTLY_ONE_RESCALING and CONSTANT_STEP_RESCALING cannot both be specified"
#	endif
#	define SEPARATE_RESCALER_ARRAY 1
#	define INT_INCR_RESCALING 0
#	define TO_EXACTLY_ONE_RESCALING 0
	const double K_RESCALING_FACTOR = 1.0e-50;
	INLINE double getLnRescalingFactor() {
		static double K_LN_RESCALING_FACTOR = 1.0;
		if (K_LN_RESCALING_FACTOR > 0.0) {
			K_LN_RESCALING_FACTOR = log(K_RESCALING_FACTOR);
		}
		return K_LN_RESCALING_FACTOR;
	}
	typedef unsigned char rescale_history_t;
#elif defined(INT_INCR_RESCALING) && INT_INCR_RESCALING
	/* In this form of rescaling, we make the ln of the rescaling factor be an
	integer. Thus, we get close to rescaling the biggest category to 1.0, but we
	are not exact (on the plus side we can store the likelihood rescale count
	as an integer).
	*/
#	if defined(TO_EXACTLY_ONE_RESCALING) && TO_EXACTLY_ONE_RESCALING
#		error "TO_EXACTLY_ONE_RESCALING and INT_INCR_RESCALING cannot both be specified"
#	endif
#	define CONSTANT_STEP_RESCALING 0
#	define TO_EXACTLY_ONE_RESCALING 0
#	define SEPARATE_RESCALER_ARRAY 1
	typedef unsigned int rescale_history_t;
#else
#	define CONSTANT_STEP_RESCALING 0
#	define INT_INCR_RESCALING 0
#	define SEPARATE_RESCALER_ARRAY 0
#	if defined(TO_EXACTLY_ONE_RESCALING)
		/* This is the default.
		When rescaling is triggered, the largest conditional likelihood is set to
		1.0
		The ln rescale factors must be stored as floating point numbers, and they
		appear after the conditional likelihoods of the state.
		*/
#		if !TO_EXACTLY_ONE_RESCALING
#			error "TO_EXACTLY_ONE_RESCALING, INT_INCR_RESCALING, or CONSTANT_STEP_RESCALING must be defined."
#		endif
#	else
#		define TO_EXACTLY_ONE_RESCALING 1
#	endif
#endif

const double K_RESCALE_THRESHOLD = 1.0e-50;
const double K_RESCALE_ERROR_THRESHOLD = 1.0e-150;
const int UNDERFLOW_ERROR_CODE = -1;




/*Internal function prototypes*/
	/*phylogenetic functions*/
	/*From MrBayes*/





INLINE void transpose_pmat_columns(double * transposed, const double ***p_mat_array, const unsigned to_state, const unsigned begin_categ, const unsigned end_categ, const unsigned nStates);
static int assign_pmat_to_leaf(LeafDataObj * leaf_d, const PMatArrayObj * pmat, const CalcContext, const CategSubsetObj);
static int do_two_tip_cla_recalc_rescale(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj* par, const CalcContext context);
static int do_two_tip_cla_recalc_no_rescale(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj* par, const CalcContext context);
INLINE static int do_two_tip_cla_recalc(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj* par, const CalcContext);
static int do_one_tip_cla_recalc_rescale(const CLAObj * fd, const double *** fmat, const LeafDataObj* sd, CLAObj* par, const CalcContext context);
static int do_one_tip_cla_recalc_no_rescale(const CLAObj * fd, const double *** fmat, const LeafDataObj* sd, CLAObj* par, const CalcContext context);
INLINE static int do_one_tip_cla_recalc(const CLAObj * fd, const PMatArrayObj *fmat, const LeafDataObj* sd, CLAObj* par, const CalcContext context);
static int do_internal_cla_rescale(const CLAObj * fd, const double *** fmat, const CLAObj* sd, const double *** smat, CLAObj* par, const CalcContext context);
static int do_internal_cla_no_rescale(const CLAObj * fd, const double *** fmat, const CLAObj* sd, const double *** smat, CLAObj* par, const CalcContext context);
static int do_add_leaf_to_cla_recalc_rescale(const LeafDataObj* fd, CLAObj* par, const CalcContext context);
static int do_add_leaf_to_cla_recalc_no_rescale(const LeafDataObj* fd, CLAObj* par, const CalcContext context);
static int do_add_leaf_to_cla_recalc(const LeafDataObj* fd, CLAObj* par, const CalcContext context);

static int do_add_internal_to_cla_rescale(const CLAObj* fd, const double *** fmat, CLAObj * par, const CalcContext context);
static int do_add_internal_to_cla_no_rescale(const CLAObj* fd, const double *** fmat, CLAObj * par, const CalcContext context);
double ln_likelihood(FullLAObj *);




void zeroLikeStructFields(LikeStructsBundle * toFree)
{
	toFree->sharedStateSetLookupStruct = 0L;
	toFree->leafData = 0L;
	toFree->clas = 0L;
	toFree->treeLike = 0L;
	toFree->asrv = 0L;
	toFree->model = 0L;
	toFree->pmats = 0L;
}

// frees all of the memory that the LikeStructsBundle object points to, but does
// NOT free the LikeStructsBundle pointer itself.
void freeLikeStructFields(LikeStructsBundle * toFree)
{
	int i;
	if (toFree == 0L)
		return;
	if (toFree->sharedStateSetLookupStruct) {
		state_set_lookup_dtor(toFree->sharedStateSetLookupStruct);
		toFree->sharedStateSetLookupStruct = 0L;
	}
	if (toFree->asrv) {
		asrv_obj_dtor(toFree->asrv);
		toFree->asrv = 0L;
	}

	if (toFree->leafData) {
		for (i = 0; i < toFree->nLeafData; ++i) {
			if (toFree->leafData[i]) {
				leaf_data_dtor(toFree->leafData[i]);
			}
		}
		free(toFree->leafData);
		toFree->leafData = 0L;
	}

	if (toFree->clas) {
		for (i = 0; i < toFree->nPartials; ++i) {
			if (toFree->clas[i]) {
				cla_dtor(toFree->clas[i]);
			}
		}
		free(toFree->clas);
		toFree->clas = 0L;
	}

	if (toFree->treeLike) {
		full_la_dtor(toFree->treeLike);
		toFree->treeLike = 0L;
	}
	if (toFree->model) {
		for (i = 0; i < toFree->nModels; ++i) {
			if (toFree->model[i]) {
				cdsctm_dtor(toFree->model[i]);
			}
		}
		free(toFree->model);
		toFree->model = 0L;
	}
	if (toFree->pmats) {
		for (i = 0; i < toFree->nPMat; ++i) {
			if (toFree->pmats[i]) {
				cpmatArrayDtor(toFree->pmats[i]);
			}
		}
		free(toFree->pmats);
		toFree->pmats = 0L;
	}
}



void printQMat(double ** qmat_obj, const unsigned nStates) {
	unsigned from_state, to_state;
	CPHYPROB_DEBUG_PRINTF("QMat:\n");
	for (from_state = 0; from_state < nStates; ++from_state) {
		for (to_state = 0; to_state < nStates; ++to_state) {
			CPHYPROB_DEBUG_PRINTF1("%f ", qmat_obj[from_state][to_state]);
		}
		CPHYPROB_DEBUG_PRINTF("\n");
	}
}


#if 0
/**
 * Adapted from MrBayes TiProbsUsingEigens in mbmath.c
 *	`pMat` and `EigValExp` are written to.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
static int prob_mat_and_derivs_from_eigensystem (
	const unsigned dim, /*the number of states*/
	double **pMat,			/**/
	double **pMatDeriv,		/**/
	double **pMatSecDeriv,	/**/
	double *EigValExp,		/**/
	const double *cijk, /*the dim*dim*dim array of temporaries used in quickly multiplying the exp(eval*branchlength) times the matrix eigenvectors and its inverse*/
	const double *eigenVals, /*array of the eigenvalues*/
	const double t, 	/*branch length (time)*/
	const double r ) /* rate */
{
	unsigned i, j, s;
	double sum, sumF, sumS;
	const double branch_len = t * r;
	const double rsq = r * r;
	assert(pMat);
	assert(*pMat);
	assert(pMatDeriv);
	assert(*pMatDeriv);
	assert(pMatSecDeriv);
	assert(*pMatSecDeriv);
	assert(cijk);
	assert(eigenVals);
	assert(t >= 0.0);
	assert(r >= 0.0);
	for (i=0; i<dim; i++)
		EigValExp[i] = exp(eigenVals[i] * branch_len);
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			sum = sumF = sumS = 0.0;
			for(s=0; s<dim; s++) {
				const double dcoeff = eigenVals[s];
				const double derivF =  dcoeff * EigValExp[s];
				const double cijk_el = *cijk++;
				sum += cijk_el * EigValExp[s];
				sumF += cijk_el * derivF;
				sumS += cijk_el * derivF * dcoeff;
			}
			sumF *= r;
			sumS *= rsq;
			pMat[i][j] = (sum < 0.0) ? 0.0 : sum;
			pMatDeriv[i][j] = sumF;
			pMatSecDeriv[i][j] = sumS;
		}
	}
	return 1;
}

#endif 

/**
 * Calculates the p->p_matrices assuming that the all use the same model with rate
 *  heteregeneity described by the asrv param.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int do_asrv_pmat_array_calc(PMatArrayObj *p, DSCTModelObj * mod, ASRVObj *asrv, double edgeLen) {
	unsigned i;
	assert(p);
	assert(mod);
	if (asrv) {
		for (i = 0; i < asrv->n; ++i) {
			p->brlenAliases[i] = asrv->val[i]*edgeLen;
			p->modelAliases[i] = mod;
		}
		return do_pmat_array_calc(p, 0, asrv->n);
	}
	p->brlenAliases[0] = edgeLen;
	p->modelAliases[0] = mod;
	return do_pmat_array_calc(p, 0, asrv->n);
}

INLINE void transpose_pmat_columns(double * transposed, const double ***p_mat_array, const unsigned to_state, const unsigned begin_categ, const unsigned end_categ, const unsigned nStates) {
	unsigned i, from_state;
	const double **pMatrix;
	for (i = begin_categ; i < end_categ; ++i) {
		pMatrix = p_mat_array[i];
		for (from_state = 0; from_state < nStates; ++from_state) {
			transposed[from_state + (i*nStates)] = pMatrix[from_state][to_state];
		}
	}
}

int copy_cla(CLAObj * dest, const CLAObj * source) {
	assert(dest);
	assert(source);
	assert(dest->cla);
	assert(source->cla);
#	if SEPARATE_RESCALER_ARRAY
		assert(dest->rescale_swap_space);
		assert(source->rescale_swap_space);
#	endif
	const unsigned subset_ind = (source->cso.n_categ_arr == 0L ? 0 : source->cso.n_categ_or_subs - 1);
	unsigned n_sites_in_last = source->cso.total_n_sites;
	if (subset_ind > 1) {
		n_sites_in_last -= source->cso.subset_char_offsets[subset_ind];
	}
	const unsigned comp_offset = (subset_ind == 0 ? 0 : source->cso.subset_cla_offsets[subset_ind]);
	const unsigned n_categ = (source->cso.n_categ_arr == 0L ? source->cso.n_categ_or_subs : source->cso.n_categ_arr[subset_ind]);
	const unsigned total_n_categ = comp_offset + n_sites_in_last*n_categ;
	memcpy(dest->cla, source->cla, total_n_categ*(source->len_per_categ)*sizeof(cla_float_t));
#	if SEPARATE_RESCALER_ARRAY
		if (source->n_rescalings) {
			memcpy(dest->rescale_swap_space, source->rescale_swap_space, (source->n_sites)*(source->n_categ)*sizeof(rescale_history_t));
			dest->n_rescalings = dest->rescale_swap_space;
		}
		else {
			dest->n_rescalings = 0L;
		}
#	endif
	dest->n_edges_since_rescaling_check = source->n_edges_since_rescaling_check;
	dest->len_per_categ = source->len_per_categ;
	dest->nStates = source->nStates;
	assert(dest->cso.n_categ_or_subs == source->cso.n_categ_or_subs);
	dest->cso.n_categ_or_subs = source->cso.n_categ_or_subs;
	dest->cso.max_n_categ = source->cso.max_n_categ;
	dest->cso.total_n_sites = source->cso.total_n_sites;
	if (source->cso.n_categ_arr) {
		assert(dest->cso.n_categ_arr);
		unsigned i; 
		for (i = 0; i < dest->cso.n_categ_or_subs; ++i) {
			dest->cso.n_categ_arr[i] = source->cso.n_categ_arr[i];
			dest->cso.subset_component_offset[i] = source->cso.subset_component_offset[i];
			dest->cso.subset_char_offsets[i] = source->cso.subset_char_offsets[i];
			dest->cso.subset_cla_offsets[i] = source->cso.subset_cla_offsets[i];
		}
	}
	return 1;
}

void printPMat(const PMatArrayObj * pmat_obj, const CategSubsetObj cso) {
	const double *** pmat_array =  (const double ***) pmat_obj->pMatrix;
	const unsigned nStates = pmat_obj->nStates;
	const unsigned subset_ind = (cso.n_categ_arr == 0L ? 0 : cso.n_categ_or_subs - 1);
	const unsigned comp_offset = (subset_ind == 0 ? 0 : cso.subset_component_offset[subset_ind]);
	const unsigned n_categ = (cso.n_categ_arr == 0L ? cso.n_categ_or_subs : cso.n_categ_arr[subset_ind]);
	const unsigned end_categ = comp_offset + n_categ;
	unsigned categ, from_state, to_state;
	CPHYPROB_DEBUG_PRINTF("PMat:\n");
	for (categ = 0; categ < end_categ; ++categ) {
		for (from_state = 0; from_state < nStates; ++from_state) {
			for (to_state = 0; to_state < nStates; ++to_state) {
				CPHYPROB_DEBUG_PRINTF1("%f ", pmat_array[categ][from_state][to_state]);
			}
			CPHYPROB_DEBUG_PRINTF("\n");
		}
		CPHYPROB_DEBUG_PRINTF("\n");		
	}
}


static int assign_pmat_to_leaf(LeafDataObj * leaf_d, const PMatArrayObj * pmat_obj, const CalcContext context, const CategSubsetObj cso) {
	const int * slook_arr;
	double ** transposed;
	double * copy_source;
	unsigned i, j, c, n_ambig_states, offset, to_state, from_state, subset_ind;
	assert(leaf_d);
	assert(pmat_obj);
 	const double *** pmat_array =  (const double ***) pmat_obj->pMatrix;
	const unsigned nStates = pmat_obj->nStates;
	const int ** slookup = (const int **) leaf_d->state_lookup;
	const unsigned n_state_sets = leaf_d->n_state_sets;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
#	if defined(PRINTING_LOTS) && PRINTING_LOTS
		printPMat(pmat_obj, cso);
#	endif
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		transposed = leaf_d->pmat_columns[subset_ind];
		for (i = 0; i < n_state_sets; ++i) {
			slook_arr = slookup[i];
			n_ambig_states = slook_arr[0];
			to_state = slook_arr[1];
			const unsigned comp_offset = (subset_ind == 0 ? 0 : cso.subset_component_offset[subset_ind]);
			const unsigned n_categ = (cso.n_categ_arr == 0L ? cso.n_categ_or_subs : cso.n_categ_arr[subset_ind]);
			const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
			const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);
			if (n_ambig_states == 1) {
				transpose_pmat_columns(transposed[to_state], pmat_array + comp_offset, to_state, categ_beg, categ_end, nStates);
			}
			else {
				/* ambiguity -- sum the ambiguous states (each member of the state
				set must have already been calculated */
				assert(to_state < i);
				double * dest = transposed[i];
				copy_source = transposed[to_state];
				for (c = categ_beg; c < categ_end; ++c) {
					offset = c*nStates;
					for (from_state = 0; from_state < nStates; ++from_state) {
						dest[from_state + offset] = copy_source[from_state + offset];
					}
				}
				for (j = 2; j <= n_ambig_states; ++j) {
					to_state = slook_arr[j];
					assert(to_state < i);
					copy_source = transposed[to_state];
					for (c = categ_beg; c < categ_end; ++c) {
						offset = c*nStates;
						for (from_state = 0; from_state < nStates; ++from_state) {
							dest[from_state + offset] += copy_source[from_state + offset];
						}
					}
				}
			}
		}
	}
	leaf_d->pmat_obj_ptr = (void *) pmat_obj;
	leaf_d->calcTimeStamp =  pmat_obj->calcTimeStamp;
	return 1;
}
#if defined(PRINTING_LOTS) && PRINTING_LOTS
	void print_cla(CLAObj * par) {
		const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
		cla_float_t *p_cla = par->cla;
		const unsigned nStates = par->nStates;
		unsigned i, j, k, subset_ind;
		for (subset_ind = 0; subset_ind < n_subsets; ++subset_ind) {
			const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
			unsigned beg_site, end_site;
			if (n_subsets > 1) {
				beg_site = par->cso.subset_char_offsets[subset_ind];
				end_site = (subset_ind == n_subsets - 1 ? par->cso.total_n_sites : par->cso.subset_char_offsets[subset_ind + 1]);
			}
			else {
				beg_site = 0;
				end_site = par->cso.total_n_sites;
			}
			for (i = beg_site; i < end_site; ++i) {
				printf("%5d ", (i+1));
				for (j = 0; j < n_categ; ++j) {
					for (k = 0; k < nStates; ++k) {
						printf("%10g ", *p_cla++);
					}
#				if TO_EXACTLY_ONE_RESCALING
					p_cla++;
#				endif
				printf("\n");
				}
			}
		}
	}
#endif
static int do_two_tip_cla_recalc_rescale(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj * par, const CalcContext context) {
	const double ** f_mat;
	const double ** s_mat;
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;	
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned subset_ind;

	double p, max_p, real_rescale;
#	if SEPARATE_RESCALER_ARRAY
#		if INT_INCR_RESCALING
			double real_rescale_ln;
#		endif
		int rescaling_done = 0;
		rescale_history_t rcount;
#	else
		double real_rescale_ln;
#	endif
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		f_mat = (const double **)fd->pmat_columns[subset_ind];
		s_mat = (const double **)sd->pmat_columns[subset_ind];
		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;

		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * f_states = (const int *) fd->ssind;
		f_states += subset_char_offsets ;
		const int * s_states = (const int *) sd->ssind;
		s_states += subset_char_offsets ;
	
		const unsigned tpm_offset = nStates*(categ_beg);
		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);

		unsigned beg_site, end_site;
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}
	
		const double * f_row;
		const double * s_row;
		unsigned site, categ, state;
#		if SEPARATE_RESCALER_ARRAY
			const unsigned subset_resc_offset = subset_cla_offsets;
			rescale_history_t * rescale_count_ptr = par->rescale_swap_space + subset_resc_offset;
#		endif
		for (site = beg_site; site < end_site; ++site) {
			f_row = f_mat[*f_states++] + tpm_offset;
			s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				max_p = 0.0;
				for (state = 0; state < nStates; ++state) {
					p = (*f_row++) * (*s_row++);
					if (p > max_p)
						max_p = p;
					p_cla[state] = p;
				}
				if (max_p < K_RESCALE_THRESHOLD) {
					if (max_p < K_RESCALE_ERROR_THRESHOLD)
						return UNDERFLOW_ERROR_CODE;
#					if SEPARATE_RESCALER_ARRAY
#						if INT_INCR_RESCALING
							rcount = 1 + (rescale_history_t)log(max_p);
							real_rescale_ln = (double) rcount;
							real_rescale = exp(real_rescale_ln);
#						elif CONSTANT_STEP_RESCALING
							real_rescale = K_RESCALING_FACTOR;
							rcount = 1;
#						else
#							error "unknown RESCALING type"
#						endif
						rescaling_done = 1;
						*rescale_count_ptr++ = rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						real_rescale_ln = log(max_p);
						real_rescale = max_p;
						p_cla[nStates] = real_rescale_ln;
#					endif
					for (state = 0; state < nStates; ++state)
						p_cla[state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = 0;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[nStates] = 0.0;
#					endif
				}
				p_cla += len_per_categ;
			}
			p_cla += post_skip;
#			if SEPARATE_RESCALER_ARRAY
				rescale_count_ptr += n_skipped_categ;
#			endif
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings  = (rescaling_done ? par->rescale_swap_space : 0L);
#	endif
	par->n_edges_since_rescaling_check = 0;
	return 1;
}
static int do_two_tip_cla_recalc_no_rescale(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj* par, const CalcContext context) {
	const double ** f_mat;
	const double ** s_mat;
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	cla_float_t tmp;
	
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned subset_ind;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		f_mat = (const double **) fd->pmat_columns[subset_ind];
		s_mat = (const double **) sd->pmat_columns[subset_ind];
		
		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;

		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * f_states = (const int *) fd->ssind;
		f_states += subset_char_offsets ;
		const int * s_states = (const int *) sd->ssind;
		s_states += subset_char_offsets ;
	

		const unsigned tpm_offset = nStates*(categ_beg);
		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);

		unsigned beg_site, end_site;
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}
	
		const double * f_row;
		const double * s_row;
		unsigned site, categ, state;

		CPHYPROB_DEBUG_PRINTF("do_two_tip_cla_recalc_no_rescale: categ cla likes:\n");
		for (site = beg_site; site < end_site; ++site) {
			CPHYPROB_DEBUG_PRINTF1("%d", site);
			f_row = f_mat[*f_states++] + tpm_offset;
			s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				for (state = 0; state < nStates; ++state) {
					tmp = (*f_row++) * (*s_row++);
					p_cla[state]  = tmp;
					CPHYPROB_DEBUG_PRINTF1("\t%f", tmp);
				}
#				if !SEPARATE_RESCALER_ARRAY
					p_cla[nStates] = 0.0;
#				endif
				p_cla += len_per_categ;
			}
			p_cla += post_skip;
			CPHYPROB_DEBUG_PRINTF("\n");
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings = 0L;
#	endif
	par->n_edges_since_rescaling_check = 2;
	return 1;
}
INLINE static int do_two_tip_cla_recalc(const LeafDataObj* fd, const LeafDataObj* sd, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(sd);
	assert(par);
	assert(fd->nStates == sd->nStates);
	assert(fd->nStates == par->nStates);
	/* 
		assert(fd->n_categ == sd->n_categ);
		assert(fd->n_categ == par->n_categ);
	*/
	if (2 < context.rescale_threshold)
		return do_two_tip_cla_recalc_no_rescale(fd, sd, par, context);
	return do_two_tip_cla_recalc_rescale(fd, sd, par, context);
}

int do_two_tip_cla(LeafDataObj* fd, const PMatArrayObj *fmat, LeafDataObj* sd, const  PMatArrayObj *smat, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(sd);
	assert(smat);
	if ((fd->pmat_obj_ptr != (void *) fmat) || (fd->calcTimeStamp != fmat->calcTimeStamp)) {
		if (!assign_pmat_to_leaf(fd, fmat, context, par->cso))
			return 0;
	}
	if ((sd->pmat_obj_ptr != (void *) smat) || (sd->calcTimeStamp != smat->calcTimeStamp)) {
		if (!assign_pmat_to_leaf(sd, smat, context, par->cso))
			return 0;
	}
	return do_two_tip_cla_recalc(fd, sd, par, context);
}

static int do_one_tip_cla_recalc_rescale(const CLAObj * fd, const double *** f_pmat_arr, const LeafDataObj* sd, CLAObj* par, const CalcContext context) {
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const double ** f_categ_pmat;
	const double * f_categ_p;
	const double ** s_mat;
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned site, categ, to_state, from_state, subset_ind, beg_site, end_site;
	cla_float_t f_des_prob;
	double p, max_p, real_rescale;
#	if SEPARATE_RESCALER_ARRAY
#		if INT_INCR_RESCALING
			double real_rescale_ln;
#		endif
		rescale_history_t rcount;
		rescale_history_t prev_rcount;
		const int has_prev_rescale = (fd->n_rescalings == 0L ? 0 : 1);
		int rescaling_done = (has_prev_rescale > 0 ? 1 : 0);
#	else
		double real_rescale_ln;
#	endif
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		s_mat = (const double **)sd->pmat_columns[subset_ind];
		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t *f_cla = (const cla_float_t*) (fd->cla + start_offset);
	
		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * s_states = (const int *) sd->ssind;
		s_states += subset_char_offsets ;


		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const unsigned tpm_offset = nStates*(categ_beg);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;

		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

#		if SEPARATE_RESCALER_ARRAY
			const unsigned subset_resc_offset = subset_cla_offsets;

			const rescale_history_t * prev_rescale_count;
			rescale_history_t * rescale_count_ptr = par->rescale_swap_space + subset_resc_offset;
			if (has_f_rescale)
				prev_rescale_count = ((const rescale_history_t *) sd->n_rescalings) + subset_resc_offset;
#		endif

		for (site = beg_site; site < end_site; ++site) {
			const double * s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				max_p = 0.0;
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					p = f_des_prob*(*s_row++);
					if (p > max_p)
						max_p = p;
					p_cla[from_state] = p;
				}
#				if SEPARATE_RESCALER_ARRAY
					prev_rcount = (has_prev_rescale ? *prev_rescale_count++ : 0);
#				endif
				if (max_p < K_RESCALE_THRESHOLD) {
					if (max_p < K_RESCALE_ERROR_THRESHOLD)
						return UNDERFLOW_ERROR_CODE;
#					if SEPARATE_RESCALER_ARRAY
#						if INT_INCR_RESCALING
							rcount = 1 + (rescale_history_t)log(max_p);
							real_rescale_ln = (double) rcount;
							real_rescale = exp(real_rescale_ln);
#						elif CONSTANT_STEP_RESCALING
							real_rescale = K_RESCALING_FACTOR;
							rcount = 1;
#						else
#							error "unknown RESCALING type"
#						endif
						rescaling_done = 1;
						*rescale_count_ptr++ = rcount + prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						real_rescale_ln = log(max_p);
						real_rescale = max_p;
						p_cla[nStates] = real_rescale_ln + f_cla[nStates];
#					endif
					for (from_state = 0; from_state < nStates; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[nStates] = f_cla[nStates];
#					endif
				}
				p_cla += len_per_categ;
				f_cla += len_per_categ;
			}
			p_cla += post_skip;
			f_cla += post_skip;
#			if SEPARATE_RESCALER_ARRAY
				rescale_count_ptr += n_skipped_categ;
				if (has_prev_rescale)
					prev_rescale_count += n_skipped_categ;
#			endif
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings  = (rescaling_done ? par->rescale_swap_space : 0L);
#	endif
	par->n_edges_since_rescaling_check = 0;
	return 1;
}
static int do_one_tip_cla_recalc_no_rescale(const CLAObj * fd, const double *** f_pmat_arr, const LeafDataObj* sd, CLAObj* par, const CalcContext context) {
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const double ** f_categ_pmat;
	const double * f_categ_p;
	const double ** s_mat;
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned site, categ, to_state, from_state, subset_ind, beg_site, end_site;
	cla_float_t f_des_prob;
	
	CPHYPROB_DEBUG_PRINTF("do_one_tip_cla_recalc_no_rescale: categ cla likes:\n");
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings = (fd->n_rescalings ? par->rescale_swap_space : 0L);
#	endif
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		s_mat = (const double **)sd->pmat_columns[subset_ind];

		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t *f_cla = (const cla_float_t*) (fd->cla + start_offset);
	
		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * s_states = (const int *) sd->ssind;
		s_states += subset_char_offsets ;


		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const unsigned tpm_offset = nStates*(categ_beg);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;

		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

		for (site = beg_site; site < end_site; ++site) {
			CPHYPROB_DEBUG_PRINTF1("%d", site);
			const double *  s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					double tmp = f_des_prob*(*s_row++);
					p_cla[from_state] = tmp;
					CPHYPROB_DEBUG_PRINTF1("\t%f", tmp);
				}
#				if ! SEPARATE_RESCALER_ARRAY
						p_cla[nStates] = f_cla[nStates];
#				endif
				p_cla += len_per_categ;
				f_cla += len_per_categ;
			}
			CPHYPROB_DEBUG_PRINTF("\n");
			p_cla += post_skip;
			f_cla += post_skip;
#		if SEPARATE_RESCALER_ARRAY
			if (has_prev_rescale) {
				const unsigned resc_offset = (par->cso.subset_cla_offsets[subset_ind]) + beg_categ;
				const rescale_history_t * prev_rescale_count = (const rescale_history_t *) fd->n_rescalings;
				prev_rescale_count += resc_offset;
				rescale_history_t * par_rescale_count = par->n_rescalings;
				par_rescale_count += resc_offset;
				unsigned i, j;
				for (i = beg_site; i < end_site; ++i) {
					for (j = categ_beg; j < categ_end; ++j) {
						par_rescale_count[j] = prev_rescale_count[j];
					}
					par_rescale_count += n_categ;
					prev_rescale_count += n_categ;
				}
			}
#		endif
		}
	}
	par->n_edges_since_rescaling_check = 2 + fd->n_edges_since_rescaling_check;
	return 1;
}
INLINE static int do_one_tip_cla_recalc(const CLAObj * fd, const PMatArrayObj *fmat, const LeafDataObj* sd, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(sd);
	assert(par);
	assert(fd->cso.total_n_sites == par->cso.total_n_sites);
	assert(fd->nStates == sd->nStates);
	assert(fd->nStates == par->nStates);
	/*
		assert(fd->n_categ == sd->n_categ);
	*/
	assert(fd->cso.n_categ_or_subs == par->cso.n_categ_or_subs);
	const double *** f_pmat_arr = (const double ***)fmat->pMatrix;
	if ((fd->n_edges_since_rescaling_check + 2) < context.rescale_threshold)
		return do_one_tip_cla_recalc_no_rescale(fd, f_pmat_arr, sd, par, context);
	return do_one_tip_cla_recalc_rescale(fd, f_pmat_arr, sd, par, context);
}

int do_one_tip_cla(const CLAObj * fd, const PMatArrayObj *fmat, LeafDataObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(sd);
	assert(smat);
	if ((sd->pmat_obj_ptr != (void *) smat) || (sd->calcTimeStamp != smat->calcTimeStamp)) {
		if (!assign_pmat_to_leaf(sd, smat, context, par->cso))
			return 0;
	}
	return do_one_tip_cla_recalc(fd, fmat, sd, par, context);
}
static int do_internal_cla_rescale(const CLAObj * fd, const double *** f_pmat_arr, const CLAObj* sd, const double *** s_pmat_arr, CLAObj* par, const CalcContext context) {
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const double ** f_categ_pmat;
	const double * f_categ_p;
	const double ** s_categ_pmat;
	const double * s_categ_p;
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	double real_rescale, p, max_p;
#	if SEPARATE_RESCALER_ARRAY
#		if INT_INCR_RESCALING
			double real_rescale_ln;
#		endif
		const int has_f_rescale = (fd->n_rescalings == 0L ? 0 : 1);
		const int has_s_rescale = (sd->n_rescalings == 0L ? 0 : 1);
		int rescaling_done = ((has_f_rescale > 0 || has_s_rescale > 0)? 1 : 0);

		rescale_history_t rcount;
		rescale_history_t prev_rcount;
#	else
		double real_rescale_ln;
#	endif
	unsigned subset_ind;
	cla_float_t f_des_prob, s_des_prob;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {

		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t *f_cla = (const cla_float_t*) (fd->cla + start_offset);
		const cla_float_t *s_cla = (const cla_float_t*) (sd->cla + start_offset);
	

		const unsigned trailing_categ = n_categ - categ_end;
#		if SEPARATE_RESCALER_ARRAY
			const unsigned n_skipped_categ = trailing_categ + categ_beg;		
#		endif
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;
		const double *** s_pmat_arr_o = s_pmat_arr + comp_offset;
		unsigned beg_site, end_site, site, categ, from_state, to_state;
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

#		if SEPARATE_RESCALER_ARRAY
			const unsigned subset_resc_offset = subset_cla_offsets;

			rescale_history_t * rescale_count_ptr = par->rescale_swap_space + subset_resc_offset;
			const rescale_history_t * f_prev_rescale_count;
			const rescale_history_t * s_prev_rescale_count;
			if (has_f_rescale)
				f_prev_rescale_count = ((const rescale_history_t *) fd->n_rescalings) + subset_resc_offset;
			if (has_s_rescale)
				s_prev_rescale_count = ((const rescale_history_t *) sd->n_rescalings) + subset_resc_offset;
#		endif

		for (site = beg_site; site < end_site; ++site) {
			for (categ = categ_beg; categ < categ_end; ++categ) {
				max_p = 0.0;
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				s_categ_pmat = s_pmat_arr_o[categ];
				s_categ_p = *s_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					s_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
						s_des_prob += s_cla[to_state]*(*s_categ_p++);
					}
					p = f_des_prob*s_des_prob;
					if (p > max_p)
						max_p = p;
					p_cla[from_state] = p;
				}
#				if SEPARATE_RESCALER_ARRAY
					prev_rcount = (has_f_rescale ? *f_prev_rescale_count++ : 0);
					prev_rcount += (has_s_rescale ? *s_prev_rescale_count++ : 0);
#				endif
				if (max_p < K_RESCALE_THRESHOLD) {
					if (max_p < K_RESCALE_ERROR_THRESHOLD)
						return UNDERFLOW_ERROR_CODE;
#					if SEPARATE_RESCALER_ARRAY
#						if INT_INCR_RESCALING
							rcount = 1 + (rescale_history_t)log(max_p);
							real_rescale_ln = (double) rcount;
							real_rescale = exp(real_rescale_ln);
#						elif CONSTANT_STEP_RESCALING
							real_rescale = K_RESCALING_FACTOR;
							rcount = 1;
#						else
#							error "unknown RESCALING type"
#						endif
						rescaling_done = 1;
						*rescale_count_ptr++ = rcount + prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						real_rescale_ln = log(max_p);
						real_rescale = max_p;
						p_cla[nStates] = real_rescale_ln + f_cla[nStates] + s_cla[nStates];
#					endif
					for (from_state = 0; from_state < nStates; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[nStates] = f_cla[nStates] + s_cla[nStates];
#					endif
				}
				p_cla += len_per_categ;
				f_cla += len_per_categ;
				s_cla += len_per_categ;
			}
			p_cla += post_skip;
			f_cla += post_skip;
			s_cla += post_skip;
#			if SEPARATE_RESCALER_ARRAY
				rescale_count_ptr += n_skipped_categ;
				if (has_f_rescale)
					f_prev_rescale_count += n_skipped_categ;
				if (has_s_rescale)
					s_prev_rescale_count += n_skipped_categ;
#			endif
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings  = (rescaling_done ? par->rescale_swap_space : 0L);
#	endif
	par->n_edges_since_rescaling_check = 0;
	return 1;
}

static int do_internal_cla_no_rescale(const CLAObj * fd, const double *** f_pmat_arr, const CLAObj* sd, const double *** s_pmat_arr, CLAObj* par, const CalcContext context) {
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const double ** f_categ_pmat;
	const double * f_categ_p;
	const double ** s_categ_pmat;
	const double * s_categ_p;
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned site, categ, to_state, from_state, subset_ind, beg_site, end_site;
	cla_float_t f_des_prob, s_des_prob;
#	if SEPARATE_RESCALER_ARRAY
		const int has_f_rescale = (fd->n_rescalings ? 1 : 0);
		const int has_s_rescale = (sd->n_rescalings ? 1 : 0);
		par->n_rescalings = (has_f_rescale == 0 && has_s_rescale == 0 ? 0L : par->rescale_swap_space);
#	endif
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {

		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t *f_cla = (const cla_float_t*) (fd->cla + start_offset);
		const cla_float_t *s_cla = (const cla_float_t*) (sd->cla + start_offset);
	

		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;
		const double *** s_pmat_arr_o = s_pmat_arr + comp_offset;

		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

		for (site = beg_site; site < end_site; ++site) {
			for (categ = categ_beg; categ < categ_end; ++categ) {
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				s_categ_pmat = s_pmat_arr_o[categ];
				s_categ_p = *s_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					s_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
						s_des_prob += s_cla[to_state]*(*s_categ_p++);
					}
					p_cla[from_state] = f_des_prob*s_des_prob;
				}
#				if ! SEPARATE_RESCALER_ARRAY
						p_cla[nStates] = f_cla[nStates] + s_cla[nStates];
#				endif
				p_cla += len_per_categ;
				f_cla += len_per_categ;
				s_cla += len_per_categ;
			}
			p_cla += post_skip;
			f_cla += post_skip;
			s_cla += post_skip;
		}
#		if SEPARATE_RESCALER_ARRAY
			if (has_f_rescale || has_s_rescale) {
				const unsigned subset_resc_offset = subset_cla_offsets;
				const rescale_history_t * s_rescale_count;
				rescale_history_t * par_rescale_count = par->rescale_swap_space + subset_resc_offset;
				unsigned i, j;
				if (has_f_rescale) {
					const rescale_history_t * f_rescale_count = ((const rescale_history_t *) fd->n_rescalings) + subset_resc_offset;
					if (has_s_rescale){
						s_rescale_count = ((const rescale_history_t *) sd->n_rescalings) + subset_resc_offset;
						for (i = beg_site; i < end_site; ++i) {
							for (j = categ_beg; j < categ_end; ++j) {
								par_rescale_count[j] = f_rescale_count[j] + s_rescale_count[j];
							}
							f_rescale_count += n_categ;
							s_rescale_count += n_categ;
							prev_rescale_count += n_categ;
						}
					}
					else {
						for (i = beg_site; i < end_site; ++i) {
							for (j = categ_beg; j < categ_end; ++j) {
								par_rescale_count[j] = f_rescale_count[j];
							}
							f_rescale_count += n_categ;
							prev_rescale_count += n_categ;
						}
					}
				}
				else if (has_s_rescale){
					s_rescale_count = ((const rescale_history_t *) sd->n_rescalings) + subset_resc_offset;
					for (i = beg_site; i < end_site; ++i) {
						for (j = categ_beg; j < categ_end; ++j) {
							par_rescale_count[j] = s_rescale_count[j];
						}
						s_rescale_count += n_categ;
						prev_rescale_count += n_categ;
					}
				}
			}
#		endif
	}
	par->n_edges_since_rescaling_check = 2 + fd->n_edges_since_rescaling_check + sd->n_edges_since_rescaling_check;
	return 1;
}

int do_internal_cla(const CLAObj * fd, const PMatArrayObj *fmat, const CLAObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(sd);
	assert(smat);
	assert(par);
	assert(fd->cso.total_n_sites == sd->cso.total_n_sites);
	assert(fd->cso.total_n_sites == par->cso.total_n_sites);
	assert(fd->nStates == sd->nStates);
	assert(fd->nStates == par->nStates);
	assert(fd->cso.n_categ_or_subs == sd->cso.n_categ_or_subs);
	assert(fd->cso.n_categ_or_subs == par->cso.n_categ_or_subs);
	const double *** f_pmat_arr = (const double ***)fmat->pMatrix;
	const double *** s_pmat_arr = (const double ***)smat->pMatrix;
	if ((fd->n_edges_since_rescaling_check + sd->n_edges_since_rescaling_check + 2) < context.rescale_threshold)
		return do_internal_cla_no_rescale(fd, f_pmat_arr, sd, s_pmat_arr, par, context);
	return do_internal_cla_rescale(fd, f_pmat_arr, sd, s_pmat_arr, par, context);
}

static int do_add_leaf_to_cla_recalc_rescale(const LeafDataObj* fd, CLAObj* par, const CalcContext context) {
	const double ** f_mat;
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	const double * f_row;
	double max_p, real_rescale;
#	if SEPARATE_RESCALER_ARRAY
#		if INT_INCR_RESCALING
			double real_rescale_ln;
#		endif
		int rescaling_done = (par->n_rescalings == 0L ? 0 : 1);
		rescale_history_t rcount;
#	else
		double real_rescale_ln;
#	endif
	unsigned subset_ind, beg_site, end_site;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		 f_mat = (const double **)fd->pmat_columns[subset_ind];
		
		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;

		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * f_states = (const int *) fd->ssind;
		f_states += subset_char_offsets ;
	

		const unsigned tpm_offset = nStates*(categ_beg);
		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned n_skipped_categ = trailing_categ + categ_beg;
		const unsigned post_skip = len_per_categ*(n_skipped_categ);
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

#		if SEPARATE_RESCALER_ARRAY
			const unsigned resc_offset = (par->cso.subset_cla_offsets[subset_ind]) + beg_categ;
			rescale_history_t * rescale_count_ptr = par->rescale_swap_space + resc_offset;
			if (has_f_rescale) {
				f_prev_rescale_count = (const rescale_history_t *) fd->n_rescalings;
				f_prev_rescale_count += resc_offset;
			}
#		endif

		unsigned state, site, categ;
		for (site = beg_site; site < end_site; ++site) {
			f_row = f_mat[*f_states++];
			f_row += tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				max_p = 0.0;
				for (state = 0; state < nStates; ++state) {
					p_cla[state] *= *f_row++;
					if (p_cla[state] > max_p)
						max_p = p_cla[state];
				}
				if (max_p < K_RESCALE_THRESHOLD) {
					if (max_p < K_RESCALE_ERROR_THRESHOLD)
						return UNDERFLOW_ERROR_CODE;
#					if SEPARATE_RESCALER_ARRAY
#						if INT_INCR_RESCALING
							rcount = 1 + (rescale_history_t)log(max_p);
							real_rescale_ln = (double) rcount;
							real_rescale = exp(real_rescale_ln);
#						elif CONSTANT_STEP_RESCALING
							real_rescale = K_RESCALING_FACTOR;
							rcount = 1;
#						else
#							error "unknown RESCALING type"
#						endif
						rescaling_done = 1;
						*rescale_count_ptr++ += rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						real_rescale_ln = log(max_p);
						real_rescale = max_p;
						p_cla[nStates] += real_rescale_ln;
#					endif
					for (state = 0; state < nStates; ++state)
						p_cla[state] /= real_rescale;
				}
				p_cla += len_per_categ;
			}
			p_cla += post_skip;
#			if SEPARATE_RESCALER_ARRAY
				rescale_count_ptr += n_skipped_categ;
#			endif
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings  = (rescaling_done ? par->rescale_swap_space : 0L);
#	endif
	par->n_edges_since_rescaling_check = 0;
	return 1;
}


static int do_add_leaf_to_cla_recalc_no_rescale(const LeafDataObj* fd, CLAObj* par, const CalcContext context) {
	const double ** f_mat;
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned subset_ind;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {
		f_mat = (const double **)fd->pmat_columns[subset_ind];
		
		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;

		const unsigned subset_char_offsets = (subset_ind == 0 ? 0 : char_offsets[subset_ind]);
		const int * f_states = (const int *) fd->ssind;
		f_states += subset_char_offsets ;
	
		const unsigned tpm_offset = nStates*(categ_beg);
		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);

		unsigned beg_site, end_site;
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}
	
		const double * f_row;
		unsigned site, categ, state;
		for (site = beg_site; site < end_site; ++site) {
			f_row = f_mat[*f_states++];
			f_row += tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				for (state = 0; state < nStates; ++state) {
					*p_cla++ *= (*f_row++) ;
				}
			}
			p_cla += post_skip;
		}
	}
	par->n_edges_since_rescaling_check += 1;
	return 1;
}

static int do_add_leaf_to_cla_recalc(const LeafDataObj* fd, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(par);
	assert(fd->nStates == par->nStates);
	/*
		assert(fd->n_categ == par->n_categ);
	*/
	if (2 < context.rescale_threshold)
		return do_add_leaf_to_cla_recalc_rescale(fd, par, context);
	return do_add_leaf_to_cla_recalc_no_rescale(fd, par, context);
}

int do_add_leaf_to_cla(LeafDataObj* fd, const PMatArrayObj *fmat, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(par);
	if ((fd->pmat_obj_ptr != (void *) fmat) || (fd->calcTimeStamp != fmat->calcTimeStamp)) {
		if (!assign_pmat_to_leaf(fd, fmat, context, par->cso))
			return 0;
	}
	return do_add_leaf_to_cla_recalc(fd, par, context);
}

static int do_add_internal_to_cla_rescale(const CLAObj * fd, const double *** f_pmat_arr, CLAObj* par, const CalcContext context) {
	const double ** f_categ_pmat;
	const double * f_categ_p;
	
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	double p, max_p, real_rescale;
#	if SEPARATE_RESCALER_ARRAY
#		if INT_INCR_RESCALING
			double real_rescale_ln;
#		endif
		rescale_history_t rcount;
		rescale_history_t prev_rcount;
		const rescale_history_t * f_prev_rescale_count = (const rescale_history_t *) fd->n_rescalings;
		const int has_f_rescale = (fd->n_rescalings == 0L ? 0 : 1);
		int rescaling_done = (has_f_rescale > 0 ? 1 : 0);
#	else
		double real_rescale_ln;
#	endif
	unsigned site, categ, from_state, to_state, subset_ind, beg_site, end_site;
	cla_float_t f_des_prob;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {

		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t * f_cla = (const cla_float_t*) (fd->cla + start_offset);
	

		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned n_skipped_categ = trailing_categ + categ_beg;
		const unsigned post_skip = len_per_categ*(n_skipped_categ);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;

		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

#		if SEPARATE_RESCALER_ARRAY
			const unsigned resc_offset = (par->cso.subset_cla_offsets[subset_ind]) + beg_categ;
			rescale_history_t * rescale_count_ptr = par->rescale_swap_space + resc_offset;
			if (has_f_rescale) {
				f_prev_rescale_count += (const rescale_history_t *) fd->n_rescalings;
				f_prev_rescale_count += resc_offset;
			}
#		endif
		for (site = beg_site; site < end_site; ++site) {
			for (categ = categ_beg; categ < categ_end; ++categ) {
				max_p = 0.0;
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					p = f_des_prob*p_cla[from_state];
					if (p > max_p)
						max_p = p;
					p_cla[from_state] = p;
				}
#				if SEPARATE_RESCALER_ARRAY
					prev_rcount = (has_f_rescale ? *f_prev_rescale_count++ : 0);
#				endif
				if (max_p < K_RESCALE_THRESHOLD) {
					if (max_p < K_RESCALE_ERROR_THRESHOLD)
						return UNDERFLOW_ERROR_CODE;
#					if SEPARATE_RESCALER_ARRAY
#						if INT_INCR_RESCALING
							rcount = 1 + (rescale_history_t)log(max_p);
							real_rescale_ln = (double) rcount;
							real_rescale = exp(real_rescale_ln);
#						elif CONSTANT_STEP_RESCALING
							real_rescale = K_RESCALING_FACTOR;
							rcount = 1;
#						else
#							error "unknown RESCALING type"
#						endif
						rescaling_done = 1;
						*rescale_count_ptr++ += rcount + prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						real_rescale_ln = log(max_p);
						real_rescale = max_p;
						p_cla[nStates] += real_rescale_ln + f_cla[nStates];
#					endif
					for (from_state = 0; from_state < nStates; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ += prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[nStates] += f_cla[nStates];
#					endif
				}
				p_cla += len_per_categ;
				f_cla += len_per_categ;
			}
			p_cla += post_skip;
			f_cla += post_skip;
#			if SEPARATE_RESCALER_ARRAY
				rescale_count_ptr += n_skipped_categ; 
				if (has_f_rescale)
					f_prev_rescale_count += n_skipped_categ;
#			endif
		}
	}
#	if SEPARATE_RESCALER_ARRAY
		par->n_rescalings  = (rescaling_done ? par->rescale_swap_space : 0L);
#	endif
	par->n_edges_since_rescaling_check = 0;
	return 1;
}


static int do_add_internal_to_cla_no_rescale(const CLAObj * fd, const double *** f_pmat_arr, CLAObj* par, const CalcContext context) {
	const double ** f_categ_pmat;
	const double * f_categ_p;
	
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned nStates = par->nStates;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned site, categ, to_state, from_state, subset_ind, beg_site, end_site;
	cla_float_t  f_des_prob;
	for (subset_ind = subset_beg; subset_ind < subset_end; ++subset_ind) {

		const unsigned n_categ = (par->cso.n_categ_arr == 0L ? par->cso.n_categ_or_subs : par->cso.n_categ_arr[subset_ind]);
		const unsigned categ_beg = (context.categ_beg > n_categ ? n_categ : context.categ_beg);
		const unsigned categ_end = (context.categ_end > n_categ ? n_categ : context.categ_end);

		const unsigned subset_cla_offsets = (subset_ind == 0 ? 0 : par->cso.subset_cla_offsets[subset_ind]);
		const unsigned start_offset = len_per_categ*(subset_cla_offsets + categ_beg);
		cla_float_t * p_cla = par->cla + start_offset;
		const cla_float_t *f_cla = (const cla_float_t*) (fd->cla + start_offset);
	

		const unsigned trailing_categ = n_categ - categ_end;
		const unsigned post_skip = len_per_categ*(trailing_categ + categ_beg);
		const unsigned comp_offset = (subset_ind == 0 ? 0 : par->cso.subset_component_offset[subset_ind]);
		const double *** f_pmat_arr_o = f_pmat_arr + comp_offset;
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}

		for (site = beg_site; site < end_site; ++site) {
			for (categ = categ_beg; categ < categ_end; ++categ) {
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < nStates; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < nStates; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					p_cla[from_state] *= f_des_prob;
				}
#			if ! SEPARATE_RESCALER_ARRAY
						p_cla[nStates] += f_cla[nStates];
#			endif
				p_cla += len_per_categ;
				f_cla += len_per_categ;
			}
			p_cla += post_skip;
			f_cla += post_skip;
		}
#		if SEPARATE_RESCALER_ARRAY
			if (fd->n_rescalings) {
				const unsigned subset_resc_offset = subset_cla_offsets;
				const rescale_history_t * prev_rescale_count = ((const rescale_history_t *) fd->n_rescalings) + subset_resc_offset;
				rescale_history_t * par_rescale_count =  par->rescale_swap_space + subset_resc_offset;
				par->n_rescalings =  par->rescale_swap_space;
				unsigned i, j;
				for (i = beg_site; i < end_site; ++i) {
					for (j = categ_beg; j < categ_end; ++j) {
						par_rescale_count[j] += prev_rescale_count[j];
					}
					par_rescale_count += n_categ;
					prev_rescale_count += n_categ;
				}
			}
#		endif
	}
	par->n_edges_since_rescaling_check = 1 + fd->n_edges_since_rescaling_check;
	return 1;
}

int do_add_internal_to_cla(const CLAObj* fd, const  PMatArrayObj *fmat, CLAObj* par,  const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(par);
	assert(fd->cso.total_n_sites == par->cso.total_n_sites);
	assert(fd->nStates == par->nStates);
	assert(fd->cso.max_n_categ == par->cso.max_n_categ);
	const double *** f_pmat_arr = (const double ***)fmat->pMatrix;
	if ((fd->n_edges_since_rescaling_check + par->n_edges_since_rescaling_check + 1) < context.rescale_threshold)
		return do_add_internal_to_cla_rescale(fd, f_pmat_arr, par, context);
	return do_add_internal_to_cla_no_rescale(fd, f_pmat_arr, par, context);
}
double internal_edge_ln_likelihood(const CLAObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj *full_la, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(par);
	assert(full_la);
	CLAObj* full_cla = full_la->full_cla;
	assert(full_cla);
	CalcContext c = context;
	c.rescale_threshold = UINT_MAX; // no point in rescaling the last multiplication.
	copy_cla(full_cla, par);
	do_add_internal_to_cla(fd, fmat, full_cla, c);
	return partitioned_ln_likelihood(full_la, 0L);
}
double terminal_edge_ln_likelihood(LeafDataObj* fd, const PMatArrayObj *fmat, const CLAObj* par, FullLAObj * full_la, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(par);
	assert(full_la);
	CLAObj* full_cla = full_la->full_cla;
	assert(full_cla);
	CalcContext c = context;
	c.rescale_threshold = UINT_MAX; // no point in rescaling the last multiplication.
	copy_cla(full_cla, par);
	do_add_leaf_to_cla(fd, fmat, full_cla, c);
	return partitioned_ln_likelihood(full_la, 0L);
}
double ln_likelihood(FullLAObj *full_la) {
	return partitioned_ln_likelihood(full_la, 0L);
}

double partitioned_ln_likelihood(FullLAObj *full_la, double * subsetLnL) {
	assert(full_la);
	const unsigned total_n_sites = full_la->cso.total_n_sites;
	const unsigned nStates = full_la->nStates;
	const unsigned * char_offsets = full_la->cso.subset_char_offsets;
	const unsigned n_subsets = (full_la->cso.n_categ_arr == 0L ? 1 : full_la->cso.n_categ_or_subs);
	unsigned subset_ind, beg_site, end_site, n_categ;
	double total_ln_L = 0.0;
	double cumul_subset_ln_L = 0.0;
	unsigned max_n_categ = full_la->cso.max_n_categ;
	const double ** const sf_mat = (const double ** const) (full_la->state_categ_freqs);
	assert(sf_mat);
	const double * sf_arr;
	cla_float_t * lnL_wts = full_la->pat_lnL_wts;
	const CLAObj * const cla_obj = (const CLAObj * const ) full_la->full_cla;
	assert(cla_obj);
	const cla_float_t * cla = (const cla_float_t *) cla_obj->cla;
	assert(cla);
	double * const f_mult = full_la->scratch;
	assert(f_mult);
	double * const categ_like = f_mult + (nStates*max_n_categ);
#	if SEPARATE_RESCALER_ARRAY
		const rescale_history_t * n_rescalings = (const rescale_history_t *) cla_obj->n_rescalings;
		assert(site_rescaler);
		rescale_history_t min_rescale = 0;
		rescale_history_t max_rescale = 0;
		rescale_history_t curr_rescale = 0;
#		if CONSTANT_STEP_RESCALING
			const double kLnRescalingFactor = getLnRescalingFactor();
#		endif
#	else
		double * site_rescale_ln = categ_like + max_n_categ;
		assert(site_rescale_ln);
		const double kNoRescalingTol = -1e-10;
		double min_rescale = 0.0;
		double max_rescale = 0.0;
#	endif
	unsigned min_rescale_cat = 0;
	unsigned max_rescale_cat = 0;
	double * curr_categ_like;
	const double * curr_mult;
	unsigned site, categ, state;
	double site_ln_L, site_L, wt, cat_freq;
	/* calculate the (categ freq)*(state_freq) that is the multiplier of each
		conditional likelihood.
		We use the scratch field as the space to hold this.
	*/
	for (subset_ind = 0; subset_ind < n_subsets; ++subset_ind) {
		cumul_subset_ln_L = 0.0;
		n_categ = (full_la->cso.n_categ_arr == 0L ? full_la->cso.n_categ_or_subs : full_la->cso.n_categ_arr[subset_ind]);
		const unsigned compOffset = (full_la->cso.subset_component_offset == 0L ? 0 : full_la->cso.subset_component_offset[subset_ind]);
		for (categ = 0; categ < n_categ; ++categ) {
			sf_arr = sf_mat[compOffset + categ];
			cat_freq = sf_arr[nStates];
			for (state = 0; state < nStates; ++state) {
				f_mult[nStates*categ + state] = cat_freq*sf_arr[state];
			}
		}
		if (n_subsets > 1) {
			beg_site = char_offsets[subset_ind];
			end_site = (subset_ind == n_subsets - 1 ? total_n_sites : char_offsets[subset_ind + 1]);
		}
		else {
			beg_site = 0;
			end_site = total_n_sites;
		}
			
		CPHYPROB_DEBUG_PRINTF("categ likes, site lnL, site weight:\n");
		for (site = beg_site; site < end_site; ++site) {
			site_ln_L = 0.0;
			site_L = 0.0;
			curr_mult = (const double *) f_mult;
			curr_categ_like = categ_like;
#			if SEPARATE_RESCALER_ARRAY
				min_rescale = max_rescale =  0;
#			else
				min_rescale = 0.0;
				max_rescale =  -1e300;
#			endif
			min_rescale_cat = 0;
			max_rescale_cat = 0;
			CPHYPROB_DEBUG_PRINTF1("%d", site);
			for (categ = 0; categ < n_categ; ++categ) {
				*curr_categ_like = 0.0;
				for (state = 0; state < nStates; ++state) {
					*curr_categ_like += (*curr_mult++)*(*cla++);
				}
				CPHYPROB_DEBUG_PRINTF1("\t%g", *curr_categ_like);
				site_L += *curr_categ_like++; /*if there is rescaling, this summation to site_ln_L will be wrong, but will be overwritten below*/
#				if SEPARATE_RESCALER_ARRAY
					if (n_rescalings) {
						curr_rescale = n_rescalings[categ];
						if (categ == 0) {
							min_rescale = max_rescale = curr_rescale;
						}
						else {
							if (curr_rescale < min_rescale) {
								min_rescale = curr_rescale;
								min_rescale_cat = categ;
							}
							else if (curr_rescale > max_rescale) {
								max_rescale = curr_rescale;
								max_rescale_cat = categ;
							}
						}
					}
#				else /* end SEPARATE_RESCALER_ARRAY */
					if (*curr_mult < kNoRescalingTol) {
						site_rescale_ln[categ] = *cla;
						if (*curr_mult < min_rescale) {
							min_rescale = *cla;
							min_rescale_cat = categ;
						}
						else if (*curr_mult > max_rescale) {
							max_rescale = *cla;
							max_rescale_cat = categ;
						}
					}
					else {
						site_rescale_ln[categ] = 0.0;
						max_rescale = 0.0;
						max_rescale_cat = categ;
					}
					cla++;
#				endif /* end ! SEPARATE_RESCALER_ARRAY */
			}
#			if SEPARATE_RESCALER_ARRAY
#				if INT_INCR_RESCALING
					if (min_rescale < 0) {
						site_L = 0.0;
						for (categ = 0; categ < n_categ; ++categ) {
							site_L  += exp(n_rescalings[categ]-max_rescale) * curr_categ_like[categ];
						}
						site_ln_L = max_rescale + log(site_L);
					}
					else
						site_ln_L = log(site_L);
#				elif CONSTANT_STEP_RESCALING
					if (max_rescale > 0) {
						site_L = 0.0;
						for (categ = 0; categ < n_categ; ++categ) {
							if (n_rescalings[categ] == min_rescale)
								site_L  += curr_categ_like[categ];
							else
								site_L  += pow(K_RESCALING_FACTOR, min_rescale-n_rescalings[categ])*curr_categ_like[categ];
						}
						site_ln_L = (min_rescale*kLnRescalingFactor) + log(site_L);
					}
					else
						site_ln_L = log(site_L);
#				else
#					error "unknown RESCALING type"
#				endif
				if (n_rescalings)
					n_rescalings += n_categ;
#			else
				if (min_rescale < kNoRescalingTol) {
						site_L = 0.0;
						for (categ = 0; categ < n_categ; ++categ) {
							site_L  += exp(site_rescale_ln[categ]-max_rescale) * curr_categ_like[categ];
						}
						site_ln_L = max_rescale + log(site_L);
				}
				else {
					site_ln_L = log(site_L);
				}
#			endif
			wt = lnL_wts[1];
			lnL_wts[0] = site_ln_L;
			lnL_wts += 2;
			cumul_subset_ln_L += wt * site_ln_L;
			CPHYPROB_DEBUG_PRINTF2("\t%g\t%f\n", site_ln_L, wt);
		}
		if (subsetLnL)
			subsetLnL[subset_ind] = cumul_subset_ln_L;
		total_ln_L += cumul_subset_ln_L;
		
	}
	return total_ln_L;
}


LeafDataObj* leaf_data_ss_partitioned_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl)  {
	assert(n_sites > 0);
	assert(n_categ > 0);
	assert(ssl);
	CPHYPROB_DEBUG_PRINTF("In leaf_data_new\n");
	const unsigned len_pmat_columns_row = ssl->nStates*n_categ;
	LeafDataObj * leaf_data = PyObject_New(LeafDataObj, &leaf_data_type);
	if (leaf_data) {
		Py_INCREF(ssl);
		leaf_data->state_set_lookup = ssl;
		leaf_data->nStates = ssl->nStates;
		leaf_data->n_state_sets = ssl->n_state_sets;
		leaf_data->state_lookup = ssl->state_lookup;
		leaf_data->ssind = 0L;
		leaf_data->pmat_columns = 0L;
		leaf_data->calcTimeStamp = -1;
		leaf_data->pmat_obj_ptr = (void*)0L;
		
		leaf_data->pmat_columns = allocateDbl3DMatrix(n_sites, ssl->n_state_sets, len_pmat_columns_row);
		if (leaf_data->pmat_columns == 0L)
			goto errorExit;
		leaf_data->ssind = (int *) malloc(n_sites*sizeof(int));
		if (leaf_data->ssind == 0L)
			goto errorExit;
	}
	return leaf_data;
errorExit:
	CPHYPROB_DEBUG_PRINTF("In leaf_data_new errorExit\n");
	leaf_data_dtor(leaf_data);
	return 0L;
}

LeafDataObj* leaf_data_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl)  {
	assert(n_sites > 0);
	assert(n_categ > 0);
	assert(ssl);
	CPHYPROB_DEBUG_PRINTF("In leaf_data_new\n");
	const unsigned len_pmat_columns_row = ssl->nStates*n_categ;
	LeafDataObj * leaf_data = PyObject_New(LeafDataObj, &leaf_data_type);
	if (leaf_data) {
		Py_INCREF(ssl);
		leaf_data->state_set_lookup = ssl;
		leaf_data->nStates = ssl->nStates;
		leaf_data->n_state_sets = ssl->n_state_sets;
		leaf_data->state_lookup = ssl->state_lookup;
		leaf_data->ssind = 0L;
		leaf_data->pmat_columns = 0L;
		leaf_data->calcTimeStamp = -1;
		leaf_data->pmat_obj_ptr = (void*)0L;
		leaf_data->pmat_columns = allocateDbl3DMatrix(1, ssl->n_state_sets, len_pmat_columns_row);
		if (leaf_data->pmat_columns == 0L)
			goto errorExit;
		leaf_data->ssind = (int *) malloc(n_sites*sizeof(int));
		if (leaf_data->ssind == 0L)
			goto errorExit;
	}
	return leaf_data;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In leaf_data_new errorExit\n");
		leaf_data_dtor(leaf_data);
		return 0L;
}

/*
	Allocates a CLAObj object for an unpartitioned calculation that mixes n_categ models. 
*/
CLAObj * cla_new(unsigned n_sites, unsigned nStates,  unsigned n_categ)  {
	unsigned cla_len;
#	if SEPARATE_RESCALER_ARRAY
		unsigned rescale_len;
#	endif
	assert(n_sites > 0);
	assert(nStates > 0);
	assert(n_categ > 0);
	CPHYPROB_DEBUG_PRINTF("In cla_new\n");
	CLAObj * cla = PyObject_New(CLAObj, &cla_type);
	if (cla) {
		cla->cso.total_n_sites = n_sites;
		cla->nStates = nStates;
		cla->cso.n_categ_or_subs = n_categ;
		cla->cso.max_n_categ = n_categ;
		cla->cso.n_categ_arr = 0L;
		cla->cso.subset_component_offset = 0L;
		cla->cso.subset_char_offsets = 0L;
		cla->cso.subset_cla_offsets = 0L;
#		if SEPARATE_RESCALER_ARRAY
			cla->len_per_categ = n_sites;
			rescale_len = n_sites*n_categ;
#		else
			cla->len_per_categ = nStates + 1;/*The +1 is for the rescale ln*/
#		endif
		cla_len = n_sites*n_categ*(cla->len_per_categ);
		cla->cla = 0L;
#		if SEPARATE_RESCALER_ARRAY
			cla->n_rescalings = 0L;
			cla->rescale_swap_space = 0L;
			cla->rescale_swap_space = (rescale_history_t *) malloc(rescale_len*sizeof(rescale_history_t));
			if (cla->rescale_swap_space == 0L)
				goto errorExit;
#		endif /*SEPARATE_RESCALER_ARRAY*/
		cla->cla = (cla_float_t *) malloc(cla_len*sizeof(cla_float_t));
		if (cla->cla == 0L)
			goto errorExit;
	}
	return cla;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In cla_new errorExit\n");
		cla_dtor(cla);
		return 0L;
}


/* Allocates the field of a CategSubsetObj object for a model in which each site gets its own model 
	(and there are no mixtures used).
	returns 0 if there is an error.
*/
int new_cso_ss_partitioned(CategSubsetObj * cso, unsigned n_sites, unsigned nStates) {
	unsigned i;
	cso->n_categ_or_subs = n_sites;
	cso->max_n_categ = 1;
	cso->n_categ_arr = 0L;
	cso->subset_component_offset = 0L;
	cso->subset_char_offsets = 0L;
	cso->subset_cla_offsets = 0L;
	cso->n_categ_arr = (unsigned *) malloc(n_sites*sizeof(unsigned));
	if (cso->n_categ_arr == 0L)
		return 0;
	cso->subset_component_offset = (unsigned *) malloc(n_sites*sizeof(unsigned));
	if (cso->subset_component_offset == 0L)
		return 0;
	cso->subset_char_offsets = (unsigned *) malloc(n_sites*sizeof(unsigned));
	if (cso->subset_char_offsets == 0L)
		return 0;
	cso->subset_cla_offsets = (unsigned *) malloc(n_sites*sizeof(unsigned));
	if (cso->subset_cla_offsets == 0L)
		return 0;
	for (i = 0; i < n_sites; ++i) {
		cso->n_categ_arr[i] = 1;
		cso->subset_component_offset[i] = i;
		cso->subset_char_offsets[i] = i;
		cso->subset_cla_offsets[i] = i;
	}
	
	return 1;
}

void cso_dtor(CategSubsetObj * cso) {
	if (cso == 0L)
		return;
	if (cso->n_categ_arr)
		free(cso->n_categ_arr);
	if (cso->subset_component_offset)
		free(cso->subset_component_offset);
	if (cso->subset_char_offsets)
		free(cso->subset_char_offsets);
	if (cso->subset_cla_offsets)
		free(cso->subset_cla_offsets);
	
}


/* Allocates a claObj object for a model in which each site gets its own model 
	(and there are no mixtures used).
*/
CLAObj* cla_ss_partitioned_new(unsigned n_sites, unsigned nStates) {
	unsigned cla_len;
	unsigned n_categ = 1;
#	if SEPARATE_RESCALER_ARRAY
		unsigned rescale_len;
#	endif
	assert(n_sites > 0);
	assert(nStates > 0);
	assert(n_categ > 0);
	CPHYPROB_DEBUG_PRINTF("In cla_ss_partitioned_new\n");
	CLAObj * cla = PyObject_New(CLAObj, &cla_type);
	if (cla) {
		cla->cso.total_n_sites = n_sites;
		cla->nStates = nStates;
		if (new_cso_ss_partitioned(&(cla->cso), n_sites, nStates) == 0)
			goto errorExit;
#		if SEPARATE_RESCALER_ARRAY
			cla->len_per_categ = n_sites;
			rescale_len = n_sites*n_categ;
#		else
			cla->len_per_categ = nStates + 1;/*The +1 is for the rescale ln*/
#		endif
		cla_len = n_sites*n_categ*(cla->len_per_categ);
		cla->cla = 0L;
#		if SEPARATE_RESCALER_ARRAY
			cla->n_rescalings = 0L;
			cla->rescale_swap_space = 0L;
			cla->rescale_swap_space = (rescale_history_t *) malloc(rescale_len*sizeof(rescale_history_t));
			if (cla->rescale_swap_space == 0L)
				goto errorExit;
#		endif /*SEPARATE_RESCALER_ARRAY*/
		cla->cla = (cla_float_t *) malloc(cla_len*sizeof(cla_float_t));
		if (cla->cla == 0L)
			goto errorExit;
	}
	return cla;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In cla_ss_partitioned_new errorExit\n");
		cla_dtor(cla);
		return 0L;
}

/* Allocates a FullLAObj object for a model in which each site gets its own model 
	(and there are no mixtures used).
*/
FullLAObj* full_la_ss_partitioned_new(unsigned n_sites, unsigned nStates) {
	CPHYPROB_DEBUG_PRINTF("In full_la_ss_partitioned_new\n");
	unsigned i;
	
	const unsigned lnL_wt_len = 2*(n_sites);
	FullLAObj * full_la = PyObject_New(FullLAObj, &full_la_type);
	if (full_la) {
		full_la->cso.total_n_sites = n_sites;
		full_la->nStates = nStates;
		full_la->full_cla  = 0L;
		full_la->pat_lnL_wts  = 0L;
		full_la->state_categ_freqs  = 0L;
		if (new_cso_ss_partitioned(&(full_la->cso), n_sites, nStates) == 0)
			goto errorExit;
		full_la->full_cla = cla_ss_partitioned_new(n_sites, nStates);
		if (full_la->full_cla == 0L)
			goto errorExit;
		full_la->pat_lnL_wts  = (cla_float_t *) malloc(lnL_wt_len*sizeof(cla_float_t));
		if (full_la->pat_lnL_wts == 0L)
			goto errorExit;
		for (i = 0; i < n_sites; ++i)
			full_la->pat_lnL_wts[1 + (2*i)] = 1.0;
		full_la->scratch = (double *) malloc(2*nStates*sizeof(double));
		if (full_la->scratch == 0L)
			goto errorExit;
		full_la->state_categ_freqs = allocateDblMatrix(n_sites, nStates + 1);
		if (full_la->state_categ_freqs == 0L)
			goto errorExit;
	}
	return full_la;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In full_la_ss_partitioned_new errorExit\n");
		full_la_dtor(full_la);
		return 0L;
}

/*
	Allocates a FullLAObj object for an unpartitioned calculation that mixes n_categ models. 
*/
FullLAObj * full_la_new(unsigned n_sites, unsigned nStates,  unsigned n_categ)  {
	CPHYPROB_DEBUG_PRINTF("In full_la_new\n");
	unsigned i;
	const unsigned lnL_wt_len = 2*(n_sites);
	FullLAObj * full_la = PyObject_New(FullLAObj, &full_la_type);
	if (full_la) {
		full_la->cso.total_n_sites = n_sites;
		full_la->nStates = nStates;
		full_la->cso.n_categ_or_subs = n_categ;
		full_la->cso.max_n_categ = n_categ;
		full_la->full_cla  = 0L;
		full_la->pat_lnL_wts  = 0L;
		full_la->state_categ_freqs  = 0L;
		full_la->cso.n_categ_arr = 0L;
		full_la->cso.subset_component_offset = 0L;
		full_la->cso.subset_char_offsets = 0L;
		full_la->cso.subset_cla_offsets = 0L;
		full_la->full_cla = cla_new(n_sites, nStates,  n_categ);
		if (full_la->full_cla == 0L)
			goto errorExit;
		full_la->pat_lnL_wts  = (cla_float_t *) malloc(lnL_wt_len*sizeof(cla_float_t));
		if (full_la->pat_lnL_wts == 0L)
			goto errorExit;
		for (i = 0; i < n_sites; ++i)
			full_la->pat_lnL_wts[1 + (2*i)] = 1.0;
		full_la->scratch = (double *) malloc(2*nStates*n_categ*sizeof(double));
		if (full_la->scratch == 0L)
			goto errorExit;
		full_la->state_categ_freqs = allocateDblMatrix(n_categ, nStates + 1);
		if (full_la->state_categ_freqs == 0L)
			goto errorExit;
	}
	return full_la;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In full_la_new errorExit\n");
		full_la_dtor(full_la);
		return 0L;
}

DSCTModelObj* dsctModelNew(unsigned dim)  {
	assert(dim > 1);
	CPHYPROB_DEBUG_PRINTF("In cdsctm_ctor\n");
	unsigned arr_len = (3 + dim*dim) * dim;
	DSCTModelObj * dsct_model = PyObject_New(DSCTModelObj, &dsct_model_type);
	if (dsct_model) {
		dsct_model->dim = dim;
		dsct_model->eigenCalcIsDirty = 1;
		dsct_model->mat_del_handle = 0L;
		dsct_model->arr_del_handle = 0L;
		dsct_model->qMat = 0L;
		dsct_model->eigenVectors = 0L;
		dsct_model->invEigenVectors = 0L;
		dsct_model->workMat = 0L;
		dsct_model->eigenValues = 0L;
		dsct_model->imEigenValues = 0L;
		dsct_model->cijk = 0L;
		dsct_model->dWork = 0L;
		dsct_model->iWork = 0L;
		dsct_model->mat_del_handle = allocateDbl3DMatrix(4, dim, dim);
		if (dsct_model->mat_del_handle == 0L)
			goto errorExit;
		dsct_model->arr_del_handle = (double *)malloc(arr_len*sizeof(double));
		if (dsct_model->arr_del_handle == 0L)
			goto errorExit;
		dsct_model->iWork = (int *)malloc(dim*sizeof(int));
		if (dsct_model->iWork == 0L)
			goto errorExit;
		dsct_model->qMat = dsct_model->mat_del_handle[0];
		dsct_model->eigenVectors = dsct_model->mat_del_handle[1];
		dsct_model->invEigenVectors = dsct_model->mat_del_handle[2];
		dsct_model->workMat = dsct_model->mat_del_handle[3];
		dsct_model->eigenValues = dsct_model->arr_del_handle;
		dsct_model->imEigenValues = dsct_model->eigenValues + dim;
		dsct_model->dWork = dsct_model->imEigenValues + dim;
		dsct_model->cijk = dsct_model->dWork + dim;
	}
	return dsct_model;
	errorExit:
		CPHYPROB_DEBUG_PRINTF("In cdsctm_ctor errorExit\n");
		cdsctm_dtor(dsct_model);
		return 0L;
}



void leaf_data_dtor(LeafDataObj * leaf_data) {
	CPHYPROB_DEBUG_PRINTF("In leaf_data_dtor\n");
	if (leaf_data == 0L)
		return;
	if (leaf_data->state_set_lookup) {
		Py_DECREF(leaf_data->state_set_lookup);
	}
	else if (leaf_data->state_lookup) {
		if (leaf_data->state_lookup[0])
			free(leaf_data->state_lookup[0]);
		free(leaf_data->state_lookup);
	}
	if (leaf_data->ssind)
		free(leaf_data->ssind);
	freeDbl3DMatrix(leaf_data->pmat_columns);
	PyObject_Del(leaf_data);
}

void cla_dtor(CLAObj * cla) {
	CPHYPROB_DEBUG_PRINTF("In cla_dtor\n");
	if (cla == 0L)
		return;
#	if SEPARATE_RESCALER_ARRAY
		if (cla->rescale_swap_space) {
			free(cla->rescale_swap_space);
		}
#	endif
	if (cla->cla) {
		free(cla->cla);
	}
	cso_dtor(&(cla->cso));

	PyObject_Del(cla);
}

void full_la_dtor(FullLAObj * full_la) {
	CPHYPROB_DEBUG_PRINTF("In full_la_dtor\n");
	if (full_la == 0L)
		return;
	if (full_la->full_cla) {
		cla_dtor(full_la->full_cla);
	}
	else if (full_la->pat_lnL_wts) {
		free(full_la->pat_lnL_wts);
	}
	if (full_la->scratch)
		free(full_la->scratch);
	if (full_la->state_categ_freqs) {
		freeDblMatrix(full_la->state_categ_freqs);
	}
	cso_dtor(&(full_la->cso));
	PyObject_Del(full_la);
}






int configure_context(int ibeg_subset, int iend_subset, unsigned n_subsets, int ibeg_cat, int iend_cat, int irescale_thresh, unsigned n_cat, CalcContext * context) {
	assert(context);
	if (ibeg_subset < 0)
		ibeg_subset += n_subsets;
	if (ibeg_subset < 0) {
		PyErr_SetString(PyExc_ValueError, "The beginning subset index must be 0 or greater");
		return 0L;
	}
	if (iend_subset < 1)
		iend_subset += n_subsets;
	if (iend_subset < ibeg_subset) {
		PyErr_SetString(PyExc_ValueError, "The end subset must be greater that the beginning subset index");
		return 0L;
	}
	if (iend_subset > n_subsets) {
		PyErr_SetString(PyExc_ValueError, "The end subsets cannot be greater that the number of subsets in the cla");
		return 0L;
	}
	if (ibeg_cat < 0)
		ibeg_cat += n_cat;
	if (ibeg_cat < 0) {
		PyErr_SetString(PyExc_ValueError, "The beginning category index must be 0 or greater");
		return 0L;
	}
	if (iend_cat < 1)
		iend_cat += n_cat;
	if (iend_cat < ibeg_cat) {
		PyErr_SetString(PyExc_ValueError, "The end category must be greater that the beginning category index");
		return 0L;
	}
	if (iend_cat > n_cat) {
		PyErr_SetString(PyExc_ValueError, "The end category cannot be greater that the number of categories in the cla");
		return 0L;
	}
	if (irescale_thresh < 1)
		irescale_thresh = 1;
	context->subset_beg = (unsigned) ibeg_subset;
	context->subset_end = (unsigned) iend_subset;
	context->categ_beg = (unsigned) ibeg_cat;
	context->categ_end = (unsigned) iend_cat;
	context->rescale_threshold = (unsigned) irescale_thresh;
	return 1;
}


#if !defined(STANDALONE_C_PROG)
#else  /* STANDALONE_C_PROG */
	int main(int argc, char *argv[]) {
		int i,useMean;
		double alpha, beta;
		unsigned ncat;
		double * f;
		double * r;
		double total = 0.0;
		if (argc < 4) {
			printf("Expecting 3 arguments: alpha n_cat use_mean\nwhere < use_mean> is 0 to use the median, and any other integer to use the mean.\n");
			return 1;
		}
		alpha = atof(argv[1]);
		beta = 1/alpha;
		ncat = atoi(argv[2]);
		useMean = atoi(argv[3]);
		f = malloc(ncat*sizeof(double));
		r = malloc(ncat*sizeof(double));
		if (f != 0L && r != 0L) {
			DiscreteGamma(f, r, alpha, beta, ncat, useMean);
			for (i = 0; i < ncat; ++i) {
				total += r[i];
				printf("%d %f\n", i, r[i]);
			}
			printf("total %f\n",total);
			free(f);
			free(r);
			return 0;
		}
		return 0;
	}
#endif /*STANDALONE_C_PROG*/








#endif // defined(ALL_IN_ONE)

#endif // 0
