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
/**
 * Pure C code
 */
#if defined(ALL_IN_ONE) || !defined(BUILDING_FOR_PYTHON)


#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#if defined(STANDALONE_C_PROG) || (defined(PRINTING_LOTS) && PRINTING_LOTS)
#	include <stdio.h>
#endif

#include "dsct_model.h"


int new_cso_ss_partitioned(CategSubsetObj * cso, unsigned n_sites, unsigned n_states);

void cso_dtor(CategSubsetObj * cso);

#if !defined(BUILDING_FOR_PYTHON)
	void PyErr_NoMemory() {}
	void PyErr_SetString(int v, const char *c)
	{
		PRINTF2("Error: code %d\nmessage: %s\n", v, c);
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


/**
 * (internal) allocates a 2-D array where the memory for the doubles is a
 * contiguous block.
 *
 *	\assert (n_rows > 0) and (n_cols > 0)
 */
double **allocateDblMatrix(unsigned n_rows, unsigned n_cols) {
	double **mat_p;
	double *arr;
	const int arr_len = n_rows*n_cols*(sizeof(double));
	unsigned i, j;
	PRINTF1("In allocateDblMatrix len = %d \n", n_rows*n_cols);
	assert(n_rows > 0);
	assert(n_cols > 0);
	mat_p = (double**)malloc(n_rows*sizeof(double*));
	if (mat_p == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	arr = (double*)malloc(arr_len);
	if (arr == 0L) {
		free(mat_p);
		PyErr_NoMemory();
		return 0L;
	}
	for (i = 0; i < n_rows; ++i) {
		mat_p[i] = arr;
		for (j = 0; j < n_cols; ++j)
			*arr++ = 0.0;
	}
	return mat_p;
}
/**
 * (internal) allocates a 2-D array where the memory for the doubles is a
 * contiguous block.
 *
 *	\assert (n_mat > 0) and (n_rows > 0) and (n_cols > 0)
 */
double ***allocateDbl3DMatrix(unsigned nm, unsigned nr, unsigned nc) {
	double ***mat_list_p;
	double **mat_walk;
	double *arr;
	const int arr_len = nm*nr*nc*(sizeof(double));
	unsigned i, j, k;
	PRINTF1("In allocateDbl3DMatrix len =%d\n", arr_len);
	assert(nm > 0);
	assert(nr > 0);
	assert(nc > 0);
	mat_list_p = (double***)malloc(nm*sizeof(double**));
	if (mat_list_p == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	mat_walk = (double**)malloc(nr*nm*sizeof(double*));
	if (mat_walk == 0L) {
		free(mat_list_p);
		PyErr_NoMemory();
		return 0L;
	}
	arr = (double*)malloc(arr_len);
	if (arr == 0L) {
		free(mat_list_p);
		free(mat_walk);
		PyErr_NoMemory();
		return 0L;
	}
	for (i = 0; i < nm; ++i) {
		mat_list_p[i] = mat_walk;
		for (j = 0; j < nr; ++j) {
			*mat_walk++ = arr;
			for (k = 0; k < nc; ++k) {
				*arr++ = 0.0;
			}
		}
	}
	return mat_list_p;
}
/**
 * (internal) frees a 2-D array allocated with allocateDblMatrix
 *
 *	\assert (n_rows > 0) and (n_cols > 0)
 */
void freeDblMatrix(double **p) {
	PRINTF("In freeDblMatrix\n");
	if (p == 0L)
		return;
	if (p[0] != 0L)
		free(p[0]);
	free(p);
}
void freeDbl3DMatrix(double ***p) {
	PRINTF("In freeDbl3DMatrix\n");
	if (p == 0L)
		return;
	if (p[0] != 0L) {
		if (p[0][0] != 0L) {
			free(p[0][0]);
		}
		free(p[0]);
	}
	free(p);
}
/* Phylogenetic and Numerical code */
/* from Ziheng Yang's PAML */
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
		   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion		if (alpha < x || x <= 1)
   (2) continued fraction	otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha) {
	int i;
	double p = alpha, g = ln_gamma_alpha;
	const double ACCURACY_TOL = 1e-10;
	const double OVERFLOW_VAL = 1e60;
	double factor, gin = 0, rn = 0, a = 0,b = 0,an = 0,dif = 0, term = 0, pn[6];
	if (x == 0)
		return 0;
	if (x < 0 || p <= 0)
		return -1;
	factor = exp(p*log(x) - x - g);
	if (x > 1 && x >= p) {
		/* continued fraction */
		a = 1 - p;
		b = a + x + 1;
		term = 0;
		pn[0] = 1;
		pn[1] = x;
		pn[2] = x + 1;
		pn[3] = x*b;
		gin = pn[2]/pn[3];
		for (;;){
			a++;
			b += 2;
			term++;
			an = a*term;
			for (i = 0; i < 2; i++)
				pn[i+4] = b*pn[i+2] - an*pn[i];
			if (pn[5] != 0) {
				rn = pn[4]/pn[5];
				dif = fabs(gin-rn);
				if ((dif <= ACCURACY_TOL) && (dif <= ACCURACY_TOL*rn))
					return 1 - factor*gin;
				gin = rn;
			}
			for (i = 0; i < 4; i++)
				pn[i]= pn[i+2];
			if (fabs(pn[4]) >= OVERFLOW_VAL) {
				for (i = 0; i < 4; i++)
					pn[i] /= OVERFLOW_VAL;
			}
		}
	}
	else {
		/* series expansion */
		gin = 1;
		term = 1;
		rn = p;
		while (term > ACCURACY_TOL) {
			rn++;
			term *= x/rn;
			gin += term;
		}
		gin *= factor/p;
	}
	return (gin);
}
long factorial(int n) {
	long f, i;
	assert(n <= 10);
	for (i = 2, f = 1; i <= (long)n; i++)
		f*= i;
	return f;
}
double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
	double f = 0, fneg = 0, z, tmp1, tmp2;
	int nx =(int)x-1;
	if((double)nx == x && nx > 0 && nx < 10)
		return log((double)factorial(nx));
	assert(x>0);
	if (x < 7) {
		f = 1;
		z = x-1;
		while (++z < 7)
			f *= z;
		x = z;
		f = -log(f);
	}
	z = 1/(x*x);
	tmp1 = fneg+ f + (x-0.5)*log(x) - x + .918938533204673;
	tmp2 = (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z +.083333333333333);
	return tmp1 + (tmp2/x);
}
/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/
double PointNormal (double prob) {
/* returns z so that Prob{x < z}= prob where x ~ N(0,1) and (1e-12)< prob < 1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)
   Newer methods:
	 Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	   normal distribution.	 37: 477-484.
	 Beasley JD & Springer SG  (1977).	Algorithm AS 111: the percentage
	   points of the normal distribution.  26: 118-121.
*/
	double a0 =-.322232431088, a1 =-1, a2 =-.342242088547, a3 =-.0204231210245;
	double a4 =-.453642210148e-4, b0 =.0993484626060, b1 =.588581570495;
	double b2 =.531103462366, b3 =.103537752850, b4 =.0038560700634;
	double y, z = 0, p = prob, p1;
	p1 = (p < 0.5 ? p : 1-p);
	if (p1 < 1e-20)
		z = 999;
	else {
		y = sqrt (log(1/(p1*p1)));
		z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	}
	return (p < 0.5 ? -z : z);
}
double PointChi2 (double prob, double v) {
	/* returns z so that Prob{x < z}= prob where x is Chi2 distributed with df = v
	returns -1 if in error.	  0.000002 < prob < 0.999998
	RATNEST FORTRAN by
	Best DJ & Roberts DE (1975) The percentage points of the
	Chi2 distribution.	Applied Statistics 24: 385-388.	 (AS91)
	Converted into C by Ziheng Yang, Oct. 1993.
	*/
	double e =.5e-6, aa =.6931471805, p = prob, g, TINY_PROB = 1e-6;
	double xx, c, ch, a = 0,q = 0,p1 = 0,p2 = 0,t = 0,x = 0,b = 0,s1,s2,s3,s4,s5,s6;
	int doL3;
	assert(v > 0);
	if (p < TINY_PROB)
		return 0.0;
	if (p > 1.0 - TINY_PROB)
		return 9999.0;
	xx = v/2;
	g = LnGamma(xx);
	c = xx-1;
	if (v < -1.24*log(p)){
		ch = pow((p*xx*exp(g+xx*aa)), 1/xx);
		if (ch-e < 0)
			return (ch);
	}
	else {
		doL3 = 1;
		if (v <= .32) {
			ch = 0.4;
			a = log(1-p);
			for (;;) {
				q = ch;
				p1 = 1+ch*(4.67+ch);
				p2 = ch*(6.73+ch*(6.66+ch));
				t = -0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
				ch -= (1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
				if (fabs(q/ch-1)-.01 <= 0) {
					doL3 = 0;
					break;
				}
			}
		}
		if (doL3) {
			x = PointNormal (p);
			p1 = 0.222222/v;
			ch = v*pow((x*sqrt(p1)+1-p1), 3.0);
			if (ch>2.2*v+6)
				ch = -2*(log(1-p)-c*log(.5*ch)+g);
		}
	}
	do {
		q = ch;
		p1 = .5*ch;
		t = IncompleteGamma (p1, xx, g);
		assert(t >= 0);
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		s1 = (210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
		s2 = (420+a*(735+a*(966+a*(1141+1278*a))))/2520;
		s3 = (210+a*(462+a*(707+932*a)))/2520;
		s4 = (252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
		s5 = (84+264*a+c*(175+606*a))/2520;
		s6 = (120+c*(346+127*c))/5040;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	} while(fabs(q/ch-1) > e);
	return ch;
}
/*
	Fills the first `n_cat` elements of `rates` with the mean rate of the
	corresponding quantile of the gamma distrib with rates alpha and gamma.
	The first element corresponds to the quantile demarcated by inverse CDF of
		0 up to 1/n_cat.
	The second goes from 1/n_cat up to 2/n_cat
	etc.
   Ziheng's comment on	DiscreteMean:
   discretization of gamma distribution with equal proportions in each
   category.
   MTH comment:
   I reworked this code from Ziheng's DiscreteMean to remove the need for a freqK
   array (which always ended up as n elements with values 1/n, but was used
   as scratch space internally).
*/
void DiscreteGammaMean (double rates[], const double alpha, const double beta, const unsigned n_cat)
{
	double lnga1;
	unsigned last_cat, i;
	double current;
	const double inv_two_beta = .5/beta;
	const double d_n_cat = (double)n_cat;
	const double cat_freq = 1.0/d_n_cat;
	double prop = cat_freq;
	double rate_sum = 0.0;
	double prev = 0.0;
	const double factor = (alpha*n_cat)/beta;
	assert(n_cat > 0);
	assert(rates != 0L);
	assert(alpha > 0.0);
	assert(beta > 0.0);
	if (n_cat < 2) {
		rates[0] = 1.0;
		return;
	}
	last_cat = n_cat - 1;
	lnga1 = LnGamma(alpha + 1.0);
	for (i = 0; i < last_cat; i++) {
		/*cutting points, Eq. 9 -- Ziheng's comment */
		current = inv_two_beta*PointChi2(prop, 2.0*alpha);
		current = IncompleteGamma(current*beta, alpha + 1, lnga1);
		rates[i] = (current - prev)*factor;
		rate_sum += rates[i];
		prev = current;
		prop +=	 cat_freq;
	}
	rates[last_cat] = (1.0 - prev)*factor;
	rate_sum += rates[last_cat];
	rate_sum /= d_n_cat;
	for (i = 0; i < n_cat; ++i)
		rates[i] /= rate_sum;
 }
/*
	Fills the first `n_cat` elements of `rates` with the mean rate of the
	corresponding quantile of the gamma distrib with rates alpha and gamma.
	The first element corresponds to the quantile demarcated by inverse CDF of
		0 up to 1/n_cat.
	The second goes from 1/n_cat up to 2/n_cat
	etc.
	Ziheng's comment on	DiscreteMean:
	discretization of gamma distribution with equal proportions in each
	category.
	MTH comment:
	I reworked this code from Ziheng's DiscreteMean to remove the need for a freqK
	array (which always ended up as n elements with values 1/n, but was used
	as scratch space internally).
*/
void DiscreteGammaMedian (double rates[], const double alpha, const double beta, unsigned n_cat)
{
	unsigned i;
	double total = 0.0;
	const double cat_freq = 1.0/((double)n_cat);
	double prop = cat_freq*0.5;
	double lnga1;
	assert(n_cat > 0);
	assert(rates != 0L);
	assert(alpha > 0.0);
	assert(beta > 0.0);
	if (n_cat < 2) {
		rates[0] = 1.0;
		return;
	}
	lnga1 = LnGamma(alpha + 1.0);
	for(i = 0; i < n_cat; ++i) {
		rates[i] = PointChi2(prop, 2.0*alpha);
		total += rates[i];
		prop +=	 cat_freq;
	}
	total /= (double)n_cat;
	for(i = 0; i < n_cat; ++i)
		rates[i] /= total;
}

void DiscreteGamma(double *f,double *r, double alpha, double beta, int ncat, int useMean) {
	const double freq = 1.0/((double)ncat);
	int i;
	for (i = 0; i < ncat; ++i)
		f[i] = freq;
	if (useMean)
		DiscreteGammaMean(r, alpha, beta, ncat);
	else
		DiscreteGammaMedian(r, alpha, beta, ncat);
}
/* end from Ziheng Yang's PAML */
/* end phylogenetic and numerical code */




/*Internal function prototypes*/
	/*phylogenetic functions*/
	/*From MrBayes*/
static int recalc_eigen_mat(DSCTModelObj *mod);
static int get_eigens(unsigned dim,double **q, double *eigen_values, double *im_eigen_values, double **eigen_vectors, double **inv_eigen_vectors, double **work_mat, double *d_work, int *i_work, unsigned *is_complex);
static void calc_c_ijk(unsigned dim, double *c_ijk, const double **eigen_vectors, const double **inv_eigen_vectors);
/*
	static int prob_mat_and_derivs_from_eigensystem (const unsigned dim, double **pMat, double **pMatDeriv, double **pMatSecDeriv, double *EigValExp, const double *cijk, const double *eigenVals, const double v, const double r);
*/
static int prob_mat_from_eigensystem (const unsigned dim, double **pMat, double *EigValExp, const double *cijk, const double *eigenVals, const double branch_len);
INLINE int prob_mat_from_clean_model(double **p, DSCTModelObj *mod, double branch_len);

INLINE void transpose_pmat_columns(double * transposed, const double ***p_mat_array, const unsigned to_state, const unsigned begin_categ, const unsigned end_categ, const unsigned n_states);
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

static const double MB_ETA_TOL = 1E-30;
static const double MB_TINY_TOL = 1.0e-20;
/*---------------------------------------------------------------------------------
|   CopyDoubleMatrices
|   Copies the contents of one matrix of doubles to another matrix.
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void CopyDoubleMatrices (int dim, const double **from, double **to) {
	int			i, j;
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			to[i][j] = from[i][j];
		}
	}
}
/*---------------------------------------------------------------------------------
|   D_sign
|   This function is called from "Hqr2".
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
double D_sign (double a, double b) {
	double		x;
	x = (a >= 0 ? a : -a);
	return (b >= 0 ? x : -x);
}
/*---------------------------------------------------------------------------------
|
|   ComplexDivision2
|
|   Returns the complex quotient of two complex numbers. It does not require that
|   the numbers be in a complex structure.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ComplexDivision2 (double ar, double ai, double br, double bi, double *cr, double *ci)
{
	double		s, ais, bis, ars, brs;
	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs*brs + bis*bis;
	*cr = (ars*brs + ais*bis) / s;
	*ci = (ais*brs - ars*bis) / s;
}
/*---------------------------------------------------------------------------------
|
|   Exchange
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void Exchange (int j, int k, int l, int m, int n, double **a, double *scale)
{
	int			i;
	double		f;
	scale[m] = (double)j;
	if (j != m) {
		for (i = 0; i <= l; i++) {
			f = a[i][j];
			a[i][j] = a[i][m];
			a[i][m] = f;
		}
		for (i = k; i < n; i++) {
			f = a[j][i];
			a[j][i] = a[m][i];
			a[m][i] = f;
		}
	}
}
/*---------------------------------------------------------------------------------
|
|   ElmHes
|
|   Given a real general matrix, this subroutine
|   reduces a submatrix situated in rows and columns
|   low through high to upper Hessenberg form by
|   stabilized elementary similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc.  if  balanc  has not been used,
|      set low = 1, high = dim.
|
|    * a contains the input matrix.
|
|   On output:
|
|    * a contains the hessenberg matrix.  The multipliers
|      which were used in the reduction are stored in the
|      remaining triangle under the hessenberg matrix.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction.
|
|   Only elements low through high are used.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ElmHes (int dim, int low, int high, double **a, int *interchanged) {
	int			i, j, m, la, mm1, kp1, mp1;
	double		x, y;
	la = high - 1;
	kp1 = low + 1;
	if (la < kp1)
		return;
	for (m = kp1; m <= la; m++) {
		mm1 = m - 1;
		x = 0.0;
		i = m;
		for (j = m; j <= high; j++) {
			if (fabs(a[j][mm1]) > fabs(x)) {
				x = a[j][mm1];
				i = j;
			}
		}
		interchanged[m] = i;
		if (i != m)  {
			/* interchange rows and columns of a */
			for (j = mm1; j < dim; j++) {
				y = a[i][j];
				a[i][j] = a[m][j];
				a[m][j] = y;
			}
			for (j = 0; j <= high; j++) {
				y = a[j][i];
				a[j][i] = a[j][m];
				a[j][m] = y;
			}
		}
		if (fabs(x) > MB_ETA_TOL) {
			mp1 = m + 1;
			for (i = mp1; i <= high; i++) {
				y = a[i][mm1];
				if (fabs(y) > MB_ETA_TOL) {
					y /= x;
					a[i][mm1] = y;
					for (j = m; j < dim; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j <= high; j++)
						a[j][m] += y * a[j][i];
				}
			}
		}
	}
}
/*---------------------------------------------------------------------------------
|
|   Balanc
|
|   This subroutine balances a real matrix and isolates
|   eigenvalues whenever possible.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * a contains the input matrix to be balanced
|
|   On output:
|
|    * a contains the balanced matrix.
|
|    * low and high are two integers such that a(i,j)
|      is equal to zero if
|         (1) i is greater than j and
|         (2) j = 1,...,low-1 or i = igh+1,...,n.
|
|    * scale contains information determining the
|      permutations and scaling factors used.
|
|   Suppose that the principal submatrix in rows pLow through pHigh
|   has been balanced, that p(j) denotes the index interchanged
|   with j during the permutation step, and that the elements
|   of the diagonal matrix used are denoted by d(i,j). Then
|      scale(j) = p(j),    for j = 1,...,pLow-1
|               = d(j,j),      j = pLow,...,pHigh
|               = p(j)         j = pHigh+1,...,dim.
|   The order in which the interchanges are made is dim to pHigh+1,
|   then 1 to pLow-1.
|
|   Note that 1 is returned for pHigh if pHigh is zero formally.
|
|   The algol procedure exc contained in balance appears in
|   balanc in line.  (Note that the algol roles of identifiers
|   k,l have been reversed.)
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   This function was converted from FORTRAN by D. L. Swofford.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void Balanc (int dim, double **a, int *low, int *high, double *scale) {
	int			i, j, k, l, m, noconv;
	double		c, f, g, r, s, b2;
	b2 = FLT_RADIX * FLT_RADIX;
	k = 0;
	l = dim - 1;
	for (j = l; j>= 0; j--) {
		for (i = 0; i <= l; i++) {
			if ((i != j) && (fabs(a[j][i]) > MB_ETA_TOL)) {
				goto next_j1;
			}
		}
		/* bug that DLS caught */
		/*m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (l < 0)
			goto leave;
		else
			j = --l;*/
		m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (--l < 0)
			goto leave;
		next_j1:
			;
	}
	for (j = k; j <= l; j++) {
		for (i = k; i <= l; i++) {
			if ((i != j) && (fabs(a[i][j]) > MB_ETA_TOL)) {
				goto next_j;
			}
		}
		m = k;
		Exchange(j, k, l, m, dim, a, scale);
		k++;
		next_j:
			;
	}
	for (i = k; i <= l; i++)
		scale[i] = 1.0;
	noconv = 1;
	while (noconv) {
		noconv = 0;
		for (i = k; i <= l; i++) {
			c = 0.0;
			r = 0.0;
			for (j = k; j <= l; j++) {
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			}
			if (fabs(c) > MB_ETA_TOL && fabs(r) > MB_ETA_TOL) {
				g = r / FLT_RADIX;
				f = 1.0;
				s = c + r;
				while (c < g) {
					f *= FLT_RADIX;
					c *= b2;
				}
				g = r * FLT_RADIX;
				while (c >= g) {
					f /= FLT_RADIX;
					c /= b2;
				}
				if ((c + r) / f < s * .95) {
					g = 1.0 / f;
					scale[i] *= f;
					noconv = 1;
					for (j = k; j < dim; j++)
						a[i][j] *= g;
					for (j = 0; j <= l; j++)
						a[j][i] *= f;
				}
			}
		}
	}
	leave:
		*low = k;
		*high = l;
}
/*---------------------------------------------------------------------------------
|
|   ElTran
|
|   This subroutine accumulates the stabilized elementary
|   similarity transformations used in the reduction of a
|   real general matrix to upper Hessenberg form by ElmHes.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If Balanc has not been used,
|      set low = 0, high = dim-1.
|
|    * a contains the multipliers which were used in the
|      reduction by  ElmHes in its lower triangle
|      below the subdiagonal.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction by ElmHes.
|      only elements low through high are used.
|
|   On output:
|
|    * z contains the transformation matrix produced in the
|      reduction by ElmHes.
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ElTran (int dim, int low, int high, double **a, int *interchanged, double **z)
{
	int			i, j, mp;
	/* initialize z to identity matrix */
	for (j = 0; j < dim; j++)
		{
		for (i = 0; i < dim; i++)
			z[i][j] = 0.0;
		z[j][j] = 1.0;
		}
	for (mp = high-1; mp>= low+1; mp--) /* there were a number of additional    */
		{                            /* variables (kl, la, m, mm, mp1) that  */
		for (i = mp+1; i <= high; i++)   /* have been eliminated here simply by  */
			z[i][mp] = a[i][mp-1];   /* initializing variables appropriately */
		i = interchanged[mp];        /* in the loops                         */
		if (i != mp) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (j = mp; j <= high; j++)
				{
				z[mp][j] = z[i][j];
				z[i][j] = 0.0;
				}
			z[i][mp] = 1.0;
			}
		}
}
/*---------------------------------------------------------------------------------
|
|   Hqr2
|
|   This subroutine finds the eigenvalues and eigenvectors
|   of a real upper Hessenberg matrix by the QR method. The
|   eigenvectors of a real general matrix can also be found
|   if ElmHes  and ElTran or OrtHes and OrTran have
|   been used to reduce this general matrix to Hessenberg form
|   and to accumulate the similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If  balanc has not been used,
|      set low = 0, high = dim-1.
|
|    * h contains the upper hessenberg matrix. Information about
|      the transformations used in the reduction to Hessenberg
|      form by  ElmHes  or OrtHes, if performed, is stored
|      in the remaining triangle under the Hessenberg matrix.
|
|   On output:
|
|    * h has been destroyed.
|
|    * wr and wi contain the real and imaginary parts,
|      respectively, of the eigenvalues. The eigenvalues
|      are unordered except that complex conjugate pairs
|      of values appear consecutively with the eigenvalue
|      having the positive imaginary part first. If an
|      error exit is made, the eigenvalues should be correct
|      for indices j,...,dim-1.
|
|    * z contains the transformation matrix produced by ElTran
|      after the reduction by ElmHes, or by OrTran after the
|      reduction by OrtHes, if performed. If the eigenvectors
|      of the Hessenberg matrix are desired, z must contain the
|      identity matrix.
|
|   Calls ComplexDivision2 for complex division.
|
|   This function returns:
|      zero       for normal return,
|      j          if the limit of 30*n iterations is exhausted
|                 while the j-th eigenvalue is being sought.
|
|   This subroutine is a translation of the ALGOL procedure HQR2,
|   Num. Math. 14, 219,231(1970) by Martin, Peters, and Wilkinson.
|   Handbook for Automatic Computation, vol. II - Linear Algebra,
|   pp. 357-391 (1971).
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int Hqr2 (int dim, int low, int high, double **h, double *wr, double *wi, double **z)
{
	int			i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;
	double		norm, p = 0.0, q = 0.0, r = 0.0, s = 0.0, t, w = 0.0, x, y = 0.0, ra, sa, vi, vr, zz = 0.0, tst1, tst2;
	norm = 0.0;
	k = 0;  /* used for array indexing. FORTRAN version: k = 1 */
	/* store roots isolated by balance, and compute matrix norm */
	for (i = 0; i < dim; i++)
		{
		for (j = k; j < dim; j++)
			norm += fabs(h[i][j]);
		k = i;
		if ((i < low) || (i > high))
			{
			wr[i] = h[i][i];
			wi[i] = 0.0;
			}
		}
	en = high;
	t = 0.0;
	itn = dim * 30;
	/* search for next eigenvalues */
	while (en >= low) /* changed from an "if(en < lo)" to eliminate a goto statement */
		{
		its = 0;
		na = en - 1;
		enm2 = na - 1;
		twoRoots = 0;
		for (;;)
			{
			for (l = en; l>low; l--) /* changed indexing, got rid of lo, ll */
				{
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
				if (fabs(s) <= MB_ETA_TOL) /* == 0.0 */
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l][l-1]);
				if (fabs(tst2 - tst1) < MB_ETA_TOL) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
			/* form shift */
			x = h[en][en];
			if (l == en) /* changed to break to remove a goto statement */
				break;
			y = h[na][na];
			w = h[en][na] * h[na][en];
			if (l == na)         /* used to return to other parts of the code */
				{
				twoRoots = 1;
				break;
				}
			if (itn == 0)
				return (en);
			/* form exceptional shift */
			if ((its == 10) || (its == 20)) /* changed to remove a goto statement */
				{
				t += x;
				for (i = low; i <= en; i++)
					h[i][i] -= x;
				s = fabs(h[en][na]) + fabs(h[na][enm2]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
				}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for (m = enm2; m>= l; m--)
				{
				/* removed m = enm2 + l - mm and above loop to remove variables */
				zz = h[m][m];
				r = x - zz;
				s = y - zz;
				p = (r * s - w) / h[m+1][m] + h[m][m+1];
				q = h[m+1][m+1] - zz - r - s;
				r = h[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break; /* changed to break to remove a goto statement */
				tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
				tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
				if (fabs(tst2 - tst1) < MB_ETA_TOL) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
			mp2 = m + 2;
			for (i = mp2; i <= en; i++)
				{
				h[i][i-2] = 0.0;
				if (i != mp2) /* changed "==" to "!=" to remove a goto statement */
					h[i][i-3] = 0.0;
				}
			/* double QR step involving rows l to en and columns m to en */
			for (k = m; k <= na; k++)
				{
				notlas = (k != na);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					{
					p = h[k][k-1];
					q = h[k+1][k-1];
					r = 0.0;
					if (notlas)
						r = h[k+2][k-1];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x < MB_ETA_TOL) /* == 0.0 */
						continue; /* changed to continue remove a goto statement */
					p /= x;
					q /= x;
					r /= x;
					}
		        /*
		        s = sqrt(p*p+q*q+r*r);
		        sgn = (p < 0)?-1:(p>0);
		        s = sgn*sqrt(p*p+q*q+r*r);
		        */
				s = D_sign(sqrt(p*p + q*q + r*r), p);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					h[k][k-1] = -s * x;
				else if (l != m) /* else if gets rid of another goto statement */
					h[k][k-1] = -h[k][k-1];
				p += s;
				x = p / s;
				y = q / s;
				zz = r / s;
				q /= p;
				r /= p;
				if (!notlas) /* changed to !notlas to remove goto statement (see **) */
					{
					/* row modification */
					for (j = k; j < dim; j++)
						{
						p = h[k][j] + q * h[k+1][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						}
				    j = k + 3;
				    if (en < j)
				    	j = en;
					/* column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						}
					/* accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						}
					}
				else /* (**) also put in else */
					{
					/* row modification */
					for (j = k; j < dim; j++)
						{
						p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						h[k+2][j] -= p * zz;
						}
					j = k + 3;
					if (en < j)
						j = en;
					/* column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						h[i][k+2] -= p * r;
						}
					/* accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						z[i][k+2] -= p * r;
						}
					}
				}
			}
		if (twoRoots)
			{
			/* two roots found */
			p = (y - x) / 2.0;
			q = p * p + w;
			zz = sqrt(fabs(q));
			h[en][en] = x + t;
			x = h[en][en];
			h[na][na] = y + t;
			if (q >= -1e-12) /* change "<" to ">=", and also change "0.0" to */
				{            /* a small number (Swofford's change)           */
				/* real pair */
				zz = p + D_sign(zz, p);
				wr[na] = x + zz;
				wr[en] = wr[na];
				if (fabs(zz) > MB_ETA_TOL) /* != 0.0 */
					wr[en] = x - w/zz;
				wi[na] = 0.0;
				wi[en] = 0.0;
				x = h[en][na];
				s = fabs(x) + fabs(zz);
				p = x / s;
				q = zz / s;
				r = sqrt(p*p + q*q);
				p /= r;
				q /= r;
				/* row modification */
				for (j = na; j < dim; j++)
					{
					zz = h[na][j];
					h[na][j] = q * zz + p * h[en][j];
					h[en][j] = q * h[en][j] - p * zz;
					}
				/* column modification */
				for (i = 0; i <= en; i++)
					{
					zz = h[i][na];
					h[i][na] = q * zz + p * h[i][en];
					h[i][en] = q * h[i][en] - p * zz;
					}
				/* accumulate transformations */
				for (i = low; i <= high; i++)
					{
					zz = z[i][na];
					z[i][na] = q * zz + p * z[i][en];
					z[i][en] = q * z[i][en] - p * zz;
					}
				}
			else
				{
				/* complex pair */
				wr[na] = x + p;
				wr[en] = x + p;
				wi[na] = zz;
				wi[en] = -zz;
				}
			en = enm2;
			}
		else
			{
			/* one root found */
			h[en][en] = x + t;
			wr[en] = h[en][en];
			wi[en] = 0.0;
			en = na;
			}
		}
	if (fabs(norm) < MB_ETA_TOL) /* == 0.0 */
		return (0); /* was a goto end of function */
	for (en = dim-1; en>= 0; en--)
		{
		/*en = n - nn - 1; and change for loop */
		p = wr[en];
		q = wi[en];
		na = en - 1;
		if (q < -1e-12)
			{
			/* last vector component chosen imaginary so that eigenvector
			   matrix is triangular */
			m = na;
			if (fabs(h[en][na]) > fabs(h[na][en]))
				{
				h[na][na] = q / h[en][na];
				h[na][en] = -(h[en][en] - p) / h[en][na];
				}
			else
				ComplexDivision2 (0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);
			h[en][na] = 0.0;
			h[en][en] = 1.0;
			enm2 = na - 1;
			if (enm2 >= 0) /* changed direction to remove goto statement */
				{
				for (i = enm2; i>= 0; i--)
					{
					w = h[i][i] - p;
					ra = 0.0;
					sa = 0.0;
					for (j = m; j <= en; j++)
						{
						ra += h[i][j] * h[j][na];
						sa += h[i][j] * h[j][en];
						}
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						r = ra;
						s = sa;
						}
					else
						{
						m = i;
						if (fabs(wi[i]) < MB_ETA_TOL) /* == 0.0 */ /* changed direction to remove goto statement */
							ComplexDivision2 (-ra, -sa, w, q, &h[i][na], &h[i][en]);
						else
							{
							/* solve complex equations */
							x = h[i][i+1];
							y = h[i+1][i];
							vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
							vi = (wr[i] - p) * 2.0 * q;
							if ((fabs(vr) < MB_ETA_TOL) && (fabs(vi) < MB_ETA_TOL))
								{
								tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
								vr = tst1;
								do	{
									vr *= .01;
									tst2 = tst1 + vr;
									}
									while (tst2 > tst1); /* made into a do/while loop */
								}
							ComplexDivision2 (x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
							if (fabs(x) > fabs(zz) + fabs(q)) /* changed direction to remove goto statement */
								{
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
								}
							else
								ComplexDivision2 (-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);
							}
						/* overflow control */
						tst1 = fabs(h[i][na]);
						tst2 = fabs(h[i][en]);
						t = (tst1 > tst2 ? tst1 : tst2);
						if (t > MB_ETA_TOL) /* t != 0.0 */
							{
							tst1 = t;
							tst2 = tst1 + 1.0 / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									{
									h[j][na] /= t;
									h[j][en] /= t;
									}
								}
							}
						}
					}
				}
			}
		else if (fabs(q)< MB_ETA_TOL)
			{
			/* real vector */
			m = en;
			h[en][en] = 1.0;
			if (na >= 0)
				{
				for (i = na; i>= 0; i--)
					{
					w = h[i][i] - p;
					r = 0.0;
					for (j = m; j <= en; j++)
						r += h[i][j] * h[j][en];
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						s = r;
						continue;  /* changed to continue to remove goto statement */
						}
					else
						{
						m = i;
						if (fabs(wi[i])< MB_ETA_TOL) /* changed to remove goto statement */
							{
							t = w;
							if (fabs(t)< MB_ETA_TOL)  /* changed to remove goto statement */
								{
								tst1 = norm;
								t = tst1;
								do	{
									t *= .01;
									tst2 = norm + t;
									}
									while (tst2 > tst1);
								}
							h[i][en] = -r / t;
							}
						else
							{
							/* solve real equations */
							x = h[i][i+1];
							y = h[i+1][i];
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
							t = (x * s - zz * r) / q;
							h[i][en] = t;
							if (fabs(x) > fabs(zz))  /* changed direction to remove goto statement */
								h[i+1][en] = (-r - w * t) / x;
							else
								h[i+1][en] = (-s - y * t) / zz;
							}
						/* overflow control */
						t = fabs(h[i][en]);
						if (t > MB_ETA_TOL)
							{
							tst1 = t;
							tst2 = tst1 + 1. / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									h[j][en] /= t;
								}
							}
						}
					}
				}
			}
		}
	for (i = 0; i < dim; i++)
		{
		if ((i < low) || (i > high)) /* changed to rid goto statement */
			{
			for (j = i; j < dim; j++)
				z[i][j] = h[i][j];
			}
		}
	/* multiply by transformation matrix to give vectors of original
	   full matrix */
	for (j = dim-1; j>= low; j--)
		{
		m = (j < high ? j : high);
		for (i = low; i <= high; i++)
			{
			zz = 0.0;
			for (k = low; k <= m; k++)
				zz += z[i][k] * h[k][j];
			z[i][j] = zz;
			}
		}
	return (0);
}
/*---------------------------------------------------------------------------------
|
|   BalBak
|
|   This subroutine forms the eigenvectors of a real general
|   matrix by back transforming those of the corresponding
|   balanced matrix determined by  balance.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by  balance
|
|    * scale contains information determining the permutations
|      and scaling factors used by balance
|
|    * m is the number of columns of z to be back transformed
|
|    * z contains the real and imaginary parts of the eigen-
|      vectors to be back transformed in its first m columns
|
|   On output:
|
|    * z contains the real and imaginary parts of the
|      transformed eigenvectors in its first m columns
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void BalBak (int dim, int low, int high, double *scale, int m, double **z)
{
	int			i, j, k, ii;
	double		s;
	if (m != 0) /* change "==" to "!=" to eliminate a goto statement */
		{
		if (high != low) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (i = low; i <= high; i++)
				{
				s = scale[i];
				for (j = 0; j < m; j++)
					z[i][j] *= s;
				}
			}
		for (ii = 0; ii < dim; ii++)
			{
			i = ii;
			if ((i < low) || (i > high)) /* was (i >= lo) && (i <= hi) but this */
				{                        /* eliminates a goto statement        */
				if (i < low)
					i = low - ii;
				k = (int)scale[i];
				if (k != i) /* change "==" to "!=" to eliminate a goto statement */
					{
					for (j = 0; j < m; j++)
						{
						s = z[i][j];
						z[i][j] = z[k][j];
						z[k][j] = s;
						}
					}
				}
			}
		}
}
/*---------------------------------------------------------------------------------
|
|   eigens_for_real_matrix
|
|   The matrix of interest is a. The ouptut is the real and imaginary parts of the
|   eigenvalues (wr and wi). z contains the real and imaginary parts of the
|   eigenvectors. iv2 and fv1 are working vectors.
|
|	Returns 0 on failure, 1 for success.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int eigens_for_real_matrix (int dim, double **a, double *wr, double *wi, double **z, int *iv1, double *fv1)
{
	static int	is1, is2;
	int			ierr;
	Balanc (dim, a, &is1, &is2, fv1);
	ElmHes (dim, is1, is2, a, iv1);
	ElTran (dim, is1, is2, a, iv1, z);
	ierr = Hqr2 (dim, is1, is2, a, wr, wi, z);
	if (ierr == 0)
		BalBak (dim, is1, is2, fv1, dim, z);
	else {
		PyErr_SetString(PyExc_RuntimeError,"Error in Hqr2: max iterations exceeded.");
		return 0;
	}
	return 1;
}
/*---------------------------------------------------------------------------------
|
|   compute_eigen_system
|
|   Calculates the eigenvalues, eigenvectors, and the inverse of the eigenvectors
|   for a matrix of real numbers.
|
|	Returns 0 for an error or complex eigen values.  If the error is ONLY the
| 		presence of complex eigenvalues then, *is_complex = 1.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int compute_eigen_system (int dim, double **a, double *v, double *vi, double **u, int *iwork, double *dwork, unsigned *is_complex) {
	int i;
	*is_complex = 0;
	if (0 == eigens_for_real_matrix (dim, a, v, vi, u, iwork, dwork))
		return 0;
	for (i = 0; i < dim; i++) {
		if (fabs(vi[i]) > MB_ETA_TOL) { /* != 0.0 */
			*is_complex = 1;
			return 0;
		}
	}
	return 1;
}
/*---------------------------------------------------------------------------------
|
|   LUBackSubstitution
|
|   Back substitute into an LU-decomposed matrix.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void LUBackSubstitution (int dim, double **a, int *i_work, double *b) {
	int			i, ip, j, ii = -1;
	double		sum;
	for (i = 0; i < dim; i++) {
		ip = i_work[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) {
			for (j = ii; j <= i-1; j++)
				sum -= a[i][j] * b[j];
		}
		else if (fabs(sum) > MB_ETA_TOL)
			ii = i;
		b[i] = sum;
	}
	for (i = dim-1; i>= 0; i--) {
		sum = b[i];
		for (j = i+1; j < dim; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}
/*---------------------------------------------------------------------------------
|
|   LUDecompose
|
|   Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int LUDecompose (int dim, double **a, double *d_work, int *i_work) {
	int			i, imax = 0, j, k;
	double		big, dum, sum, temp, d;
	d = 1.0;
	for (i = 0; i < dim; i++) {
		big = 0.0;
		for (j = 0; j < dim; j++) {
			temp = fabs(a[i][j]);
			if (temp > big)
				big = temp;
			}
		if (fabs(big) < MB_ETA_TOL) {
			PyErr_SetString(PyExc_RuntimeError,"Error in LUDecompose: largest abs value of the elements in a row is too small.");
			return 0;
		}
		d_work[i] = 1.0 / big;
	}
	for (j = 0; j < dim; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < dim; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			dum = d_work[i] * fabs(sum);
			if (dum >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < dim; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			d_work[imax] = d_work[j];
		}
		i_work[j] = imax;
		if (fabs(a[j][j]) < MB_ETA_TOL)
			a[j][j] = MB_TINY_TOL;
		if (j != dim - 1) {
			dum = 1.0 / a[j][j];
			for (i = j + 1; i < dim; i++)
				a[i][j] *= dum;
		}
	}
	return 1;
}
/*---------------------------------------------------------------------------------
|
|   InvertMatrix
|
|   Calculates aInv = a^{-1} using LU-decomposition. The input matrix a is
|   destroyed in the process. The program returns an error if the matrix is
|   singular. d_work and i_work are work vectors.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int InvertMatrix (
	int dim,
 	double **a, /*dim by dim - valid on input, but overwritten */
 	double *d_work, /* scratch array, len dim */
 	int *i_work, /* scratch array, len dim */
 	double **aInv/*dim by dim - valid on completion unless 0 is returned */
 	) {
	unsigned	i, j;
	if (LUDecompose (dim, a, d_work, i_work) == 0L)
		return 0;
	for (j = 0; j < dim; j++)
		{
		for (i = 0; i < dim; i++)
			d_work[i] = 0.0;
		d_work[j] = 1.0;
		LUBackSubstitution (dim, a, i_work, d_work);
		for (i = 0; i < dim; i++)
			aInv[i][j] = d_work[i];
		}
	return 1;
}
/*---------------------------------------------------------------------------------
|	get_eigens
|
|	\Returns 0 on error, *is_complex will be 1 if the cause of the failure is
|		complex eigenvalues.
|	Taken from MrBayes GetEigens.
|	MrBayes seems to flag any complex eigensystems with a return code and abort
|		calculations.  Thus, I removed the complex ** arguments and code in the
|		isComplex branch.
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static int get_eigens(
  unsigned dim,
  double **q,
  double *eigen_values,
  double *im_eigen_values,
  double **eigen_vectors,
  double **inv_eigen_vectors,
  double **work_mat, /*scratch  dim by dim matrix*/
  double *d_work, /*scratch of length dim*/
  int *i_work, /*scratch of length dim*/
  unsigned *is_complex) {
	int rc;
	*is_complex = 0;
	memset(d_work, 0, (size_t) (dim * sizeof(double)));
	memset(i_work, 0, (size_t) (dim * sizeof(int)));
	/* calculate eigenvalues and eigenvectors */
	rc = compute_eigen_system (dim, q, eigen_values, im_eigen_values, eigen_vectors, i_work, d_work, is_complex);
	if (rc == 0) {
		if (*is_complex == 1){
			PyErr_SetString(PyExc_ValueError,"Complex Eigenvalues found");
		}
		return 0;
	}
	CopyDoubleMatrices (dim, (const double **) eigen_vectors, work_mat);
	if (InvertMatrix (dim, work_mat, d_work, i_work, inv_eigen_vectors) == 0L)
		return 0;
	return 1;
}
/*---------------------------------------------------------------------------------
|
|   from MrBayes CalcCijk
|
|   This function precalculates the product of the eigenvectors and their
|   inverse for faster calculation of transition probabilities. The output
|   is a vector of precalculated values. The input is the eigenvectors and
|   the inverse of the eigenvector matrix.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void calc_c_ijk(unsigned dim, double *c_ijk, const double **eigen_vectors, const double **inv_eigen_vectors) {
	register int 	i, j, k;
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			for (k = 0; k < dim; k++)
			 	*c_ijk++ = eigen_vectors[i][k] * inv_eigen_vectors[k][j];
}

void printQMat(double ** qmat_obj, const unsigned n_states) {
	unsigned from_state, to_state;
	PRINTF("QMat:\n");
	for (from_state = 0; from_state < n_states; ++from_state) {
		for (to_state = 0; to_state < n_states; ++to_state) {
			PRINTF1("%f ", qmat_obj[from_state][to_state]);
		}
		PRINTF("\n");
	}
}

/**
 * Recalculates the eigensystem and temporaries.  Returns 0 if any steps fail
 *  or the complex eigenvalues are encountered.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
static int recalc_eigen_mat(DSCTModelObj *mod) {
	unsigned is_complex;
	int rc;
	
#	if defined(PRINTING_LOTS) && PRINTING_LOTS
		printQMat(mod->q_mat, mod->dim);
#	endif 	
	
	rc = get_eigens(mod->dim,
					mod->q_mat,
					mod->eigen_values,
					mod->im_eigen_values,
					mod->eigen_vectors,
					mod->inv_eigen_vectors,
					mod->work_mat,
					mod->d_work,
					mod->i_work,
					&is_complex);
	if (rc == 0) {
		if (is_complex == 1)
			PyErr_SetString(PyExc_RuntimeError, "Error: Complex eigenvalues found.");
		return 0;
	}
	calc_c_ijk(mod->dim,
			   mod->cijk,
			   (const double **)mod->eigen_vectors,
			   (const double **)mod->inv_eigen_vectors);
	mod->eigen_calc_dirty = 0;
	return 1;
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
 * Adapted from MrBayes TiProbsUsingEigens in mbmath.c
 *	`pMat` and `EigValExp` are written to.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
static int prob_mat_from_eigensystem (
	const unsigned dim, /*the number of states*/
	double **pMat,			/**/
	double *EigValExp,		/**/
	const double *cijk, /*the dim*dim*dim array of temporaries used in quickly multiplying the exp(eval*branchlength) times the matrix eigenvectors and its inverse*/
	const double *eigenVals, /*array of the eigenvalues*/
	const double branch_len ) /* branch length (rate*time) */
{
	unsigned i, j, s;
	double sum;
	assert(pMat);
	assert(*pMat);
	assert(cijk);
	assert(eigenVals);
	assert(branch_len >= 0.0);
	for (i=0; i<dim; i++)
		EigValExp[i] = exp(eigenVals[i] * branch_len);
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			sum = 0.0;
			for(s=0; s<dim; s++)
				sum += (*cijk++) * EigValExp[s];
			pMat[i][j] = (sum < 0.0) ? 0.0 : sum;
		}
	}
	return 1;
}
INLINE int prob_mat_from_clean_model(double **p, DSCTModelObj *mod, double branch_len) {
	assert(mod);
	assert(p);
	assert(*p);
	return prob_mat_from_eigensystem (mod->dim, p, mod->d_work, mod->cijk, mod->eigen_values, branch_len);
}
/**
 * Calculates a P-Matrix from mod and branch_len.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int do_pmat_calc(double **p, DSCTModelObj *mod, double branch_len) {
	if (mod == 0L) {
		PyErr_SetString(PyExc_ValueError, "Illegal NULL Model pointer");
		return 0;
	}
	if (branch_len < 0.0) {
		PyErr_SetString(PyExc_ValueError, "Illegal (negative) branch length");
		return 0;
	}
	if (p == 0L || *p == 0L) {
		PyErr_SetString(PyExc_ValueError, "Illegal NULL P-Matrix pointer");
		return 0;
	}
	PRINTF("About to check if eigensystem is calculated\n");
	if (mod->eigen_calc_dirty) {
		PRINTF("About to recalc eigensystem\n");
		if (!recalc_eigen_mat(mod))
			return 0;
	}
	PRINTF("About to calc pmats from eigensystem\n");
	if (!prob_mat_from_clean_model(p, mod, branch_len))
		return 0;
	return 1;
}
/**
 * Calculates the p->p_mat[i] for i in [start_ind, end_ind)
 * the branch length and models used come from p->model_aliases, p->brlen_aliases
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int do_pmat_array_calc(PMatArrayObj *p, unsigned start_ind, unsigned end_ind) {
	unsigned i;
	if (p == 0L) {
		PyErr_SetString(PyExc_TypeError, "Expecting a PMatArrayObj object");
		return 0;
	}
	if (p->p_mat == 0L || p->model_aliases == 0L || p->brlen_aliases == 0L) {
		PyErr_SetString(PyExc_ValueError, "Invalid PMatArrayObj object");
		return 0;
	}
	if (p->n_mat < start_ind || p->n_mat < (end_ind -1)) {
		PyErr_SetString(PyExc_IndexError, "End matrix exceeds size of PMatrix Array");
		return 0;
	}
	for (i = start_ind; i < end_ind; ++i) {
		if (!do_pmat_calc(p->p_mat[i], p->model_aliases[i], p->brlen_aliases[i])) {
			return 0;
		}
	}
	p->calc_time_stamp = next_calc_stamp();
	return 1;
}

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
			p->brlen_aliases[i] = asrv->val[i]*edgeLen;
			p->model_aliases[i] = mod;
		}
		return do_pmat_array_calc(p, 0, asrv->n);
	}
	p->brlen_aliases[0] = edgeLen;
	p->model_aliases[0] = mod;
	return do_pmat_array_calc(p, 0, asrv->n);
}

INLINE void transpose_pmat_columns(double * transposed, const double ***p_mat_array, const unsigned to_state, const unsigned begin_categ, const unsigned end_categ, const unsigned n_states) {
	unsigned i, from_state;
	const double **p_mat;
	for (i = begin_categ; i < end_categ; ++i) {
		p_mat = p_mat_array[i];
		for (from_state = 0; from_state < n_states; ++from_state) {
			transposed[from_state + (i*n_states)] = p_mat[from_state][to_state];
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
	dest->n_states = source->n_states;
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
	const double *** pmat_array =  (const double ***) pmat_obj->p_mat;
	const unsigned n_states = pmat_obj->n_states;
	const unsigned subset_ind = (cso.n_categ_arr == 0L ? 0 : cso.n_categ_or_subs - 1);
	const unsigned comp_offset = (subset_ind == 0 ? 0 : cso.subset_component_offset[subset_ind]);
	const unsigned n_categ = (cso.n_categ_arr == 0L ? cso.n_categ_or_subs : cso.n_categ_arr[subset_ind]);
	const unsigned end_categ = comp_offset + n_categ;
	unsigned categ, from_state, to_state;
	PRINTF("PMat:\n");
	for (categ = 0; categ < end_categ; ++categ) {
		for (from_state = 0; from_state < n_states; ++from_state) {
			for (to_state = 0; to_state < n_states; ++to_state) {
				PRINTF1("%f ", pmat_array[categ][from_state][to_state]);
			}
			PRINTF("\n");
		}
		PRINTF("\n");		
	}
}


static int assign_pmat_to_leaf(LeafDataObj * leaf_d, const PMatArrayObj * pmat_obj, const CalcContext context, const CategSubsetObj cso) {
	const int * slook_arr;
	double ** transposed;
	double * copy_source;
	unsigned i, j, c, n_ambig_states, offset, to_state, from_state, subset_ind;
	assert(leaf_d);
	assert(pmat_obj);
 	const double *** pmat_array =  (const double ***) pmat_obj->p_mat;
	const unsigned n_states = pmat_obj->n_states;
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
				transpose_pmat_columns(transposed[to_state], pmat_array + comp_offset, to_state, categ_beg, categ_end, n_states);
			}
			else {
				/* ambiguity -- sum the ambiguous states (each member of the state
				set must have already been calculated */
				assert(to_state < i);
				double * dest = transposed[i];
				copy_source = transposed[to_state];
				for (c = categ_beg; c < categ_end; ++c) {
					offset = c*n_states;
					for (from_state = 0; from_state < n_states; ++from_state) {
						dest[from_state + offset] = copy_source[from_state + offset];
					}
				}
				for (j = 2; j <= n_ambig_states; ++j) {
					to_state = slook_arr[j];
					assert(to_state < i);
					copy_source = transposed[to_state];
					for (c = categ_beg; c < categ_end; ++c) {
						offset = c*n_states;
						for (from_state = 0; from_state < n_states; ++from_state) {
							dest[from_state + offset] += copy_source[from_state + offset];
						}
					}
				}
			}
		}
	}
	leaf_d->pmat_obj_ptr = (void *) pmat_obj;
	leaf_d->calc_time_stamp =  pmat_obj->calc_time_stamp;
	return 1;
}
#if defined(PRINTING_LOTS) && PRINTING_LOTS
	void print_cla(CLAObj * par) {
		const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
		cla_float_t *p_cla = par->cla;
		const unsigned n_states = par->n_states;
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
					for (k = 0; k < n_states; ++k) {
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
	const unsigned n_states = par->n_states;	
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
	
		const unsigned tpm_offset = n_states*(categ_beg);
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
				for (state = 0; state < n_states; ++state) {
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
						p_cla[n_states] = real_rescale_ln;
#					endif
					for (state = 0; state < n_states; ++state)
						p_cla[state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = 0;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[n_states] = 0.0;
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
	const unsigned n_states = par->n_states;
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
	

		const unsigned tpm_offset = n_states*(categ_beg);
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

		PRINTF("do_two_tip_cla_recalc_no_rescale: categ cla likes:\n");
		for (site = beg_site; site < end_site; ++site) {
			PRINTF1("%d", site);
			f_row = f_mat[*f_states++] + tpm_offset;
			s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				for (state = 0; state < n_states; ++state) {
					tmp = (*f_row++) * (*s_row++);
					p_cla[state]  = tmp;
					PRINTF1("\t%f", tmp);
				}
#				if !SEPARATE_RESCALER_ARRAY
					p_cla[n_states] = 0.0;
#				endif
				p_cla += len_per_categ;
			}
			p_cla += post_skip;
			PRINTF("\n");
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
	assert(fd->n_states == sd->n_states);
	assert(fd->n_states == par->n_states);
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
	if ((fd->pmat_obj_ptr != (void *) fmat) || (fd->calc_time_stamp != fmat->calc_time_stamp)) {
		if (!assign_pmat_to_leaf(fd, fmat, context, par->cso))
			return 0;
	}
	if ((sd->pmat_obj_ptr != (void *) smat) || (sd->calc_time_stamp != smat->calc_time_stamp)) {
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
	const unsigned n_states = par->n_states;
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
		const unsigned tpm_offset = n_states*(categ_beg);
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
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
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
						p_cla[n_states] = real_rescale_ln + f_cla[n_states];
#					endif
					for (from_state = 0; from_state < n_states; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[n_states] = f_cla[n_states];
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
	const unsigned n_states = par->n_states;
	const unsigned subset_beg = context.subset_beg;
	const unsigned subset_end = context.subset_end;
	const unsigned * char_offsets = par->cso.subset_char_offsets;
	const unsigned len_per_categ = par->len_per_categ;
	unsigned site, categ, to_state, from_state, subset_ind, beg_site, end_site;
	cla_float_t f_des_prob;
	
	PRINTF("do_one_tip_cla_recalc_no_rescale: categ cla likes:\n");
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
		const unsigned tpm_offset = n_states*(categ_beg);
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
			PRINTF1("%d", site);
			const double *  s_row = s_mat[*s_states++] + tpm_offset;
			for (categ = categ_beg; categ < categ_end; ++categ) {
				f_categ_pmat = f_pmat_arr_o[categ];
				f_categ_p = *f_categ_pmat; /* pointer math in the next 2 loops keeps this pointing to the right elemente*/
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					double tmp = f_des_prob*(*s_row++);
					p_cla[from_state] = tmp;
					PRINTF1("\t%f", tmp);
				}
#				if ! SEPARATE_RESCALER_ARRAY
						p_cla[n_states] = f_cla[n_states];
#				endif
				p_cla += len_per_categ;
				f_cla += len_per_categ;
			}
			PRINTF("\n");
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
	assert(fd->n_states == sd->n_states);
	assert(fd->n_states == par->n_states);
	/*
		assert(fd->n_categ == sd->n_categ);
	*/
	assert(fd->cso.n_categ_or_subs == par->cso.n_categ_or_subs);
	const double *** f_pmat_arr = (const double ***)fmat->p_mat;
	if ((fd->n_edges_since_rescaling_check + 2) < context.rescale_threshold)
		return do_one_tip_cla_recalc_no_rescale(fd, f_pmat_arr, sd, par, context);
	return do_one_tip_cla_recalc_rescale(fd, f_pmat_arr, sd, par, context);
}

int do_one_tip_cla(const CLAObj * fd, const PMatArrayObj *fmat, LeafDataObj* sd, const PMatArrayObj *smat, CLAObj* par, const CalcContext context) {
	assert(fd);
	assert(fmat);
	assert(sd);
	assert(smat);
	if ((sd->pmat_obj_ptr != (void *) smat) || (sd->calc_time_stamp != smat->calc_time_stamp)) {
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
	const unsigned n_states = par->n_states;
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
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					s_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
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
						p_cla[n_states] = real_rescale_ln + f_cla[n_states] + s_cla[n_states];
#					endif
					for (from_state = 0; from_state < n_states; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ = prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[n_states] = f_cla[n_states] + s_cla[n_states];
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
	const unsigned n_states = par->n_states;
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
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					s_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
						s_des_prob += s_cla[to_state]*(*s_categ_p++);
					}
					p_cla[from_state] = f_des_prob*s_des_prob;
				}
#				if ! SEPARATE_RESCALER_ARRAY
						p_cla[n_states] = f_cla[n_states] + s_cla[n_states];
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
	assert(fd->n_states == sd->n_states);
	assert(fd->n_states == par->n_states);
	assert(fd->cso.n_categ_or_subs == sd->cso.n_categ_or_subs);
	assert(fd->cso.n_categ_or_subs == par->cso.n_categ_or_subs);
	const double *** f_pmat_arr = (const double ***)fmat->p_mat;
	const double *** s_pmat_arr = (const double ***)smat->p_mat;
	if ((fd->n_edges_since_rescaling_check + sd->n_edges_since_rescaling_check + 2) < context.rescale_threshold)
		return do_internal_cla_no_rescale(fd, f_pmat_arr, sd, s_pmat_arr, par, context);
	return do_internal_cla_rescale(fd, f_pmat_arr, sd, s_pmat_arr, par, context);
}

static int do_add_leaf_to_cla_recalc_rescale(const LeafDataObj* fd, CLAObj* par, const CalcContext context) {
	const double ** f_mat;
	const unsigned n_subsets = (par->cso.n_categ_arr == 0L ? 1 : par->cso.n_categ_or_subs);
	const unsigned total_n_sites = par->cso.total_n_sites;
	const unsigned n_states = par->n_states;
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
	

		const unsigned tpm_offset = n_states*(categ_beg);
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
				for (state = 0; state < n_states; ++state) {
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
						p_cla[n_states] += real_rescale_ln;
#					endif
					for (state = 0; state < n_states; ++state)
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
	const unsigned n_states = par->n_states;
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
	
		const unsigned tpm_offset = n_states*(categ_beg);
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
				for (state = 0; state < n_states; ++state) {
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
	assert(fd->n_states == par->n_states);
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
	if ((fd->pmat_obj_ptr != (void *) fmat) || (fd->calc_time_stamp != fmat->calc_time_stamp)) {
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
	const unsigned n_states = par->n_states;
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
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
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
						p_cla[n_states] += real_rescale_ln + f_cla[n_states];
#					endif
					for (from_state = 0; from_state < n_states; ++from_state)
						p_cla[from_state] /= real_rescale;
				}
				else {
#					if SEPARATE_RESCALER_ARRAY
						*rescale_count_ptr++ += prev_rcount;
#					else  /*SEPARATE_RESCALER_ARRAY */
						p_cla[n_states] += f_cla[n_states];
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
	const unsigned n_states = par->n_states;
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
				for (from_state = 0; from_state < n_states; ++from_state) {
					f_des_prob = 0.0;
					for (to_state = 0; to_state < n_states; ++to_state) {
						f_des_prob += f_cla[to_state]*(*f_categ_p++);
					}
					p_cla[from_state] *= f_des_prob;
				}
#			if ! SEPARATE_RESCALER_ARRAY
						p_cla[n_states] += f_cla[n_states];
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
	assert(fd->n_states == par->n_states);
	assert(fd->cso.max_n_categ == par->cso.max_n_categ);
	const double *** f_pmat_arr = (const double ***)fmat->p_mat;
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
	const unsigned n_states = full_la->n_states;
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
	double * const categ_like = f_mult + (n_states*max_n_categ);
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
			cat_freq = sf_arr[n_states];
			for (state = 0; state < n_states; ++state) {
				f_mult[n_states*categ + state] = cat_freq*sf_arr[state];
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
			
		PRINTF("categ likes, site lnL, site weight:\n");
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
			PRINTF1("%d", site);
			for (categ = 0; categ < n_categ; ++categ) {
				*curr_categ_like = 0.0;
				for (state = 0; state < n_states; ++state) {
					*curr_categ_like += (*curr_mult++)*(*cla++);
				}
				PRINTF1("\t%g", *curr_categ_like);
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
			PRINTF2("\t%g\t%f\n", site_ln_L, wt);
		}
		if (subsetLnL)
			subsetLnL[subset_ind] = cumul_subset_ln_L;
		total_ln_L += cumul_subset_ln_L;
		
	}
	return total_ln_L;
}

/* ctor ("internal" factory-function) and dtor */
StateSetLookupStruct* sslookup_new(unsigned n_states, unsigned n_state_sets, int ** p)  {
	PRINTF("In sslookup_new\n");
	StateSetLookupStruct * sslookup = PyObject_New(StateSetLookupStruct, &state_set_lookup_type);
	if (sslookup) {
		sslookup->state_lookup = p;
		sslookup->n_states = n_states;
		sslookup->n_state_sets = n_state_sets;
	}
	return sslookup;
}


LeafDataObj* leaf_data_ss_partitioned_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl)  {
	assert(n_sites > 0);
	assert(n_categ > 0);
	assert(ssl);
	PRINTF("In leaf_data_new\n");
	const unsigned len_pmat_columns_row = ssl->n_states*n_categ;
	LeafDataObj * leaf_data = PyObject_New(LeafDataObj, &leaf_data_type);
	if (leaf_data) {
		Py_INCREF(ssl);
		leaf_data->sslookup = ssl;
		leaf_data->n_states = ssl->n_states;
		leaf_data->n_state_sets = ssl->n_state_sets;
		leaf_data->state_lookup = ssl->state_lookup;
		leaf_data->ssind = 0L;
		leaf_data->pmat_columns = 0L;
		leaf_data->calc_time_stamp = -1;
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
	PRINTF("In leaf_data_new errorExit\n");
	leaf_data_dtor(leaf_data);
	return 0L;
}

LeafDataObj* leaf_data_new(unsigned n_sites, unsigned n_categ, StateSetLookupStruct * ssl)  {
	assert(n_sites > 0);
	assert(n_categ > 0);
	assert(ssl);
	PRINTF("In leaf_data_new\n");
	const unsigned len_pmat_columns_row = ssl->n_states*n_categ;
	LeafDataObj * leaf_data = PyObject_New(LeafDataObj, &leaf_data_type);
	if (leaf_data) {
		Py_INCREF(ssl);
		leaf_data->sslookup = ssl;
		leaf_data->n_states = ssl->n_states;
		leaf_data->n_state_sets = ssl->n_state_sets;
		leaf_data->state_lookup = ssl->state_lookup;
		leaf_data->ssind = 0L;
		leaf_data->pmat_columns = 0L;
		leaf_data->calc_time_stamp = -1;
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
		PRINTF("In leaf_data_new errorExit\n");
		leaf_data_dtor(leaf_data);
		return 0L;
}

/*
	Allocates a CLAObj object for an unpartitioned calculation that mixes n_categ models. 
*/
CLAObj * cla_new(unsigned n_sites, unsigned n_states,  unsigned n_categ)  {
	unsigned cla_len;
#	if SEPARATE_RESCALER_ARRAY
		unsigned rescale_len;
#	endif
	assert(n_sites > 0);
	assert(n_states > 0);
	assert(n_categ > 0);
	PRINTF("In cla_new\n");
	CLAObj * cla = PyObject_New(CLAObj, &cla_type);
	if (cla) {
		cla->cso.total_n_sites = n_sites;
		cla->n_states = n_states;
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
			cla->len_per_categ = n_states + 1;/*The +1 is for the rescale ln*/
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
		PRINTF("In cla_new errorExit\n");
		cla_dtor(cla);
		return 0L;
}


/* Allocates the field of a CategSubsetObj object for a model in which each site gets its own model 
	(and there are no mixtures used).
	returns 0 if there is an error.
*/
int new_cso_ss_partitioned(CategSubsetObj * cso, unsigned n_sites, unsigned n_states) {
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
CLAObj* cla_ss_partitioned_new(unsigned n_sites, unsigned n_states) {
	unsigned cla_len;
	unsigned n_categ = 1;
#	if SEPARATE_RESCALER_ARRAY
		unsigned rescale_len;
#	endif
	assert(n_sites > 0);
	assert(n_states > 0);
	assert(n_categ > 0);
	PRINTF("In cla_ss_partitioned_new\n");
	CLAObj * cla = PyObject_New(CLAObj, &cla_type);
	if (cla) {
		cla->cso.total_n_sites = n_sites;
		cla->n_states = n_states;
		if (new_cso_ss_partitioned(&(cla->cso), n_sites, n_states) == 0)
			goto errorExit;
#		if SEPARATE_RESCALER_ARRAY
			cla->len_per_categ = n_sites;
			rescale_len = n_sites*n_categ;
#		else
			cla->len_per_categ = n_states + 1;/*The +1 is for the rescale ln*/
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
		PRINTF("In cla_ss_partitioned_new errorExit\n");
		cla_dtor(cla);
		return 0L;
}

/* Allocates a FullLAObj object for a model in which each site gets its own model 
	(and there are no mixtures used).
*/
FullLAObj* full_la_ss_partitioned_new(unsigned n_sites, unsigned n_states) {
	PRINTF("In full_la_ss_partitioned_new\n");
	unsigned i;
	
	const unsigned lnL_wt_len = 2*(n_sites);
	FullLAObj * full_la = PyObject_New(FullLAObj, &full_la_type);
	if (full_la) {
		full_la->cso.total_n_sites = n_sites;
		full_la->n_states = n_states;
		full_la->full_cla  = 0L;
		full_la->pat_lnL_wts  = 0L;
		full_la->state_categ_freqs  = 0L;
		if (new_cso_ss_partitioned(&(full_la->cso), n_sites, n_states) == 0)
			goto errorExit;
		full_la->full_cla = cla_ss_partitioned_new(n_sites, n_states);
		if (full_la->full_cla == 0L)
			goto errorExit;
		full_la->pat_lnL_wts  = (cla_float_t *) malloc(lnL_wt_len*sizeof(cla_float_t));
		if (full_la->pat_lnL_wts == 0L)
			goto errorExit;
		for (i = 0; i < n_sites; ++i)
			full_la->pat_lnL_wts[1 + (2*i)] = 1.0;
		full_la->scratch = (double *) malloc(2*n_states*sizeof(double));
		if (full_la->scratch == 0L)
			goto errorExit;
		full_la->state_categ_freqs = allocateDblMatrix(n_sites, n_states + 1);
		if (full_la->state_categ_freqs == 0L)
			goto errorExit;
	}
	return full_la;
	errorExit:
		PRINTF("In full_la_ss_partitioned_new errorExit\n");
		full_la_dtor(full_la);
		return 0L;
}

/*
	Allocates a FullLAObj object for an unpartitioned calculation that mixes n_categ models. 
*/
FullLAObj * full_la_new(unsigned n_sites, unsigned n_states,  unsigned n_categ)  {
	PRINTF("In full_la_new\n");
	unsigned i;
	const unsigned lnL_wt_len = 2*(n_sites);
	FullLAObj * full_la = PyObject_New(FullLAObj, &full_la_type);
	if (full_la) {
		full_la->cso.total_n_sites = n_sites;
		full_la->n_states = n_states;
		full_la->cso.n_categ_or_subs = n_categ;
		full_la->cso.max_n_categ = n_categ;
		full_la->full_cla  = 0L;
		full_la->pat_lnL_wts  = 0L;
		full_la->state_categ_freqs  = 0L;
		full_la->cso.n_categ_arr = 0L;
		full_la->cso.subset_component_offset = 0L;
		full_la->cso.subset_char_offsets = 0L;
		full_la->cso.subset_cla_offsets = 0L;
		full_la->full_cla = cla_new(n_sites, n_states,  n_categ);
		if (full_la->full_cla == 0L)
			goto errorExit;
		full_la->pat_lnL_wts  = (cla_float_t *) malloc(lnL_wt_len*sizeof(cla_float_t));
		if (full_la->pat_lnL_wts == 0L)
			goto errorExit;
		for (i = 0; i < n_sites; ++i)
			full_la->pat_lnL_wts[1 + (2*i)] = 1.0;
		full_la->scratch = (double *) malloc(2*n_states*n_categ*sizeof(double));
		if (full_la->scratch == 0L)
			goto errorExit;
		full_la->state_categ_freqs = allocateDblMatrix(n_categ, n_states + 1);
		if (full_la->state_categ_freqs == 0L)
			goto errorExit;
	}
	return full_la;
	errorExit:
		PRINTF("In full_la_new errorExit\n");
		full_la_dtor(full_la);
		return 0L;
}
DSCTModelObj* dsct_model_new(unsigned dim)  {
	assert(dim > 1);
	PRINTF("In cdsctm_ctor\n");
	unsigned arr_len = (3 + dim*dim) * dim;
	DSCTModelObj * dsct_model = PyObject_New(DSCTModelObj, &dsct_model_type);
	if (dsct_model) {
		dsct_model->dim = dim;
		dsct_model->eigen_calc_dirty = 1;
		dsct_model->mat_del_handle = 0L;
		dsct_model->arr_del_handle = 0L;
		dsct_model->q_mat = 0L;
		dsct_model->eigen_vectors = 0L;
		dsct_model->inv_eigen_vectors = 0L;
		dsct_model->work_mat = 0L;
		dsct_model->eigen_values = 0L;
		dsct_model->im_eigen_values = 0L;
		dsct_model->cijk = 0L;
		dsct_model->d_work = 0L;
		dsct_model->i_work = 0L;
		dsct_model->mat_del_handle = allocateDbl3DMatrix(4, dim, dim);
		if (dsct_model->mat_del_handle == 0L)
			goto errorExit;
		dsct_model->arr_del_handle = (double *)malloc(arr_len*sizeof(double));
		if (dsct_model->arr_del_handle == 0L)
			goto errorExit;
		dsct_model->i_work = (int *)malloc(dim*sizeof(int));
		if (dsct_model->i_work == 0L)
			goto errorExit;
		dsct_model->q_mat = dsct_model->mat_del_handle[0];
		dsct_model->eigen_vectors = dsct_model->mat_del_handle[1];
		dsct_model->inv_eigen_vectors = dsct_model->mat_del_handle[2];
		dsct_model->work_mat = dsct_model->mat_del_handle[3];
		dsct_model->eigen_values = dsct_model->arr_del_handle;
		dsct_model->im_eigen_values = dsct_model->eigen_values + dim;
		dsct_model->d_work = dsct_model->im_eigen_values + dim;
		dsct_model->cijk = dsct_model->d_work + dim;
	}
	return dsct_model;
	errorExit:
		PRINTF("In cdsctm_ctor errorExit\n");
		cdsctm_dtor(dsct_model);
		return 0L;
}
void sslookup_dtor(StateSetLookupStruct* ssl_struct) {
	PRINTF("In sslookup_dtor\n");
	if (ssl_struct == 0L)
		return;
	if (ssl_struct->state_lookup) {
		if (ssl_struct->state_lookup[0])
			free(ssl_struct->state_lookup[0]);
		free(ssl_struct->state_lookup);
	}
	PyObject_Del(ssl_struct);
	PRINTF("Leaving sslookup_dtor\n");
}
void leaf_data_dtor(LeafDataObj * leaf_data) {
	PRINTF("In leaf_data_dtor\n");
	if (leaf_data == 0L)
		return;
	if (leaf_data->sslookup) {
		Py_DECREF(leaf_data->sslookup);
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
	PRINTF("In cla_dtor\n");
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
	PRINTF("In full_la_dtor\n");
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

void cdsctm_dtor(DSCTModelObj* dsct_model) {
	PRINTF("In cdsctm_dtor\n");
	if (dsct_model == 0L)
		return;
	if (dsct_model->mat_del_handle)
		freeDbl3DMatrix(dsct_model->mat_del_handle);
	if (dsct_model->arr_del_handle)
		free(dsct_model->arr_del_handle);
	if (dsct_model->i_work)
		free(dsct_model->i_work);
	PyObject_Del(dsct_model);
}



PMatArrayObj* cpmat_array_new(unsigned n_mat, unsigned n_states)  {
	assert(n_mat >= 1);
	PRINTF("In cpmat_array_new\n");
	PMatArrayObj * pmat_array_obj = PyObject_New(PMatArrayObj, &pmat_array_type);
	if (pmat_array_obj) {
		pmat_array_obj->n_mat = n_mat;
		pmat_array_obj->n_states = n_states;
		pmat_array_obj->model_aliases = 0L;
		pmat_array_obj->brlen_aliases = 0L;
		pmat_array_obj->p_mat = allocateDbl3DMatrix(n_mat, n_states, n_states);
		pmat_array_obj->calc_time_stamp = 0;
		PRINTF("PMat alloc\n");
		if (pmat_array_obj->p_mat == 0L)
			goto errorExit;
		pmat_array_obj->model_aliases = (DSCTModelObj **)malloc(n_mat*sizeof(DSCTModelObj *));
		pmat_array_obj->brlen_aliases = (double*)malloc(n_mat*sizeof(double));
		pmat_array_obj->first_el = **pmat_array_obj->p_mat;
	}
	return pmat_array_obj;
	errorExit:
		PRINTF("In cpmat_array_new errorExit\n");
		cpmat_array_dtor(pmat_array_obj);
		return 0L;
}

int next_calc_stamp() {
	static int gCalcStamp = 0; /*used like a time stamp to flag when pmats were calculated */
	++gCalcStamp;
	if (gCalcStamp < 0)
		gCalcStamp = 0;
	return gCalcStamp;
}



void cpmat_array_dtor(PMatArrayObj* pmat_array_obj) {
	PRINTF("In pmat_array_dtor\n");
	if (pmat_array_obj == 0L)
		return;
	if (pmat_array_obj->p_mat)
		freeDbl3DMatrix(pmat_array_obj->p_mat);
	if (pmat_array_obj->model_aliases)
		free(pmat_array_obj->model_aliases);
	if (pmat_array_obj->brlen_aliases)
		free(pmat_array_obj->brlen_aliases);
	PyObject_Del(pmat_array_obj);
}

void internal_asrv_set_shape(ASRVObj *asrh, double val) {
	double beta;
	asrh->param = val;
	beta = 1.0/val;
	DiscreteGamma(asrh->freq, asrh->val, val, beta, asrh->n, asrh->style);
}


ASRVObj* asrv_obj_new(unsigned dim, int style, double param)  {
	const unsigned arr_len = dim*sizeof(double);
	assert(dim > 1);
	PRINTF("In asrv_obj_new\n");
	ASRVObj * asrh = PyObject_New(ASRVObj, &asrv_type);
	if (asrh) {
		asrh->n = dim;
		asrh->val = 0L;
		asrh->freq = 0L;
		asrh->style = style;
		asrh->param = param;
		if (dim > 0) {
			asrh->val = (double*)malloc(arr_len*sizeof(double));
			if (asrh->val == 0L) {
				PyErr_NoMemory();
				goto errorExit;
			}
			asrh->freq = (double*)malloc(arr_len*sizeof(double));
			if (asrh->val == 0L) {
				goto errorExit;
			}
		internal_asrv_set_shape(asrh, param);
		}
	}
	return asrh;
	errorExit:
		PRINTF("In asrv_obj_new errorExit\n");
		asrv_obj_dtor(asrh);
		return 0L;
}

void asrv_obj_dtor(ASRVObj * asrh) {
	PRINTF("In asrh dtor\n");
	if (asrh == 0L)
		return;
	if (asrh->val != 0L)
		free(asrh->val);
	if (asrh->freq != 0L)
		free(asrh->freq);
	PyObject_Del(asrh);
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
