#if ! defined(NON_BEAGLE_IMPL_H)
#define NON_BEAGLE_IMPL_H

#ifdef __cplusplus
extern "C" 
{
#endif

/* All functions that return an int return 0 for failure and non-zero for success */



int prob_mat_from_eigensystem (
		const unsigned dim,      /* IN the number of states*/
		double **pMat,			 /* OUT Filled on output  dim x dim transition probability matrix*/
		double *workspace,		 /* SCRATCH array of length `dim` */
		const double *cijk, 	 /* IN the dim*dim*dim array of temporaries used in quickly multiplying
									   the exp(eval*branchlength) times the matrix eigenvectors 
									   and its inverse. Call calc_c_ijk first to fill this */
		const double *eigenVals, /* IN array of the eigenvalues*/
		const double branchLen); /* IN length of branch - rate x time duration */



void calc_c_ijk(
	unsigned dim,  /* IN number of states */
	double *c_ijk, /* OUT must be dim*dim*dim array.  Will be filled with temporaries
						used in prob_mat_from_eigensystem*/
	const double **eigenVectors, /* IN dim*dim matrix of eigenvectors */
	const double **invEigenVectors); /* IN dim*dim inverse of matrix of eigenvectors */


#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif
