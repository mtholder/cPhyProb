#include <assert.h>
#include <math.h>
#include "non_beagle_impl.h"
/**
 * Adapted from MrBayes TiProbsUsingEigens in mbmath.c
 *	`pMat` and `EigValExp` are written to.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int prob_mat_from_eigensystem (
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
void calc_c_ijk(unsigned dim, double *c_ijk, const double **eigenVectors, const double **invEigenVectors) {
	register int 	i, j, k;
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			for (k = 0; k < dim; k++)
			 	*c_ijk++ = eigenVectors[i][k] * invEigenVectors[k][j];
}

