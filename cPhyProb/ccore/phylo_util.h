#if ! defined(PHYLO_UTIL_H)
#define PHYLO_UTIL_H

#ifdef __cplusplus
extern "C" 
{
#endif

double **allocateDblMatrix(unsigned n_rows, unsigned n_cols);
double ***allocateDbl3DMatrix(unsigned nm, unsigned nr, unsigned nc);
void freeDblMatrix(double **p);
void freeDbl3DMatrix(double ***p);
void DiscreteGamma(double *f,double *r, double alpha, double beta, int ncat, int useMean);

int get_eigens(unsigned dim,double **q, double *eigenValues, double *imEigenValues, double **eigenVectors, double **invEigenVectors, double **workMat, double *dWork, int *iWork, unsigned *is_complex);


#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif
