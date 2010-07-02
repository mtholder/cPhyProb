#include "asrv.h"
#include "phylo_util.h"

void internal_asrv_set_shape(ASRVObj *asrh, double val) {
	CPHYPROB_DEBUG_PRINTF("In internal_asrv_set_shape\n");
	double beta;
	asrh->param = val;
	beta = 1.0/val;
	DiscreteGamma(asrh->freq, asrh->val, val, beta, asrh->n, asrh->style);
}



void asrv_obj_dtor(ASRVObj * asrh) {
	CPHYPROB_DEBUG_PRINTF("In asrh dtor\n");
	if (asrh == 0L)
		return;
	if (asrh->val != 0L)
		free(asrh->val);
	if (asrh->freq != 0L)
		free(asrh->freq);
	PyObject_Del(asrh);
}


