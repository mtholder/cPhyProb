#include "cphyprob_defs.h"
#include "cphyprob.h"

#if defined(USE_BEAGLE_LIB)
	int getProbMatIndex(struct LikeCalculatorInstance * b, int *indices, int numIndices) {
		return 1; /* @TEMPORARY */
	}
#else
	double *** getProbMatPtrAlias(struct LikeCalculatorInstance * b, const int * indices, int numIndices) {
		int i = 0;
		assert(b);
		assert(indices);
		for (; i < numIndices; ++i) {
			const int index = indices[i];
			b->probMatAliasArray[i] = b->probMatArray[index];
		}
		return b->probMatAliasArray; /* NOT thread safe assumes that only will thread will ask for this alias buffer at a time */
	}

#endif /* USE_BEAGLE_LIB*/
