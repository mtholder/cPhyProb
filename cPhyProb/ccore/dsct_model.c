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

#if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
#	include "libhmsbeagle/beagle.h"
#endif

#include "dsct_model.h"
#include "phylo_util.h"
#if ! defined (USE_BEAGLE_LIB)
#	include "non_beagle_impl.h"
#endif

#if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
#	include "libhmsbeagle/beagle.h"
	INLINE void BEAGLE_ERROR_MSG(int rc, const char * funcName) {
			if (rc != BEAGLE_SUCCESS) {
				CPHYPROB_DEBUG_PRINTF2("Beagle error code %d in function \"%s\"", rc, funcName);
			}
		}
#endif




int do_asrv_pmat_array_calc(PMatArrayObj *p, DSCTModelObj * mod, ASRVObj *asrv, double edgeLen);
int calc_pmat_array(PMatArrayObj *p, unsigned start_ind, unsigned end_ind);
int calc_pmat(double **p, DSCTModelObj *mod, double branch_len);



static int next_calc_stamp(void);

static int recalc_eigen_mat(DSCTModelObj *mod);

INLINE int prob_mat_from_clean_model(double **p, DSCTModelObj *mod, double branch_len);

int next_calc_stamp(void) {
	static int gCalcStamp = 0; /*used like a time stamp to flag when pmats were calculated */
	++gCalcStamp;
	if (gCalcStamp < 0)
		gCalcStamp = 0;
	return gCalcStamp;
}


EigenSolutionStruct * getEigenSolutionStruct(DSCTModelObj * mod, int essIndex) {
	assert(0);
	/* @TEMPORARY */
	static EigenSolutionStruct s;
	return &s;
}

EigenCalcScratchpad * getEigenCalcScratchpad(DSCTModelObj * mod, int essIndex) {
	assert(0);
	/* @TEMPORARY */
	static EigenCalcScratchpad s;
	return &s;
}

static void printQMat(double ** qmat_obj, const unsigned nStates) {
	unsigned from_state, to_state;
	CPHYPROB_DEBUG_PRINTF("QMat:\n");
	for (from_state = 0; from_state < nStates; ++from_state) {
		for (to_state = 0; to_state < nStates; ++to_state) {
			CPHYPROB_DEBUG_PRINTF1("%f ", qmat_obj[from_state][to_state]);
		}
		CPHYPROB_DEBUG_PRINTF("\n");
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
		printQMat(mod->qMat, mod->dim);
#	endif

	int essIndex = 0 ; /* @TEMPORARY!!! */
	EigenSolutionStruct * eigenSolutionStruct = getEigenSolutionStruct(mod, essIndex);
	EigenCalcScratchpad * eigenCalcScratchpad = getEigenCalcScratchpad(mod, essIndex);

	rc = get_eigens(mod->dim,
					mod->qMat,
					eigenSolutionStruct->eigenValues,
					eigenSolutionStruct->imEigenValues,
					eigenSolutionStruct->eigenVectors,
					eigenSolutionStruct->invEigenVectors,
					eigenCalcScratchpad->workMat,
					eigenCalcScratchpad->dWork,
					eigenCalcScratchpad->iWork,
					&is_complex);
	if (rc == 0) {
		if (is_complex == 1)
			PyErr_SetString(PyExc_RuntimeError, "Error: Complex eigenvalues found.");
		return 0;
	}
#	if ! defined(USE_BEAGLE_LIB)
		calc_c_ijk(mod->dim, eigenSolutionStruct->cijk,
			       (const double **)eigenSolutionStruct->eigenVectors,
			       (const double **)eigenSolutionStruct->invEigenVectors);
#	endif
	mod->eigenCalcIsDirty = 0;
	return 1;
}



INLINE int prob_mat_from_clean_model(double **p, DSCTModelObj *mod, double branch_len) {
	assert(mod);
	assert(p);
	assert(*p);
	struct LikeCalculatorInstance * likeCalcInstanceAlias = mod->likeCalcInstanceAlias;
	if (!likeCalcInstanceValid(likeCalcInstanceAlias))
		return 0;

	int essIndex = 0 ; /* @TEMPORARY!!! */
	EigenSolutionStruct * eigenSolutionStruct = getEigenSolutionStruct(mod, essIndex);

#	if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
		assert(eigenSolutionStruct->beagleEigenBufferIndex >= 0);
		int probMatIndex;
		getProbMatIndex(likeCalcInstanceAlias, &probMatIndex, 1);
		int rc = beagleUpdateTransitionMatrices(likeCalcInstanceAlias->beagleInstanceIndex,
									   eigenSolutionStruct->beagleEigenBufferIndex, /* index of eigen solution in beagle */
									   &probMatIndex, /* index of prob mat in beagle */
									   0L, /* no first derivs */
									   0L, /* no second derivs */
									   &branch_len, /* array of branch lengths */
									   1); /* only one p mat being updated this way */
		BEAGLE_ERROR_MSG(rc, "beagleUpdateTransitionMatrices");
		if (rc == BEAGLE_SUCCESS) {
			rc = beagleGetTransitionMatrix(likeCalcInstanceAlias->beagleInstanceIndex,
									       eigenSolutionStruct->beagleEigenBufferIndex,
									       p[0]);
			BEAGLE_ERROR_MSG(rc, "beagleGetTransitionMatrix");
		}
		return (rc == BEAGLE_SUCCESS ? 1 : 0);
#	else
		EigenCalcScratchpad * eigenCalcScratchpad = getEigenCalcScratchpad(mod, essIndex);

		return prob_mat_from_eigensystem (mod->dim,
										  p,
										  eigenCalcScratchpad->dWork,
										  eigenSolutionStruct->cijk,
										  eigenSolutionStruct->eigenValues,
										  branch_len);
#	endif
}









/**
 * Calculates a P-Matrix from mod and branch_len.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int calc_pmat(double **p, DSCTModelObj *mod, double branch_len) {
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
	CPHYPROB_DEBUG_PRINTF("About to check if eigensystem is calculated\n");
	if (mod->eigenCalcIsDirty) {
		CPHYPROB_DEBUG_PRINTF("About to recalc eigensystem\n");
		if (!recalc_eigen_mat(mod))
			return 0;
	}
	CPHYPROB_DEBUG_PRINTF("About to calc pmats from eigensystem\n");
	if (!prob_mat_from_clean_model(p, mod, branch_len)) {
		PyErr_SetString(PyExc_ValueError, "prob_mat_from_clean_model failed");
		return 0;
	}
	return 1;
}


/**
 * Calculates a P-Matrix from mod and branch_len.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int calcPmatInBuffers(DSCTModelObj *mod, const int * indices, const double * branchLengths, int numMats) {
	int i;
	if (mod == 0L) {
		PyErr_SetString(PyExc_ValueError, "Illegal NULL Model pointer");
		return 0;
	}
	if (indices == 0L) {
		PyErr_SetString(PyExc_ValueError, "Illegal NULL indices pointer");
		return 0;
	}
	if (branchLengths == 0L) {
		PyErr_SetString(PyExc_ValueError, "Illegal NULL branchLengths pointer");
		return 0;
	}
	if (numMats < 0) {
		PyErr_SetString(PyExc_ValueError, "number of probability matrices to calculate cannot be negative");
		return 0;
	}
	for (i = 0; i < numMats; ++i) {
		if (branchLengths[i] < 0.0) {
			PyErr_SetString(PyExc_ValueError, "Illegal (negative) branch length");
			return 0;
		}
		if (indices[i] < 0) {
			PyErr_SetString(PyExc_ValueError, "Illegal P-Matrix index < 0");
			return 0;
		}
	}
	CPHYPROB_DEBUG_PRINTF("About to check if eigensystem is calculated\n");
	if (mod->eigenCalcIsDirty) {
		CPHYPROB_DEBUG_PRINTF("About to recalc eigensystem\n");
		if (!recalc_eigen_mat(mod))
			return 0;
	}
	struct LikeCalculatorInstance * likeCalcInstanceAlias = mod->likeCalcInstanceAlias;
	if (!likeCalcInstanceValid(likeCalcInstanceAlias)) {
		PyErr_SetString(PyExc_ValueError, "likeCalcInstanceAlias from model is invalid");
		return 0;
	}

	int essIndex = 0 ; /* @TEMPORARY!!! */
	EigenSolutionStruct * eigenSolutionStruct = getEigenSolutionStruct(mod, essIndex);

#	if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
		assert(eigenSolutionStruct->beagleEigenBufferIndex >= 0);
		int rc = beagleUpdateTransitionMatrices(likeCalcInstanceAlias->beagleInstanceIndex,
									   eigenSolutionStruct->beagleEigenBufferIndex, /* index of eigen solution in beagle */
									   indices, /* index of prob mat in beagle */
									   0L, /* no first derivs */
									   0L, /* no second derivs */
									   branchLengths, /* array of branch lengths */
									   numMats); /* only one p mat being updated this way */
		BEAGLE_ERROR_MSG(rc, "beagleUpdateTransitionMatrices");
		return (rc == BEAGLE_SUCCESS ? 1 : 0);
#	else
		double ***pMatArray = getProbMatPtrAlias(likeCalcInstanceAlias, indices, numMats);
		if (!pMatArray) {
			PyErr_SetString(PyExc_ValueError, "number of probability matrices exceeds working space");
			return 0;
		}
		EigenCalcScratchpad * eigenCalcScratchpad = eigenSolutionStruct->scratchPad;

		for (i = 0; i < numMats; ++i) {
			double ** p = pMatArray[i];
			int rc = prob_mat_from_eigensystem (mod->dim,
												p,
												eigenCalcScratchpad->dWork,
												eigenSolutionStruct->cijk,
												eigenSolutionStruct->eigenValues,
												branchLengths[i]);
			if (!rc) {
				PyErr_SetString(PyExc_ValueError, "prob_mat_from_eigensystem failed");
				return 0;
			}
		}
#	endif
	return 1;
}



void cpmatArrayDtor(PMatArrayObj* pmat_array_obj) {
	CPHYPROB_DEBUG_PRINTF("In pmat_array_dtor\n");
	if (pmat_array_obj == 0L)
		return;
	if (pmat_array_obj->pMatrix)
		freeDbl3DMatrix(pmat_array_obj->pMatrix);
	if (pmat_array_obj->modelAliases)
		free(pmat_array_obj->modelAliases);
	if (pmat_array_obj->brlenAliases)
		free(pmat_array_obj->brlenAliases);
	PyObject_Del(pmat_array_obj);
}




/**
 * Calculates the p->pMatrix[i] for i in [start_ind, end_ind)
 * the branch length and models used come from p->modelAliases, p->brlenAliases
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int calc_pmat_array(PMatArrayObj *p, unsigned start_ind, unsigned end_ind) {
	unsigned i;
	if (p == 0L) {
		PyErr_SetString(PyExc_TypeError, "Expecting a PMatArrayObj object");
		return 0;
	}
	if (p->pMatrix == 0L || p->modelAliases == 0L || p->brlenAliases == 0L) {
		PyErr_SetString(PyExc_ValueError, "Invalid PMatArrayObj object");
		return 0;
	}
	if (p->nMatrices < start_ind || p->nMatrices < (end_ind -1)) {
		PyErr_SetString(PyExc_IndexError, "End matrix exceeds size of PMatrix Array");
		return 0;
	}
	for (i = start_ind; i < end_ind; ++i) {
		if (!calc_pmat(p->pMatrix[i], p->modelAliases[i], p->brlenAliases[i])) {
			return 0;
		}
	}
	p->calcTimeStamp = next_calc_stamp();
	return 1;
}



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
