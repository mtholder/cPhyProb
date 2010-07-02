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


struct LikeCalculatorInstance;

typedef struct {
		double ** workMat; /*dim by dim*/
		double * dWork;  /*len dim*/
		int * iWork; /*len dim*/	
} EigenCalcScratchpad;

typedef struct {
#	if defined(USE_BEAGLE_LIB)
		int beagleEigenBufferIndex;
#	else
		double * cijk; /*len dim^3 */ /* alias */
#	endif
	double * eigenValues; /* alias */
	double * imEigenValues; /* alias */
	double ** eigenVectors;
	double ** invEigenVectors;
	EigenCalcScratchpad * scratchPad;
	unsigned int dim;
	double *** matDelHandle; 
	double * arrDelHandle; 
} EigenSolutionStruct;


/* Define the type that corresponds to a model (holds the Q-Matrix and other
	temporary fields used to speed up calculation).
*/
typedef struct {
	PyObject_HEAD
	unsigned dim;
	double **qMat;
	struct LikeCalculatorInstance * likeCalcInstanceAlias;
	int eigenBufferIndex; /* index where this eigen system is currently stored */
	

	int eigenCalcIsDirty; /*0 if the eigen calculations are up-to-date*/
} DSCTModelObj;
/* Define the type that corresponds to a model (holds the Q-Matrix and other
	temporary fields used to speed up calculation).
*/

EigenSolutionStruct * getEigenSolutionStruct(DSCTModelObj *, int essIndex);
EigenCalcScratchpad * getEigenCalcScratchpad(DSCTModelObj *, int essIndex);



typedef struct {
	PyObject_HEAD
	unsigned nMatrices;
	unsigned nStates;
	double *** pMatrix; /* [categ_ind][from_state][to_state] */
	double * firstEl; /*alias to the first element*/
	DSCTModelObj ** modelAliases; /*aliases to the model used for this pmat (may not be stable over repeated invocations -- these are used as scratch) */
	double * brlenAliases; /*copies of the branch length used for this pmat (may not be stable over repeated invocations -- these are used as scratch) */
	int calcTimeStamp;
} PMatArrayObj;
#if defined(BUILDING_FOR_PYTHON)
	/* Forward declare the python type object*/
	staticforward PyTypeObject pmat_array_type;
#endif

PMatArrayObj* cpmatArrayNew(unsigned nMatrices, unsigned nStates);
void cpmatArrayDtor(PMatArrayObj* pmat_array_obj);


/* BEAGLE style */
/* Calculate */
int calcPmatInBuffers(DSCTModelObj *mod, const int * indices, const double *  branchLengths, int numMats);
/* Retreive */
int getPmatInBuffers(DSCTModelObj *mod, const int * indices, double ***pMatPtrs, int numMats);


EigenSolutionStruct * eigenSolutionStructNew(unsigned dim);
EigenCalcScratchpad * eigenSolutionScratchpadNew(unsigned dim);
void eigenSolutionStructDtor(EigenSolutionStruct * p);
void eigenSolutionScratchpadDtor(EigenCalcScratchpad *);


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
