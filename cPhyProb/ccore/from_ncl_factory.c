//	Copyright (C) 2008 Mark T. Holder
//
//	chimne_ssweep is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	chimne_ssweep is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details. 
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

#include "dsct_model.h"
#include "from_ncl_factory.h"
#if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
#	include "libhmsbeagle/beagle.h"
#endif

#include "ncl/nxscdiscretematrix.h"




LikeStructsBundle newPartitionedLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nPartials, 
  unsigned nPMats,
  unsigned nSubsets)
{
	unsigned arr_len = 2*matrix.nStates; 
	unsigned i, j;
	int ** p = 0L;
	LikeStructsBundle toReturn;
	zeroLikeStructFields(&toReturn);
	toReturn.nModels = nSubsets;

	if (!(matrix.matrix && matrix.stateList))
		goto errorExit;
		

	/*Allocate p - a copy of the stateList*/
	/* each of the "fundamental" states needs to spaces in the array" */
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		arr_len += 1; //add one for the # of states element
		unsigned pos = matrix.stateListPos[i] ;
		arr_len += matrix.stateList[pos];
	}

	/*allocate a ragged two-D array with the memory contiguous.*/
	p = (int**)malloc(matrix.nObservedStateSets*sizeof(int*));
	if (p == 0L)
		goto errorExit;	
	p[0] = (int*)malloc(arr_len*sizeof(int));
	if (p[0] == 0L) 
		goto errorExit;

	int * curr_p = p[0];
	for (i = 0 ; i < matrix.nStates; ++i) {
		p[i] = curr_p;
		p[i][0] = 1;
		p[i][1] = i;
		curr_p += 2;
	}
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		p[i] = curr_p;
		unsigned pos = matrix.stateListPos[i];
		p[i][0] = matrix.stateList[pos];
		curr_p++;
		unsigned j;
		for (j = 1; j <= p[i][0]; ++j) {
			*curr_p++ = matrix.stateList[pos + j];
		}
	}
	
	toReturn.sharedStateSetLookupStruct = sslookup_new(matrix.nStates, matrix.nObservedStateSets, p);
	if (toReturn.sharedStateSetLookupStruct == 0L)
		goto errorExit;
	p = 0L; /* p has given the pointer to toReturn.sharedStateSetLookupStruct so we set it to NULL to avoid double deletion */
	
	
	
	/* allocate the leafData and zero leafData pointer array*/
	toReturn.nLeafData = matrix.nTax;
	toReturn.leafData = (LeafDataObj **)malloc(matrix.nTax*sizeof(LeafDataObj*));	
	if (toReturn.leafData == 0L)
		goto errorExit;
	for (i = 0; i < matrix.nTax; ++i)
		toReturn.leafData[i] = NULL;
	/* allocate each leafData object and fill its data array with the same state codes as were used in the NCL NxsCDiscreteStateSet */
	for (i = 0; i < matrix.nTax; ++i) {
		toReturn.leafData[i] = leaf_data_ss_partitioned_new(matrix.nChar, 1, toReturn.sharedStateSetLookupStruct);
		if (!(toReturn.leafData[i]))
			goto errorExit;
		int * dest = toReturn.leafData[i]->ssind;
		assert(dest);
		NxsCDiscreteStateSet * sourceRow = matrix.matrix[i];
		if (!sourceRow)
			goto errorExit;
		for (j = 0; j < matrix.nChar; ++j)
			*dest++ = (int)(*sourceRow++);
	}
	
	toReturn.nPartials = nPartials;
	toReturn.clas = (CLAObj **)malloc(nPartials*sizeof(CLAObj *));	
	if (toReturn.clas == 0L)
		goto errorExit;
	for (i = 0; i < nPartials; ++i)
		toReturn.clas[i] = NULL;
	for (i = 0; i < nPartials; ++i) {
		toReturn.clas[i] = cla_ss_partitioned_new(matrix.nChar, matrix.nStates);
		if (!(toReturn.clas[i]))
			goto errorExit;
	}
	
	toReturn.treeLike = full_la_ss_partitioned_new(matrix.nChar, matrix.nStates);
	if (toReturn.treeLike == 0L)
		goto errorExit;
	
	for (i = 0; i < nSubsets; ++i) {
		(*toReturn.treeLike).state_categ_freqs[i][matrix.nStates] = 1.0;
	}

	toReturn.asrv = 0L;

	toReturn.model = (DSCTModelObj **)malloc(nSubsets*sizeof(DSCTModelObj *));	
	if (!toReturn.model)
		goto errorExit;
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = 0L;
	}
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = dsctModelNew(matrix.nStates);
		if (!toReturn.model[i])
			goto errorExit;
	}

	toReturn.nPMat = nPMats;
	toReturn.pmats = (PMatArrayObj **)malloc(nPMats*sizeof(PMatArrayObj *));	
	if (toReturn.pmats == 0L)
		goto errorExit;
	for (i = 0; i < nPMats; ++i)
		toReturn.pmats[i] = NULL;
	for (i = 0; i < nPMats; ++i) {
		toReturn.pmats[i] = cpmatArrayNew(nSubsets, matrix.nStates);
		if (!(toReturn.pmats[i]))
			goto errorExit;
	}

	return toReturn;

	errorExit:
		PRINTF("In newLikeStructs errorExit\n");
		if ((toReturn.sharedStateSetLookupStruct == 0L) && p) {
			free(p[0]);
			free(p);
		}
		freeLikeStructFields(&toReturn);
		return toReturn;
}

LikeStructsBundle newLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nPartials, 
  unsigned nPMats,
  unsigned nRates)
{
	unsigned arr_len = 2*matrix.nStates; 
	unsigned i, j;
	double rateCatFreq;
	int ** p = 0L;
	LikeStructsBundle toReturn;
	zeroLikeStructFields(&toReturn);
	toReturn.nModels = 1;

	if (!(matrix.matrix && matrix.stateList))
		goto errorExit;
		

	/*Allocate p - a copy of the stateList*/
	/* each of the "fundamental" states needs to spaces in the array" */
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		arr_len += 1; //add one for the # of states element
		unsigned pos = matrix.stateListPos[i] ;
		arr_len += matrix.stateList[pos];
	}

	/*allocate a ragged two-D array with the memory contiguous.*/
	p = (int**)malloc(matrix.nObservedStateSets*sizeof(int*));
	if (p == 0L)
		goto errorExit;	
	p[0] = (int*)malloc(arr_len*sizeof(int));
	if (p[0] == 0L) 
		goto errorExit;

	int * curr_p = p[0];
	for (i = 0 ; i < matrix.nStates; ++i) {
		p[i] = curr_p;
		p[i][0] = 1;
		p[i][1] = i;
		curr_p += 2;
	}
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		p[i] = curr_p;
		unsigned pos = matrix.stateListPos[i];
		p[i][0] = matrix.stateList[pos];
		curr_p++;
		unsigned j;
		for (j = 1; j <= p[i][0]; ++j) {
			*curr_p++ = matrix.stateList[pos + j];
		}
	}
	
	toReturn.sharedStateSetLookupStruct = sslookup_new(matrix.nStates, matrix.nObservedStateSets, p);
	if (toReturn.sharedStateSetLookupStruct == 0L)
		goto errorExit;
	p = 0L; /* p has given the pointer to toReturn.sharedStateSetLookupStruct so we set it to NULL to avoid double deletion */
	
	
	
	/* allocate the leafData and zero leafData pointer array*/
	toReturn.nLeafData = matrix.nTax;
	toReturn.leafData = (LeafDataObj **)malloc(matrix.nTax*sizeof(LeafDataObj*));	
	if (toReturn.leafData == 0L)
		goto errorExit;
	for (i = 0; i < matrix.nTax; ++i)
		toReturn.leafData[i] = NULL;
	/* allocate each leafData object and fill its data array with the same state codes as were used in the NCL NxsCDiscreteStateSet */
	for (i = 0; i < matrix.nTax; ++i) {
		toReturn.leafData[i] = leaf_data_new(matrix.nChar, nRates, toReturn.sharedStateSetLookupStruct);
		if (!(toReturn.leafData[i]))
			goto errorExit;
		int * dest = toReturn.leafData[i]->ssind;
		assert(dest);
		NxsCDiscreteStateSet * sourceRow = matrix.matrix[i];
		if (!sourceRow)
			goto errorExit;
		for (j = 0; j < matrix.nChar; ++j)
			*dest++ = (int)(*sourceRow++);
	}
	
	toReturn.nPartials = nPartials;
	toReturn.clas = (CLAObj **)malloc(nPartials*sizeof(CLAObj *));	
	if (toReturn.clas == 0L)
		goto errorExit;
	for (i = 0; i < nPartials; ++i)
		toReturn.clas[i] = NULL;
	for (i = 0; i < nPartials; ++i) {
		toReturn.clas[i] = cla_new(matrix.nChar, matrix.nStates, nRates);
		if (!(toReturn.clas[i]))
			goto errorExit;
	}
	
	toReturn.treeLike = full_la_new(matrix.nChar, matrix.nStates, nRates);
	if (toReturn.treeLike == 0L)
		goto errorExit;
	
	rateCatFreq = 1.0/((double) nRates);
	for (i = 0; i < nRates; ++i) {
		(*toReturn.treeLike).state_categ_freqs[i][matrix.nStates] = rateCatFreq;
	}

	toReturn.asrv = asrv_obj_new(nRates, 1, 0.5);
	if (toReturn.asrv == 0L)
		goto errorExit;

	toReturn.model = (DSCTModelObj **)malloc(1*sizeof(DSCTModelObj *));	
	if (!toReturn.model)
		goto errorExit;
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = 0L;
	}
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = dsctModelNew(matrix.nStates);
		if (!toReturn.model[i])
			goto errorExit;
	}

	toReturn.nPMat = nPMats;
	toReturn.pmats = (PMatArrayObj **)malloc(nPMats*sizeof(PMatArrayObj *));	
	if (toReturn.pmats == 0L)
		goto errorExit;
	for (i = 0; i < nPMats; ++i)
		toReturn.pmats[i] = NULL;
	for (i = 0; i < nPMats; ++i) {
		toReturn.pmats[i] = cpmatArrayNew(nRates, matrix.nStates);
		if (!(toReturn.pmats[i]))
			goto errorExit;
	}

	return toReturn;

	errorExit:
		PRINTF("In newLikeStructs errorExit\n");
		if ((toReturn.sharedStateSetLookupStruct == 0L) && p) {
			free(p[0]);
			free(p);
		}
		freeLikeStructFields(&toReturn);
		return toReturn;
}



