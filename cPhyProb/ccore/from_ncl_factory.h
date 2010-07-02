//	Copyright (C) 2009 Mark T. Holder
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

#if !defined(FROM_NCL_FACTORY_H)
#define FROM_NCL_FACTORY_H


#include "dsct_model.h"
#include "ncl/nxscdiscretematrix.h"


#ifdef __cplusplus
extern "C" 
{
#endif

/* 
	Assumes that the NxsCDiscreteMatrix outlives the LikeStructsBundle (some
	pointers may be aliases.
*/
LikeStructsBundle newLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nCLAs, 
  unsigned nPMats, 
  unsigned nRates
  );

LikeStructsBundle newPartitionedLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nPartials, 
  unsigned nPMats,
  unsigned nSubsets);


#ifdef __cplusplus
}
#endif

#endif
