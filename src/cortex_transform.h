/**
 * @brief decomposition into channels, visual masking, phase uncertainty
 * 
 * This file is a part of HDR VDP.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2003,2004 Rafal Mantiuk and Grzegorz Krawczyk
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ---------------------------------------------------------------------- 
 *
 * @author Rafal Mantiuk, <mantiuk@mpi-inf.mpg.de>
 *
 * $Id: cortex_transform.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef CORTEX_TRANSFORM_H
#define CORTEX_TRANSFORM_H

#include "fftw_array2d.h"
#include "fftutils.h"

class DetectionMechanism
{
protected:
  float psychometricFunctionSlope;
public:
  DetectionMechanism( float psychometricFunctionSlope ); 
  
  virtual void process( BidomainArray2D *target, BidomainArray2D *mask, pfs::Array2D *probabilityMap ) = 0;
  
};

class WatsonDM: public DetectionMechanism
{
protected:
  float maskingSlope;
  float normalizationFactor;
  bool normalizeContrast;
  bool phaseUncertainty;
public:
  WatsonDM( float psychometricFunctionSlope, float maskingSlope,
    bool normalizeContrast, float normalizationFactor, bool phaseUncertainty );
  
  void process( BidomainArray2D *target, BidomainArray2D *mask, pfs::Array2D *probabilityMap );
  
};

class SingleChannelDM: public DetectionMechanism
{
public:
  SingleChannelDM( float psychometricFunctionSlope );
  
  void process( BidomainArray2D *target, BidomainArray2D *mask, pfs::Array2D *probabilityMap );
  
};


#endif

