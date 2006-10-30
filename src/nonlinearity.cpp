/**
 * @brief photoreceptor nonlinearity
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
 * $Id: nonlinearity.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <math.h>

#include "nonlinearity.h"
#include "csf_psi_nonlinearity.h"
#include "array2d_algorithm.h"

void JNDScaledNonlinearity::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const pfs::Array2D *aIn = in->getSpatial();
  pfs::Array2D *aOut = out->setSpatial();
  const int size = aIn->getRows()*aIn->getCols();
  
  for( int i = 0; i < size; i++ ) {
    (*aOut)(i) = getPsiValue( (*aIn)(i) ) / (peakSensitivity / 0.01f);
  }
  
}


void VDPNonlinearity::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const pfs::Array2D *aIn = in->getSpatial();
  pfs::Array2D *aOut = out->setSpatial();
  const int size = aIn->getRows()*aIn->getCols();

  // Renormalize to get back from normalized response to luminance
  // The paper does not say clearly if such step is required, but it
  // was done like this in Karol's implementation
  const float displayMaxLum = 100; // 100 cd/m^2

  // This variable needs to be calibrated to give reasonable results
  float renormalize = displayMaxLum;

  // Several alternatices:
  //float maxY = max( aIn );
  //float renormalize = maxY / (maxY / (maxY + powf( 12.6f*maxY, 0.63f )));
  //float renormalize = displayMaxLum/logAverage(aIn)*(displayMaxLum + powf( 12.6f*displayMaxLum, 0.63f ));
  
  for( int i = 0; i < size; i++ ) {
    const float y = (*aIn)(i);
    (*aOut)(i) = renormalize * y / (y + 12.6f*powf( y, 0.63f ));    
  }
  
}

