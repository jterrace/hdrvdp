/**
 * @brief convert luminance to JND-scaled space
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
 * $Id: csf_psi_nonlinearity.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <assert.h>
#include <math.h>

float csf_psi_mapping[] = {
#include "csf_psi_mapping.h"
};

static const int lum_lookup_size = sizeof( csf_psi_mapping ) / sizeof( float );

/**
 * Converts luiminance y in cd/m^2 to psi mapping.
 */
float getPsiValue( float y )
{
  // 5 units in psi space equals 1JDN under the peak sensitivity of 1%
  const float f = 5;       // To normalize values to 1JND (1% peak) 

  //Clamp if beyond boundaries
  if( y < csf_psi_mapping[0] ) return 0;
  if( y > csf_psi_mapping[lum_lookup_size-1] ) return (lum_lookup_size-1)/f;
  
  //binary search for best fitting luminance
  int l, r;
  l = 0;
  r = lum_lookup_size-1;
  while( l+1 < r ) {
    int m = (l+r)/2;
    if( y < csf_psi_mapping[m] ) r = m;
    else l = m;
  }
  assert( r - l == 1 );
  assert( (y-csf_psi_mapping[l])/(csf_psi_mapping[r]-csf_psi_mapping[l]) <= 1 );
  
    
  //Interpolate to get more acurate result  
  return ((float)l + (y-csf_psi_mapping[l])/(csf_psi_mapping[r]-csf_psi_mapping[l]) )/f;
  
}

/*
void applyPsiAmplitudeNonLinearity( float **image, int width, int height )
{ 
  for( int y = 1; y <= height; y++ )
    for( int x = 1; x <= width; x++ ) {
      image[y][x] = getPsiValue( image[y][x] );      
    }
}

*/
