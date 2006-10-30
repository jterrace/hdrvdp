/**
 * @brief operations on 2D arrays
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
 * $Id: array2d_algorithm.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <math.h>

#include "array2d_algorithm.h"

float average( const pfs::Array2D *array )
{
  float sum = 0;
  const int size = array->getCols()*array->getRows();
  for( int i = 0; i < size; i++ ) {
    sum += (*array)(i);
  }
  return sum/(float)size;
}

float logAverage( const pfs::Array2D *array )
{
  float sum = 0;
  const int size = array->getCols()*array->getRows();
  for( int i = 0; i < size; i++ ) {
    assert( (*array)(i) > 0 );
    sum += logf((*array)(i));
  }
  return expf( sum/(float)size );
}

float max(  const pfs::Array2D *array )
{
  float maxV = -999999999;
  const int size = array->getCols()*array->getRows();
  for( int i = 0; i < size; i++ ) {
    if( (*array)(i) > maxV ) maxV = (*array)(i);
  }
  return maxV;
}
  
void normalizeToSum( pfs::Array2D *array )
{
  float sum = 0;
  const int size = array->getCols()*array->getRows();
  for( int i = 0; i < size; i++ ) {
    sum += (*array)(i);
  }
  for( int i = 0; i < size; i++ ) {
    (*array)(i) /= sum;
  }
}

void normalizeToMax( pfs::Array2D *array )
{
  float max = -99999;
  const int size = array->getCols()*array->getRows();
  for( int i = 0; i < size; i++ ) {
    if( max < (*array)(i) ) max = (*array)(i);
  }
  for( int i = 0; i < size; i++ ) {
    (*array)(i) /= max;
  }
}

/**
 * Clamp values in the array to positive integers.  Negative, zero,
 * NaN and -Inf are set to the smallest positive value. +Inf is set to
 * the largest value, which is not +Inf.
 */
void clampToPositive( const pfs::Array2D *from, pfs::Array2D *to )
{
  float minPositive = 1e12;
  float maxNonInf = -1;  
  const int size = from->getCols()*from->getRows();
  for( int i = 0; i < size; i++ ) {
    if( !isfinite((*from)(i)) )
      continue;
    if( (*from)(i) > 0 && (*from)(i) < minPositive )
      minPositive = (*from)(i);
    if( (*from)(i) > maxNonInf )
      maxNonInf = (*from)(i);
  }

  for( int i = 0; i < size; i++ ) {
    if( isinf( (*from)(i) ) == -1 || isnanf( (*from)(i) ) || (*from)(i) <= 0 )
      (*to)(i) = minPositive;
    else if( isinf( (*from)(i) ) == 1 )
      (*to)(i) = maxNonInf;
    else
      (*to)(i) = (*from)(i);
  }
  
}
