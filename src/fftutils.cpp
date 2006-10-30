/**
 * @brief wrapper for FFT
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
 * $Id: fftutils.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include "fftutils.h"

#include <string.h>
#include <assert.h>


BidomainArray2D::BidomainArray2D( int cols, int rows ) : cols( cols ), rows( rows )
{
  spatialValid = frequencyValid = false;
  spatial = NULL;
  frequency = NULL;
}

BidomainArray2D& BidomainArray2D::operator=(BidomainArray2D &bda)
{
  assert( rows == bda.rows && cols == bda.cols );
  spatialValid = bda.spatialValid;
  if( spatialValid ) {
    if( spatial == NULL )
      spatial = new FFTWArray2D( cols, rows );
    memcpy( spatial->getData(), bda.spatial->getData(),
      sizeof( float ) * cols*rows );
  }

  frequencyValid = bda.frequencyValid;
  if( frequencyValid ) {
    if( frequency == NULL ) 
      frequency = new FFTWComplexArray( cols, rows );
    memcpy( frequency->getData(), bda.frequency->getData(),
      sizeof( fftwf_complex ) * frequency->getFFTCols()*frequency->getFFTRows() );
  }
  
  return *this;
}


BidomainArray2D::~BidomainArray2D()
{
  delete frequency;
  delete spatial;
}

const FFTWArray2D* BidomainArray2D::getSpatial()
{
  assert( spatialValid || frequencyValid );

  if( spatialValid ) {
    assert( spatial != NULL );
    return spatial;
  }

  // IFFT required
  if( spatial == NULL ) 
    spatial = new FFTWArray2D( cols, rows );

  //Potential pitfall: c2r tranformation damages input array,
  //so it is necessary to make a copy (or invalidate it)
  FFTWComplexArray temp( cols, rows );

  // FFTW_MEASURE could return different FFT routine for each
  // run. This in turn could give slightly different results for
  // different runs.
  fftwf_plan inverseFFT = fftwf_plan_dft_c2r_2d( rows, cols,
    temp.getData(), spatial->getData(), FFTW_ESTIMATE );
  
  copyArray( frequency, &temp );
  
  fftwf_execute(inverseFFT);
  fftwf_destroy_plan(inverseFFT);
  normalizeFFTW( spatial );
  spatialValid = true;

  return spatial;
}

const FFTWComplexArray* BidomainArray2D::getFrequency()
{
  assert( spatialValid || frequencyValid );

  if( frequencyValid ) {
    assert( frequency != NULL );
    return frequency;
  }

  // FFT required
  if( frequency == NULL ) 
    frequency = new FFTWComplexArray( cols, rows );

  fftwf_plan forwardFFT = fftwf_plan_dft_r2c_2d( rows, cols,
    spatial->getData(), frequency->getData(), FFTW_ESTIMATE ); // MEASURE would damage the data
  
  fftwf_execute(forwardFFT);
  fftwf_destroy_plan(forwardFFT);
  frequencyValid = true;

  return frequency;

}

FFTWArray2D* BidomainArray2D::setSpatial()
{
  if( spatial == NULL ) 
    spatial = new FFTWArray2D( cols, rows );

  spatialValid = true;
  frequencyValid = false;
  
  return spatial;
}
  
FFTWComplexArray* BidomainArray2D::setFrequency()
{
  if( frequency == NULL ) 
    frequency = new FFTWComplexArray( cols, rows );

  spatialValid = false;
  frequencyValid = true;
  
  return frequency;
}

//--------------------------------------------------
// Filtering
//--------------------------------------------------

void filterFFTW( const fftwf_complex *input, fftwf_complex *output, int cols, int rows, const pfs::Array2D *filter )
{
  const int fftCols = filter->getCols();
  const int fftRows = filter->getRows();

  assert( (fftCols-1)*2 == (cols & ~1) );
  assert( (fftRows-1)*2 == (rows & ~1) || fftRows == rows );  

  bool negativeFrequency = (fftRows == rows);

  if( negativeFrequency ) {
    
    // A simple multiplication
    for( int x = 0; x < fftCols; x++ )
      for( int y = 0; y < fftRows; y++ ) { 
        const float filterV = (*filter)( x, y );
        const int ind1 = x + y*fftCols;
        output[ind1][0] = input[ind1][0]*filterV;
        output[ind1][1] = input[ind1][1]*filterV;
      }
    
  } else {

    // Filter fy=0 & DC
    for( int x = 0; x < fftCols; x++ ) {
      const float filterV = (*filter)( x, 0 );
      output[x][0] = input[x][0]*filterV;
      output[x][1] = input[x][1]*filterV;
    }
      
    // Filter fx=0:1/N-1, fy=1:1/N-1
    for( int x = 0; x < fftCols; x++ )
      for( int y = 1; y < fftRows-!(rows&1); y++ ) { // one less if # rows even
        const float filterV = (*filter)( x, y );
        const int ind1 = x + y*fftCols, ind2 = x + (rows-y)*fftCols;
        output[ind1][0] = input[ind1][0]*filterV;
        output[ind1][1] = input[ind1][1]*filterV;
        output[ind2][0] = input[ind2][0]*filterV;
        output[ind2][1] = input[ind2][1]*filterV;
      }

    // Filter fy=fc=1/2 (Nyquist limit)
    if( !(rows & 1) )             // If # rows even
      for( int x = 0; x < fftCols; x++ ) {
        const int y = fftRows-1;
        const float filterV = (*filter)( x, y );
        const int ind = x+y*fftCols;
        output[ind][0] = input[ind][0]*filterV;
        output[ind][1] = input[ind][1]*filterV;
      }
  
  }
  
}

void normalizeFFTW( pfs::Array2D *array )
{
  const int size = array->getRows()*array->getCols();
  float scale = (float)size;
  for( int i = 0; i < size; i++ ) {
    (*array)(i) /= scale;
  }  
}

/**
 * Perform multiplication: z = x * y. z can be the same as x or y.
 */
void multiplyArray( FFTWComplexArray *z, const FFTWComplexArray *x,
  const FFTWComplexArray *y )
{
  const fftwf_complex *ax = x->getData();
  const fftwf_complex *ay = y->getData();
  fftwf_complex *az = z->getData();

  const int size = x->getFFTCols()*x->getFFTRows();
  for( int i = 0; i < size; i++ ) {
    az[i][0] = ax[i][0]*ay[i][0] - ax[i][1]*ay[i][1];
    az[i][1] = ax[i][1]*ay[i][0] + ax[i][0]*ay[i][1];
  }
  
}
