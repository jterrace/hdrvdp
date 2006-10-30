/**
 * @brief Array2D subclass adjusted to work with fftw library
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
 * $Id: fftw_array2d.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef FFTW_ARRAY2D_H
#define FFTW_ARRAY2D_H

#include <array2d.h>
#include <fftw3.h>

#include <assert.h>

/**
 * Implementation of Array2D that uses fftw library to ensure proper
 * data alignment.  FFT on such array should be faster than on
 * a non-aligned array.
 */
class FFTWArray2D : public pfs::Array2D
{
  float *data;
  int cols, rows;
    
public:

  FFTWArray2D( int cols, int rows ) : cols( cols ), rows( rows )
    {
      data = (float*)fftwf_malloc(sizeof(float) * cols * rows );
      assert( data != NULL );
    }
    
  ~FFTWArray2D()
    {
      fftwf_free( data );
    }

  inline int getCols() const { return cols; }
  inline int getRows() const { return rows; }

  inline float& operator()( int x, int y ) {
    assert( x < cols );
    assert( y < rows );
    return data[ x+y*cols ];
  }
  inline const float& operator()( int x, int y ) const { return data[ x+y*cols ]; }

  inline float& operator()( int rowMajorIndex ) { return data[rowMajorIndex]; }
  inline const float& operator()( int rowMajorIndex ) const { return data[rowMajorIndex]; }

  inline float *getData() 
    {
      return data;
    }
};


/**
 * Fourier representation of a 2D image of real values
 */
class FFTWComplexArray
{
  fftwf_complex *data;
  int cols, rows;
    
public:

  FFTWComplexArray( int cols, int rows ) : cols( cols ), rows( rows )
    {
      data = (fftwf_complex*)fftwf_malloc( sizeof(fftwf_complex) * (cols/2+1) * rows );
      assert( data != NULL );
    }
    
  ~FFTWComplexArray()
    {
      fftwf_free( data );
    }

  inline int getCols() const { return cols; }
  inline int getRows() const { return rows; }

  inline int getFFTCols() const { return cols/2+1; }
  inline int getFFTRows() const { return rows; }
  
  inline fftwf_complex *getData() 
    {
      return data;
    }  

  inline const fftwf_complex *getData() const
    {
      return data;
    }  

};


/**
 * Copy data from one Array2D to another. Dimensions of the arrays must be the same.
 */
inline void copyArray(const FFTWComplexArray *from, FFTWComplexArray *to)
{
  assert( from->getRows() == to->getRows() );
  assert( from->getCols() == to->getCols() );
  
  const int elements = from->getFFTRows()*from->getFFTCols();
  const fftwf_complex *fromA = from->getData();
  fftwf_complex *toA = to->getData();
  for( int i = 0; i < elements; i++ ) {
    toA[i][0] = fromA[i][0];
    toA[i][1] = fromA[i][1];
  }
  
}


#endif
