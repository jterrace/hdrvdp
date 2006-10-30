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
 * $Id: fftutils.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef FFTUTILS_H
#define FFTUTILS_H

#include <fftw3.h>
#include <array2d.h>

#include "fftw_array2d.h"

/**
 * Filter 2d array of fftwf_complex with filter stored in Array2D object.
 * Note that the filter must be of the size [cols/2+1; rows/2+1].
 * The input may be the same as the output.
 */
void filterFFTW( const fftwf_complex *input, fftwf_complex *output, int cols, int rows, const pfs::Array2D *filter );

void normalizeFFTW( pfs::Array2D *array );

/**
 * Perform multiplication: z = x * y. z can be the same as x or
 * y. Multiplication is done using complex number arithmetic.
 */
void multiplyArray( FFTWComplexArray *z, const FFTWComplexArray *x,
  const FFTWComplexArray *y );

/* /\** */
/*  * Sets all elements of the array to a complex number (re, im). */
/*  *\/ */
/* inline void setArray(FFTWComplexArray *array, const float re, const float im ) */
/* { */
/*   const int elements = array->getFFTRows()*array->getFFTCols(); */
/*   for( int i = 0; i < elements; i++ ) { */
/*     (*array)(i)[0] = re; */
/*     (*array)(i)[1] = im; */
/*   } */
/* } */
  


/**
 * This class can store 2D arrays (images) in both spatial and
 * frequency domain. It can then transparently convert from one
 * representation to the other, depending on what particular
 * application requires.
 */
class BidomainArray2D
{
  int cols, rows;
  bool spatialValid;
  bool frequencyValid;

  FFTWArray2D *spatial;
  FFTWComplexArray *frequency;
  
public:
  BidomainArray2D( int cols, int rows );
  ~BidomainArray2D();
  
  const FFTWArray2D* getSpatial();
  const FFTWComplexArray* getFrequency();

  inline int getCols() const { return cols; }
  inline int getRows() const { return rows; }

  FFTWArray2D* setSpatial();
  FFTWComplexArray* setFrequency();

  BidomainArray2D& operator=(BidomainArray2D &bda);
  
};




#endif

