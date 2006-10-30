/**
 * @brief filter image using CSF
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
 * $Id: csf_filter.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <math.h>
#include <iostream>

#include "csf_filter.h"
#include "CSF.h"

#include "dump_image_aspect.h"

MultiAdaptationCSF::MultiAdaptationCSF( int cols, int rows,
  const ViewingConditions &viewCond, const pfs::Array2D *otf ) 
{
  
  // TODO: Find min, max and decide what adaptation levels should be used
//   float amMin = 9999999;
//   float amMax = -9999999;
//   const int size = in->getRows()*in->getCols();
//   for( int i = 0; i < size; i++ ) {
//     const float v = (*adaptationMap)(i);
//     if( amMin > v ) amMin = v;
//     else if( amMax < v ) amMax = v;
//   }
  
//   static const float templateAdaptationLevels[] = 
//     {
//       0.0001, 0.01, 0.1, 1, 10, 100
//     };
//   for( int i = 0; i < sizeof( templateAdaptationLevels ) / sizeof( float ); i++ ) {
//     if( templateAdaptationLevels[i+1]
//   }

  static const float templateAdaptationLevels[] = 
    {
      0.0001, 0.01, 0.1, 1, 10, 100
    };
  const int templateAdaptationLevelsCount = sizeof( templateAdaptationLevels ) / sizeof( float );
  
  for( int i = 0; i < templateAdaptationLevelsCount; i++ ) 
    adaptationLevels[i] = templateAdaptationLevels[i];
  adaptationLevelsCount = templateAdaptationLevelsCount;

  filters = new pfs::Array2DImpl*[adaptationLevelsCount];
  
  // Prepare CSF filters
  for( int i = 0; i < adaptationLevelsCount; i++ ) { // For each adaptation level
    std::cerr << "Creating CSF filter: ";
    // Build CSF filter for the viewing conditions
    filters[i] = new pfs::Array2DImpl( cols/2+1, rows/2+1 );

    float meanObserverDistance = (viewCond.minDistance + viewCond.maxDistance)/2;
    float pixelsPerDeg = viewCond.getPixelsPerDegree( meanObserverDistance );
    float imgSizeVD = ((float)viewCond.xResolution / pixelsPerDeg) *
      ((float)viewCond.yResolution / pixelsPerDeg);
    createCSFFilter( filters[i], adaptationLevels[i], imgSizeVD,
      pixelsPerDeg, meanObserverDistance, CSF_DALY_NORMALIZED );

    if( otf != NULL ) {
      // Remove OTF part of CSF if OTF is given
      assert( otf->getRows() == filters[i]->getRows() &&
        otf->getCols() == filters[i]->getCols() );
      
      const int pixels = otf->getRows() * otf->getCols();
      for( int j = 0; j < pixels; j++ ) {
        (*filters[i])(j) /= (*otf)(j);
        if( (*filters[i])(j) > 1 )
          (*filters[i])(j) = 1;
      }
    }
    
  }  
}

MultiAdaptationCSF::~MultiAdaptationCSF()
{
  for( int i = 0; i < adaptationLevelsCount; i++ ) {
    delete filters[i];
  }
  delete[] filters;
}


void MultiAdaptationCSF::process( BidomainArray2D *in, BidomainArray2D *out,
  BidomainArray2D *adaptationMap )
{
  const int cols = in->getCols(), rows = in->getRows();

  assert( cols == adaptationMap->getCols() );
  assert( rows == adaptationMap->getRows() );
  
  const FFTWComplexArray *freqOriginal = in->getFrequency(); 
  FFTWComplexArray freqFiltered( cols, rows );
  FFTWArray2D spatialTemp( cols, rows );
  
  fftwf_plan inverseFFT = fftwf_plan_dft_c2r_2d( rows, cols,
    freqFiltered.getData(), spatialTemp.getData(), FFTW_ESTIMATE ); // MEASURE would damage the data

  //NOT compatible with new Cygwin version of gcc.
  //pfs::Array2DImpl **filteredImage = new (pfs::Array2DImpl*)[adaptationLevelsCount]; 

  // Results of filtering in spatial domain are stored there
  pfs::Array2DImpl **filteredImage = new pfs::Array2DImpl*[adaptationLevelsCount]; 
  
  
  for( int i = 0; i < adaptationLevelsCount; i++ ) { // For each adaptation level
    
    filterFFTW( freqOriginal->getData(), freqFiltered.getData(), cols, rows, filters[i] );
      
//    dumpPFS( "fft_image.pfs", freqFiltered, cols/2+1, rows, "Y" );

    fftwf_execute(inverseFFT);

    // Copy to filteredImage and normalize
    filteredImage[i] = new pfs::Array2DImpl( cols, rows );
    for( int pix = 0; pix < cols*rows; pix++ )
      (*filteredImage[i])(pix) = spatialTemp(pix)/(cols*rows);

//     // Some debug info
//     char buf[100];
//     sprintf( buf, "csf_filtered_%g.pfs", adaptationLevels[i] );
//     dumpPFS( buf, filteredImage[i], "Y" );

    std::cerr << ".";
    
  }

  std::cerr << "\n";

  const pfs::Array2D *adaptationMapArray = adaptationMap->getSpatial();
  
  pfs::Array2D *outA = out->setSpatial(); // output array
  // Linear intepolation between adaptation levels
  {
    int ind = 0;
    for( int ind = 0; ind < rows*cols; ind++ ) {
        float adapt = (*adaptationMapArray)( ind );
            
        if( adapt < adaptationLevels[0] )
          (*outA)(ind) = (*filteredImage[0])(ind);
        else if( adapt > adaptationLevels[adaptationLevelsCount-1] )
          (*outA)(ind) = (*filteredImage[adaptationLevelsCount-1])(ind);
        else {            // interpolate
          int l;
          for( l = 1; l < adaptationLevelsCount; l++ )
            if(adapt <= adaptationLevels[l]) break;
          assert( l > 0 && l < adaptationLevelsCount );
              
          (*outA)(ind) = (*filteredImage[l-1])(ind) +
            ((*filteredImage[l])(ind)-(*filteredImage[l-1])(ind))*
            (adapt-adaptationLevels[l-1])/(adaptationLevels[l]-adaptationLevels[l-1]);
        }
      }          
  }
//   dumpPFS( "after_csf.pfs", in, "Y" );  
  
  // Clean up
  for( int i = 0; i < adaptationLevelsCount; i++ )
    delete filteredImage[i];      
  delete[] filteredImage;

  fftwf_destroy_plan(inverseFFT);
}

float ViewingConditions::getPixelsPerDegree( float observerDistance ) const
{
  return 2*tan( 0.5 / 180. * M_PI )*observerDistance * (float)xResolution / displayWidth;  
}


//----------------------------------------------------------------
// CSF Filter used in the original Daly's VDP.
//----------------------------------------------------------------

VDPCSF::VDPCSF( int cols, int rows, const ViewingConditions &viewCond, const float Y_adapt )
{
  // Build CSF filter for the viewing conditions
  filter = new pfs::Array2DImpl( cols/2+1, rows/2+1 );

  float meanObserverDistance = (viewCond.minDistance + viewCond.maxDistance)/2;
  float pixelsPerDeg = viewCond.getPixelsPerDegree( meanObserverDistance );
  float imgSizeVD = ((float)viewCond.xResolution / pixelsPerDeg) *
    ((float)viewCond.yResolution / pixelsPerDeg);
  createCSFFilter( filter, Y_adapt, imgSizeVD,
    pixelsPerDeg, meanObserverDistance, CSF_DALY );

  dumpImageAspect->dump( "filter_csf.pfs", filter, "Y" );
  
}

VDPCSF::~VDPCSF()
{
  delete filter;
}

void VDPCSF::process( BidomainArray2D *in, BidomainArray2D *out, BidomainArray2D *adaptationMap )
{
  const int cols = in->getCols(), rows = in->getRows();
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency()->getData() if (in == out)
  filterFFTW( input->getData(), out->setFrequency()->getData(), cols, rows, filter );
}
