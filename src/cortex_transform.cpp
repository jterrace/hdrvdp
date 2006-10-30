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
 * $Id: cortex_transform.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <math.h>
#include <iostream>

#include "cortex_transform.h"
#include "dump_image_aspect.h"
#include "array2d_algorithm.h"

#include <assert.h>

inline float max( float a, float b )
{
  return a > b ? a : b;
}

inline float min( float a, float b )
{
  return a < b ? a : b;
}

void buildDomAndMesaFilter( pfs::Array2D *domFilter, pfs::Array2D *mesaFilter,
  int band );
void buildFanFilter( pfs::Array2D *filter, int orientation );

void computeThresholdElevation( const pfs::Array2D *contrast, pfs::Array2D *thesholdEl,
  float maskingSlope );
void computeProbabilityOfDetection( const pfs::Array2D *contrast, const pfs::Array2D *threshold,
  pfs::Array2D *probabilityMap, const float psychometricFunctionSlope );

DetectionMechanism::DetectionMechanism( float psychometricFunctionSlope ):
  psychometricFunctionSlope( psychometricFunctionSlope )
{
}  

// ================= Watson decomposition DM ===================

WatsonDM::WatsonDM( float psychometricFunctionSlope, float maskingSlope,
  bool normalizeContrast, float normalizationFactor, bool phaseUncertainty ) :
  DetectionMechanism( psychometricFunctionSlope ),
  maskingSlope( maskingSlope ),
  normalizeContrast( normalizeContrast ),
  normalizationFactor( normalizationFactor ),
  phaseUncertainty( phaseUncertainty )
{
}

void WatsonDM::process( BidomainArray2D *target, BidomainArray2D *mask, pfs::Array2D *probabilityMap )
{
  const int cols = target->getCols(), rows = target->getRows();

  assert( mask->getCols() == cols );
  assert( mask->getRows() == rows );  

  FFTWComplexArray freqFiltered( cols, rows );

  FFTWArray2D targetFiltered( cols, rows );
  FFTWArray2D maskFiltered( cols, rows );  

  fftwf_plan targetInverseFFT = fftwf_plan_dft_c2r_2d( rows, cols,
    freqFiltered.getData(), targetFiltered.getData(), FFTW_MEASURE );      
  fftwf_plan maskInverseFFT = fftwf_plan_dft_c2r_2d( rows, cols,
    freqFiltered.getData(), maskFiltered.getData(), FFTW_MEASURE );
  

  pfs::Array2DImpl channelProbability( cols, rows );
  pfs::Array2DImpl probabilityMapSign( cols, rows );

  pfs::setArray( &probabilityMapSign, 0 ); // probabilityMapSign(x,y) = 0;
  pfs::setArray( probabilityMap, 1 ); // probabilityMap(x,y) = 1;

  for( int band = 6; band >= 1; band-- ) {
    std::cerr << "Processing band " << band << " ";

    pfs::Array2DImpl mesaFilter( cols/2+1, rows );
    pfs::Array2DImpl domFilter( cols/2+1, rows );
    buildDomAndMesaFilter( &domFilter, &mesaFilter, band );
    
    dumpImageAspect->dump( "dom_filter_cortex_%d.pfs", &domFilter, "Y", band );
    dumpImageAspect->dump( "mesa_filter_cortex_%d.pfs", &mesaFilter, "Y", band );
    
    for( int orientation = 1; orientation <= (band == 6 ? 1 : 6); orientation++ ) {
      // only one orientation for low pass
      
      std::cerr << ".";

      pfs::Array2DImpl fanFilter( cols/2+1, rows ); // Include negative freq in the filter
      pfs::Array2D *filter;
      if( band == 6 )           // Base band
        filter = &domFilter;
      else {        
        buildFanFilter( &fanFilter, orientation );
        pfs::multiplyArray( &fanFilter, &fanFilter, &domFilter );
        filter = &fanFilter;
      }      

      filterFFTW( target->getFrequency()->getData(), freqFiltered.getData(),
        cols, rows, filter );
      fftwf_execute( targetInverseFFT ); // targetFiltered = IFFT( freqFiltered )
      normalizeFFTW( &targetFiltered );

      filterFFTW( mask->getFrequency()->getData(), freqFiltered.getData(),
        cols, rows, filter );
      fftwf_execute( maskInverseFFT ); // maskFiltered = IFFT( freqFiltered )
      normalizeFFTW( &maskFiltered );

      // in Daly's VDP contrast should be normalized by mean of the base band
      // This however does not give stable solution (average can be even <0)
      // Therefore log mean luminance is used
//       if( normalizeByMeanBase && band == 6 ) {
//         baseMean = (average( &targetFiltered ) + average( &maskFiltered ))/2.f;
//         fprintf( stderr, "Base mean: %g\n", baseMean );
//       }

      // Dump debug images
      dumpImageAspect->dump( "mask_cortex_%d_%d.pfs", &maskFiltered, "Y", band, orientation );
      dumpImageAspect->dump( "target_cortex_%d_%d.pfs", &targetFiltered, "Y", band, orientation );
      dumpImageAspect->dump( "filter_cortex_%d_%d.pfs", filter, "Y", band, orientation );

      // Compute difference of mask and target (contrast) 
      pfs::Array2DImpl diff( cols, rows );
      const int elements = maskFiltered.getRows()*maskFiltered.getCols();
      for( int i = 0; i < elements; i++ ) {
        float maskV = maskFiltered(i), targetV = targetFiltered(i);
        if( normalizeContrast ) {
          diff(i) = (targetV - maskV)/normalizationFactor;
        } else {
          diff(i) = targetV - maskV;
        }        
      }
      dumpImageAspect->dump( "difference_cortex_%d_%d.pfs", &diff, "Y", band, orientation );
      
      // Compute threshold elevation (mutual masking and phase uncertainty)
      pfs::Array2D *thresholdElevation = &maskFiltered; // Memory saving
      {
        BidomainArray2D maskTE( cols, rows ), targetTE( cols, rows );
        // In Daly's paper masking slope ranged from 0.7 for the base band
        // to 1.0 for middle frequencies (data taken from the patent)
        float usedMaskingSlope;
        if( maskingSlope != -1 ) usedMaskingSlope = maskingSlope;
        else {
          if( band == 6 ) usedMaskingSlope = 0.7; // base band
          else if( band == 5 || band == 4 ) usedMaskingSlope = 0.8;
          else usedMaskingSlope = 1;
        }
        computeThresholdElevation( &targetFiltered, targetTE.setSpatial(),
          usedMaskingSlope );
        computeThresholdElevation( &maskFiltered, maskTE.setSpatial(),
          usedMaskingSlope );
        //Implementation of "phase uncertainty" - described only in the patent
        if( phaseUncertainty ) {
          const fftwf_complex *source = targetTE.getFrequency()->getData();
          filterFFTW( source, targetTE.setFrequency()->getData(),
            cols, rows, &mesaFilter );
          source = maskTE.getFrequency()->getData();
          filterFFTW( source, maskTE.setFrequency()->getData(),
            cols, rows, &mesaFilter );
        }
        // Find minimum (mutual masking)
        const pfs::Array2D *targetTEA = targetTE.getSpatial();
        const pfs::Array2D *maskTEA = maskTE.getSpatial();
        for( int i = 0; i < elements; i++ ) {
          (*thresholdElevation)(i) = max( 1, min( (*targetTEA)(i), (*maskTEA)(i) ) );
//           (*thresholdElevation)(i) = max( 1, min( fabs((*targetTEA)(i)), fabs((*maskTEA)(i)) ) );
        }
      }      
      dumpImageAspect->dump( "te_cortex_%d_%d.pfs", thresholdElevation, "Y", band, orientation );      
      
      // Compute probability of detection (psychometric function)
      computeProbabilityOfDetection( &diff, thresholdElevation, &channelProbability, psychometricFunctionSlope );
      
      dumpImageAspect->dump( "probability_cortex_%d_%d.pfs", &channelProbability, "Y", band, orientation );
      
      {
        // Probability summation
        const int elements = probabilityMap->getRows()*probabilityMap->getCols();
        for( int i = 0; i < elements; i++ ) {
          (*probabilityMap)(i) *= (1.0f - fabs(channelProbability(i)));
          probabilityMapSign(i) += fabs(channelProbability(i));
	}
      }
      
    }
    std::cerr << "\n";    
  }

  fftwf_destroy_plan(targetInverseFFT);
  fftwf_destroy_plan(maskInverseFFT);
  
  {
    // The final step of probability summation
    // TODO: include sign
    const int elements = probabilityMap->getRows()*probabilityMap->getCols();
    for( int i = 0; i < elements; i++ )
      (*probabilityMap)(i) = 1.f - (*probabilityMap)(i);
  }
  
}

// ================= Single Channel DM ===================

SingleChannelDM::SingleChannelDM( float psychometricFunctionSlope ) : DetectionMechanism( psychometricFunctionSlope )
{
}


void SingleChannelDM::process( BidomainArray2D *target, BidomainArray2D *mask, pfs::Array2D *probabilityMap )
{
  const int cols = target->getCols(), rows = target->getRows();

  assert( mask->getCols() == cols );
  assert( mask->getRows() == rows );  

  const FFTWArray2D *targetSpatial = target->getSpatial();
  const FFTWArray2D *maskSpatial = mask->getSpatial();
  
  pfs::Array2DImpl diff( cols, rows );
  const int elements = rows*cols;
  for( int i = 0; i < elements; i++ )
    diff(i) = (*targetSpatial)(i) - (*maskSpatial)(i);

  pfs::Array2DImpl threshold( cols, rows );
  pfs::setArray( &threshold, 1 );
   
  computeProbabilityOfDetection( &diff, &threshold, probabilityMap, psychometricFunctionSlope );
}


void computeThresholdElevation( const pfs::Array2D *contrast, pfs::Array2D *thresholdEl,
  float maskingSlope )
{
  // Masking variables
  const float W = 6.0;
  const float Q = 0.7;
  const float b = 4.0;
  const float s = maskingSlope;  // valid range [0.65-1.0] 
  const float k1 = powf(W, 1.0 - 1.0 / (1.0 - Q));
  const float k2 = powf(W, 1.0 / (1.0 - Q));
  
  const int elements = contrast->getRows()*contrast->getCols();
  for( int i = 0; i < elements; i++ ) {
//    float sign = (*contrast)(i) < 0 ? -1 : 1;
    (*thresholdEl)(i) = powf(1.0 + powf(k1 * powf(k2 * fabs( (*contrast)(i) ), s), b), 1.0/b);
  }  
}

void computeProbabilityOfDetection( const pfs::Array2D *contrast, const pfs::Array2D *threshold,
  pfs::Array2D *probabilityMap, const float psychometricFunctionSlope )
{

  const int elements = contrast->getRows()*contrast->getCols();
  for( int i = 0; i < elements; i++ ) {
    float sign = (*contrast)(i) < 0 ? -1 : 1;
    (*probabilityMap)(i) = sign *
      (1.0 - expf( -pow( fabs((*contrast)(i))/(*threshold)(i), psychometricFunctionSlope ) ));
  }
}


#define FILTER_TRASITION_FREQENCY    0.66666666
#define EPSILON          0.00001

void buildDomAndMesaFilter( pfs::Array2D *domFilter, pfs::Array2D *mesaFilter,
  int band ) {

  assert( domFilter->getRows() == mesaFilter->getRows() );
  assert( domFilter->getCols() == mesaFilter->getCols() );
  
  int filterWidth=domFilter->getCols();
  int filterHeight=domFilter->getRows() / 2 + 1;

  // === Build mesa filter ===
  
  float ro2[7], tw[7];

  ro2[0] = 1.0;
  tw[0] = ro2[0] * FILTER_TRASITION_FREQENCY;

  const int pyramid_height = 6;
  
  for( int k = 1; k < pyramid_height; k ++)
  {
    ro2[k] = ro2[k-1] / 2.0;
    tw[k] = ro2[k] * FILTER_TRASITION_FREQENCY;
  }
  ro2[pyramid_height] = ro2[pyramid_height-1];
  tw[pyramid_height] = ro2[pyramid_height] * FILTER_TRASITION_FREQENCY;
    
  // Mesa filters parameters
  float ro1_min = ro2[band-1] - 0.5 * tw[band-1];
  float ro1_max = ro2[band-1] + 0.5 * tw[band-1];
  float arg1 = M_PI * (0.5 * tw[band-1] - ro2[band-1]) / tw[band-1];
  float fct1 = M_PI / tw[band-1];

  float ro2_min = ro2[band] - 0.5 * tw[band];
  float ro2_max = ro2[band] + 0.5 * tw[band];
  float arg2 = M_PI * (0.5 * tw[band] - ro2[band]) / tw[band];
  float fct2 = M_PI / tw[band];


  // Base filter parameters
  int k_base = pyramid_height - 1;
  float sigma = (ro2[k_base] + 0.5 * tw[k_base]) / 3.0;
  sigma *= 2.0 * sigma;
  float ro_min = ro2[k_base] + 0.5 * tw[k_base];
/* 
   printf("ro1[%f,%f], ro2[%f,%f], ro_base[%f]\n",ro1_min,ro1_max,ro2_min,ro2_max,ro_min);
*/
  ro_min *= ro_min;

  int x, y;
  float x_norm = 0.5 / filterWidth, y_norm = 0.5 / filterHeight;
  float dx, dy;
  float ro;
  float mesa, mesa1, mesa2;
  for( y = 0, dy = 0.0; y < filterHeight; y++, dy += y_norm ) {
    for( x = 0, dx = 0.0; x < filterWidth; x++, dx += x_norm ) {
//         if (x > y)  // use symmetry to speed up filter calculations
//           continue;
      
      ro = (float)(dx * dx + dy * dy);
      if (band >= k_base)       // Base band
      {
        if (ro < ro_min)
        {
          mesa2 = (float)exp(-ro / sigma);
        }
        else
        {
          mesa2 = 0.0;
        }
        if (band == pyramid_height)   // Base filter only
        {
          (*domFilter)(x,y) = (*mesaFilter)(x,y) = (*domFilter)(x,y) = mesa2;
          if( y != 0 ) (*domFilter)(x,domFilter->getRows()-y) = (*mesaFilter)(x,mesaFilter->getRows()-y) = mesa2;
          continue;
        }
      }

      ro = sqrtf(ro);
      if (ro < ro1_min)
      {
        mesa1 = 1.0;
      }
      else if (ro > ro1_max)
      {
        mesa1 = 0.0;
      }
      else
      {
        mesa1 = 0.5 * (1.0 + cosf(fct1 * ro + arg1));
      }

      if (band < k_base)
      {
        if (ro < ro2_min)
        {
          mesa2 = 1.0;
        }
        else if (ro > ro2_max)
        {
          mesa2 = 0.0;
        }
        else
        {
          mesa2 = 0.5 * (1.0 + cosf(fct2 * ro + arg2));
        }
      }

//       if (just_low_pass)
// 	mesa = mesa2;
//       else
      mesa = mesa1 - mesa2;


      if (mesa < 0.0)
      {
        if (mesa < -EPSILON)
          fprintf( stderr, "DOM filter negative: %f = %f - %f\n", mesa, mesa1, mesa2);

        mesa = 0.0;
      }

      (*domFilter)(x,y) = mesa;
      (*mesaFilter)(x,y) = mesa1;
      if( y != 0 ) {
        (*domFilter)(x,domFilter->getRows()-y) = mesa;
        (*mesaFilter)(x,mesaFilter->getRows()-y) = mesa1;
      }
      
    }
  }
}

void buildFanFilter( pfs::Array2D *filter, int orientation )
{ 
  
  float theta_tw = M_PI / 6;
  float theta_l = (orientation - 1) * theta_tw - M_PI / 2.0;
  int filterWidth=filter->getCols();
  int filterHeight=filter->getRows() / 2 + 1;
  float x_norm = 0.5 / filterWidth, y_norm = 0.5 / filterHeight;

  // Positive freq
  int x, y;
  float dx, dy;
  for( y = 0, dy = 0.0; y < filterHeight; y++, dy += y_norm ) {
    for( x = 0, dx = 0.0; x < filterWidth; x++, dx += x_norm ) {

      if( x == 0 && y == 0 ) {
        (*filter)(x,y) = 0.0;
        continue;
      }
      
      float theta = atan2f( dy, dx );
      float delta_theta = fabs( theta + theta_l );
      
      if( delta_theta <= theta_tw ) {
	(*filter)(x,y) = 0.5 * (1.0 + cosf(M_PI * delta_theta / theta_tw));
      } else
	(*filter)(x,y) = 0.0;      
    }
  }

  // Negative freq
  for( y = filter->getRows()-1, dy = 0.0; y >= filterHeight; y--, dy += y_norm ) {
    for( x = 0, dx = 0.0; x < filterWidth; x++, dx += x_norm ) {

      if( x == 0 ) {
        (*filter)(x,y) = 0.0;
        continue;
      }
      
      float theta = atan2f( -dy, dx );
      float delta_theta = fabs( theta + theta_l );
      if( delta_theta >= M_PI/2 ) delta_theta = M_PI - delta_theta;
      
      if( delta_theta <= theta_tw ) {
	(*filter)(x,y) = 0.5 * (1.0 + cosf(M_PI * delta_theta / theta_tw));
      } else
	(*filter)(x,y) = 0.0;
      
      
    }
  }
  
}
