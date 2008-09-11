/**
 * @brief Optical Tranfer Function for predicting ocular light scatter
 * and veiling glare
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
 * $Id: otf.cpp,v 1.6 2008/09/11 22:20:13 rafm Exp $
 */

#include <iostream>
#include <math.h>

#include "otf.h"
#include "array2d_algorithm.h"

#include "dump_image_aspect.h"

inline float max( float a, float b )
{
  return a > b ? a : b;
}

class Func1D
{
public:
  virtual float func( float x ) = 0;
};

void createVisDegDigitalFilterBox( Func1D *func, pfs::Array2D *filter,
  float pixelsPerDeg );


void createVisDegDigitalFilter( float(*func)(float), pfs::Array2D *filter,
  float pixelsPerDeg );

void createVisDegDigitalFilter( Func1D *func, pfs::Array2D *filter,
  float pixelsPerDeg );

void multiplyAndAddArray(pfs::Array2D *z, const pfs::Array2D *x, const float f);



/**
 * returns pupil diameter in mm for the adaptation level of Y cd/m^2
 * from:
 * P. Moon and D. Spencer, "On the Stiles-Crawford effect," J. Opt. Soc.
 * Am.  34, 319-329 (1944), Equation 3
 */
float getPupilDiameter( float Y )
{
  return (4.9 - 3 * tanhf( 0.4 * (log10f( Y * M_PI ) - 0.5) ))*1e-3;
}


//-------------------------------------------------
// Optical Tranfer Function from
// 
// Normann, R. A., Baxter, B. S., Ravindra, H., and Anderton, P. J.
// The Photoreceptor Contributions to Contrast Sensitivity:
// Applications in Radiological Diagnosis.
// IEEE Trans. on Systems, Man and Cybernetics.  SMC-13(5):944-953; 1983.
//-------------------------------------------------

class NormannBaxterOTFFunc : public Func1D
{
public:    
  float func( float angle )
  {
    const float A = 1.13, B = 2.35, C = 1.42, D = 9.85, E = 9.35e-6, F=7.8;
    angle *= 60;                  // To get minutes
    return angle <= 7 ?
      A*expf(-B*B*angle*angle) + C/(D+powf(angle,3.3)) :
      E*powf(F+angle, -2.5);
  }
};

NormannBaxterOTF::NormannBaxterOTF( const ViewingConditions &viewCond, int cols, int rows ) : filter( cols, rows )
{
  pfs::Array2D *filterSpatial = filter.setSpatial();
  const float pixelsPerDeg = viewCond.getPixelsPerDegree(
    (viewCond.maxDistance + viewCond.minDistance)/2. );
  
  std::cerr << "Creating OTF filter..";

  std::cerr << "(Pixels per deg: " << pixelsPerDeg << ")";

  NormannBaxterOTFFunc normannBaxterOTFFunc;
  createVisDegDigitalFilter( &normannBaxterOTFFunc, filterSpatial, pixelsPerDeg );

  dumpImageAspect->dump( "filter_otf.pfs", filterSpatial, "Y" );
  dumpImageAspect->dumpFrequency( "filter_otf_fft.pfs", &filter );
  
  std::cerr << ".\n";
  
}


void NormannBaxterOTF::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency() if (in == out)  
  multiplyArray( out->setFrequency(), input, filter.getFrequency() );
}




//-------------------------------------------------
// Optical Tranfer Function from
//
// Westheimer (1986)
// used in Sarnoff Visual Discrimination Model
//-------------------------------------------------


class WestheimerOTFFunc : public Func1D
{
public:    
  float func( float angle )
  {
    angle *= 60;                  // To get minutes
    return 0.952*exp( -2.59*powf( fabsf(angle), 1.36 ) ) +
      0.048*exp( -2.43* powf( fabsf(angle), 1.74 ) );
  }
};

WestheimerOTF::WestheimerOTF( const ViewingConditions &viewCond, int cols, int rows ) : filter( cols, rows )
{
  pfs::Array2D *filterSpatial = filter.setSpatial();
  const float pixelsPerDeg = viewCond.getPixelsPerDegree(
    (viewCond.maxDistance + viewCond.minDistance)/2. );
  
  std::cerr << "Creating OTF filter..";
  
  WestheimerOTFFunc normannBaxterOTFFunc;
  createVisDegDigitalFilter( &normannBaxterOTFFunc, filterSpatial, pixelsPerDeg );

  dumpImageAspect->dump( "filter_otf.pfs", filterSpatial, "Y" );
  dumpImageAspect->dumpFrequency( "filter_otf_fft.pfs", &filter );
  
  std::cerr << ".\n";
  
}


void WestheimerOTF::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency() if (in == out)  
  multiplyArray( out->setFrequency(), input, filter.getFrequency() );
}


//-------------------------------------------------
// Optical Tranfer Function from
//
// @inproceedings{218466,
//  author = {Greg Spencer and Peter Shirley and Kurt Zimmerman and Donald P. Greenberg},
//  title = {Physically-based glare effects for digital images},
//  booktitle = {Proceedings of the 22nd annual conference on Computer graphics and interactive techniques},
//  year = {1995},
//  isbn = {0-89791-701-4},
//  pages = {325--334},
//  doi = {http://doi.acm.org/10.1145/218380.218466},
//  publisher = {ACM Press},
//  }
//  
//-------------------------------------------------


static float spencerPSF0( float theta )
{
  return 2610000*expf( -(theta/0.002)*(theta/0.002) );
}

static float spencerPSF1( float theta )
{
  return 20.91/powf(theta+0.02, 3);
}

static float spencerPSF2( float theta )
{
  return 72.37/((theta+0.02)*(theta+0.02));
}

// Lenticular halo
static float spencerPSF3( float theta, float lambda )
{
  return 436.9 * (568/lambda)*expf( -19.75*(theta-3*lambda/568)*(theta-3*lambda/568) );  
}

// Lenticular halo for white light 
inline static float spencerPSF3_Y( float theta )
{
  return 0.2126 * spencerPSF3( theta, 558 ) + // R
    0.7152 * spencerPSF3( theta, 531 ) +      // G
    0.0722 * spencerPSF3( theta, 419 );       // B
}

class SpencerFunc : public Func1D
{
  const float w0, w1, w2, w3;
public:
  SpencerFunc( float w0, float w1, float w2, float w3 ):
    w0( w0 ), w1( w1 ), w2( w2 ), w3( w3 )
  {
  }
    
  float func( float angle )
  {
    if( w3 != 0 )
      return w0*spencerPSF0( angle ) +
        w1*spencerPSF1( angle ) +
        w2*spencerPSF2( angle ) +
        w3*spencerPSF3_Y( angle );
    else 
      return w0*spencerPSF0( angle ) +
        w1*spencerPSF1( angle ) +
        w2*spencerPSF2( angle );
  }
};

SpencerOTF::SpencerOTF( const ViewingConditions &viewCond, int cols, int rows, float adaptationLuminance ) : filter( cols, rows )
{
  const float pixelsPerDeg = viewCond.getPixelsPerDegree( (viewCond.maxDistance +
                                                            viewCond.minDistance)/2. );
  std::cerr << "Creating OTF filter ";

  const float SCOTOPIC = 0.01, PHOTOPIC = 3;

  float w[4];
  
  const float scotopicW[4] = { 0.282, 0.478, 0.207, 0.033 };
//  const float mesopicW[4] = { 0.368, 0.478, 0.138, 0.016 };
  const float photopicW[4] = { 0.384, 0.478, 0.138, 0 };

  if( adaptationLuminance > PHOTOPIC ) {
    for( int i = 0; i < 4; i++ )
      w[i] = photopicW[i];
    std::cerr << "for photopic conditions..";
  } else if( adaptationLuminance <= SCOTOPIC ) {
    for( int i = 0; i < 4; i++ )
      w[i] = scotopicW[i];
    std::cerr << "for scotopic conditions..";
  } else {
    float ad = (log( adaptationLuminance ) - log( SCOTOPIC )) /
      (log(PHOTOPIC)-log(SCOTOPIC));
    for( int i = 0; i < 4; i++ )
      w[i] = scotopicW[i]*(1-ad) + photopicW[i]*ad;
    std::cerr << "for mesopic conditions..";    
  }

  pfs::Array2DImpl f( cols, rows );
  pfs::Array2D *filterSpatial = filter.setSpatial();
  pfs::setArray( filterSpatial, 0 );
  
  createVisDegDigitalFilter( spencerPSF0, &f, pixelsPerDeg );
  normalizeToSum( &f );
  multiplyAndAddArray( filterSpatial, &f, w[0] );

  createVisDegDigitalFilter( spencerPSF1, &f, pixelsPerDeg );
  normalizeToSum( &f );
  multiplyAndAddArray( filterSpatial, &f, w[1] );

  createVisDegDigitalFilter( spencerPSF2, &f, pixelsPerDeg );
  normalizeToSum( &f );
  multiplyAndAddArray( filterSpatial, &f, w[2] );

  createVisDegDigitalFilter( spencerPSF3_Y, &f, pixelsPerDeg );
  normalizeToSum( &f );
  multiplyAndAddArray( filterSpatial, &f, w[3] );
  
  normalizeToSum( filterSpatial );
  
  dumpImageAspect->dump( "filter_otf.pfs", filterSpatial, "Y" );
  dumpImageAspect->dumpFrequency( "filter_otf_fft.pfs", &filter );  
  std::cerr << ".\n";
}

void SpencerOTF::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency() if (in == out)  
  multiplyArray( out->setFrequency(), input, filter.getFrequency() );
}

class Func2D
{
public:
  virtual float func( float x, float y ) = 0;
};

class VisDegFunc2D : public Func2D
{
  Func1D *innerFunc;
  const float pixelsPerDeg;
public:
  VisDegFunc2D( Func1D *innerFunc, float pixelsPerDeg ):
    innerFunc( innerFunc ), pixelsPerDeg( pixelsPerDeg )
  {
  }
  
  virtual float func( float x, float y )
  {
    float angle = sqrtf( x*x + y*y ) / pixelsPerDeg; // Given in degrees
    return innerFunc->func( angle );
  }
  
};


float trapzd2d(Func2D *func, float x1, float x2,
  float y1, float y2, int n);

/**
 * Create a digital filter from a function of visual angle. Sample
 * points near low angles more densely.
 * @param func Function returning filter value for given visual
 * angle in vidual degrees
 */
void createVisDegDigitalFilterBox( Func1D *func, pfs::Array2D *filter,
  float pixelsPerDeg )
{
  const int cols = filter->getCols(), rows = filter->getRows();

  const float sampleMax = 10, sampleMin = 2,
    maxR = sqrtf( cols*cols/4 + rows*rows/4 );

  VisDegFunc2D visDegFunc2D( func, pixelsPerDeg );
  
  for( int x = 0; x < cols/2; x++ )
    for( int y = 0; y < rows/2; y++ ) {
      float delta = 1/(sampleMin + sampleMax -
        (sqrtf( x*x + y*y )/maxR)*(sampleMax-sampleMin));
      
      float v = trapzd2d( &visDegFunc2D, (float)x-0.5, (float)x+0.5,
        (float)y-0.5, (float)y+0.5, (x<15 && y<15) ? 8 : 4);
      
      (*filter)(x,y) = v;
      if( x != 0 ) (*filter)(cols-x,y) = v;
      if( y != 0 ) (*filter)(x,rows-y) = v;      
      if( x != 0 && y != 0 )(*filter)(cols-x,rows-y) = v;
    }
  
  normalizeToSum( filter );
  
}

/*
 * Trapezoidal integration from
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * + some extensions to 2D 
 */

/**
 * This routine computes the nth stage of refinement of an extended
 * trapezoidal rule. func is input as a pointer to the function to be
 * integrated between limits a and b, also input. When called with
 * n=1, the routine returns the crudest estimate of b a
 * f(x)dx. Subsequent calls with n=2,3,... (in that sequential order)
 * will improve the accuracy by adding 2n-2 additional interior
 * points.
 */
float trapzd(Func2D *func, float y, float a, float b, int n)
{
  float x,tnm,sum,del;
  static float s;
  int it,j;
  if (n == 1) {
    return (s=0.5*(b-a)*(func->func(a, y)+func->func(b, y)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm; //This is the spacing of the points to be added.
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += func->func(x, y);
    s=0.5*(s+(b-a)*sum/tnm);  //This replaces s by its refined value.
    return s;
  }
}


class TrapezFunc2D : public Func2D
{
  Func2D *innerFunc;  
  const float x1, x2;
public:
  TrapezFunc2D( Func2D *innerFunc, float x1, float x2 ) :
    innerFunc( innerFunc ), x1( x1 ), x2( x2 )
  {
  }
  virtual float func( float y, float n )
  {
    return trapzd( innerFunc, y, x1, x2, (int)n );
  }
};
  

/**
 * 2D trapezoidal integration
 */
float trapzd2d(Func2D *func, float x1, float x2,
  float y1, float y2, int n)
{
  TrapezFunc2D trapezFunc2D( func, x1, x2 );

  return trapzd( &trapezFunc2D, n, y1, y2, n );
}


void createVisDegDigitalFilter( float(*func)(float), pfs::Array2D *filter,
  float pixelsPerDeg )
{
  const int cols = filter->getCols(), rows = filter->getRows();

  for( int x = 0; x < cols/2; x++ )
    for( int y = 0; y < rows/2; y++ ) {
      float delta = sqrtf( x*x + y*y )/pixelsPerDeg;
      float v = func( delta );
      
      (*filter)(x,y) = v;
      if( x != 0 ) (*filter)(cols-x,y) = v;
      if( y != 0 ) (*filter)(x,rows-y) = v;      
      if( x != 0 && y != 0 )(*filter)(cols-x,rows-y) = v;
    }
  
  normalizeToSum( filter );
  
}

void createVisDegDigitalFilter( Func1D *func, pfs::Array2D *filter,
  float pixelsPerDeg )
{
  const int cols = filter->getCols(), rows = filter->getRows();

  for( int x = 0; x < cols/2; x++ )
    for( int y = 0; y < rows/2; y++ ) {
      float delta = sqrtf( x*x + y*y )/pixelsPerDeg;
      float v = func->func( delta );
      
      (*filter)(x,y) = v;
      if( x != 0 ) (*filter)(cols-x,y) = v;
      if( y != 0 ) (*filter)(x,rows-y) = v;      
      if( x != 0 && y != 0 )(*filter)(cols-x,rows-y) = v;
    }
  
  normalizeToSum( filter );
  
}


/**
 * z = z + x*f
 */
void multiplyAndAddArray( pfs::Array2D *z, const pfs::Array2D *x, const float f )
{
  const int elements = x->getRows()*x->getCols();
  for( int i = 0; i < elements; i++ )
    (*z)(i) += (*x)(i) * f;
}


//-------------------------------------------------
// Optical Tranfer Function from
//
// Marimont D.H. and Wandell B.A., 1994
// Matching  colour images; the effects of axial chromatic aberration
// Journal of the Optical Society of America A, 11 (no 12) 2113-3122 
//  
//-------------------------------------------------

struct SpectralSensitivity
{
  float lambda, weight;
};

SpectralSensitivity luminanceEfficiencyFunction[] = {
     380, 0.0002000,
     385, 0.0003960,
     390, 0.0008000,
     395, 0.0015500,
     400, 0.0028000,
     405, 0.0046600,
     410, 0.0074000,
     415, 0.0118000,
     420, 0.0175000,
     425, 0.0227000,
     430, 0.0273000,
     435, 0.0326000,
     440, 0.0379000,
     445, 0.0424000,
     450, 0.0468000,
     455, 0.0521000,
     460, 0.0600000,
     465, 0.0739000,
     470, 0.0910000,
     475, 0.1130000,
     480, 0.1390000,
     485, 0.1690000,
     490, 0.2080000,
     495, 0.2590000,
     500, 0.3230000,
     505, 0.4073000,
     510, 0.5030000,
     515, 0.6082000,
     520, 0.7100000,
     525, 0.7932000,
     530, 0.8620000,
     535, 0.9148500,
     540, 0.9540000,
     545, 0.9803000,
     550, 0.9949500,
     555, 1.0000000,
     560, 0.9950000,
     565, 0.9786000,
     570, 0.9520000,
     575, 0.9154000,
     580, 0.8700000,
     585, 0.8163000,
     590, 0.7570000,
     595, 0.6949000,
     600, 0.6310000,
     605, 0.5668000,
     610, 0.5030000,
     615, 0.4412000,
     620, 0.3810000,
     625, 0.3200000,
     630, 0.2650000,
     635, 0.2170000,
     640, 0.1750000,
     645, 0.1382000,
     650, 0.1070000,
     655, 0.0816000,
     660, 0.0610000,
     665, 0.0446000,
     670, 0.0320000,
     675, 0.0232000,
     680, 0.0170000,
     685, 0.0119000,
     690, 0.0082100,
     695, 0.0057200,
     700, 0.0041000,
     705, 0.0029300,
     710, 0.0020900,
     715, 0.0014800,
     720, 0.0010500,
     725, 0.0007400,
     730, 0.0005200,
     735, 0.0003610,
     740, 0.0002490,
     745, 0.0001720,
     750, 0.0001200,
     755, 0.0000848,
     760, 0.0000600,
     765, 0.0000424,
     770, 0.0000300,
     775, 0.0000212,
     780, 0.0000150,
};

MarimontOTF::MarimontOTF( const ViewingConditions &viewCond,
  int cols, int rows, float adaptationLuminance )
{
  const float n_p = 1.336;
  const float f_p = 22.2888e-3;
  const float D0 = n_p/f_p;

  const float c = 3.434e3;

  const float q1 = 1.7312;
  const float q2 = 0.63346;
  const float q3 = 0.21410;

  filter = new pfs::Array2DImpl( cols/2+1, rows/2+1 );

  
  int filterWidth=filter->getCols();
  int filterHeight=filter->getRows();
  
  float meanObserverDistance = (viewCond.minDistance + viewCond.maxDistance)/2;
  float pix_per_deg = viewCond.getPixelsPerDegree( meanObserverDistance );

  float x_norm = 0.5 / (filterWidth-1), y_norm = 0.5 / (filterHeight-1);

  float dx, dy;
  int x, y;

  float p = getPupilDiameter( adaptationLuminance )/2; // Pupil radius in meters

  std::cerr << "Creating OTF filter..";

  std::cerr << "(pupil diameter " << p*1e3*2 << "mm for adaptation lum " << adaptationLuminance << "cd/m^2)";  

  for (y = 0, dy = 0.0; y < filterHeight; y ++, dy += y_norm)
  {
    for (x = 0, dx = 0.0; x < filterWidth; x ++, dx += x_norm)
    {
      float v = pix_per_deg * sqrtf((float)(dx * dx + dy * dy));
      if( v == 0 ) v = 0.1;

      const int lefSize = sizeof( luminanceEfficiencyFunction ) / sizeof( SpectralSensitivity );
      float H = 0;
      for( int l = 0; l < lefSize; l++ ) {
        
        const float lambda = luminanceEfficiencyFunction[l].lambda*1e-9;
    
        const float s = c*(lambda / (D0*p)) * v;
        
        const float D_lambda = q1 - (q2/(lambda*1e6-q3));
        const float w20 = p*p/2 * D0 * D_lambda / (D0 + D_lambda);
        const float alpha = 4*M_PI/lambda * w20 * fabsf( s );
        const float py_max = sqrtf( (1-(s/2))*(1-(s/2)) );

        const float dpy = py_max/10;
        float h = 0;
        for( float py = 0; py <= py_max; py += dpy ) {  
          h += 4/(M_PI*alpha)*sinf( alpha*(sqrtf(1 - py*py) - fabsf(s)/2) )*dpy;
        }
        H += h*luminanceEfficiencyFunction[l].weight;
      }      

      (*filter)(x,y) = H * (0.3481+0.6519*exp(-0.1212*v));
      
    }
  }

  (*filter)(0,0) = max( (*filter)(1,0), (*filter)(0,1) );

  normalizeToMax( filter );

  dumpImageAspect->dump( "filter_otf.pfs", filter, "Y" );
  
  std::cerr << ".\n";  
   
}

MarimontOTF::~MarimontOTF()
{
  delete filter;
}


void MarimontOTF::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const int cols = in->getCols(), rows = in->getRows();
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency()->getData() if (in == out)
  filterFFTW( input->getData(), out->setFrequency()->getData(), cols, rows, filter );
}


//-------------------------------------------------
// Optical Tranfer Function from
//
// Deeley, R.J., Drasdo, N., & Charman, W. N. (1991)
// A simple parametric model of the human ocular modulation transfer function.
// Ophthalmology and Physiological Optics, 11, 91-93
//
//-------------------------------------------------

DeeleyOTF::DeeleyOTF( const ViewingConditions &viewCond, int cols, int rows, float adaptationLuminance )
{

  filter = new pfs::Array2DImpl( cols/2+1, rows/2+1 );

  int filterWidth=filter->getCols();
  int filterHeight=filter->getRows();
  
  float meanObserverDistance = (viewCond.minDistance + viewCond.maxDistance)/2;
  float pix_per_deg = viewCond.getPixelsPerDegree( meanObserverDistance );

  float x_norm = 0.5 / filterWidth, y_norm = 0.5 / filterHeight;

  float dx, dy;
  int x, y;

  float p = getPupilDiameter( adaptationLuminance );
  float p_mm = p * 1e3;

  std::cerr << "Creating OTF filter..";


    std::cerr << "(pupil diameter " << p_mm << "mm for adaptation lum " << adaptationLuminance << "cd/m^2)";  

  for (y = 0, dy = 0.0; y < filterHeight; y ++, dy += y_norm)
  {
    for (x = 0, dx = 0.0; x < filterWidth; x ++, dx += x_norm)
    {
      float v = pix_per_deg * sqrtf((float)(dx * dx + dy * dy));
      if( v == 0 ) v = 0.1;
      
      (*filter)(x,y) = expf( -pow(v/ (20.9 - 2.1*p_mm), 1.3-0.07*p_mm) );
      
    }
  }

  (*filter)(0,0) = max( (*filter)(1,0), (*filter)(0,1) );

  normalizeToMax( filter );

  dumpImageAspect->dump( "filter_otf.pfs", filter, "Y" );
  
  std::cerr << ".\n";
}

DeeleyOTF::~DeeleyOTF()
{
  delete filter;
}


void DeeleyOTF::process( BidomainArray2D *in, BidomainArray2D *out )
{
  const int cols = in->getCols(), rows = in->getRows();
  const FFTWComplexArray *input = in->getFrequency(); // This must be executed before out->setFrequency()->getData() if (in == out)
  filterFFTW( input->getData(), out->setFrequency()->getData(), cols, rows, filter );
}
