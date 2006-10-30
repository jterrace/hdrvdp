/**
 * @brief visualization of the HDR VDP result 
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
 * $Id: vdpvis.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <pfs.h>
#include <getopt.h>

#define FIFTY_THRESHOLD       0.50
#define FIFTY1_THRESHOLD      0.625
#define PROBABILITY_THRESHOLD 0.75
#define ALMOST_SURE_THRESHOLD 0.95


#define PROG_NAME "vdpvis"

void colorCode( const pfs::Array2D *vdp, const pfs::Array2D *L, pfs::Array2D *R, pfs::Array2D *G, pfs::Array2D *B );

void computePercentageOfPixels( const pfs::Array2D *vdp, const float prob1, float &procPixelsProb1,
  const float prob2, float &procPixelsProb2 );


class QuietException 
{
};

static void errorCheck( bool condition, char *string )
{
  if( !condition ) {
    throw pfs::Exception( string );
  }
}

void printHelp()
{
  fprintf( stderr, PROG_NAME " [--summary=file.csv] [--verbose] [--help]\n"
    "See man page for more information.\n" );
}

void visualizeVDP(int argc, char **argv)
{
  bool verbose = false;
  const char *outputFile = NULL;

  static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { "summary", required_argument, NULL, 's' },
    { NULL, 0, NULL, 0 }
  };

  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long(argc, argv, "s:hv", cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
      break;
    case 's':
      outputFile = optarg;      
      break;      
    case '?':
      throw QuietException();
    case ':':
      throw QuietException();
    }
  }


  pfs::DOMIO pfsio;
    
  // Load frame with VDP output from stdin
  pfs::Frame *frame = pfsio.readFrame( stdin );
  errorCheck( frame != NULL, "Frame not found in pfs stream" );
  pfs::Array2D *vdp = frame->getChannel( "VDP" );
  errorCheck( vdp != NULL, "There is no vdp output in the input frame" );
  
  pfs::Channel *Y = frame->getChannel( "Y" );
  errorCheck( Y != NULL, "There is no luminance channel in the input frame" );

  float percent75, percent95;
  computePercentageOfPixels( vdp, 0.75f, percent75, 0.95f, percent95 );
  fprintf( stderr, "Probability of detection P>75%% for %2.4f%% of pixels\n", percent75*100.f );
  fprintf( stderr, "Probability of detection P>95%% for %2.4f%% of pixels\n", percent95*100.f );
  if( outputFile != NULL ) {
    std::ofstream file( outputFile );
    file << percent75 << ", " << percent95 << "\n";
    file.close();    
  }
  
  pfs::Channel *L = Y;          // Lightness
  {                             // Tone map image and compute lightness
    const int size = Y->getCols()*Y->getRows();

    // TODO: Take also variance into account to modify 0.63
    float count = 0;
    float logMean = 0;    
    for( int i = 0; i < size; i++ )
      if( (*Y)(i) > 0.0001 ) { 
        logMean += log( (*Y)(i) );
        count++;
      }
    logMean = exp( logMean / count );    
    
    for( int i = 0; i < size; i++ ) {
      (*L)(i) = 1 - 1/(1+powf( (*Y)(i)/logMean, 0.63 ));
    }
  }
  

  pfs::Channel *X, *Z;
  X = frame->createChannel( "X" );
  Z = frame->createChannel( "Z" );
  
  pfs::Channel *R = X, *G = Y, *B = Z;

  colorCode( vdp, L, R, G, B );
  
  pfs::transformColorSpace( pfs::CS_RGB, R, G, B, pfs::CS_XYZ, X, Y, Z );

  frame->getTags()->setString( "LUMINANCE", "DISPLAY" ); // LDR image
  
  pfsio.writeFrame( frame, stdout );

  pfsio.freeFrame( frame );

  
  
}


int main( int argc, char* argv[] )
{
  try {
    visualizeVDP( argc, argv );
  }
  catch( pfs::Exception ex ) {
    fprintf( stderr, PROG_NAME " error: %s\n", ex.getMessage() );
    return EXIT_FAILURE;
  }        
  catch( QuietException  ex ) {
    return EXIT_FAILURE;
  }        
  return EXIT_SUCCESS;
}

/**
 * Count pixels that have higher probability than prob1 and prob2
 */

void computePercentageOfPixels( const pfs::Array2D *vdp, const float prob1, float &procPixelsProb1, const float prob2, float &procPixelsProb2 )
{
  const int size = vdp->getCols()*vdp->getRows();

  procPixelsProb1 = 0;
  procPixelsProb2 = 0;
  
  for( int i = 0; i < size; i++ ) {
    const float prob = fabs( (*vdp)(i) );
    if( prob >= prob1 ) procPixelsProb1++;
    if( prob >= prob2 ) procPixelsProb2++;
  }

  procPixelsProb1 /= (float)size;
  procPixelsProb2 /= (float)size;
}


void colorCode( const pfs::Array2D *vdp, const pfs::Array2D *L, pfs::Array2D *R, pfs::Array2D *G, pfs::Array2D *B )
{
  const int size = vdp->getCols()*vdp->getRows();
  
  // Color coding
  for( int i = 0; i < size; i++ ) {
    const float prob = fabs( (*vdp)(i) );
    const float l = (*L)(i) * 0.7 + 0.2;    // Lightness (scaled to
                                            // avoid desaturated
                                            // colors for low and high
                                            // lightness

      // My simple coloring
//     if( prob < FIFTY_THRESHOLD ) {
//       (*R)(i) = 0.f;
//       (*G)(i) = 1.f;
//       (*B)(i) = 0.f;
//     } else if( prob < PROBABILITY_THRESHOLD ) {
//       (*R)(i) = 1.f;
//       (*G)(i) = 1.f;
//       (*B)(i) = 0.f;
//     } else if( prob < ALMOST_SURE_THRESHOLD ) {
//       (*R)(i) = 1.f;
//       (*G)(i) = 0.f;
//       (*B)(i) = 0.f;
//     } else {
//       (*R)(i) = 1.f;
//       (*G)(i) = 0.5f;
//       (*B)(i) = 0.5f;
//     }


      // Karol's PAPUGA coloring
      if( prob < PROBABILITY_THRESHOLD ) {
        if( prob > FIFTY_THRESHOLD ) {
          if( prob > FIFTY1_THRESHOLD ) {
            float t = (prob - FIFTY1_THRESHOLD) / (0.5 * (PROBABILITY_THRESHOLD - FIFTY_THRESHOLD));
            (*R)(i) = 1.0;
            (*G)(i) = 1.0 - t;
            (*B)(i) = 0.0;
          } else {
            float t = (prob - FIFTY_THRESHOLD) / (0.5 * (PROBABILITY_THRESHOLD - FIFTY_THRESHOLD));
            (*R)(i) = t;
            (*G)(i) = 1.0;
            (*B)(i) = 0.0;
          }
        } else {
          if( prob < FIFTY_THRESHOLD*0.5 ) {
            (*R)(i) = 0.8f;
            (*G)(i) = 0.8f;
            (*B)(i) = 0.8f;
          } else {
            float t = prob / FIFTY_THRESHOLD;
            (*R)(i) = 0.0;
            (*G)(i) = t;
            (*B)(i) = 0.0;
          }
        }
      } else if( prob > ALMOST_SURE_THRESHOLD ) {
        float t = (prob - ALMOST_SURE_THRESHOLD) / (1.0 - ALMOST_SURE_THRESHOLD);
        (*R)(i) = 1.f;
        (*G)(i) = 0.f;
        (*B)(i) = 0.5f * t;
      } else {
        (*R)(i) = 1.f;
        (*G)(i) = 0.f;
        (*B)(i) = 0.f;
      }    

      // Karol's no PAPUGA coloring        
//     if( prob < PROBABILITY_THRESHOLD )
//     {
//       float t = prob / PROBABILITY_THRESHOLD;
//       t -= 0.5;
//       (*R)(i) = 3.f * t;
//       (*G)(i) = 3.f * t;
//       (*B)(i) = 0.f;
//     }
//     else if( prob > ALMOST_SURE_THRESHOLD ) {
//       float t = (prob - ALMOST_SURE_THRESHOLD) / (1.0 - ALMOST_SURE_THRESHOLD);
//       (*R)(i) = 1.f;
//       (*G)(i) = 0.f;
//       (*B)(i) = 0.5f * t;
//     } else {
//       (*R)(i) = 1.f;
//       (*G)(i) = 0.f;
//       (*B)(i) = 0.f;
//     }    

      // Influence of lightness
       (*R)(i) *= l;
       (*G)(i) *= l;
       (*B)(i) *= l;
    
    
  }
}
