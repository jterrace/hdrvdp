/**
 * @brief Blur images in PFS stream
 *
 * @author Rafal Mantiuk, <mantiuk@mpi-sb.mpg.de>
 *
 * $Id: pfsblur.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include <pfs.h>

#define PROG_NAME "pfsblur"


class QuietException 
{
};

void blur( pfs::Array2D *array, const float sigma, const int support );


void printHelp()
{
  fprintf( stderr, PROG_NAME " [--sigma <val>] [--verbose] [--help]\n"
    "See man page for more information.\n" );
}

void blurFrames( int argc, char* argv[] )
{
  pfs::DOMIO pfsio;

  float sigma = 2;
  int support = -1;
  bool verbose = false;

  static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { "sigma", required_argument, NULL, 's' },
    { NULL, 0, NULL, 0 }
  };

  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long (argc, argv, "s:", cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
      break;
    case 's':
      sigma = (float)strtod( optarg, NULL );
      break;
    case '?':
      throw QuietException();
    case ':':
      throw QuietException();
    }
  }

  support = (int)(sigma+1);

  if( verbose ) fprintf( stderr, "Bluring channels with gaussian of sigma %g\n", sigma );
   
  while( true ) {
    pfs::Frame *frame = pfsio.readFrame( stdin );
    if( frame == NULL ) break; // No more frames
        
    pfs::Channel *X, *Y, *Z;
    frame->getXYZChannels( X, Y, Z );

    if( X != NULL ) {           // Color, XYZ
      blur( X, sigma, support );
      blur( Y, sigma, support );
      blur( Z, sigma, support );
    } else if( (Y = frame->getChannel( "Y" )) != NULL ) {
      // Luminance only
      blur( Y, sigma, support );
    } else
      throw pfs::Exception( "Missing X, Y, Z channels in the PFS stream" );        

    pfsio.writeFrame( frame, stdout );
    pfsio.freeFrame( frame );        
  }
}



void blur( pfs::Array2D *array, const float sigma, const int support )
{
  const int fullSupport = support*2 +1;
  const int rows = array->getRows();
  const int cols = array->getCols();
  
  float *kernel = new float[ fullSupport ];

  // Generate kernel
  float sum = 0;
  for( int i = 0; i < fullSupport; i++ ) {
    float x = (float)(i-support);
    kernel[i] = expf( -x*x / (2.*sigma*sigma) );
    sum += kernel[i];
  }
    
  // Process rows
  for( int r = 0; r < rows; r++ ) {
    for( int c = 0; c < cols; c++ ) {
      float sum = 0, k = 0;
      for( int i = 0; i < fullSupport; i++ ) {
        int ic = c + i - support;
        if( ic < 0 || ic >= cols ) continue;
        sum += (*array)(ic,r)*kernel[i];
        k += kernel[i];
      }
      (*array)(c,r) = sum / k;
    }
  }      

  // Process columns
  for( int c = 0; c < cols; c++ ) {
    for( int r = 0; r < rows; r++ ) {
      float sum = 0, k = 0;
      for( int i = 0; i < fullSupport; i++ ) {
        int ir = r + i - support;
        if( ir < 0 || ir >= rows ) continue;
        sum += (*array)(c,ir)*kernel[i];
        k += kernel[i];
      }
      (*array)(c,r) = sum / k;
    }
  }      
  
  delete[] kernel;
}



int main( int argc, char* argv[] )
{
  try {
    blurFrames( argc, argv );
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
