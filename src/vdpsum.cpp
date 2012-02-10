/**
 * @brief sums the probabilities of VDP output
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
 * $Id: vdpvis.cpp,v 1.2 2008/08/21 17:48:00 rafm Exp $
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <pfs.h>
#include <getopt.h>

#define PROG_NAME "vdpsum"

double computeSumOfPixels( const pfs::Array2D *vdp );

class QuietException
{
};

static void errorCheck( bool condition, const char *string )
{
  if( !condition ) {
    throw pfs::Exception( string );
  }
}

void printHelp()
{
  fprintf( stderr, PROG_NAME " [--verbose] [--help]\n"
    "See man page for more information.\n" );
}

void sumVDP(int argc, char **argv)
{
  bool verbose = false;
  const char *outputFile = NULL;

  static struct option cmdLineOptions[] = {
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };

  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long(argc, argv, ":hv", cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
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

  pfs::ChannelIteratorPtr it( frame->getChannelIterator() );
  while( it->hasNext() ) {
      pfs::Channel *ch = it->getNext();
      fprintf(stderr, "found channel: %s\n", ch->getName());
  }

  pfs::Array2D *vdp = frame->getChannel( "VDP" );
  errorCheck( vdp != NULL, "There is no vdp output in the input frame" );

  pfs::Channel *Y = frame->getChannel( "Y" );
  errorCheck( Y != NULL, "There is no luminance channel in the input frame" );

  double sum = computeSumOfPixels( vdp );
  fprintf( stderr, "probability sum: %g\n", sum);

}


int main( int argc, char* argv[] )
{
  try {
    sumVDP( argc, argv );
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
 * Sum probability of each pixel
 */

double computeSumOfPixels( const pfs::Array2D *vdp )
{
  const int size = vdp->getCols()*vdp->getRows();

  double sum = 0;

  for( int i = 0; i < size; i++ ) {
    const float prob = fabs( (*vdp)(i) );
    sum += prob;
  }

  return sum;
}
