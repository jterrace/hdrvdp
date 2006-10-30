/**
 * @brief dump partial result as a pfs image
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
 * $Id: dump_image_aspect.cpp,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#include <stdarg.h>
#include <stdio.h>
#include <fnmatch.h>
#include <math.h>

#include "dump_image_aspect.h"

#include <pfs.h>

DumpImageAspect *dumpImageAspect;

void initializeDumpImageAspect()
{
  dumpImageAspect = new DumpImageAspect();
}


static void dumpArray2D( const pfs::Array2D *data,
  const char *fileName, const char *channelName )
{
  pfs::DOMIO pfsio;

  pfs::Frame *frame = pfsio.createFrame( data->getCols(), data->getRows() );
  pfs::Channel *channel = frame->createChannel( channelName );
  pfs::copyArray( data, channel );

  FILE *out = fopen( fileName, "wb" );
  pfsio.writeFrame( frame, out );
  fclose( out );
  pfsio.freeFrame( frame );
  
}
    
void DumpImageAspect::dump( const char *fileName, const pfs::Array2D *data, const char *channelName, ... )
{
  va_list v_args;
  va_start( v_args, channelName );
  
  char fileNameBuf[256];
  vsprintf( fileNameBuf, fileName, v_args );

//  fprintf( stderr, "Dumping image: %s\n", fileNameBuf );

  PatternList::iterator it;
  for( it = patterns.begin(); it != patterns.end(); it++ )
    if( fnmatch( *it, fileNameBuf, 0 ) == 0 ) break;

  // No pattern match - do not dump
  if( it == patterns.end() ) return;

  dumpArray2D( data, fileNameBuf, channelName );
  
  va_end( v_args );

    
}

void DumpImageAspect::dumpFrequency( const char *fileName, BidomainArray2D *data, ... )
{
  va_list v_args;
  va_start( v_args, data );
  
  char fileNameBuf[256];
  vsprintf( fileNameBuf, fileName, v_args );

//  fprintf( stderr, "Dumping image: %s\n", fileNameBuf );

  PatternList::iterator it;
  for( it = patterns.begin(); it != patterns.end(); it++ )
    if( fnmatch( *it, fileNameBuf, 0 ) == 0 ) break;

  // No pattern match - do not dump
  if( it == patterns.end() ) return;

  pfs::DOMIO pfsio;

  const FFTWComplexArray *complexArray = data->getFrequency();
  pfs::Frame *frame = pfsio.createFrame( complexArray->getFFTCols(),
    complexArray->getFFTRows() );
  pfs::Channel *channel = frame->createChannel( "abs" );
  const int size = complexArray->getFFTCols()*complexArray->getFFTRows();
  const fftwf_complex *d = complexArray->getData();
  for( int i = 0; i < size; i++ )
    (*channel)(i) = sqrt( d[i][0]*d[i][0] + d[i][1]*d[i][1] );

  channel = frame->createChannel( "angle" );
  for( int i = 0; i < size; i++ )
    (*channel)(i) = atan( d[i][1] / d[i][0] );
  
  FILE *out = fopen( fileNameBuf, "wb" );
  pfsio.writeFrame( frame, out );
  fclose( out );
  pfsio.freeFrame( frame );
  
  va_end( v_args );
    
}

void DumpImageAspect::dump( const char *fileName, BidomainArray2D *data, const char *channelName, ... )
{
  va_list v_args;
  va_start( v_args, channelName );

  char fileNameBuf[256];
  vsprintf( fileNameBuf, fileName, v_args );

  PatternList::iterator it;
  for( it = patterns.begin(); it != patterns.end(); it++ )
    if( fnmatch( *it, fileNameBuf, 0 ) == 0 ) break;

  // No pattern match - do not dump
  if( it == patterns.end() ) return;
  
  const pfs::Array2D *array = data->getSpatial();

  dumpArray2D( array, fileNameBuf, channelName );
  
  
  va_end( v_args );
}



void DumpImageAspect::addIncludePattern( char *pattern )
{
  patterns.push_front( pattern );
}

