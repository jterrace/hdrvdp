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
 * $Id: dump_image_aspect.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef DUMP_IMAGE_ASPECT_H
#define DUMP_IMAGE_ASPECT_H

#include <array2d.h>
#include <list>

#include "fftutils.h"

class DumpImageAspect 
{

  typedef std::list<const char*> PatternList;
  PatternList patterns;
  
public:
  void dump( const char *fileName, const pfs::Array2D *data, const char *channelName, ... );
  void dump( const char *fileName, BidomainArray2D *data, const char *channelName, ... );
  void dumpFrequency( const char *fileName, BidomainArray2D *data, ... );
  
  /**
   * Dump all images that match the pattern
   */
  void addIncludePattern( char *pattern );
  
  
};

extern DumpImageAspect *dumpImageAspect;

void initializeDumpImageAspect();


#endif
