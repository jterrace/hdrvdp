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
 * $Id: csf_filter.h,v 1.2 2010/12/25 16:56:06 rafm Exp $
 */

#ifndef CSF_FILTER_H
#define CSF_FILTER_H

#include "fftw_array2d.h"
#include "fftutils.h"

/**
 * Base class for all CSF filters
 */
class CSFFilter 
{
public:
  virtual void process( BidomainArray2D *in, BidomainArray2D *out,
	  BidomainArray2D *adaptationMap ) = 0;

  virtual const char *getName() const = 0;  
};


/**
 * This class defines viewing conditions for a display (distance,
 * resolution, display size, etc.)
 */
class ViewingConditions
{
public:
  
  int xResolution, yResolution;
  float displayWidth, displayHeight;
  float minDistance, maxDistance;
  float pixelsPerDegree;
    
  ViewingConditions( int xResolution, int yResolution,
    float displayWidth, float displayHeight,
    float minDistance, float maxDistance ) :
      xResolution( xResolution ), yResolution( yResolution ),
      displayWidth( displayWidth ), displayHeight( displayHeight ),
      minDistance( minDistance ), maxDistance( maxDistance ),
      pixelsPerDegree( -1 )
    {
    }

  ViewingConditions( float pixelsPerDegree,
    float minDistance, float maxDistance ) :
      xResolution( -1 ), yResolution( -1 ),
	displayWidth( -1 ), displayHeight( -1 ),
	minDistance( minDistance ), maxDistance( maxDistance ),
	pixelsPerDegree( pixelsPerDegree )
    {
    }

  float getPixelsPerDegree( float observerDistance ) const;
    
};


#define MAX_ADAPTATION_LEVELS 6

/**
 * CSF Filtering that accounts for local adaptation of the eye.
 * Used in HDRVDP.
 */
class MultiAdaptationCSF : public CSFFilter
{
  float adaptationLevels[MAX_ADAPTATION_LEVELS];
  int adaptationLevelsCount;
  pfs::Array2DImpl **filters;
  
public:
  MultiAdaptationCSF( int cols, int rows, const ViewingConditions &viewCond,
    const pfs::Array2D *otf = NULL );
  
  void process( BidomainArray2D *in, BidomainArray2D *out,
    BidomainArray2D *adaptationMap  );

  const char *getName() const
    {
      return "Local-adaptation CSF";
    }
  
  
  ~MultiAdaptationCSF();
  
};

/**
 * CSF Filter used in the original Daly's VDP.
 * It assumes a single level of luminance adaptation.
 */
class VDPCSF : public CSFFilter
{
  pfs::Array2DImpl *filter;
  
public:
  VDPCSF( int cols, int rows, const ViewingConditions &viewCond, const float Y_adapt );
  
  void process( BidomainArray2D *in, BidomainArray2D *out,
    BidomainArray2D *adaptationMap );

  const char *getName() const
    {
      return "VDP original CSF filter";
    }
  
  
  ~VDPCSF();
  
};


#endif
