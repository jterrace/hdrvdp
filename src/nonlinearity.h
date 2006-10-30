/**
 * @brief photoreceptor nonlinearity
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
 * $Id: nonlinearity.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef NONLINEARITY_H
#define NONLINEARITY_H

#include <iostream>
#include <array2d.h>

#include "fftutils.h"

/**
 * Photoreceptor nonlinearity
 */
class PhotoreceptorNonlin 
{
public:
  virtual void process( BidomainArray2D *in, BidomainArray2D *out ) = 0;
  virtual const char *getName() = 0;
};

class JNDScaledNonlinearity : public PhotoreceptorNonlin
{
  float peakSensitivity;
public:
  JNDScaledNonlinearity( float peakSensitivity ) : peakSensitivity(peakSensitivity)
    {
    }
  
  void process( BidomainArray2D *in, BidomainArray2D *out );
  
  const char *getName() 
    {
      return "JND Scaling";
    }
  
};

/**
 * This is a simplified version of a Normann & Baxter model of
 * photoreceptor used in the original Daly's VDP.
 */
class VDPNonlinearity : public PhotoreceptorNonlin
{
public:
  void process( BidomainArray2D *in, BidomainArray2D *out );
  
  const char *getName() 
    {
      return "Simplified Normann & Baxter";
    }
  
};

#endif
