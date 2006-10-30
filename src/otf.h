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
 * $Id: otf.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef OTF_H
#define OTF_H

#include "fftutils.h"
#include "csf_filter.h"

class OTF
{
public:
  virtual void process( BidomainArray2D *in, BidomainArray2D *out ) = 0;
  virtual const char *getName() = 0;
  virtual const pfs::Array2D *getFilter()
    {
      return NULL;
    }  

  virtual ~OTF() 
    {
    }
};


class NormannBaxterOTF : public OTF
{
  BidomainArray2D filter;
public:
  NormannBaxterOTF( const ViewingConditions &viewCond, int cols, int rows );
  void process( BidomainArray2D *in, BidomainArray2D *out );
  
  const char *getName() 
    {
      return "Normann & Baxter's OTF model";
    }
  
};

class WestheimerOTF : public OTF
{
  BidomainArray2D filter;
public:
  WestheimerOTF( const ViewingConditions &viewCond, int cols, int rows );
  void process( BidomainArray2D *in, BidomainArray2D *out );
  
  const char *getName() 
    {
      return "Westheimer's OTF model (used in Sarnoff Visual Discrimination Model)";
    }
  
};

class SpencerOTF : public OTF
{
  BidomainArray2D filter;
public:
  SpencerOTF( const ViewingConditions &viewCond, int cols, int rows, float adaptationLuminance );
  void process( BidomainArray2D *in, BidomainArray2D *out );
  
  const char *getName() 
    {
      return "Spencer'95 OTF model";
    }
  
};


class MarimontOTF : public OTF
{
  pfs::Array2D *filter;
public:
  MarimontOTF( const ViewingConditions &viewCond, int cols, int rows, float adaptationLuminance );
  void process( BidomainArray2D *in, BidomainArray2D *out );
  ~MarimontOTF();
  
  const char *getName() 
    {
      return "Marimont and Wandell OTF model";
    }
  const pfs::Array2D *getFilter()
    {
      return filter;
    }  
};


class DeeleyOTF : public OTF
{
  pfs::Array2D *filter;
public:
  DeeleyOTF( const ViewingConditions &viewCond, int cols, int rows, float adaptationLuminance );
  void process( BidomainArray2D *in, BidomainArray2D *out );
  ~DeeleyOTF();
  
  const char *getName() 
    {
      return "Deeley OTF model";
    }

  const pfs::Array2D *getFilter()
    {
      return filter;
    }
};

#endif
