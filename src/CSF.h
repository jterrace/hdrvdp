/**
 * @brief CSF from the original Daly's VDP
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
 * $Id: CSF.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef CSF_H
#define CSF_H

#define CSF_NULL              0
#define CSF_DALY              1
#define CSF_MANNOS            2
#define CSF_MARTIN            3
#define CSF_NULL_KEEP_DALY_AMPL 4
#define CSF_SPATIO_VELOCITY   5
#define CSF_DALY_NORMALIZED   6
#define CSF_DALY_MULTIADAPTATION 7

#include <array2d.h>


void createCSFFilter( pfs::Array2D *filter,
  float adaptationLuminance,
  float imageVDSize, float pixelsPerDegree, float observerDistance,
  int CSF_type );



#endif
