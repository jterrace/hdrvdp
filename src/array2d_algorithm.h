/**
 * @brief operations on 2D arrays
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
 * $Id: array2d_algorithm.h,v 1.1 2006/10/30 19:07:48 rafm Exp $
 */

#ifndef ARRAY2D_ALGORITHM
#define ARRAY2D_ALGORITHM

#include <array2d.h>

float average( const pfs::Array2D *array );
float logAverage( const pfs::Array2D *array );
float max(  const pfs::Array2D *array );
void normalizeToSum( pfs::Array2D *array );
void normalizeToMax( pfs::Array2D *array );

/**
 * Clamp values in the array to positive integers.  Negative, zero,
 * NaN and -Inf are set to the smallest positive value. +Inf is set to
 * the largest value, which is not +Inf.
 */
void clampToPositive( const pfs::Array2D *from, pfs::Array2D *to );



#endif
