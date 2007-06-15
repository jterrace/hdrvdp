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
 * $Id: CSF.cpp,v 1.2 2007/06/15 15:19:15 rafm Exp $
 */

#include <math.h>
#include <iostream>

#include "CSF.h"

void createCSFFilter( pfs::Array2D *filter,
  float adaptationLuminance,
  float imageVDSize, float pixelsPerDegree, float observerDistance,
  int CSF_type )
{

  // Sensitivity calibration constant (from Daly's paper: 250)
  const float P = 250.f;
  
// Frequency domain
  int filterWidth=filter->getCols();
  int filterHeight=filter->getRows();
  float eps = 0.9;
  float i_pow2 = imageVDSize;
  float l = adaptationLuminance;
  float A = 0.801 * powf(1 + 0.7 / l, -0.20);
  float B = 0.3 * powf(1 + 100 / l, 0.15);

  std::cerr << "l_adapt: " << l << " i_pow2: " << i_pow2 << " pix_per_deg: " << pixelsPerDegree << "\n";

  float r_a = 0.856 * powf(observerDistance, 0.14);
  float e = 0.0;     // eccentricity in visual degrees
  float r_e = 1.0 / (1.0 + 0.24 * e);
  float r_a_r_e = r_a * r_e;

  float pix_per_deg = pixelsPerDegree;
  
  int x, y;
  float x_norm = 0.5 / filterWidth, y_norm = 0.5 / filterHeight;
  float dx, dy;
  float theta, r_theta;
  float ro;
  float S1, S2;
  for (y = 0, dy = 0.0; y < filterHeight; y ++, dy += y_norm)
  {
    for (x = 0, dx = 0.0; x < filterWidth; x ++, dx += x_norm)
    {
      // TODO: ro should depend on image size - cycles in the image
      ro = pix_per_deg * sqrtf((float)(dx * dx + dy * dy));

      switch (CSF_type)
      {
      case CSF_DALY:
      case CSF_DALY_NORMALIZED:
      case CSF_DALY_MULTIADAPTATION:
      {
        float B1 = B*eps*ro;
        if (B1 > 50.0  ||  (x == 0  &&  y == 0))
          S1 = 0.0;
        else
          S1 = powf(powf(3.23 * powf(ro*ro * i_pow2, -0.3), 5.0) + 1.0, -0.2)*
            A*eps*ro * expf(-B1) * sqrtf(1 + 0.06*expf(B1));
      }

      if (dx == 0.0)
      {
        if (dy == 0.0)
          theta = 0.0;
        else
          theta = 0.5 * M_PI;
      }
      else
        theta = atan2f(dy, dx);
      
      {  
        const float ob=0.78;
        r_theta = (1-ob)/2 * cosf(4.0 * theta) + (1+ob)/2;
        ro = ro / (r_a_r_e * r_theta);
      }
      

      {
        float B1 = B*eps*ro;
        if (B1 > 50.0  ||  (x == 0  &&  y == 0))
          S2 = 0.0;
        else
          S2 = powf(powf(3.23 * powf(ro*ro * i_pow2, -0.3), 5.0) + 1.0, -0.2)*
            A*eps*ro * expf(-B1) * sqrtf(1 + 0.06*expf(B1));
      }
      
      (*filter)(x,y) = (S1 > S2 ? S2 : S1) * P;
      break;

      case CSF_MANNOS:
        (*filter)(x,y) = P * (0.0499 + 0.2964 * ro) * expf(-powf(0.114 * ro, 1.1));
        break;
      case CSF_MARTIN:
        if (ro > 0.0)
        {
          float X = log10f(ro);
          float Y = log10f(adaptationLuminance);
          float X2 = X*X;
          float Y2 = Y*Y;
          float Y3 = Y*Y2;
          float logCSF = 2.094 + 0.6019*X - 0.9730*X2 + 0.2218*Y + 0.6828*X*Y -
            0.3656*X2*Y - 0.10258*Y2 + 0.07854*X*Y2 - 0.05234*X2*Y2 -
            0.01557*Y3 - 0.02197*X*Y3 + 0.03920*X2*Y3;

          (*filter)(x,y) = powf(10.0, logCSF);
        }
        break;
      }
    }
  }
  // Estimating DC filter component
  //(*filter)(0,0)=(*filter)(0,1);
  (*filter)(0,0)= 0; 

  if( CSF_type == CSF_DALY_NORMALIZED || CSF_type == CSF_DALY_MULTIADAPTATION) {
    // Normalize CSF, so that maximum sensitivity == 1

    float maxV = -99999;
    for (y = 0; y < filterHeight; y++) {
      for (x = 0; x < filterWidth; x++) {
        if( (*filter)(x,y) > maxV ) maxV = (*filter)(x,y);
      }
    }

    for( int i = 0; i < filterWidth*filterHeight; i++ )
        (*filter)(i) /= maxV;    
  }
  
}
