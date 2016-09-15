/*
 * Copyright (c) 1999-2009 Delft University of Technology, The Netherlands
 *
 * This file is part of Doris, the Delft o-o radar interferometric software.
 *
 * Doris program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Doris is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 */
/****************************************************************
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/referencephase.hh,v $
 * $Revision: 3.8 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Routines for reference phase computation
 * - flatearth
 * - radarcodedem
 * - interpolation routines delaunay?
 ****************************************************************/


#ifndef REFERENCEPHASE_H
#define REFERENCEPHASE_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // typedefs
#include "readinput.hh"                 // input structs
#include "orbitbk.hh"                   // my orbit class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class




// ====== Prototypes ======
// ______ Compute reference phase ______
void flatearth(
        const input_comprefpha  &flatearthinput,
        const input_ell         &ellips,
        const slcimage  &master,
        const slcimage  &slave,
        const productinfo       &interferogram,
        orbit                   &masterorbit,
        orbit                   &slaveorbit);

// ______ DEM assisted coregistration ______
void demassist(
        const input_gen         &generalinput,
        const input_ell         &ellips,
        const input_demassist   &demassistinput,
        const slcimage          &master,
        const slcimage          &slave,
        orbit                   &masterorbit,
        orbit                   &slaveorbit);

// ______ Compute reference phase ______
void radarcodedem(
        const input_gen         &generalinput,
        const input_ell         &ellips,
        const input_comprefdem  &comprefdeminput,
        const slcimage          &master,
        const slcimage          &slave,
        const productinfo       &interferogram,
        orbit                   &masterorbit,
        orbit                   &slaveorbit);

// ______ Get indices of corners DEM ________________
void getcorners(
                const real8            &l0,
                const real8            &lN,
                const real8            &p0,
                const real8            &pN,
                const real8            &extralat,
                const real8            &extralong,
                const real8            &lat0,
                const real8            &long0,
                const real8            &DEMdeltalat,
                const real8            &DEMdeltalong,
                const int32            &Nlatpixels,
                const int32            &Nlongpixels,
                const input_ell        &ellips,
                const slcimage         &master,
                orbit                  &masterorbit,
                real8            &phimin,
                real8            &phimax,
                real8            &lambdamin,
                real8            &lambdamax,
                int32            &indexphi0DEM,
                int32            &indexphiNDEM,
                int32            &indexlambda0DEM,
                int32            &indexlambdaNDEM);


// _______ griddatalinear ______________________
void griddatalinear(
                    const matrix<real8> &x_in,          // line values
                    const matrix<real8> &y_in,          // pixel values
                    const matrix<real8> &z_in,          // z. hei, amp in dem coords
                    const real8         &x_min,
                    const real8         &x_max,
                    const real8         &y_min,         // pixel extend min
                    const real8         &y_max,         // pixel extend max
                    const int32         &x_inc,         // increment based on mlL
                    const int32         &y_inc,
                    const real8         &r_az_ratio,
                    const real8         &offset,
                    const real8         &NODATA,
                    matrix<real8>       &grd
                    );

// _______ pointintriangle ______________________
int pointintriangle(
                    double *xt, 
                    double *yt,
                    double x ,
                    double y);


#define rint(x) (floor((x)+0.5))
#define irint(x) ((int)rint(x))
#define x_to_i(x,x0,dx,off,nx) (irint(((((x) - (x0)) / (dx)) - (off))))
#define y_to_j(y,y0,dy,off,ny) (irint(((((y) - (y0)) / (dy)) - (off))))
#define i_to_x(i,x0,x1,dx,off,nx) ((x0) + (i) * (dx) + (off))
#define j_to_y(j,y0,y1,dy,off,ny) ((y0) + (j) * (dy) + (off))

//  // ______ Return true if P(l,p) is inside triangle P1P2P3 ______
//  // ______ newtriangle is true if last call was for same P1P2P3 ______
//  bool pointintriangle(
//          real4 l,
//          real4 p,
//          real4 l1,
//          real4 p1,
//          real4 l2,
//          real4 p2,
//          real4 l3,
//          real4 p3,
//          bool newtriangle);

// ______ Return true if P(l,p) is inside triangle P1P2P3 ______
// ______ newtriangle is true if last call was for same P1P2P3 ______
// ______ if not real8 then totally wrong... (or eps = AP1P2P3/1e6) ______
bool pointintriangle(
        real8 l,
        real8 p,
        real8 l1,
        real8 p1,
        real8 l2,
        real8 p2,
        real8 l3,
        real8 p3,
        bool newtriangle);




#endif // REFERENCEPHASE_H



