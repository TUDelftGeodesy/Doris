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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/filtering.hh,v $
 * $Revision: 3.9 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Declaration of routines for filtering 
 ****************************************************************/


#ifndef FILTERING_H                     // guard
#define FILTERING_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // typedefs
#include "readinput.hh"                 // input structs
#include "matrixbk.hh"                  // my matrix class
#include "slcimage.hh"                  // my slc image class
#include "orbitbk.hh"                   // my orbit class
#include "productinfo.hh"               // my 'products' class




// ______ Adaptive range filtering, estimate local fringe frequency _____
void rangefilter(
        const input_gen       &generalinput,
        const slcimage        &master,
        const slcimage        &slave,
        const productinfo     &interferogram,
        const input_filtrange &filtrangeinput);

// ______ Iteratively called by rangefilter ______
void rfilterblock(
        matrix<complr4>       &MASTER,          // updated
        matrix<complr4>       &SLAVE,           // updated
        int32                  nlmean,
        real8                  SNRthreshold,
        real8                  RSR,             // in MHz?
        real8                  RBW,             // in MHz?
        real8                  hammingalpha,    
        int32                  oversamplefactor,   // pow2
        bool                   docorrectcorrel,
        real8                 &meanSNR,            // returned
        real8                 &percentnotfiltered);// returned

// ______ Range filtering based on orbit _____
void rangefiltporbits(
        const input_gen       &generalinput,
        const input_filtrange &filtrangeinput,
        const input_ell       &ellips,
        const slcimage        &master,
        const slcimage        &slave,
              orbit           &masterorbit,
              orbit           &slaveorbit
        );

// ====== Adaptive phase filtering, weight spectrum, goldstein ======
void phasefilter(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput,
	const productinfo     &coherence);

// ______ goldstein filter per buffer ______
matrix<complr4> goldstein(
        const matrix<complr4> &CINT,
        const real4            ALPHA,
        const int32            OVERLAP,
        const matrix<real4>   &smoothkernel);

// ______ modified goldstein filter per buffer ______
matrix<complr4> modgoldstein(
        const matrix<complr4> &CINT,
        const matrix<real4>   &COH,
        const int32            OVERLAP,
        const matrix<real4>   &smoothkernel);

// ______ Smooth by spatial circular averaging over 2N+1 box ______
matrix<real4> smooth(
        const matrix<real4> &A,
        int32 halblocksize);

// ______ Smooth by FFT's, KERNEL2D once computed for ______
// ______ recursive calls to smooth routine, save time ______
matrix<real4> smooth(
        const matrix<real4>   &A,
        const matrix<complr4> &KERNEL2D);

// ====== Spatial convolution filter per buffer ======
void spatialphasefilt(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput);

// ______ Spatial convolution filter per block _____
matrix<complr4> convbuffer(
        const matrix<complr4> &CINT,
        const matrix<complr4> &KERNEL2D,
        const int32            OVERLAP);

// ====== phase filtering, spectral ======
void phasefilterspectral(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput);

// ______ spectral filter per buffer ______
matrix<complr4> spectralfilt(
        const matrix<complr4> &CINT,
        const matrix<real4>   &KERNEL2D,
        const int32            OVERLAP);

// ====== azimuth filter ======
void azimuthfilter(
        const input_gen       &generalinput,
        const input_filtazi   &filtaziinput,
        slcimage        &master,
        slcimage        &slave);

matrix<complr4> blockazifilt(
        const matrix<complr4> &SLCIMAGE,
        const slcimage        &master,          // PRF, BW, fd0
        const slcimage        &slave,           // PRF, BW, fd0
        const real8            HAMMING);



#endif // FILTERING_H



