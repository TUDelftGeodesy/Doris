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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/coregistration.hh,v $
 * $Revision: 3.13 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Routines for coregistration.
 ****************************************************************/


#ifndef COREGISTRATION_H
#define COREGISTRATION_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // typedefs
#include "readinput.hh"                 // input structs
#include "orbitbk.hh"                   // my orbit class, includes my matrix class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class
#include "bk_baseline.hh"               // my 'baseline' class



// ______ Prototypes ______
#ifdef __USE_VECLIB_LIBRARY__
// ______This typedef is needed for output of veclib routine cdot?______
typedef struct {real4 re,im;} STUPID_cr4;

extern "C" { void  cgemv (char*, int32*, int32*, complr4*, complr4*,  int32*,
                         complr4*, int32*,
                         complr4*, complr4*, int32*, int32); }

extern "C" { STUPID_cr4 cdotu (int32*, complr4*, int32*, complr4*, int32*); }
#endif



// ====== Prototypes ======
// ______ Compute coarse coregistration ______
void coarseporbit(
        const input_ell         &ell,
        const slcimage  &master,
        const slcimage  &slave,
        orbit           &masterorbit,
        orbit           &slaveorbit,
        const BASELINE  &baseline);


// ______ Coarse coregistration ______
void coarsecorrel(
        const input_coarsecorr  &input, 
        const slcimage  &minfo,
        const slcimage  &sinfo);


// ______ Corr by fft ______
void coarsecorrelfft(
        const input_coarsecorr  &input, 
        const slcimage  &minfo,
        const slcimage  &sinfo);

// ______ Sim. Amplitude coregistration (magspace) ______
void mtiming_correl(
        // const input_coarsecorr  &input,
        const input_mtiming  &input,
        const slcimage       &minfo,      // master
        const productinfo    &sinfo);     // simulated amplitude

// ______ Sim. Amplitude coregistration (magfft) ______
void mtiming_correlfft(
        const input_mtiming &input,
        const slcimage      &minfo,
        const productinfo   &sinfo);

// ______ Corr by fft ______
// superseded by coherencefft with factor=1
//real4 corrfft( 
//      const matrix<real4>     &magnitudeMaster,
//      const matrix<real4>     &magnitudeMask,
//      real4                   &offsetL,
//      real4                   &offsetP);


// ______ Distribute nW windows over win ______
//matrix<uint> distributepoints(
// [FvL] for correct folding of points outside overlap window when inserted by file
matrix<int> distributepoints(
        real4                     numberofpoints,
        const window             &win);


// ______ Estimate offset based on consistency ______
void getoffset(
        const matrix<real4>     &Result,// not sorted on output
        int32                   &offsetLines,
        int32                   &offsetPixels);

// ______ Estimate offset based on consistency ______
void getmodeoffset(
        const matrix<real4>     &Result,// not sorted on output
        int32                   &offsetLines,
        int32                   &offsetPixels);


// ______ Fine coregistration ______
//void finecoreg(
//        const input_fine        &fineinput,
 //       const slcimage  &minfo,
 //       const slcimage  &sinfo);

void finecoreg(
        const input_fine &fineinput,
        const slcimage   &minfo,
        const slcimage   &sinfo,
        const input_ell &ell,
        orbit           &masterorbit,  // cannot be const for spline
        orbit           &slaveorbit,   // cannot be const for spline
        const BASELINE  &baseline);

// ______ Correlation with FFT MCC______
real4 coherencefft(
        const matrix<complr4>   &Master,
        const matrix<complr4>   &Mask,
        const int32 factor,              // ovs factor (1 for not)
        const int32 AccL,                // search window to oversample
        const int32 AccP,                // search window to oversample
        real4                   &offsetL,
        real4                   &offsetP);


// ______ Correlation with FFT ______
real4 crosscorrelate(
        const matrix<complr4>   &Master,
        const matrix<complr4>   &Mask,
        const int32 factor,              // ovs factor (1 for not)
        const int32 AccL,                // search window to oversample
        const int32 AccP,                // search window to oversample
        real4                   &offsetL,
        real4                   &offsetP);

// ______ Correlation with FFT ______
real4 intensity(
        const matrix<complr4>   &Master,
        const matrix<complr4>   &Mask,
        const int32 factor,              // ovs factor (1 for not)
        const int32 AccL,                // search window to oversample
        const int32 AccP,                // search window to oversample
        real4                   &offsetL,
        real4                   &offsetP);


// ______ Correlation in space domain ______
real4 coherencespace(
        const input_fine        &fineinput, 
        const matrix<complr4>   &Master,
        const matrix<complr4>   &Mask,
        real4                   &offsetL,
        real4                   &offsetP);


// ______ Compute coregistration parameters ______
//      const window            &originalmaster,
void coregpm(
        const slcimage          &master,
    const slcimage      &slave, //[FvL]
        const char              *i_resfile,
        const input_coregpm     &coregpminput,
    const int16             &demassist); //[FvL]
        //const uint             oversamplingsfactor);


// ______ Read observations from file ______
matrix<real4> getofffile(
        const char              *file,
        real4                    threshold);


// ______ Resample slave ______
void resample(
        const input_gen         &generalinput,
        const input_resample    &resampleinput,
        const slcimage          &master,
        const slcimage          &slave,
        const matrix<real8>     &cpmL,
        const matrix<real8>     &cpmP,
        const int16             &demassist,
        const matrix<real8>     &minMaxL,//[MCC]
        const matrix<real8>     &minMaxP);//[MCC]

// ______ Compute master-slave timing error ______
void ms_timing_error(
        const slcimage          &master,
        const char              *i_resfile,
        const input_reltiming   &timinginput,
        int32                   &coarse_orbit_offsetL,
        int32                   &coarse_orbit_offsetP);

// ______Interpolation kernals______
matrix<real4> cc4(const  matrix<real4> &x);
matrix<real4> cc6(const  matrix<real4> &x);
matrix<real4> ts6(const  matrix<real4> &x);
matrix<real4> ts8(const  matrix<real4> &x);
matrix<real4> ts16(const matrix<real4> &x);
matrix<real4> rect(const matrix<real4> &x);
matrix<real4> tri(const  matrix<real4> &x);
// ___ knab: oversampling factor of signal CHI, number of points N ___
matrix<real4> knab(const matrix<real4> &x, const real4 CHI, const int32 N);
// ___ raised cosine: oversampling factor of signal CHI, number of points N ___
matrix<real4> rc_kernel(const matrix<real4> &x, const real4 CHI, const int32 N);

#endif // COREGISTRATION_H

