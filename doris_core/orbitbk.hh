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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/orbitbk.hh,v $    *
 * $Revision: 3.13 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of orbit class.                               *
 * - orbit interpolation.                                       * 
 * - baseline estimation.                                       *
 * - utility dumping etc.                                       *
 *                                                              *
 * Compilation with: g++ -D__TESTMAIN__ orbit.cc -o testorbit   *
 * creates a standalone executable for testing of the functions *
 * of the orbit class. Please see below: "main program".        *
 ****************************************************************/


#ifndef ORBITBK_H               // guard
#define ORBITBK_H


using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrixbk.hh"          // my matrix class, only 2b included once...
#include "slcimage.hh"          // my slc image class
#include "constants.hh"         // global constants
#include "readinput.hh"         // input structs




// ====== Some inline functions ======
// make these private..., but doesn't seem to work ok that way, so do it like this
// ______ The setofequation for t2xyz ______
inline real8 eq1_doppler(const cn &velocity, const cn &dsat_P)
  {return velocity.in(dsat_P);}
inline real8 eq2_range(const cn &dsat_P, const real8 &rangetime)
  {return dsat_P.in(dsat_P)-sqr(SOL*rangetime);}
inline real8 eq3_ellipsoid(const cn &P, const real8 &semimajora, const real8 &semiminorb)
  {return ((sqr(P.x)+sqr(P.y))/sqr(semimajora))+sqr(P.z/semiminorb)-1.0;}

// ______ The partial to Pxyz of the 3 eq. for t2xyz ______
// they are best evaluated within calling function:
// dE1/dx  dE1/dy  dE1/dz   vel.x                ~y  ~z
// dE2/dx  dE2/dy  dE2/dz = 2*(pos.x-sat.x)      ~y  ~z
// dE3/dx  dE3/dy  dE3/dz   (2*pos.x)/sqr(ell.a) ~y  (2*pos.z)/sqr(ell.b);

// ______ The partial to time of dopplereq. for xyz2t ______
inline real8 eq1_doppler_dt(const cn &dsat_P, const cn &velocity, const cn &accerelation)
  {return accerelation.in(dsat_P)-sqr(velocity.x)-sqr(velocity.y)-sqr(velocity.z);}




// ====== The orbit class itself ======
class orbit
  {
  // ====== Private data, functions ======
  private:
    // ______ Data ______
    //char system[80];                  // WGS84, GRS80, Bessel, etc.
    int16 interp_method;                // 0: cubic spline, x: polyfit(x)
    int32 numberofpoints;               // ...
    int16 orbvector_type;              // 30: positional vectors only  31: position + velocity vectors [MA]
    int32 klo;                          // index in timevector correct
    int32 khi;                          // +part piecewise polynomial
    matrix<real8> time;                 // vector (numberofpoints,1)
    matrix<real8> data_x;               // vector (numberofpoints,1)
    matrix<real8> data_y;               // vector (numberofpoints,1)
    matrix<real8> data_z;               // vector (numberofpoints,1)
    matrix<real8> coef_x;               // vector (numberofpoints,1)
    matrix<real8> coef_y;               // vector (numberofpoints,1)
    matrix<real8> coef_z;               // vector (numberofpoints,1)
    matrix<real8> data_xv;               // vector (numberofpoints,1)
    matrix<real8> data_yv;               // vector (numberofpoints,1)
    matrix<real8> data_zv;               // vector (numberofpoints,1)
    matrix<real8> coef_xv;               // vector (numberofpoints,1)
    matrix<real8> coef_yv;               // vector (numberofpoints,1)
    matrix<real8> coef_zv;               // vector (numberofpoints,1)

    // ______ Functions ______
    void computecoefficients();         // pp coeffs.
    //void computecoefficientsV();         // pp coeffs.
    void getklokhi(real8 t);            // piecewize polynomial indices


  // ====== Public data, functions ======
  public:
    //default method is polynomial degree=npoints, but max degree=5.
    //orbit() {interp_method=ORB_DEFAULT;numberofpoints=0;}// constructor
    orbit() {interp_method=ORB_DEFAULT;numberofpoints=0;orbvector_type=ORB_PRM_POS;}// constructor
    int32 npoints()     {return numberofpoints;}
    void set_interp_method(int16 m){interp_method=m;}
    void set_orbvector_type(int16 t){orbvector_type=t;}

    // ______ Read from file&store t,x,y,z; compute coefficients ______ 
    // ______ function initialize should be called first! ______
    void initialize (const char *file); // do all...
    bool is_initialized() {return (numberofpoints!=0) ? true : false;}// constructor
    void showdata   ();                 // debugging...

    // ______ Modification from baseline error polynomial [HB] ______
    matrix<real8> modify(const matrix<real8> &coeff,
			 orbit               &referenceorbit,
			 const real8         tshift,
			 const real8         tmin,
			 const real8         tmax);
        
    // ====== Interpolation ======
    //cn getxyz       (real8 time) const; // not const, klo/hi updated
    cn getxyz       (real8 time);
    cn getxyzdot    (real8 time);
    cn getxyzddot   (real8 time);

    /*
    // ====== Conversion between coordinate systems ======
    // ______ radar coordinate line/pixel to xyz on ellipsoid ______
    friend int32 lp2xyz(
        real8            line,
        real8            pixel,
        const input_ell &ell,
        const slcimage  &image,
        orbit           &orb,
        cn              &returnpos,
        int32            MAXITER=10,            // [.] defaults
        real8            CRITERPOS=1e-6);       // [m]

    // ______ xyz cartesian on ellipsoid to orbital coord. ______
    friend int32 xyz2orb(
        cn              &returnpossat,
        const slcimage  &image,
        orbit           &orb,
        const cn        &pointonellips,
        int32            MAXITER=10,            // [.] defaults
        real8            CRITERTIM=1e-10);      // [s]

    // ______ xyz cartesian on ellipsoid to azimuth/range time ______
    friend int32 xyz2t(
        real8           &returntazi,            // azimuth
        real8           &returntran,            // and range time
        const slcimage  &image,
        orbit           &orb,
        const cn        &pos,
        int32            MAXITER=10,            // [.] defaults
        real8            CRITERTIM=1e-10);      // [s]

    // ______ convert xyz-ellipsoid to line/pixel ______
    friend int32 xyz2lp(
        real8           &returnline,
        real8           &returnpixel,
        const slcimage  &image,
        orbit           &orb,
        const cn        &pos,
        int32            MAXITER=10,             // defaults
        real8            CRITERTIM=1e-10);       // seconds

    // ______ Convert ellipsoid to radar coordinates ______
    friend int32 ell2lp(
        real8           &returnline,
        real8           &returnpixel,
        const input_ell &ell,
        const slcimage  &image,
        orbit           &orb,
        real8            phi,
        real8            lambda,
        real8            height,
        int32            MAXITER=10,
        real8            CRITERTIM=1e-10);  // seconds

    // ______ Convert radar coordinates to ellipsoidal coordinates ______
    friend int32 lp2ell(
        real8            line,
        real8            pixel,
        const input_ell &ell,
        const slcimage  &image,
        orbit           &orb,
        real8           &returnphi,
        real8           &returnlambda,
        real8           &returnheight,
        int32            MAXITER=10,
        real8            CRITERPOS=1e-6);   // meter
    */
    // ====== Information/debugging ======
    // ______ dump computed coeffs, interpolated orbit, etc. ______
    void dumporbit(
        const input_pr_orbits &inputorb,
        const int16      ID);

    // ______ compute baseline on a grid ______
    friend void compbaseline(
        const input_gen &generalinput,
        const slcimage  &master,
        const slcimage  &slave,
        orbit           &masterorbit,
        orbit           &slaveorbit);

  }; // END class orbit

int32 lp2xyz(
    real8            line,
    real8            pixel,
    const input_ell &ell,
    const slcimage  &image,
    orbit           &orb,
    cn              &returnpos,
    int32            MAXITER=10,            // [.] defaults
    real8            CRITERPOS=1e-6);       // [m]

// ______ xyz cartesian on ellipsoid to orbital coord. ______
int32 xyz2orb(
    cn              &returnpossat,
    const slcimage  &image,
    orbit           &orb,
    const cn        &pointonellips,
    int32            MAXITER=10,            // [.] defaults
    real8            CRITERTIM=1e-10);      // [s]

// ______ xyz cartesian on ellipsoid to azimuth/range time ______
int32 xyz2t(
    real8           &returntazi,            // azimuth
    real8           &returntran,            // and range time
    const slcimage  &image,
    orbit           &orb,
    const cn        &pos,
    int32            MAXITER=10,            // [.] defaults
    real8            CRITERTIM=1e-10);      // [s]

// ______ convert xyz-ellipsoid to line/pixel ______
int32 xyz2lp(
    real8           &returnline,
    real8           &returnpixel,
    const slcimage  &image,
    orbit           &orb,
    const cn        &pos,
    int32            MAXITER=10,             // defaults
    real8            CRITERTIM=1e-10);       // seconds

// ______ Convert ellipsoid to radar coordinates ______
int32 ell2lp(
    real8           &returnline,
    real8           &returnpixel,
    const input_ell &ell,
    const slcimage  &image,
    orbit           &orb,
    real8            phi,
    real8            lambda,
    real8            height,
    int32            MAXITER=10,
    real8            CRITERTIM=1e-10);  // seconds

// ______ Convert radar coordinates to ellipsoidal coordinates ______
int32 lp2ell(
    real8            line,
    real8            pixel,
    const input_ell &ell,
    const slcimage  &image,
    orbit           &orb,
    real8           &returnphi,
    real8           &returnlambda,
    real8           &returnheight,
    int32            MAXITER=10,
    real8            CRITERPOS=1e-6);   // meter

#endif // ORBITBK_H guard


