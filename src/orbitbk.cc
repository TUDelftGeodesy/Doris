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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/orbitbk.cc,v $    *
 * $Revision: 3.20 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of orbit class.                               *
 * - orbit interpolation (spline and polynomial).               * 
 * - baseline estimation.                                       *
 * - utility dumping etc.                                       *
 *                                                              *
 * Compilation with: make testorbit                             *
 *  g++ -D__TESTMAIN__ orbit.cc -o testorbit                    *
 * creates a standalone executable for testing of the functions *
 * of the orbit class. Please see below: "main program".        *
 ****************************************************************/


#include "orbitbk.hh"                   // declarations, matrix class
#include "constants.hh"                 // global constants
#include "ioroutines.hh"                // error messages
#include "utilities.hh"                 // solve33
#include "exceptions.hh"                 // my exceptions class
#include "estorbit.hh"                  // probably only temporary [HB]

#include <cstdio>                       // some compilers, remove function
#include <strstream>                    // for memory stream
#include <iomanip>                      // setw


// ====== Global variables, declared extern in constants.h ======
// === Declared extern in constants.h, set in main to be able to use them ===
// === in readinput.cc they are set via the input file ======================
// [MA] copied from processor.cc to build orbit as a standalone binary
#ifdef __TESTMAIN__
bk_messages TRACE;
bk_messages DEBUG;
bk_messages INFO;
bk_messages PROGRESS;
bk_messages WARNING;
bk_messages ERROR;
bk_messages matDEBUG;
bk_messages matERROR;
#endif

// ====== Prototypes in this file ======
matrix<real8> splineinterpol(
        const matrix<real8> &time,
        const matrix<real8> &data);
matrix<real8> polyfit(
        const matrix<real8> &time,
        const matrix<real8> &y, 
        const int32 DEGREE);



/****************************************************************
 *    orbit::initialize                                         *
 * Fills data vectors, computes coefficients for interpolation  *
 * Expects time sorted data in file after string:               *
 * "NUMBER_OF_DATAPOINTS:   3"                                  *
 *  t1.xxx x1.xxx y1.xxx z1.xxx                                 *
 *  t2.xxx x2.xxx y2.xxx z2.xxx                                 *
 *  t3.xxx x3.xxx y3.xxx z3.xxx                                 *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 #%// BK 17-Jul-2000: changed to class                          *
 ****************************************************************/
void orbit::initialize(const char* file)
  {
  TRACE_FUNCTION("orbit::initialize (BK 17-Jul-2000)")
  char           dummyline[4*ONE27];
  char           word[EIGHTY];
  bool           foundsection = false;
  register int32 i;

  // ______ Open file ______
  ifstream infile(file, ios::in);
  bk_assert(infile,file,__FILE__,__LINE__);
  numberofpoints = 0;

  // ====== Search file for data section ======
  while (infile)
    {
    infile >> word;
    if (strcmp("NUMBER_OF_DATAPOINTS:",word))           // no pattern match.
      {
      infile.getline(dummyline,4*ONE27,'\n');           // goto next line.
      }
    else                                                // in data section
      {
      foundsection=true;
      infile >> numberofpoints;
      klo=0;                            // initial guess
      khi=1;                            // initial guess, place in constructor
      time.resize(numberofpoints,1);                                    //
      data_x.resize(numberofpoints,1);                          //
      data_y.resize(numberofpoints,1);                          //
      data_z.resize(numberofpoints,1);                          // 
      data_xv.resize(numberofpoints,1);                          //
      data_yv.resize(numberofpoints,1);                          //
      data_zv.resize(numberofpoints,1);                          // 

cerr << "[orbit.cc] " << "orbvector_type: " << orbvector_type << " vs velo prm fixed: " << ORB_PRM_VEL << endl;

      // ______ Actually read data ______
      for (i=0;i<numberofpoints;i++)
        {
        infile.getline(dummyline,ONE27,'\n');           // goto next record (start data)
          infile >> time(i,0)
                 >> data_x(i,0)
                 >> data_y(i,0)
                 >> data_z(i,0);
        }
      } // else not foundsection
    } // file

  // ______ Tidy up ______
  infile.close();
  if (!foundsection)
    {
    WARNING << "string: \"NUMBER_OF_DATAPOINTS:\" not found (i.e., orbit) in file: "
         << file << ".";
    WARNING.print();
    }
  INFO << numberofpoints << " datapoints (t,x,y,z) read from: \""
       << file << "\"";
  INFO.print();
 

  // ______ Compute interpolation coefficients ______
  DEBUG << "value of interp_method: " << interp_method;
  DEBUG.print();
  if (foundsection)
    {
    // ______ Check if time sorted data ______
    for (i=0; i<numberofpoints-1; ++i)
      {
      const real8 h = time(i+1,0) - time(i,0);  // delta_t, might not be const.
      if (h < EPS)                              // ~= 0.)
        {
        PRINT_ERROR("orbit time axis: require distinct, time sorted data")
        throw(input_error);// could simply warn as alternative
        }
      }
    // ___ set default degree if not spline or set explicitly ___
    if (interp_method==ORB_DEFAULT)
      {
      INFO.print("Setting default orbit interpolation method.");
      // use interpolation, since datapoints are result of orbit 
      // propogator, but don't be stupid.
      // use some redundancy to be able to detect outlier!
      interp_method = (numberofpoints>6) ? 5 : numberofpoints-2;// 1 redundant for testing!
      // for radarsat however, do use interpolation, seems to work
      // even when there are warnings?
      //if (numberofpoints==15 && (time(1,0)-time(0,0))>479.9)// RSAT?
      if ((time(1,0)-time(0,0))>479.9 && (time(1,0)-time(0,0))<481.1)// RSAT?
        {
        WARNING.print("Assuming RADARSAT: use highest polyfit recommended.");
        interp_method = numberofpoints-1;// even polyfit(14) seems to work..
        }
      }
    // ___ call function that computes coefficients ___
    computecoefficients();              // cubic spline or polynomial
    PROGRESS.print("Orbit: interpolation coefficients computed.");
    }
  #ifdef __DEBUG
  showdata();
  const real8 t_tmp = time(0,0);
  DEBUG.width(22);
  DEBUG.precision(20);
  DEBUG.rewind();
  DEBUG << "check: computing orbit for t between t0 and t1 = " << t_tmp;
  DEBUG.print();
  cn pos_tmp = getxyz(t_tmp);
  cn vel_tmp = getxyzdot(t_tmp);
  cn acc_tmp = getxyzddot(t_tmp);
  DEBUG << "interp. pos: " << pos_tmp.x << " " << pos_tmp.y << " " << pos_tmp.z;
  DEBUG.print();
  DEBUG << "interp. vel: " << vel_tmp.x << " " << vel_tmp.y << " " << vel_tmp.z;
  DEBUG.print();
  DEBUG << "interp. acc: " << acc_tmp.x << " " << acc_tmp.y << " " << acc_tmp.z;
  DEBUG.print();
  #endif
  } // END initialize



/****************************************************************
 *    Compute coefficients for piecewise polynomial.            *
 * time axis and x,y,z should be filled already.                *
 * call splineinterpol routine.                                 *
 #%// BK 18-Jul-2000                                            *
 ****************************************************************/
void orbit::computecoefficients()
  {
  TRACE_FUNCTION("orbit::computecoefficients (BK 18-Jul-2000)")
  DEBUG << "value of interp_method: " << interp_method;
  DEBUG.print();
  switch (interp_method)// either spline method or degree of polynomial
    {
    case ORB_SPLINE:
      INFO.print("Computing coefficients for orbit spline interpolation");
      coef_x = splineinterpol(time,data_x);
      coef_y = splineinterpol(time,data_y);
      coef_z = splineinterpol(time,data_z);
      break;
    default:
      INFO << "Computing coefficients for orbit polyfit degree: "
           << interp_method;
      INFO.print();
      if (interp_method < 0)
        {
        PRINT_ERROR("degree of polynomial < 0")
        throw(input_error);
        }
      coef_x = polyfit(time,data_x,interp_method);// method==degree
      coef_y = polyfit(time,data_y,interp_method);// method==degree
      coef_z = polyfit(time,data_z,interp_method);// method==degree
    }
  } // END computecoefficients



/****************************************************************
 *    private getklokhi                                         *
 *                                                              *
 * Returns klo/khi for spline interpolation routines.           *
 * klo is index in time axis before point t, khi is index after * 
 * old values are checked first because lot of successive       *
 * interpolation between same t0 t1 is expected.                *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 #%// BK 18-Jul-2000: orbit class                               *
 ****************************************************************/
void orbit::getklokhi(real8 t)
  {
  #ifdef __DEBUG
  // ___ Called a lot, within debug ___
  TRACE_FUNCTION("orbit::getinterval (BK 18-Jul-2000)")
  #endif

  // ______ Check if last interval still applies, init. to [0;1] ______
  if (time(klo,0) <= t && time(khi,0) >= t) 
    return;

  // _____ Else compute correct interval ______
  int32 k;
  klo = 0;
  khi = numberofpoints - 1;
  while(khi-klo>1)
    {
    k = (khi+klo) >> 1;                         // bisection, divide by two
    if (time(k,0) > t) 
      khi = k;
    else
      klo = k;
    }
  } // END getinterval



/****************************************************************
 * public getxyz                                                *
 *                                                              *
 * Returns (natural) cubic spline interpolant.                  *
 * based on numerical recipes routine: splint pp.113            *
 * input:                                                       *
 *  - matrix time and position info                             *
 *  - matrix from splineinterpol                                *
 *  - t to be evaluated                                         *
 * output:                                                      *
 *  - interpolated value at t                                   *
 *                                                              *
 * use routine getsatpos, later we will make an orbit class...  *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 *    Bert Kampes, 02-Jun-1999 streamlined, const variables etc.*
 ****************************************************************/
cn orbit::getxyz(
        real8 t) // const // not, klo,khi
  {
#ifdef __DEBUG
  TRACE_FUNCTION("orbit::getxyz (BK 02-Jun-1999)")
  if (t < time(0,0) || t > time(numberofpoints-1,0))
    {
    WARNING << "interpolation at: " << t << " is outside interval time axis: ("
         << time(0,0) << ", " << time(numberofpoints-1,0) << ").";
    WARNING.print();
    }
#endif

  cn satpos;
  if (interp_method==ORB_SPLINE)
    {
    // ______ Compute correct interval ______
    getklokhi(t);                               // return correct interval
    const real8 h = time(khi,0) - time(klo,0);  // delta_t, might not be const.
    const real8 a = (time(khi,0) - t) / h;
    const real8 b = 1. - a;
    // ______ Evaluate function in x,y,z ______
    satpos.x = a*data_x(klo,0)+b*data_x(khi,0) +
      (((sqr(a)*a-a)*coef_x(klo,0)+(sqr(b)*b-b)*coef_x(khi,0))*sqr(h))/6.; 
    satpos.y = a*data_y(klo,0)+b*data_y(khi,0) +
      (((sqr(a)*a-a)*coef_y(klo,0)+(sqr(b)*b-b)*coef_y(khi,0))*sqr(h))/6.; 
    satpos.z = a*data_z(klo,0)+b*data_z(khi,0) +
      (((sqr(a)*a-a)*coef_z(klo,0)+(sqr(b)*b-b)*coef_z(khi,0))*sqr(h))/6.; 
    }
  // ______ Orbit interpolator is simple polynomial ______
  // ______ unfortunately, i normalize data here each time.
  // ______ we could normalize input t, but then also t_azi1
  // ______ and maybe more consequences for code.  see howlong this takes
  // ______ we normalize here only by a shift to avoid rescaling xdot and xddot
  // ______ (it seems else vdot=vdot*.5(min+max) for some reason)
  else
    {
    const real8 t_tmp = (t-time(numberofpoints/2,0))/real8(10.0);
    satpos.x = polyval1d(t_tmp, coef_x);
    satpos.y = polyval1d(t_tmp, coef_y);
    satpos.z = polyval1d(t_tmp, coef_z);
    }
  return satpos;
  } // END getxyz



/****************************************************************
 *    getxyzdot                                                 *
 *                                                              *
 * Returns 1st derivative natural cubic spline interpolant.     *
 * Numerical recipes splint pp113                               *
 * input:                                                       *
 *  - t to be evaluated                                         *
 * output:                                                      *
 *  - interpolated value at t for x,y,z                         *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 *    Bert Kampes, 02-Jun-1999 streamlined, const variables etc.*
 *    Bert Kampes, 06-Jul-2000 bugfix dt!=1                     *
 ****************************************************************/
cn orbit::getxyzdot(
        real8 t)
  {
#ifdef __DEBUG
  TRACE_FUNCTION("orbit::getxyzddot (BK 06-Jul-2000)")
  if (t < time(0,0) || t > time(numberofpoints-1,0))
    {
    WARNING << "interpolation at: " << t << " is outside interval time axis: ("
         << time(0,0) << ", " << time(numberofpoints-1,0) << ").";
    WARNING.print();
    }
#endif

  cn satvel;
  // ______ Orbit interpolator is splines piecewise polynomial ______
  if (interp_method==ORB_SPLINE)
    {
    // ______ Compute correct interval ______
    getklokhi(t);                                       // return correct interval
    const real8 h = time(khi,0) - time(klo,0);  // delta_t, might not be const.
    const real8 a = (time(khi,0) - t) / h;
    const real8 b = 1. - a;
    // ______ Evaluate 1st derivative of function in x,y,z ______
    satvel.x = ((data_x(khi,0)-data_x(klo,0)) / h) +
      (h * (((1-(3*sqr(a)))*coef_x(klo,0)) + (((3*sqr(b))-1)*coef_x(khi,0)))) / 6.;
    satvel.y = ((data_y(khi,0)-data_y(klo,0)) / h) +
      (h * (((1-(3*sqr(a)))*coef_y(klo,0)) + (((3*sqr(b))-1)*coef_y(khi,0)))) / 6.;
    satvel.z = ((data_z(khi,0)-data_z(klo,0)) / h) +
      (h * (((1-(3*sqr(a)))*coef_z(klo,0)) + (((3*sqr(b))-1)*coef_z(khi,0)))) / 6.;
    }
  // ______ Orbit interpolator is simple polynomial ______
  else
    {
    const int32 DEGREE = coef_x.lines()-1;
    satvel.x = coef_x(1,0);// a_1*t^0 + 2a_2*t^1 + 3a_3*t^2 + ...
    satvel.y = coef_y(1,0);
    satvel.z = coef_z(1,0);
    const real8 t_tmp = (t-time(numberofpoints/2,0))/real8(10.0);
    for (int32 i=2; i<=DEGREE; ++i)
      {
      real8 powt = real8(i)*pow(t_tmp,real8(i-1));
      satvel.x += coef_x(i,0)*powt;
      satvel.y += coef_y(i,0)*powt;
      satvel.z += coef_z(i,0)*powt;
      }
    satvel.x /= real8(10.0);// normalization
    satvel.y /= real8(10.0);// normalization
    satvel.z /= real8(10.0);// normalization
    }
  return satvel;
  } // END getxyzdot



/****************************************************************
 *    getxyzddot                                                *
 *                                                              *
 * Returns 2nd derivative natural cubic spline interpolant.     *
 * Numerical recipes splint pp113                               *
 * input:                                                       *
 *  - t to be evaluated                                         *
 * output:                                                      *
 *  - interpolated value at t for x,y,z                         *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 *    Bert Kampes, 02-Jun-1999 streamlined, const variables etc.*
 *    Bert Kampes, 06-Jul-2000 bugfix dt!=1                     *
 ****************************************************************/
cn orbit::getxyzddot(
        real8 t)
  {
#ifdef __DEBUG
  TRACE_FUNCTION("orbit::getxyzddot (BK 06-Jul-2000)")
  if (t < time(0,0) || t > time(numberofpoints-1,0))
    {
    WARNING << "interpolation at: " << t << " is outside interval time axis: ("
         << time(0,0) << ", " << time(numberofpoints-1,0) << ").";
    WARNING.print();
    }
#endif

  cn satacc;
  // ______ Orbit interpolator is splines piecewise polynomial ______
  if (interp_method==ORB_SPLINE)
    {
    // ______ Compute correct interval ______
    getklokhi(t);                                       // return correct interval
    const real8 h = time(khi,0) - time(klo,0);  // delta_t, might not be const.
    const real8 a = (time(khi,0) - t) / h;
    const real8 b = 1. - a;

    // ______ Evaluate 2nd derivative of function in x,y,z ______
    satacc.x = a*coef_x(klo,0) + b*coef_x(khi,0);
    satacc.y = a*coef_y(klo,0) + b*coef_y(khi,0);
    satacc.z = a*coef_z(klo,0) + b*coef_z(khi,0);
    }
  // ______ Orbit interpolator is simple polynomial ______
  else
    {
    // 2a_2 + 2*3a_3*t^1 + 3*4a_4*t^2...
    const int32 DEGREE = coef_x.lines()-1;
    satacc.x = 0.0;
    satacc.y = 0.0;
    satacc.z = 0.0;
    const real8 t_tmp = (t-time(numberofpoints/2,0))/real8(10.0);
    for (int32 i=2; i<=DEGREE; ++i)
      {
      real8 powt = real8((i-1)*i)*pow(t_tmp,real8(i-2));
      satacc.x += coef_x(i,0)*powt;
      satacc.y += coef_y(i,0)*powt;
      satacc.z += coef_z(i,0)*powt;
      }
    satacc.x /= real8(100.0);//normalization
    satacc.y /= real8(100.0);//normalization
    satacc.z /= real8(100.0);//normalization
    }
  return satacc;
  } // END getxyzddot


// ====== Friends ======
/****************************************************************
 *    lp2xyz                                                    *
 *                                                              *
 * converts line,pixels to x,y,z on (which) ellipsoid?          *
 *                                                              *
 * input:                                                       *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default ~ 10-6 m ?)                              *
 * output:                                                      *
 *  - cn (x,y,z) (returnpos is filled)                          *
 *  - number of iterations is returned                          *
 *                                                              *
 *    Bert Kampes, 04-Jan-1999                                  *
 ****************************************************************/
int32 lp2xyz(
        real8            line,
        real8            pixel,
        const input_ell &ell,
        const slcimage  &image,
        orbit           &orb,
        cn              &returnpos,
        int32            MAXITER,
        real8            CRITERPOS)
  {
  TRACE_FUNCTION("lp2xyz (BK 04-Jan-1999)")

  // ______ Convert lp to time ______
  const real8 aztime = image.line2ta(line);
  const real8 ratime = image.pix2tr(pixel);

  const cn possat = orb.getxyz(aztime);
  const cn velsat = orb.getxyzdot(aztime);

  // ______ Set up system of equations and iterate ______
  returnpos.x  = image.approxcentreoriginal.x;  // iteration 0
  returnpos.y  = image.approxcentreoriginal.y;  // iteration 0
  returnpos.z  = image.approxcentreoriginal.z;  // iteration 0


// ______ Save some alloc, init and dealloc time by declaring static (15%?) ______
  static matrix<real8> solxyz(3,1);             // solution of system
  static matrix<real8> equationset(3,1);        // observations
  static matrix<real8> partialsxyz(3,3);        // partials to position due to nonlinear

  register int32 iter;
  for (iter=0; iter<=MAXITER; ++iter)
    {

    // ______ Update equations + solve system ______
    const cn dsat_P  =  returnpos.min(possat);  // vector satellite, P on ellipsoid
    equationset(0,0) = -eq1_doppler(velsat, dsat_P);
    equationset(1,0) = -eq2_range(dsat_P, ratime);
    equationset(2,0) = -eq3_ellipsoid(returnpos,ell.a,ell.b);
    partialsxyz(0,0) =  velsat.x;
    partialsxyz(0,1) =  velsat.y;
    partialsxyz(0,2) =  velsat.z;
    partialsxyz(1,0) =  2*dsat_P.x;
    partialsxyz(1,1) =  2*dsat_P.y;
    partialsxyz(1,2) =  2*dsat_P.z;
    partialsxyz(2,0) = (2*returnpos.x)/sqr(ell.a);
    partialsxyz(2,1) = (2*returnpos.y)/sqr(ell.a);
    partialsxyz(2,2) = (2*returnpos.z)/sqr(ell.b);

    // ______ Solve system ______
    solve33(solxyz, equationset,partialsxyz);

    // ______Update solution______
    returnpos.x += solxyz(0,0);                         // update approx. value
    returnpos.y += solxyz(1,0);                         // update approx. value
    returnpos.z += solxyz(2,0);                         // update approx. value

    // ______Check convergence______
    if (abs(solxyz(0,0)) < CRITERPOS &&                         // dx
        abs(solxyz(1,0)) < CRITERPOS &&                         // dy
        abs(solxyz(2,0)) < CRITERPOS   )                        // dz
      break; // converged
    }

  // ______ Check number of iterations ______
  if (iter >= MAXITER)
    {
    WARNING << "line, pix -> x,y,z: maximum iterations (" << MAXITER << ") reached. "
         << "Criterium (m): "<< CRITERPOS
         << "dx,dy,dz=" << solxyz(0,0) << ", " << solxyz(1,0) << ", " << solxyz(2,0);
    WARNING.print();
    }

  // ______ (Converged) Result is in returnpos ______
  return iter;
  } // END lp2xyz


/****************************************************************
 *    xyz2orb                                                   *
 *                                                              *
 * converts xyz (which) ellipsoid? to orbital coordinates       *
 * ellips, xyz is in system of orbit ephemerides                * 
 *                                                              *
 * input:                                                       *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default = 10-15s)                                *
 * output:                                                      *
 *  - x,y,z                                                     *
 *  - number of iterations                                      *
 *                                                              *
 #%// BK 22-Sep-2000                                            *
 ****************************************************************/
int32 xyz2orb(
        cn              &possat,
        const slcimage  &image,
        orbit           &orb,       // non const, khi/klo
        const cn        &pointonellips,
        int32            MAXITER,   // defaults
        real8            CRITERTIM) // seconds
  {
  TRACE_FUNCTION("xyz2orb (BK 22-Sep-2000)");
  // ______ Initial value azimuth time ______
  real8 sol = 0.0;
  register int32 iter;
  real8 tazi = image.line2ta(.5*(image.currentwindow.linehi-image.currentwindow.linelo));
  for(iter=0 ;iter<=MAXITER; ++iter)           // break at convergence
    {
    // ______ Update equations ______
    possat = orb.getxyz(tazi);
    const cn velsat = orb.getxyzdot(tazi);
    const cn accsat = orb.getxyzddot(tazi);
    const cn delta  = pointonellips.min(possat);

    // ______ Update solution ______
    sol   = -eq1_doppler(velsat, delta) /
             eq1_doppler_dt(delta, velsat, accsat);
    tazi +=  sol;

    // ______ Check convergence ______
    if (abs(sol) < CRITERTIM)                   // dta
      break;
    }

  // ______ Check number of iterations _____
  if (iter >= MAXITER)
    {
    WARNING << "x,y,z -> line, pix: maximum iterations (" << MAXITER << ") reached. "
         << "Criterium (s):" << CRITERTIM
         << "dta (s)=" << sol;
    WARNING.print();
    }

  // ====== Compute range time ======
  // ______ Update equations ______
  possat = orb.getxyz(tazi);
  return iter;
  } // END xyz2orb




/****************************************************************
 *    xyz2t                                                     *
 *                                                              *
 * converts xyz (which) ellipsoid? to azimuth time with         *
 *  zero doppler equation, use slant range eq. for range time   *
 #%// BK 18-Jul-2000                                            *
 * ellips, xyz is in system of orbit ephemerides                * 
 *                                                              *
 * input:                                                       *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default = 10-15s)                                *
 * output:                                                      *
 *  - (updated tazi,tran)                                       *
 *  - number of iterations                                      *
 *                                                              *
 *    Bert Kampes, 04-Jan-1999                                  *
 ****************************************************************/
int32 xyz2t(
        real8           &tazi,      // azimuth
        real8           &tran,      // and range time
        const slcimage  &image,
        orbit           &orb,       // non const, khi/klo
        const cn        &pos,
        int32            MAXITER,   // defaults
        real8            CRITERTIM) // seconds
  {
  TRACE_FUNCTION("xyz2t (BK 04-Jan-1999)")

  // ______ Compute initial value ______
  tazi = image.line2ta(0.5*(image.currentwindow.linehi-image.currentwindow.linelo));
  // _____ Start _______
  real8 sol = 0.0;
  register int32 iter;
  for(iter=0 ;iter<=MAXITER; ++iter)           // break at convergence
    {
    // ______ Update equations ______
    const cn possat = orb.getxyz(tazi);
    const cn velsat = orb.getxyzdot(tazi);
    const cn accsat = orb.getxyzddot(tazi);
    const cn delta  = pos.min(possat);

    // ______ Update solution ______
    sol   = -eq1_doppler(velsat, delta) /
             eq1_doppler_dt(delta, velsat, accsat);
    tazi +=  sol;

    // ______ Check convergence ______
    if (abs(sol) < CRITERTIM)                   // dta
      break;
    }

  // ______ Check number of iterations _____
  if (iter >= MAXITER)
    {
    WARNING << "x,y,z -> line, pix: maximum iterations (" << MAXITER << ") reached. "
         << "Criterium (s):" << CRITERTIM
         << "dta (s)=" << sol;
    WARNING.print();
    }

  // ====== Compute range time ======
  // ______ Update equations ______
  const cn possat = orb.getxyz(tazi);
  const cn delta  = pos.min(possat);
  tran      = delta.norm() / SOL;

  return iter;
  } // END xyz2t



/****************************************************************
 *    xyz2lp                                                    *
 *                                                              *
 * converts xyz (which) ellipsoid? to line,pixels               *
 *                                                              *
 * input:                                                       *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default ~ 10-15s)                                *
 * output:                                                      *
 *  - cn (x,y,z) (returnpos is filled)                          *
 *  - number of iterations                                      *
 *                                                              *
 *    Bert Kampes, 04-Jan-1999                                  *
 ****************************************************************/
int32 xyz2lp(
        real8           &returnline,
        real8           &returnpixel,
        const slcimage &image,
        orbit           &orb,
        const cn        &pos,               // point at ground
        int32            MAXITER,
        real8            CRITERTIM)
  {
  TRACE_FUNCTION("xyz2lp (BK 04-Jan-1999)");
  real8 tazi;
  real8 tran;

  // ______ Compute tazi, tran ______
  int32 iter = xyz2t(tazi,tran,image,orb,pos,MAXITER,CRITERTIM);

  // ______ Convert time to pixel ______
  // ______ (Converged) Result is in returnline/pixel ______
  returnline  = image.ta2line(tazi);
  returnpixel = image.tr2pix (tran);
  return iter;
  } // END xyz2lp



/****************************************************************
 *    ell2lp                                                    *
 *                                                              *
 * converts ell to x,y,z, to l,p                                *
 *                                                              *
 * input:                                                       *
 *  - ??                                                        *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default ~ 10-6 m ?)                              *
 * output:                                                      *
 *  - phi,lambda,height (updated)                               *
 *  - number of iterations                                      *
 *                                                              *
 *    Bert Kampes, 27-Jan-1999                                  *
 ****************************************************************/
int32 ell2lp(
        real8           &returnline,
        real8           &returnpixel,
        const input_ell &ellips,
        const slcimage &image,
        orbit           &orb,
        real8            phi,
        real8            lambda,
        real8            height,
        int32            MAXITER,
        real8            CRITERTIM)
  {
  TRACE_FUNCTION("ell2lp (BK 27-Jan-1999)")
  // ______ Transform ell2xyz ______
  cn returnpos = ellips.ell2xyz(phi, lambda, height);
  DEBUG << "tmp result phi,lambda,h --> x,y,z: " 
        << phi << ", " << lambda << ", " << height << " --> "
        << returnpos.x << " " << returnpos.y << " " << returnpos.z;
  DEBUG.print();
  // ______ Transform xyz2lp ______
  int32 n_iter = xyz2lp(returnline, returnpixel,
            image, orb, returnpos,
            MAXITER, CRITERTIM);
  return n_iter;                            // number of iterations
  } // END ell2lp



/****************************************************************
 *    lp2ell                                                    *
 *                                                              *
 * converts line,pixels to x,y,z to geodetic coordinates ell    *
 *                                                              *
 * input:                                                       *
 *  - line,pixel, etc.                                          *
 *  - MAXITER (default = 10)                                    *
 *  - CRITER  (default ~ 10-6 m ?)                              *
 * output:                                                      *
 *  - phi,lambda,height (updated)                               *
 *  - number of iterations                                      *
 *                                                              *
 *    Bert Kampes, 27-Jan-1999                                  *
 ****************************************************************/
int32 lp2ell(
        real8            line,
        real8            pixel,
        const input_ell &ellips,
        const slcimage  &image,
        orbit           &orb,
        real8           &returnphi,
        real8           &returnlambda,
        real8           &returnheight,
        int32            MAXITER,
        real8            CRITERPOS)
  {
  TRACE_FUNCTION("lp2ell (BK: 27-Jan-1999)");
  // ______ Transform lp2xyz ______
  cn returnpos;
  int32 iter = lp2xyz(line, pixel, ellips, image, orb,
                       returnpos, MAXITER, CRITERPOS);
  // ______ Transform xyz2ell, compute phi,lambda,h ______
  ellips.xyz2ell(returnpos, returnphi, returnlambda, returnheight);
  return iter;
  } // END lp2ell



/****************************************************************
 *    dumporbit                                                 *
 *                                                              *
 * rescale data to interval [-1,1]                              *
 * data -> X(data-min)/(max-min) ; data+SS                      *
 * X=1--1==EE-SS; SS=-1, EE=1                                   *
 *                                                              *
 * input:                                                       *
 * output:                                                      *
 *    Bert Kampes, 03-Jul-2000                                  *
 ****************************************************************/
void orbit::dumporbit(
        const input_pr_orbits &inputorb,
        const int16 ID)
  {
  TRACE_FUNCTION("orbit::dumporbit (BK 03-Jul-2000)")
  // ___ prevent error if called after readfiles, but no orbit section ___
  if (numberofpoints==0)
    {
    INFO.print("Exiting dumporbit, no orbit data available.");
    return;
    }
  // ====== Compute positions for center pixel ======
  // ______ *dumporbit initialized to -1.0 in readinput.c ______
  real8 dt=0.;
  char ofile[EIGHTY];
  if (ID==MASTERID)
    {
    if (inputorb.dumpmasterorbit<0) return;
    dt = inputorb.dumpmasterorbit;
    strcpy(ofile,"masterorbit.dat");
    PROGRESS.print("Dumping master orbit.");
    }
  else if (ID==SLAVEID)
    {
    dt = inputorb.dumpslaveorbit;
    if (inputorb.dumpslaveorbit<0) return;
    strcpy(ofile,"slaveorbit.dat");
    PROGRESS.print("Dumping slave orbit.");
    }
  else
    {
    PRINT_ERROR("Panic: not possible.")
    throw(unhandled_case_error);
    }
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;
  INFO << "dumporbits: MAXITER: "   << MAXITER   << "; "
                   << "CRITERPOS: " << CRITERPOS << " m; "
                   << "CRITERTIM: " << CRITERTIM << " s";
  INFO.print();

  //  ______ Evaluate polynomial orbit for t1:dt:tN ______
  int32 outputlines = 1 + int32((time(numberofpoints-1,0)-time(0,0)) / dt);
  real8 tazi = time(0,0);
  ofstream fo(ofile,ios::out | ios::trunc);
  matassert(fo,ofile,__FILE__,__LINE__);
  fo.precision(3);
  fo.width(11);
  // this is ok, but test if not: fo.setf(ios::left);
  fo.setf(ios::fixed);
  fo.setf(ios::showpoint);
  for (register int32 i=0; i<outputlines; ++i)
    {
    const cn position = getxyz(tazi);
    const cn velocity = getxyzdot(tazi);
    const cn accerela = getxyzddot(tazi);
    fo << tazi << " "
       << position.x << " " << position.y << " " << position.z << " "
       << velocity.x << " " << velocity.y << " " << velocity.z << " "
       << accerela.x << " " << accerela.y << " " << accerela.z << "\n";
    tazi += dt;
    }
  fo.close();


  // ______ dump coeff. as well for testing ... ______
  #ifdef __DEBUG
  if (ID==MASTERID)
    {
    DEBUG.print("dumping files m_t, m_x, m_y, m_z, m_cx, m_cy, m_cz for spline interpolation.");
    dumpasc("m_t",time);   
    dumpasc("m_x",data_x);  dumpasc("m_y",data_y);  dumpasc("m_z",data_z);
    dumpasc("m_cx",coef_x); dumpasc("m_cy",coef_y); dumpasc("m_cz",coef_z);
    }
  else
    {
    DEBUG.print("dumping files s_t, s_x, s_y, s_z, s_cx, s_cy, s_cz for spline interpolation.");
    dumpasc("s_t",time);
    dumpasc("s_x",data_x);  dumpasc("s_y",data_y);  dumpasc("s_z",data_z);
    dumpasc("s_cx",coef_x); dumpasc("s_cy",coef_y); dumpasc("s_cz",coef_z);
    }
  #endif
  } // END dumporbit




// ====== helper functions ======
/****************************************************************
 *    splineinterpol                                            *
 *                                                              *
 * Based on Numerical recipes spline pp.113                     *
 * time <-> x                                                   *
 * data <-> y                                                   *
 * pp   <-> y2                                                  *
 * rhs  <-> u                                                   *
 *                                                              *
 * Compute coefficients of piecewize polynomials for            *
 *  cubic splines. Use natural ones if __NATURALSPLINE__ is     *
 *  defined (by default internal in this routine.)              *
 * Else use boundary condition from nrc pp113.                  *
 * (fix first derivative.)                                      *
 * There does not seem to be a big difference between both      *
 * boundary conditions, particularly for middle segment.        *
 *  A pp = rhs; L Lt pp = rhs; (?)                              *
 *  pptmp=Lt pp); L pptmp = rhs; Lt pp = pptmp; (?)             *
 * input:                                                       *
 *  - matrix by getdata with time and position info             *
 * output:                                                      *
 *  - matrix with 2nd order derivatives                         *
 *    (input for interp. routines)                              *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 #%// BK 17-Jul-2000: added nonnatural                          *
 ****************************************************************/
matrix<real8> splineinterpol(
        const matrix<real8> &time,
        const matrix<real8> &data)
  {
  TRACE_FUNCTION("splineinterpol (BK 17-Jul-2000)")
  if (time.pixels() != 1 || data.pixels() != 1)
    {
    PRINT_ERROR("code 901: splineinterpol: wrong input.");
    throw(input_error);
    }
  if (time.lines() != data.lines())                             // number of points
    {
    PRINT_ERROR("code 901: splineinterpol: require same size vectors.");
    throw(input_error);
    }

  const int32   N = time.lines();                       // number of points
  matrix<real8> rhs(N-1,1);
  matrix<real8> pp(N,1);                        // init to 0.

  // #define __NATURALSPLINE__ // in Makefile or here if so desired.
#define __NATURALSPLINE__
  // ====== Set boundary condition first datapoint ======
  #ifdef __NATURALSPLINE__
    pp(0,0)  = 0.0;                     // natural spline boundary condition
  #else
    // ______ yp1 is first der. at point 0 ______
    // estimate it by dy/dx: nearly linear for sat. orbit...
    const real8 yp1 = (data(1,0)-data(0,0)) / (time(1,0)-time(0,0));
    pp(0,0)  = -0.5;
    rhs(0,0) = (3.0 / (time(1,0)-time(0,0))) *
               ((data(1,0)-data(0,0)) / (time(1,0)-time(0,0)) - yp1); 
  #endif

  // ====== Decomposition loop ======
  register int32 i;
  for (i=1; i<=N-2; ++i)
    {
    const int32 ip1 = i + 1;
    const int32 im1 = i - 1;
    const real8 sig = (time(i,0)-time(im1,0)) / (time(ip1,0)-time(im1,0));
    const real8 p   = sig*pp(im1,0) + 2.0;
    //
    pp(i,0)   = (sig-1.0) / p;
    rhs(i,0)  =   (data(ip1,0)-data(i,0)) / (time(ip1,0)-time(i,0))
                - (data(i,0)-data(im1,0)) / (time(i,0)-time(im1,0));
    rhs(i,0)  = (6.*rhs(i,0)/(time(ip1,0)-time(im1,0))-sig*rhs(im1,0))/p;
    }

  // ====== Set boundary condition last datapoint ======
  #ifdef __NATURALSPLINE__
    pp(N-1,0) = 0.0;            // natural spline oundary condition
  #else
    // ______ ypN is first derivative at point N-1 ______
    const real8 ypN = (data(N-1,0)-data(N-2,0)) / (time(N-1,0)-time(N-2,0));
    const real8 qn  = 0.5;
    const real8 un  = (3.0 / (time(N-1,0) - time(N-2,0))) *
                      (ypN - (data(N-1,0) - data(N-2,0)) / (time(N-1,0) - time(N-2,0)));
    pp(N-1,0) = (un - qn*rhs(N-2,0)) / (qn*pp(N-2,0) + 1.0);
  #endif

  // ====== Backsub loop ======
  for (i=N-2; i>=0; --i)
    pp(i,0) = pp(i,0)*pp(i+1,0) + rhs(i,0);
    // pp(i,0) *= pp(i+1,0) + rhs(i,0); // ??

  return pp;            // contains second derivatives at datapoints
  } // END splineinterpol



/****************************************************************
 *    polyfit                                                   *
 *                                                              *
 * Compute coefficients of x=a0+a1*t+a2*t^2+a3*t^3 polynomial   *
 * for orbit interpolation.  Do this to facilitate a method     *
 * in case only a few datapoints are given.                     *
 * Data t is normalized approximately [-x,x], then polynomial   *
 * coefficients are computed.  For poly_val this is repeated    *
 * see getxyz, etc.                                             *
 *                                                              *
 * input:                                                       *
 *  - matrix by getdata with time and position info             *
 * output:                                                      *
 *  - matrix with coeff.                                        *
 *    (input for interp. routines)                              *
 *                                                              *
 *    Bert Kampes, 16-Jun-2003                                  *
 ****************************************************************/
matrix<real8> polyfit(
        const matrix<real8> &time,
        const matrix<real8> &y, 
        const int32 DEGREE)
  {
  TRACE_FUNCTION("polyfit (BK 16-Jun-2003)")
  if (time.pixels() != 1 || y.pixels() != 1)
    {
    PRINT_ERROR("code 902: polyfit: wrong input.");
    throw(input_error);
    }
  if (time.lines() != y.lines())                // number of points
    {
    PRINT_ERROR("code 902: polyfit: require same size vectors.");
    throw(input_error);
    }
 
  // ______ Normalize t for numerical reasons ______
  const int32 Npoints = time.lines();   // number of points
  DEBUG.print("Normalizing t axis for least squares fit");
  matrix<real8> t=(time-time(Npoints/2,0))/real8(10.0);

  // ______ Check redundancy ______
  const int32 Nunk = DEGREE+1;
  DEBUG << "Degree of orbit interpolating polynomial: " << DEGREE;
  DEBUG.print();
  DEBUG << "Number of unknowns: " << Nunk;
  DEBUG.print();
  DEBUG << "Number of data points (orbit): " << Npoints;
  DEBUG.print();
  if (Npoints < Nunk)
    {
    PRINT_ERROR("polyfit: Number of points is smaller than parameters solved for.")
    throw(input_error);
    }

  // ______ Set up system of equations to solve coeff. ______
  DEBUG.print("Setting up lin. system of equations");
  matrix<real8> A(Npoints,Nunk);// designmatrix
  for (int32 j=0; j<=DEGREE; j++)
    {
    matrix<real8> t_tmp = t;
    t_tmp.mypow(real8(j));
    A.setcolumn(j,t_tmp);
    }
  DEBUG.print("Solving lin. system of equations with cholesky.");
  matrix<real8> N      = matTxmat(A,A);
  matrix<real8> rhs    = matTxmat(A,y);
  matrix<real8> Qx_hat = N;
  choles(Qx_hat);               // Cholesky factorisation normalmatrix
  solvechol(Qx_hat,rhs);        // Estimate of unknowns in rhs
  invertchol(Qx_hat);           // Covariance matrix
  // ______Test inverse______
  for (uint i=0; i<Qx_hat.lines(); i++)
    for (uint j=0; j<i; j++)
      Qx_hat(j,i) = Qx_hat(i,j);// repair matrix
  const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
  DEBUG << "polyfit orbit: max(abs(N*inv(N)-I)) = " << maxdev;
  DEBUG.print();
  // ___ report max error... (seems sometimes this can be extremely large) ___
  if (maxdev > 1e-6) WARNING.print("polyfit orbit interpolation instable!");
  matrix<real8> y_hat   = A * rhs;
  matrix<real8> e_hat   = y - y_hat;
  if (max(abs(e_hat)) > 0.02)// 0.05 is already 1 wavelength! (?)
    {
    WARNING << "Max. approximation error at datapoints (x,y,or z?): " << max(abs(e_hat)) << "m";
    WARNING.print();
    }
  else
    {
    INFO << "Max. approximation error at datapoints (x,y,or z?): " << max(abs(e_hat)) << "m";
    INFO.print();
    }
  // ___ report computations in detail ___
  DEBUG.width(16);
  DEBUG.precision(15);
  DEBUG.print("REPORTING POLYFIT LEAST SQUARES ERRORS");
  DEBUG.print("      time              y                yhat               ehat");
  for (int32 i=0; i<Npoints; i++)
    {
    DEBUG << setw(16) << setprecision(15) << time(i,0) << "   " 
          << y(i,0) << "   " << y_hat(i,0) << "   " << e_hat(i,0);
    DEBUG.print();
    }
  // ___ seems sometimes dt is wrong?? --->log it. ___
  for (int32 i=0; i<Npoints-1; i++)
    {
    // ___ check if dt is constant, not necessary for me, but may ___
    // ___ signal error in header data of SLC image ___
    real8 dt   = time(i+1,0)-time(i,0);
    DEBUG << "dt between point " << i+1 << " and " << i << "= " << dt;
    DEBUG.print();
    if(abs(dt-(time(1,0)-time(0,0))) > 0.001)// 1ms i will allow...
      WARNING.print("Orbit: data does not have equidistant time interval?");
    }
  DEBUG.reset();

  // ___ return the coefficients wrt. normalized t ___
  return rhs;
  } // END polyfit



/****************************************************************
 *    orbit::showdata                                           *
 #%// BK 17-Jul-2000                                            *
 ****************************************************************/
void orbit::showdata()
  {
  TRACE_FUNCTION("orbit::showdata (BK 17-Jul-2000)")
  DEBUG.print("Time of orbit ephemerides:");
  time.showdata();
  DEBUG.print("Orbit ephemerides x:");
  data_x.showdata();
  DEBUG.print("Orbit ephemerides y:");
  data_y.showdata();
  DEBUG.print("Orbit ephemerides z:");
  data_z.showdata();
  DEBUG.print("Estimated coefficients x(t)");
  coef_x.showdata();
  DEBUG.print("Estimated coefficients y(t)");
  coef_y.showdata();
  DEBUG.print("Estimated coefficients z(t)");
  coef_z.showdata();
  } // END showdata




#ifdef __TESTMAIN__
/****************************************************************
 *    test program for orbit routine                            *
 #%// BK 18-Jul-2000                                            *
 ****************************************************************/
int32 main()
  {
  orbit masterorbit;
  orbit slaveorbit;

  // file required named "master.out" with string
  // "NUMBER_OF_DATAPOINTS:" and ephemerides required.
  masterorbit.initialize("master.out");
  cerr << "\nShowing master orbit: "
       << masterorbit.npoints() << " points.\n";
  masterorbit.showdata();  // see matrixbk::showdata()

  // file required named "slave.out" with string
  // "NUMBER_OF_DATAPOINTS:" and ephemerides required.
  slaveorbit.initialize("slave.out");
  cerr << "\nShowing slave orbit: "
       << slaveorbit.npoints() << " points.\n";
  slaveorbit.showdata();

  // start interpolation tests
  real8 t = 38002.341;
  cn pos  = masterorbit.getxyz(t);
  cn vel  = masterorbit.getxyzdot(t);
  cn acc  = masterorbit.getxyzddot(t);
  cerr << "pos(" << t << "): " << pos.x << ", " << pos.y << ", " << pos.z << endl;
  cerr << "vel(" << t << "): " << vel.x << ", " << vel.y << ", " << vel.z << endl;
  cerr << "acc(" << t << "): " << acc.x << ", " << acc.y << ", " << acc.z << endl;

  input_pr_orbits orbitinput;
  masterorbit.dumporbit(orbitinput,MASTERID);

  slcimage masterinfo;
  slcimage slaveinfo;
  input_gen generalinput;
  generalinput.dumpbaselineL = 6;
  generalinput.dumpbaselineP = 4;

  //cerr << "\ncompbaseline\n";
  //compbaseline(generalinput,masterinfo,slaveinfo,masterorbit,slaveorbit);

  cerr << "\nOK!\n";
  return 0;
  } // END main
#endif // TESTMAIN



/****************************************************************
 *    orbit::modify                                             *
 *                                                              *
 * Modify orbit from baseline error polynomial. Return updated  *
 * datapoints.                                                  *
 *                                                              *
 * coeff:          baseline error polynomial                    *
 * referenceorbit: defines Frenet frame of unit vectors         *
 * tshift:         azimuth time shift w. r. t. reference orbit  *
 * [tmin,tmax]:    normalisation interval for polynomial coeff. *
 *                                                              *
 # Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
matrix<real8> orbit::modify(const matrix<real8> &coeff,
			    orbit               &referenceorbit,
			    const real8         tshift,
			    const real8         tmin,
			    const real8         tmax)
{
  TRACE_FUNCTION("orbit::modify (HB 17-Jul-2000)")

  // ====== update data points ======
  const uint degree = coeff.lines()/2-1;
  for (uint i=0; i<numberofpoints; i++)
    {
      real8 tref;
      cn ea, er, ex, update;

      tref = time(i,0) + tshift;
      getCoordinateFrame(ea,er,ex,tref,referenceorbit);
      tref = normalize(tref,tmin,tmax);

      // ______ update positions ______
      update = ex * polyval1d(tref, coeff.getdata(window(0,degree,0,0)))
	+ er * polyval1d(tref, coeff.getdata(window(degree+1,2*degree+1,0,0)));
      data_x(i,0) += update.x;
      data_y(i,0) += update.y;
      data_z(i,0) += update.z;

      // ______ update velocities ______
      if (orbvector_type==ORB_PRM_VEL)
	{
	  update = ex * polyval1d(tref, coeff.getdata(window(0,degree,0,0)),1)
	    + er * polyval1d(tref, coeff.getdata(window(degree+1,2*degree+1,0,0)),1);
	  cerr << data_xv(i,0) << " " << update.x << " " << data_yv(i,0) << " " << update.y << " " << " " << data_zv << " " << update.z << endl;
	  data_xv(i,0) += update.x;
	  data_yv(i,0) += update.y;
	  data_zv(i,0) += update.z;
	}
    }

  // ====== re-compute polynomial coefficients ======
  computecoefficients();

  // ====== return datapoints for output ======
  matrix<real8> datapoints;
  if (orbvector_type==ORB_PRM_POS)
    datapoints = matrix<real8>(numberofpoints,4);
  else if (orbvector_type==ORB_PRM_VEL)
    {
      datapoints = matrix<real8>(numberofpoints,7);
      datapoints.setdata(0,4,data_xv);
      datapoints.setdata(0,5,data_yv);
      datapoints.setdata(0,6,data_zv);
    }
  datapoints.setdata(0,0,time);
  datapoints.setdata(0,1,data_x);
  datapoints.setdata(0,2,data_y);
  datapoints.setdata(0,3,data_z);
  return datapoints;
} // END orbit::modify  
