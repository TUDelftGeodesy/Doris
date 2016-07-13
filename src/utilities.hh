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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/utilities.hh,v $
 * $Revision: 3.18 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * - Mathematical utilities.                                    *
 * - fast trigonometric functions sin/cos/atan2                 *
 ****************************************************************/


#ifndef UTILITIES_H                     // guard
#define UTILITIES_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // typedefs
#include "slcimage.hh"                  // my slc image class
#include "orbitbk.hh"                   // my orbit class 
#include "ioroutines.hh"                // ERROR
#include "exceptions.hh"                // exceptions class
#include <cmath>                        // floor, sqr


// === prototype requires below, so put first/.
// ______ Solve system of 3 equations ______
void solve33(
        matrix<real8> &Result,
        const matrix<real8> &rhs,
        const matrix<real8> &partials);


// ______ Solve set of 2 equations ______
matrix<real8> solve22(
        const matrix<real8> &rhs,
        const matrix<real8> &partials);



// ====== Inline ======
inline real8 deg2rad(const real8 &x){return x * PI / 180.0;}
inline real4 deg2rad(const real4 &x){return x * PI / 180.0;}
inline real8 rad2deg(const real8 &x){return x * 180.0 / PI;}
inline real4 rad2deg(const real4 &x){return x * 180.0 / PI;}

// ______ For resampling output ci2 casting ______
// ______ No clipping check? max is 2^15=-32768:32767
inline compli16 cr4toci2(const complr4 &x)
  {return compli16((x.real() > 0.0) ? short(x.real()+0.5) : short(x.real()-0.5),
                   (x.imag() > 0.0) ? short(x.imag()+0.5) : short(x.imag()-0.5));}


// ====== Inline ======
inline bool iseven(const int16 &w)              {return (w+1)%2;}
inline bool iseven(const int32 &w)              {return (w+1)%2;}
inline bool iseven(const int64 &w)              {return (w+1)%2;}
inline bool iseven(const uint  &w)              {return (w+1)%2;}
inline bool iseven(const uint64  &w)            {return (w+1)%2;}
inline bool isodd (const int16 &w)              {return w%2;}
inline bool isodd (const int32 &w)              {return w%2;}
inline bool isodd (const int64 &w)              {return w%2;}
inline bool isodd (const uint  &w)              {return w%2;}
inline bool isodd (const uint64 &w)             {return w%2;}

// --- added 1 as power of 2 bk 8/05
inline bool ispower2 (const uint &w)
  {if (w==1   || w==2    || w==4    || w==8    || w==16   ||
       w==32  || w==64   || w==128  || w==256  ||
       w==512 || w==1024 || w==2048 || w==4096   )
     return true;
   return false;}

inline int32 Ncoeffs(const int32 &degree)
  {return int32(0.5*(sqr(degree+1)+degree+1));}

inline int32 degree(const int32 &Ncoeffs)
  { return int32(0.5*(-1 + int32(sqrt(real4(1+8*Ncoeffs)))))-1;}

inline real4 remainder(const real4 &number, const real4 &divisor)       
  {return number-floor(number/divisor)*divisor;}

inline real8 remainder(const real8 &number, const real8 &divisor)       
  {return number-floor(number/divisor)*divisor;}

inline real4 sinc(const real4 &x)       
  {return ((x==0) ? 1 : sin(PI*x)/(PI*x));}

//____ Raffaele Nutricato 
inline real8 sinc(const real8 &x)       
  {return ((x==0) ? 1 : sin(PI*x)/(PI*x));}

inline real4 rect(const real4 &x)       
  {real4 ans = 0.0;
   if (x<0.5 && x>-0.5) ans = 1;
   else if (x==0.5 || x==-0.5) ans = 0.5;
   return ans;}

inline real4 tri(const real4 &x)        
  {real4 ans = 0.0;
   if (x<1.0 && x>-1.0)
     {(x<0) ? ans=1+x : ans=1-x;}
   return ans;}

// --- positve numbers ? --- // Bert Kampes, 24-Aug-2005
inline real8 onedecimal(const real8 &x)
  {return real8(int32(x*10+0.5))/10.0;}
inline real4 onedecimal(const real4 &x)
  {return real4(int32(x*10+0.5))/10.0;}



/****************************************************************
 * myrect                                                       *
 * rect window, lying vector                                    *
 * r(i) = 1 if abs(i)<=.5                                       *
 * used in filtering, see curlander, swabisch, geudtner         *
 *    Bert Kampes, 31-Mar-2000                                  *
 ****************************************************************/
inline matrix<real4> myrect(const matrix<real4> &X)
  {
  TRACE_FUNCTION("myrect (BK 31-Mar-2000)")
  #ifdef __DEBUG
    if (X.lines()!=1) 
      {
      PRINT_ERROR("myrect: only lying vectors.")
      throw(argument_error);
      }
  #endif
  matrix<real4> Res(1,X.pixels());      // init2zero
  real4 *pntX = X[0];
  real4 *pntR = Res[0];
  for (register int32 i=0; i<int32(X.pixels()); ++i)
    {
    //if (abs(*pntX++)<=0.5)
    //  *pntR=1.;
    //pntR++;
    //changed by FvL (for g++/gcc > 4.0):
    if (abs(*pntX)<=0.5)
      *pntR=1.;
    pntR++;
    *pntX++;
    }
  // BK 01-Nov-2000:  better: ??
  //for (real4 *pntX=&X[0][0]; pntX<=&X[X.lines()-1][X.pixels()-1]; ++pntX)
  //  *pntR++ = (abs(*pntX)<=0.5) ? 1 : 0;
  // rather something like: R(find(abs(X)<.5))=1;
  return Res;
  }


/****************************************************************
 * myhamming                                                    *
 * hamming window, lying vector                                 *
 * w = (a + (1.-a).*cos((2.*pi/fs).*fr)) .* myrect(fr./Br);     *
 * scale/shift filter by g(x)=f((x-xo)/s)                       *
 * alpha==1 yields a myrect window                              *
 *    Bert Kampes, 31-Mar-2000                                  *
 ****************************************************************/
inline matrix<real4> myhamming(
        const matrix<real4> &fr,
        real8 RBW,
        real8 RSR,
        real8 alpha)
  {
  TRACE_FUNCTION("myhamming (BK 31-Mar-200)")
  #ifdef __DEBUG
    if (fr.lines()!=1)
      {
      PRINT_ERROR("myhamming: only lying vectors.")
      throw(argument_error);
      }
    if (alpha<0.0 || alpha>1.0)
      {
      PRINT_ERROR("myhamming: !alpha e{0..1}.")
      throw(argument_error);
      }
    if (RBW>RSR)
      {
      PRINT_ERROR("myhamming: RBW>RSR.")
      throw(argument_error);
      }
  #endif
  matrix<real4> Res(1,fr.pixels());
  real4 *pntfr = fr[0];
  real4 *pntR  = Res[0];
  for (register int32 i=0; i<int32(fr.pixels()); ++i)
    {
    if (abs((*pntfr)/RBW)<0.5)          // rect window
      *pntR = alpha+(1.0-alpha)*cos((2*PI/RSR)*(*pntfr));
    pntfr++;
    pntR++;
    }
  return Res;
  }


/****************************************************************
 *    interpbilinear                                            *
 * bilinear interpolation with NORMALIZED data on grid.         *
 * f(row,col) = a00 + a01*c + a10*r +a11*r*c                    *
 * P1(0,0); P2(0,1); P3(1,0); P4(1,1); P[0:1,0:1];              *
 *                                                              *
 *         cols->                                               *
 *      r P1    P2                                              *
 *      o    .P                                                 *
 *      w P3    P4                                              *
 *      s                                                       *
 *                                                              *
 * input:                                                       *
 *  - normalized coordinates of point (row,col).                *
 *  - functional values at 4 corners.                           *
 * output:                                                      *
 *  - interpolated value                                        *
 *                                                              *
 *    Bert Kampes, 19-Jan-2000                                  *
 #%// BK 25-Sep-2000
 ****************************************************************/
inline real8 interpbilinear(
        const real8 r,
        const real8 c,
        const real8 P1,
        const real8 P2,
        const real8 P3,
        const real8 P4)         
  {
  #ifdef __DEBUG
  TRACE_FUNCTION("interpbilinear() (BK 25-Sep-2000)")
  #endif
  return P1 + (P2-P1)*c + (P3-P1)*r + (P4-P3-P2+P1)*r*c;
  }


/****************************************************************
 *    interp_trilinear                                          *
 * trilinear interpolation using relative coordinates.          *
 *                                                              *
 * ====== Points A,B,C and do trigrid lin. interpolation at "+" *
 * ====== "+", origin of local coordinate system                *
 *                                                              *
 *     A                                                        *
 *                                                              *
 *      +                                                       *
 *  B                                                           *
 *                                                              *
 *              C                                               *
 *                                                              *
 * ====== Then, cn A.z = f(x,y) = a00+a10*A.x+a01*A.y;          *
 * ======       cn B.z = f(x,y) = a00+a10*B.x+a01*B.y;          *
 * ======       cn C.z = f(x,y) = a00+a10*C.x+a01*C.y;          *
 * ====== solve for a00; the Z value at origin, i.e., it!       *
 *                                                              *
 * input:                                                       *
 *  - normalized coordinates of point (row,col).                *
 *  - cn A,B,C with relative coordinates and functional values  *
 * output:                                                      *
 *  - interpolated value                                        *
 *                                                              *
 * Bert Kampes, 07-Apr-2005                                     *
 ****************************************************************/
inline real8 interp_trilinear(
        const cn &A,
        const cn &B,
        const cn &C)
  {
  #ifdef __DEBUG
  TRACE_FUNCTION("interp_trilinear() (BK 07-Apr-2005)")
  if (A.x==B.x && A.y==B.y)
    {
    WARNING.print("interp_trilinear: A==B");
    return (A.z+B.z+C.z)/3.0;
    }
  if (A.x==C.x && A.y==C.y)
    {
    WARNING.print("interp_trilinear: A==C");
    return (A.z+B.z+C.z)/3.0;
    }
  if (B.x==C.x && B.y==C.y)
    {
    WARNING.print("interp_trilinear: B==C");
    return (A.z+B.z+C.z)/3.0;
    }
  // --- check if points are on line or are all same... ---
  // ??? how
  // --- maybe simply add eps to A.x, A.y and subtract it from B.x, B.y?
  // --- since how else can we guarantee solve33 won't crash?
  // we assume the caller routines guarantees this (as referencephase.cc does)
  // for speed
  #endif

  // --- Vandermonde matrix ---
  //  const real8 EPS = 0.00001;// assume points are pixels
  //  matrix<real8> DESIGNMAT(3,3); 
  //  DESIGNMAT(0,0)=1.0;  DESIGNMAT(0,1)=A.x + EPS;  DESIGNMAT(0,2) = A.y + EPS;
  //  DESIGNMAT(1,0)=1.0;  DESIGNMAT(1,1)=B.x - EPS;  DESIGNMAT(1,2) = B.y - EPS;
  //  DESIGNMAT(2,0)=1.0;  DESIGNMAT(2,1)=C.x;        DESIGNMAT(2,2) = C.y;
  //  // --- Right hand side ---
  //  matrix<real8> RHS(3,1); 
  //  RHS(0,0)=A.z;
  //  RHS(1,0)=B.z;
  //  RHS(2,0)=C.z;
  //  // --- Solve system using library for all coefficients ---
  //  matrix<real8> RESULT(3,1);// returned
  //  solve33(RESULT,RHS,DESIGNMAT);// return result = a00,a10,a01
  //  const real8 a00 = RESULT(0,0);
  //  return a00;

  // since first the last element is solved, turn it around and do it locally:
  //* ====== Then, cn A.z = f(x,y) = [ya,xa,1][a01]
  //* ======       cn B.z = f(x,y) = [yb,xb,1][a10]
  //* ======       cn C.z = f(x,y) = [yc,xc,1][a00] <--- solve this on by LU
  //* ====== solve for a00; the Z value at origin, i.e., it!
  // ______ LU decomposition ______
  const real8 EPS = 0.00000001;// assume points are pixels (zero division if A.y==0.0)
  const real8 L10 = (B.y-EPS)/(A.y+EPS);// DESIGN(1,0)/DESIGN(0,0);
  const real8 L20 =  C.y/(A.y+EPS);// DESIGN(2,0)/DESIGN(0,0);
  const real8 U11 = (B.x-EPS)-L10*(A.x+EPS);// DESIGN(1,1)-L10*DESIGN(0,1);
  // L21 had a bug in placing of brackets before v3.16
  // Bert Kampes, 25-Aug-2005
  const real8 L21 = (C.x-(A.x+EPS)*L20)/U11;// (DESIGN(2,1)-(DESIGN(0,1)*L20))/U11;
  const real8 U12 =  1.0-L10;// DESIGN(1,2)-L10*DESIGN(0,2);
  const real8 U22 =  1.0-L20-L21*U12;// DESIGN(2,2)-L20*DESIGN(0,2)-L21*U12;

  // ______ Solution: forward substitution ______
  const real8 b0  =  A.z;// rhs(0,0);
  const real8 b1  =  B.z-b0*L10;// rhs(1,0)-b0*L10;
  const real8 b2  =  C.z-b0*L20-b1*L21;// rhs(2,0)-b0*L20-b1*L21;

  // ______ Solution: backwards substitution ______
  const real8 a00 =  b2/U22;//RESULT(2,0) = b2/U22;
  //RESULT(1,0) = (b1-U12*RESULT(2,0))/U11;
  //RESULT(0,0) = (b0-DESIGN(0,1)*RESULT(1,0)-DESIGN(0,2)*RESULT(2,0))/DESIGN(0,0);
  return a00;
  }



/****************************************************************
 *    ProjectPointOnLine                                        *
 * Computes Projection of a 3D point M                          *
 * on to a vector in 3D. The vector is directed                 *
 * from  P to S                                                 *
 * this can be useed in baseline computation                    *
 #%// KKM 09-March-2005                                         *
 ****************************************************************/
inline cn ProjectPointOnLine(
        cn M,
        cn P,
        cn S)
  {
  TRACE_FUNCTION("ProjectPointOnLine() (KKM 09-Mar-2005)")
  cn U            = P.min(S);// direction cosines of PS
  const real8 U_2 = U.in(U);
  U               = U.scale(1.0/sqrt(U_2));
  // --- Umat is the matrix representation of U to compute U*U' ---
  matrix<real4> Umat(3,1);
  Umat(0,0)       = U.x;
  Umat(1,0      ) = U.y;
  Umat(2,0)       = U.z;
  const matrix<real4> UU = matxmatT(Umat,Umat);// [3x3]
  // --- X is returned, base of perpendicular from M to PS
  matrix<real4> UmatT(1,3);
  UmatT(0,0)       = M.x-S.x;
  UmatT(0,1)       = M.y-S.y;
  UmatT(0,2)       = M.z-S.z;
  UmatT            = UmatT*UU;// [1x3]
  // --- returned cn, base of perpendicular from M to PS ---
  return cn(UmatT(0,0)+S.x, UmatT(0,1)+S.y, UmatT(0,2)+S.z);
  }



/****************************************************************
 * matrix<real4> I = indgen(real4(10))                          *
 * matrix<int32> I = indgen(int32(10)), etc.                    *
 * create lying matrix [0,1,2,..N-1]                            *
 * Bert Kampes, 08-Apr-2005                                     *
 ****************************************************************/
template <class Type>
inline 
matrix<Type> indgen(
        const Type N) // use Type of N to create appropriate matrix
  {
  TRACE_FUNCTION("indgen() (BK 08-Apr-2005)")
  const uint dim = uint(N+0.5);// round positive number
  matrix<Type> e(1,dim); 
  for(uint i=0; i<dim; ++i) e(0,i)=Type(i);
  return e;
  }


/****************************************************************
 * matrix<real4> I = linspace(x0,xn,N)                          *
 * create lying matrix [x0:xn] equally space N samples          *
 * Bert Kampes, 08-Apr-2005                                     *
 ****************************************************************/
template <class Type>
inline 
matrix<Type> linspace(
        const Type x0,
        const Type xN,
        const int32 N) 
  {
  TRACE_FUNCTION("linspace() (BK 08-Apr-2005)")
  return x0+((indgen(Type(N)))*((xN-x0)/Type(N-1.0)));// use Type(N) to get ok.
  }


/****************************************************************
 * matrix<real4> E = ones(10,10)                                *
 * create matrix filled with ones.                              *
 * Bert Kampes, 08-Apr-2005                                     *
 ****************************************************************/
template <class Type>
inline
matrix<Type>  ones(
        const Type Ny,
        const Type Nx)// use Type of N to create appropriate matrix 
  {
  TRACE_FUNCTION("ones() (BK 08-Apr-2005)")
  const uint dimX = uint(Nx+0.5);// round positive number
  const uint dimY = uint(Ny+0.5);// round positive number
  matrix<Type> E(dimY,dimX); 
  E.setdata(Type(1));
  return E;
  }


/****************************************************************
 * matrix<real4> I = eye(10)                                    *
 * create identity matrix with ones on main diagonal.           *
 * Bert Kampes, 08-Apr-2005                                     *
 ****************************************************************/
template <class Type>
inline
matrix<Type>  eye(
        const Type N)// use Type of N to create appropriate matrix
  {
  TRACE_FUNCTION("eye() (BK 08-Apr-2005)")
  const uint dim = uint(N+0.5);// round positive number
  matrix<Type> I(dim,dim);
  for(uint i=0; i<dim; ++i) I(i,i)=Type(1);
  return I;
  }



// ====== Prototypes ============================================
// ______ Call getorb for precise orbits ______
void getorb(
        const slcimage &image,
        const input_pr_orbits &inputorb,
        int16 FILEID);


// ______ Converts output to desired format ______
void convertgetorbout(
        int16 FILEID,
        const char* file="scratchorbit");

// ______ Return next power of 2 ______
uint nextpow2 (
        real8 w);


// ______ Evaluate 2d polynomial on grid ______
template <class Type>
matrix<Type> polyval(                        // MA overloading func. for more precision like real8
        const matrix<real4> &x,
        const matrix<real4> &y,
        const matrix<real8> &coeff,
        int32 degree = -1);


// ______ Evaluate 2d polynomial ______
real8 polyval(
        real8 x,
        real8 y,
        const matrix<real8> &coeff);


// ______ Evaluate 2d polynomial ______
real8 polyval(
        real8 x,
        real8 y,
        const matrix<real8> &coeff,
        int32 degree);


// ______ Evaluate 1d polynomial ______
real8 polyval1d(
        real8 x,
        const matrix<real8> &coeff);
//        int32 degree);


// ______ Evaluate derivative of 1d polynomial [HB] ______
real8 polyval1d(
        real8 x,
        const matrix<real8> &coeff,
        uint deriv);


// ______ Normalize a matrix data -> [-2,2] ______
void normalize(
        matrix<real4> &data,
        real8 min,
        real8 max
        );


// ______ Normalize a matrix data -> [-2,2] ______
void normalize(
        matrix<real8> &data,
        real8 min,
        real8 max
        );


// ______ Normalize data from [min,max] to [-2,2] ______
inline real4 normalize(real4 data, real8 min, real8 max)
  {return real4((data-(0.5*(min+max)))/(0.25*(max-min)));}


// ______ Normalize data from [min,max] to [-2,2] ______
inline real8 normalize(real8 data, real8 min, real8 max)
  {return (data-(0.5*(min+max)))/(0.25*(max-min));}


// ______ Use this to compute Bperp for a few points, slow! ______
void BBparBperp(
        real8 &B, real8 &Bpar, real8 &Berp,     // returned
        const cn Master, const cn Point, const cn Slave);

// ______ Sign not computed, slow! ______
void BBhBv(
        real8 &B, real8 &Bh, real8 &Bv,         // returned
        const cn Master, const cn Slave);

// ______ Temporal baseline in days returned ______
// ______ (positive for master earlier than slave) ______
int32 Btemp(
        const char *utc_master, const char *utc_slave);

// ______ Baseline parametrizations returned ______
void BalphaBhBvBparBperpTheta(
        real8 &B, real8 &alpha, real8 &Bh, real8 &Bv,
        real8 &Bpar, real8 &Bperp, real8 &theta,
        const cn M, const cn P, const cn S);



// ______ Shift azimuth spectrum from fDC to zero, and vv. ______
void shiftazispectrum(
        matrix<complr4> &data,  // slcdata sin space domain
        const slcimage  &slave, // polynomial fDC
        const real4      shift);// abs(shift==firstpix in slave system)
                                // if shift>0 then shift from zero to fDC

// ______ Dump tiepoint INFO ______
void tiepoint(
        const input_gen     &generalinput,
        const slcimage      &master,
        const slcimage      &slave,
              orbit         &masterorbit,
              orbit         &slaveorbit,
        const input_ell     &ellips);

// ______ Line,Pixel offsets2timing ______
void offsets2timing(
        const slcimage        &slcinfo, 
        const int32           &offsetL,
        const int32           &offsetP,
        real8                 &azTIME,
        real8                 &rgTIME);

// _________  Write_KML ___________ [BO]
void write_kml(
        const input_ell        &ellips,
        const window           &input_dbow,
        const slcimage         &master,
        orbit                  &masterorbit,
	const char		       filename[4*ONE27]);

// ====== Fast trigonometry by LookUpTables ========================
// ______ Standard functions can be switched on using "-D NO_FASTTRIG"
// ______ Bert Kampes, 29-Sep-2005
#ifdef NO_FASTTRIG
//#warning "NO_FASTTRIG used" // [MA] informative
  // inline real4 fast_sin(const real4 &x)                    {return sin(x);}
  // inline real4 fast_min_sin(const real4 &x)                {return -sin(x);}
  // inline real4 fast_cos(const real4 &x)                    {return cos(x);}
  // inline real4 fast_atan(const real4 &x)                   {return atan(x);}
  // inline real4 fast_atan2(const real4 &y, const real4 &x)  {return atan2(y,x);}
  // inline real4 fast_arg(const complex<real4> &z)           {return arg(z);}
  #define fast_sin     sin             // [MA] this is better definition
  #define fast_min_sin -sin
  #define fast_cos     cos
  #define fast_atan    atan
  #define fast_atan2   atan2
  #define fast_arg     arg
#else
// ______ Global variables defined in utilities.cc ______
extern const uint32 LUT_LENGHT;          // 65536
extern const uint16 LUT_ATAN_MAX;        // 128
extern const real8 LUT_DX_SIN;           // idx=int(x*LUT_DX_SIN)
extern const real8 LUT_DX_ATAN;          // ..
extern real8 LUT_SIN[];                  // look-up table [0:dx:2pi-dx]
extern real8 LUT_ATAN[];                 // look-up table [-128:dx:128-dx]


// --- Initialize LUTs ---
// --- LUT for 0:dx; dx:2dx; Ndx:2pi-dx, etc; eval at start bin  ---
// --- This function should be called only once at start program ---
// Bert Kampes, 04-Oct-2005
inline void init_fasttrig()
  {
  for(int32 i=0; i<LUT_LENGHT; ++i)
    {
    LUT_SIN[i]  = sin( real8(i) / LUT_DX_SIN );         // eval at i  // MA type is upgrade to real8
    LUT_ATAN[i] = atan(real8(i-LUT_LENGHT/2)/LUT_DX_ATAN);
    }
  }// init_fasttrig


// --- fast_sin()
// --- fast for pos and neg x; not best for given LUT_LENGTH.
// --- Bert Kampes, 04-Oct-2005
inline real8 fast_sin(const real8 &x)
  {
   //return LUT_SIN[ int32(x*LUT_DX_SIN) & (LUT_LENGHT-1) ] ;               // original line of code
                                                                            // [MA]; binary operation on the rhs garantees to wrap the index to the size of LUT_SIN.
                                                                            // [MA] use unsigned int rather than int, we don't want rewrapping to negative values.
   return LUT_SIN [ uint64(x*LUT_DX_SIN) & (LUT_LENGHT-1) ];                // 1 soln. 
   //return LUT_SIN [ uint32((fmod(x,2*PI)*LUT_DX_SIN)) & (LUT_LENGHT-1) ]; // 2 soln. great but computational intensive due to fmod
   //return LUT_SIN [ uint16(x*LUT_DX_SIN) ];                               // 3 bad soln. max positive value 65535 but types are different won't wrap as expected. 
   }


// --- fast_min_sin(x)=-sin(x)=sin(x+pi)
// --- Bert Kampes, 04-Oct-2005. revisited in 20090602
inline real8 fast_min_sin(const real8 &x)
  { return LUT_SIN[ ( uint64(x*LUT_DX_SIN) + LUT_LENGHT/2 ) & (LUT_LENGHT-1) ]; }


// --- fast cosine using LUT
// --- Bert Kampes, 04-Oct-2005
inline real8 fast_cos(const real8 &x)
  {return LUT_SIN[ ( uint64(x*LUT_DX_SIN) + LUT_LENGHT/4 ) & (LUT_LENGHT-1) ]; }  


// --- fast_atan() using LUT
// --- Bert Kampes, 04-Oct-2005
inline real8 fast_atan(const real8 &x)
  {
  if (x<-LUT_ATAN_MAX)
    return -PI/2;
  else if (x>LUT_ATAN_MAX)
    return  PI/2;
  else
    return LUT_ATAN[uint64(x*LUT_DX_ATAN)+LUT_LENGHT/2];
  }


// --- fast_atan2() fast four-quadrant atan2(y,x) by LUT ---
// --- Bert Kampes, 04-Oct-2005
inline real8 fast_atan2(const real8 &y, const real8 &x)
  {
  if (x==0.0) return (y==0.0) ? 0.0 : (y<0.0) ? -PI/2 : PI/2;
  if (x<0.0)  return (y<0.0) ? fast_atan(y/x)-PI : fast_atan(y/x)+PI;
  else        return fast_atan(y/x);
  }


// --- fast_arg (angle of complex number) ---
// --- Bert Kampes, 04-Oct-2005
inline real8 fast_arg(const complex<real8> &z)
  {
  return fast_atan2(imag(z),real(z));
  }
#endif // NO_FASTTRIG






//  /****************************************************************
//   *    pol2xyz                                                       *
//   *                                                          *
//   * Converts         lat,lon,hei(above surface) (sphere)             * 
//   *         to       x,y,z (cartesian)                               *
//   * uses constants in constants.h (pi, radius)                       *
//   *                                                          *
//   * lat is positive for North, lon is negative for West.             *
//   *                                                          *
//   *    Bert Kampes, 21-Dec-1998                                      *
//   ****************************************************************/
//  void pol2xyz(
//      cn &xyz,
//      real8 phi,
//      real8 lon,
//      real8 height) //default =0)
//    {
//    TRACE_FUNCTION("pol2xyz (BK 21-Dec-1998)");
//    real8 Rph = RADIUS + height;
//    xyz.x = Rph *cos(phi)*cos(lon);
//    xyz.y = Rph *cos(phi)*sin(lon);
//    xyz.z = Rph *sin(phi);
//    }  // END pol2xyz
//  
//  
//  /****************************************************************
//   *    xyz2pol                                                       *
//   *                                                          *
//   * Converts x,y,z (cartesian)                                       *
//   *           to lat,lon,hei(above surface) (sphere)         * 
//   * uses constants in constants.h (pi, radius)                       *
//   *                                                          *
//   * lat is positive for North, lon is negative for West.             *
//   *                                                          *
//   *    Bert Kampes, 21-Dec-1998                                      *
//   ****************************************************************/
//  void xyz2pol(
//      const cn &xyz,
//      real8    &phi,
//      real8    &lambda,
//      real8    &height)
//    {
//    TRACE_FUNCTION("xyz2pol (BK 21-Dec-1998)");
//    real8 r = sqrt(sqr(xyz.x)+sqr(xyz.y));
//    phi     = atan2(xyz.z,r);
//    lambda  = atan2(xyz.y,xyz.x);
//    height  = (r/cos(phi)) - RADIUS;
//    }  // END xyz2pol
//  


#endif // UTILITIES_H



