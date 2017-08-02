/*
 * Copyright (c) 1999-2012 Delft University of Technology, The Netherlands
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
 */
/****************************************************************
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/constants.hh,v $  *
 * $Revision: 3.33 $                                            *
 * $Date: 2005/10/18 14:15:26 $                                 *
 * $Author: TUDelft $                                           *
 *                                                              *
 * Definitions of some global (external) structs +functions,    *
 * constants                                                    *
 * define: SWNAME SWVERSION                                     *
 * __GplusplusCOMPILER__ then implement T sqr(T), cause no ansi c++ *
 * pragma: aCC things (removed since v3.8)                      * 
 *                                                              *
 * CLASSES IN THIS FILE:                                        *
 *   - ellips:  struct + functions                              *
 *   - cn:      3D coordinates + functions                      *
 *   - window:  in matrices/files + functions                   *
 ****************************************************************/


#ifndef CONSTANTS_H
#define CONSTANTS_H

using namespace std;                    // BK 29-Mar-2003, new compiler?
                                        // TODO: SLiu avoid namespace declaration in the header file.

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "bk_messages.hh"
#include "exceptions.hh" // [MA]
#include "cstring"              // strcpy for exceptions

//#if defined __GNUC__ && __GNUC__ >= 2
//#else
//#endif
// g++ -v -c processor.cc gives (ao):
// -D__GNUC__=2 -D__GNUG__=2 -D__GNUC_MINOR__=95
#include <cmath>                // for pi (atan)
#include <complex>              // default known complex type
#include <iostream>             // cout etc.
#include <strstream>            // for memory stream



// ====== Globals for messages to stdout ======
extern bk_messages TRACE;
extern bk_messages DEBUG;
extern bk_messages INFO;
extern bk_messages PROGRESS;
extern bk_messages WARNING;
extern bk_messages ERROR;
// _____ Macros for abbreviation ______
// ___ Note that useage could be like PRINT_ERROR(WARNING.get_str()) ___
// ___ (normally WARNING.rewind() required after this, but exited) ___
// ___                                PRINT_ERROR(".") ___
// ___ w/o the ";"  (since in an if/else construction without {}{}
// ___ it will otherwise crash ___
// ___ use TRACE_FUNCTION ONLY SIMPLE WITH CHAR ARRAY "dadsa" ___
#define TRACE_FUNCTION(s) {TRACE.reset();\
  TRACE<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<s<<ends;\
  TRACE.print();}//

#define PRINT_ERROR(s) {ERROR.terminate(); char cp_s[256]; strcpy(cp_s,s); ERROR.reset();\
  ERROR<<"["<<__FILE__<<"["<<__LINE__<<"]]: "<<cp_s<<ends; ERROR.print();}//


// ====== Define time macro =======  // [MA]
#ifndef __TIMESTAMP__
#define __TIMESTAMP__ __DATE__ " " __TIME__
#endif

// ====== Define name/version ======
#define SWNAME    "Doris (Delft o-o Radar Interferometric Software)"
//#define SWVERSION "Version 2.0"       // 4-Nov-99
//#define SWVERSION "Version 2.2"       // 23-Feb-00 DEM subtraction
//#define SWVERSION "Version 2.3"       // 05-Apr-00 rangefilter
//#define SWVERSION "Version 2.4"       // 22-Jun-00 makefile for g++ all versions
//#define SWVERSION "Version 2.5"       // 07-Jul-00 linux x86 byteorder
//#define SWVERSION "Version 2.5.1"     // 19-Jul-00 orbit class
//#define SWVERSION "Version 2.6"       // 26-Oct-00 dinsar
//#define SWVERSION "Version 3.0"       // ??-???-00 filtphase, filtazi, filtrange2
//#define SWVERSION "Version 3.1"       // 04-May-01 DEM radarcode formats, general impr.
//#define SWVERSION "Version 3.3"       // 06-FEb-02 BEEP, PREVIEW,
                                        // Cygwin, COMPREFDEM (april).
//#define SWVERSION "Version 3.4 (01-NOV-2002)" // azishift resample change
                                                // ...snaphu unwrap
//#define SWVERSION "Version 3.5 (03-APR-2003)" // compatible with g++ v3.2 c++ standard
//#define SWVERSION "Version 3.6 (15-JUN-2003)" // envisat interface.
//#define SWVERSION "Version 3.7 (02-SEP-2003)" // auto coregpm. example scripts geocoding
//#define SWVERSION "Version 3.8 (22-SEP-2003)" // added fftw lib, quicklook option of runscript
//#define SWVERSION "Version 3.9 (16-OKT-2003)" // mainly release under GPL
//#define SWVERSION "Version 3.10 (16-JAN-2004)" // ovs, resample, ...
//#define SWVERSION "Version 3.11 (13-MAY-2004)" // envisat bugfix, height card
//#define SWVERSION "Version 3.12 (13-AUG-2004)" // mainly radarsat added
//#define SWVERSION "Version 3.13 (01-OCT-2004)" // fleafixes
//#define SWVERSION "Version 3.14 (07-MAR-2005)" // fixes fdc computation; .h -> .hh
//#define SWVERSION "Version 3.15 (14-JUL-2005)" // baseline class, exceptions
//#define SWVERSION "Version 3.16 (01-SEP-2005)" // OVS az.; small bugs; coreg better
//#define SWVERSION "Version 3.17 (05-NOV-2005)" // code reorg. g++-4.0 compiler
//#define SWVERSION "Version 3.18test3 (11-JUL-2006)" // 
//#define SWVERSION "Version 3.18test5 (26-FEB-2007)" // alos, bug.fix, rsmp.dbow.geo
//#define SWVERSION "Version 3.19test8 (30-JUL-2007)" // bug.fix, dem.interpolator.fix                                                           
//#define SWVERSION "Version 3.20 (02-JAN-2008)" // bug.fix, coregistration               
//#define SWVERSION "Version 3.94 (22-AUG-2007)" // Freek: DEM assistant coregistration [LG] + timing error   
//#define SWVERSION "Version 3.95 (21-AUG-2008)" //   PM & MA & GJ: bug fixes , compiler support updates, counter bug, > 4GB filesupport in 32-bit
//#define SWVERSION "Version 3.96 (06-SEP-2008)" // merge v3.20 and v3.95
//#define SWVERSION "Version 3.97 (29-OCT-2008)"  // Don, MA, Fvl: DEMA bug fixed
//#define SWVERSION "Version 3.98 (06-Nov-2008)" // Don coherence est. w/ removal topo. slope, BO & FvL & MA: simamp + master_timing, MA: several additional fixes
//#define SWVERSION "Version 3.99 (23-Dec-2008)" // wrapping up updates
//#define SWVERSION "Version  4.01 (24-Dec-2008)" //  release candidate v4.01 rc3, ALOS and RADARSAT1 fixes, some improvement in SIMAMP, casting fixes for compilers
//#define SWVERSION "Version  4.01_rc8 (21-Mar-2009)" //  release candidate v4.01 rc8 
//#define SWVERSION "Version  4.01_rc9 (27-Mar-2009)" //  release candidate v4.01 rc9 
//#define SWVERSION "Version  4.01_rc10 (13-Apr-2009)" // release candidate v4.01 rc10 PM TERRASAR-X sensor
//#define SWVERSION "Version  4.01 (14-Apr-2009)" //  release v4.01 finally! :-)
//#define SWVERSION "Version  4.02 (08-Jun-2009)" //  release v4.02, fixes on comprefdem buffer when multilooked [HB]; [MA] lapack warnings fix and lib versions, fix on sinus lut and NO_FASTTRIG
//#define SWVERSION "Version  4.03-beta (14-Jul-2009)" //  release v4.03, 
//#define SWVERSION "Version  4.03-beta2 (22-Aug-2009)" //  release v4.03, rf_nlmean fix, mte_distwinfile fix.
//#define SWVERSION "Version  4.03-beta3 (06-Sep-2009)" //  fix readres at oversample slave
//#define SWVERSION "Version  4.03-beta4 (06-Oct-2009)" //  added Radarsat-2 reader
//#define SWVERSION "Version  4.03-beta5 (08-11-2009)" //  added Radarsat-2 reader
//#define SWVERSION "Version  4.03-beta7 (04-01-2010)" //  fix productfill for long pathnames, MCC Envisat alternating polarization update.
//#define SWVERSION "Version  4.03-beta8 (09-04-2010)" "  Build: " __TIMESTAMP__ //  [MA] add SAM: demhei_lp and theta_lp
//#define SWVERSION "version  4.03-beta9 (15-04-2010)" "\n\t\t     build \t" __TIMESTAMP__ //  [MA] crosscorrelate: drop correlation >1 affects mtiming, coarsecorr(magfft), filename update for DEMA, CRD, screen info update
//#define SWVERSION "version  4.03-beta10 (23-06-2010)" "\n\t\t     build \t" __TIMESTAMP__ //  [MA] TSX heading info, [PD] CSK reader, SLIU important comments
//#define SWVERSION "version  4.03-beta11 (16-09-2010)" "\n\t\t     build \t" __TIMESTAMP__ //  [BO] GAMMA SLC, KML Generation, TSX dumpheader bugfix.
//#define SWVERSION "version  4.03-beta12 (11-11-2010)" "\n\t\t     build \t" __TIMESTAMP__ //  [MA] Read only the parameters of the step to run, not all ex: for simamp is done, see products.cc
                                                                                          //   coarsewindow size updated, simamp multilook factor is fixed to 1 due to Doris convention. defined memory_max in constants.hh, read necessart inputs for processing not all see simamp part in readinput.cc
//#define SWVERSION "version  4.03-beta13 (01-12-2010)" "\n\t\t     build \t" __TIMESTAMP__ // [MA] fix to Read only the parameters issues on simamp, skip loop1 update. 
//#define SWVERSION "version  4.04-beta1 (06-12-2010)" "\n\t\t     build \t" __TIMESTAMP__ // start to merge two branches 
//#define SWVERSION "version  4.04-beta2 (05-01-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [MA] ERS_N1 support Doris part
//#define SWVERSION "version  4.04-beta3 (07-02-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [MA] ERS_N1 support fix
//#define SWVERSION "version  4.04-beta4 (10-03-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [MA] master printout crop numlines+numpixels
//#define SWVERSION "version  4.05-beta1 (28-03-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [HB] ESTORBIT module
//#define SWVERSION "version  4.06-beta1 (23-10-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [BA] Modified Goldstein Filter
//#define SWVERSION "version  4.06-beta2 (28-12-2011)" "\n\t\t     build \t" __TIMESTAMP__ // [MA] path length fix for demassist at productfill and the rest.
//#define SWVERSION "version  4.06-beta3 (31-10-2013)" "\n\t\t     build \t" __TIMESTAMP__ // [FvL] path length fix in readcoeff function (ioroutines.cc).
//#define SWVERSION "version  4.0.7 (23-07-2014)" "\n\t\t     build \t" __TIMESTAMP__ // [FvL] removed unwanted automatic removal of pre-calculated reference phase in INTERFERO step
#define SWVERSION "version  4.0.8 (04-09-2014)" "\n\t\t     build \t" __TIMESTAMP__ // [FvL] new version based on svn trunk with still untested spotlight developments.

// ====== Typedefs for portability ======
typedef short int           int16;    // 16 bits --> 2 bytes.  It has a range of -32768 to 32767. [ from -2^15 to (2^15 - 1) ]  [MA]
typedef unsigned short int  uint16;   // 16 bits --> 2 bytes.  It has a range of 0 to 65535.      [ 2^16 - 1]
typedef int                 int32;    // 32 bits --> 4 bytes.  It has a range of -2,147,483,648 to 2,147,483,647, [ -(2^31) to (2^31 - 1) ]
typedef unsigned int        uint32;   // 32 bits --> 4 bytes.  It has a range of 0 to 4,294,967,295. [ 2^32 - 1 ]
typedef unsigned int        uint;     // sizeof (int) = sizeof (long) on 32bit data model; sizeof (int) < sizeof (long) on 64 bit.  
typedef long long           int64;    // 64 bits --> 8 bytes.  It has a range of -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807 . [ -(2^63) to 2^63 - 1]
typedef unsigned long long  uint64;   // 64 bits --> 8 bytes.  It has a range of  0 to 18,446,744,073,709,551,615. [ 2^64 - 1 ] 
typedef float               real4;    // 7 digits  x.1xxxx6e(-126 and +127)
typedef double              real8;    // 16 digits x.1xxxxxxxxxxxxx15(-1022  and +1023)
typedef complex<int16>      compli16;
typedef complex<int32>      compli32;
typedef complex<real4>      complr4;
typedef complex<real8>      complr8;
#if defined (__GNUC__) && __GNUC__ >= 4
typedef long double         real16;   // 34 digits x.1xxxxxxxxxxxxx33(-16382 and +16383) // quadruple precision  [MA]
//# else
# endif  /* GCC.  */




// ====== function: Type sqr(Type) not in stdlib math.h! ======
// ====== Possibly also for other compilers that have errors this define will fix it ===
// note though that this is not save-casted
#ifdef __GNUC__
inline int16  sqr(const int16 &x)  {return(x*x);}
inline uint16 sqr(const uint16 &x)  {return(x*x);}
inline int32  sqr(const int32 &x)  {return(x*x);}
inline int64  sqr(const int64 &x)  {return(x*x);}
inline uint64 sqr(const uint64 &x) {return(x*x);}
inline uint   sqr(const uint  &x)  {return(x*x);}
inline real4  sqr(const real4 &x)  {return(x*x);}
inline real8  sqr(const real8 &x)  {return(x*x);}
#endif


// For windows, it seems max,min is already defined, but with different arguments.
// I cannot believe this, but Jia changed all call`s to min/max to _min _max (?)
// He also changed the names of DEBUG to DEBUG1 and INFO to INFO1, etc.
// in processor.cc:
// Bert Kampes, 24-Aug-2005
//#ifdef WIN32
//#define max  _MAX
//#define min  _MIN
//#endif




// ====== Global constants ======
// ______ General constants ______
const real8     SOL             = 299792458.0;  // speed of light in m/s
const real8     EPS             = 1e-13;        // small number
const int32     NaN             = -999;         // Not a Number
const int32     Inf             =  99999;       // Infinity
#ifdef M_PI
const real8     PI              = M_PI;// conform to standard value in math.h
#else
const real8     PI              = real8(4)*atan(real8(1));
#endif
//const int16           EIGHTY          = 80;
const int16     EIGHTY          = 200;// 80 seems to short for ASAR filenames
const int16     ONE27           = 127;

// doris memory for buffer operations
const uint32    MEMORY_DEF      = 500;          // default memory (~500MB), unit Mbits [MA]
const uint32    MEMORY_MAX      = 12000;        // maximum memory (~12GB) 

// ______ Handy for software, IDs have to be enumerated (used in loop) ______
const int16     LOGID           = 1;            // general identifier for log
const int16     MASTERID        = 2;            // general identifier for master
const int16     SLAVEID         = 3;            // general identifier for slave
const int16     INTERFID        = 4;            // general identifier for interferogram
const int16     FORMATCI2       = 11;           // formatflag
const int16     FORMATCR4       = 12;           // formatflag
const int16     FORMATR4        = 13;           // formatflag
const int16     FORMATI2        = 14;           // formatflag
const int16     FORMATI2_BIGENDIAN = 17;        // formatflag
const int16     FORMATHGT       = 15;           // formatflag unwrapped band interleaved ampl(1)/phase(2)
const int16     FORMATR8        = 16;           // formatflag
const int16     FORMATI4        = 18;           // formatflag // TODO MA revisit to order and out to enum
const int16     FORMATCR8       = 19;           // formatflag

// ______ product specifier ______
const int16     SLC_ERS         = 1;            // ERS
const int16     SLC_ERS_N1     = 11;            // ERS in N1 format
const int16     SLC_ASAR    = 2;                // ENVISAT
const int16     SLC_ASAR_AP_HH = 21;            // ENVISAT AP HH
const int16     SLC_ASAR_AP_VV = 22;            // ENVISAT AP VV

const int16     SLC_S1A     = 30;               //Sentinel-1 A (TOPS mode))
const int16     SLC_RSAT    = 3;                // RadarSAT (and Atlantis processor??)
const int16     SLC_JERS    = 4;                // JERS (ceos?)
// for ALOS: [PM]
const int16     SLC_ALOS     = 5;               // ALOS (ceos)   //  [MA] TODO what about different polarizations: single 5 --> 51; dual == 52 etc.
const int16     SLC_TSX      = 6;               // TSX stripmap
// for Radarsat-2: [MA]
const int16     SLC_RS2      = 70;              // RS2 [Default for Radarsat-2 SLC]  [MA]
const int16     SLC_RS2_QUAD = 71;              // RS2  71 --> QUADPOL
// for CSK: [PD]
const int16     SLC_CSK      = 80;              // CSK stripmap SLC  
const int16     SLC_CSK_POL  = 81;              // CSK POL template
const int16	SLC_GAMMA    = 9;		// GAMMA SLC. Only for readfiles step. [BO]

// ______ processor specifier ______
const int16     SARPR_VMP       = 11;            // ESA PAFs SAR processor SLC
const int16     SARPR_ATL       = 12;            // Atlantis created SLC
const int16     SARPR_TUD       = 13;            // TUDelft processor created SLC
// Modified by LG for reading ALOS Fine
const int16     SARPR_JAX       = 14;            // JAXA  processor created SLC
const int16     SARPR_TSX       = 15;            // TSX  processor created COSAR
const int16	SARPR_GAM	= 16;		 // GAMMA processor created SLC [BO]
const int16     SARPR_RS2       = 17;            // RS2 processor created ...
const int16     SARPR_CSK       = 18;            // CSK processor created H5

// ====== Unique method selectors for interpolation ======
// ====== For ASAR (only 5 points), spline is not good ======
const int16     ORB_SPLINE      = -11; // else degree of polynomial
const int16     ORB_DEFAULT     = -12; // [degree of polynomial=numpoints-1, max=5]
const int16     ORB_PRM_POS     = 30;            // POSITION default [MA] [TODO]
const int16     ORB_PRM_VEL     = 31;            // POSIVELO         [MA]


// ====== ELLIPSOID CLASS ======
// ______flatting is different: grs80: 1/298.257222101
// ______flatting is different: wgs84: 1/298.257223563
const   double  WGS84_A  = 6378137.000;           // semimajor axis wgs84
const   double  WGS84_B  = 6356752.3142;          // semiminor axis wgs84

//const   double  GRS80_A  = 6378137.000;           // semimajor axis grs80
//const   double  GRS80_B  = 6356752.3141;          // semiminor axis grs80


// ____ used in matric.cc: convert_type ____  [GJ] && [MA] 2009
inline int getformat( const int16   x ) { return FORMATI2;  } // cannot detect endianness !!
inline int getformat( const int32   x ) { return FORMATI4;  } // cannot detect endianness !!
inline int getformat( const real4   x ) { return FORMATR4;  }
inline int getformat( const real8   x ) { return FORMATR8;  }
inline int getformat( const complr4 x ) { return FORMATCR4; }
inline int getformat( const complr8 x ) { return FORMATCR8; }


// ====== class cn for coordinates + functions =======================
class cn                                // coordinates
  {
  public:
  real8         x,                      // 
                y,                      //
                z;                      //
  // ______ Constructors ______
  cn()
    {x=0.0; y=0.0; z=0.0;}
  cn(real8 Px, real8 Py, real8 Pz)
    {x=Px; y=Py; z=Pz;}
  // ______ Copy constructor ______
  cn(const cn& P)
    {x=P.x; y=P.y; z=P.z;}
  // ______ Destructor ______
  ~cn()
    {;}// nothing to destruct that isn't destructed automatically
  // ______ Operators ______
  inline cn& operator = (const cn &P)// assignment operator
    {if (this != &P) {x=P.x; y=P.y; z=P.z;} return *this;}
  inline cn& operator += (const cn &P)
    {x+=P.x; y+=P.y; z+=P.z; return *this;}
  inline cn& operator -= (const cn &P)
    {x-=P.x; y-=P.y; z-=P.z; return *this;}
  inline cn& operator *= (const real8 s)// multiply by scalar
    {x*=s; y*=s; z*=s; return *this;}
  inline cn& operator *= (const cn &P)
    {x*=P.x; y*=P.y; z*=P.z; return *this;}
  inline cn& operator /= (const real8 s)// divide by scalar
    {x/=s; y/=s; z/=s; return *this;}
  inline cn& operator /= (const cn &P)
    {x/=P.x; y/=P.y; z/=P.z; return *this;}
  inline cn operator + () // unary plus
    {return cn(x,y,z);}
  inline cn operator + (const cn &P) // binary plus
    {return cn(x,y,z)+=P;}
  inline cn operator - ()// unary min
    {return cn(-x, -y, -z);}
  inline cn operator - (const cn &P)// binary min
    {return cn(x,y,z)-=P;}
  inline cn operator * (const real8 s) // Q=P*s
    {return cn(x,y,z)*=s;}
  inline cn operator * (const cn &P) // Q=P*Q
    {return cn(x,y,z)*=P;}
  inline cn operator / (const real8 s) // Q=P/s
    {return cn(x,y,z)/=s;}
  inline cn operator / (const cn &P) // Q=P/Q
    {return cn(x,y,z)/=P;}
  inline bool operator == (const cn &P) const // isequal operator
    {return (P.x==x && P.y==y && P.z==z) ? true : false;}
  inline bool operator != (const cn &P) const // isnotequal operator
    {return (P.x==x && P.y==y && P.z==z) ? false : true;}
  // ______ Public function in struct ______
  inline real8 in(cn P) const           // scalar product: r=P.in(Q)
    {return x*P.x+y*P.y+z*P.z;}
  inline cn out(cn P) const             // cross product cn r=P.out(Q)
    {return cn(y*P.z-z*P.y, z*P.x-x*P.z, x*P.y-y*P.x);}
  inline real8 dist(cn P) const         // distance: d=P.dist(Q)
    {return sqrt(sqr(x-P.x)+sqr(y-P.y)+sqr(z-P.z));}
  // ___ WIN32: Jia defined following, but I don;t know why:
  // cn _min(cn P) const            // cn r=P.min(Q);         //  Jia YouLiang
  //   {return cn(x-P.x, y-P.y, z-P.z);}
  inline cn min(cn P) const             // cn r=P.min(Q)
    {return cn(x-P.x, y-P.y, z-P.z);}
  inline real8 norm2() const            // n=P.norm2()
    {return sqr(x)+sqr(y)+sqr(z);}
  inline real8 norm() const             // n=P.norm()
    {return sqrt(sqr(x)+sqr(y)+sqr(z));}
  inline cn normalize() const           // cn R=P.normalize()
    {return cn(x,y,z)/norm();}
    //{const real8 n=norm(); return cn(x/n, y/n, z/n);}
  inline real8 angle(cn A) const        // angle=A.angle(B); //0,pi;
    {return acos(in(A)/(norm()*A.norm()));}
  inline cn scale(real8 val) const // scaling of a point by a scalar KKM-05
    {return cn(x*val, y*val, z*val);}
  inline void showdata() const     // displays a point by KKM-05 Mar 9, 2005
    {DEBUG << "cn.x=" << x << "; cn.y=" << y << "; cn.z=" << z; DEBUG.print();}
  // --- Test program for coordinate class ----------------------------
  // --- This test is executed in inittest() --------------------------
  inline void test()
    {
    // constructors
    cn X(1,2,3);
    DEBUG << "cn X(1,2,3): "; X.showdata();
    cn Y;
    DEBUG << "cn Y:          "; Y.showdata();
    Y.x=4; Y.y=5; Y.z=6;
    DEBUG << "Y.x=4;Y.y=5;Y.z=6: "; Y.showdata();
    cn Z=Y;
    DEBUG << "Z=Y:           "; Z.showdata();
    // operators
    Z*=X;
    DEBUG << "Z*=X:          "; Z.showdata();
    Z/=X;
    DEBUG << "Z/=X:          "; Z.showdata();
    Z+=X;
    DEBUG << "Z+=X:          "; Z.showdata();
    Z-=X;
    DEBUG << "Z-=X:          "; Z.showdata();
    // operators
    DEBUG << "X+Y:           ";  (X+Y).showdata();
    DEBUG << "X-Y:           ";  (X-Y).showdata();
    DEBUG << "X/Y:           ";  (X/Y).showdata();
    DEBUG << "X*Y:           ";  (X*Y).showdata();
    DEBUG << "+X:            ";   (+X).showdata();
    DEBUG << "-X:            ";   (-X).showdata();
    DEBUG << "X==Y:          " << (X==Y); DEBUG.print();
    DEBUG << "X!=Y:          " << (X!=Y); DEBUG.print();
    // functions
    DEBUG << "X.min(Y):      "; (X.min(Y)).showdata();
    DEBUG << "X.normalize(): "; (X.normalize()).showdata();
    DEBUG << "X.scale(5):    "; (X.scale(5)).showdata();
    DEBUG << "X.out(Y):      "; (X.out(Y)).showdata();
    DEBUG << "X.dist(Y):     " << X.dist(Y); DEBUG.print();
    DEBUG << "X.norm():      " << X.norm(); DEBUG.print();
    DEBUG << "X.norm2():     " << X.norm2(); DEBUG.print();
    DEBUG << "X.angle(Y):    " << X.angle(Y); DEBUG.print();
    DEBUG << "X.in(Y):       " << X.in(Y); DEBUG.print();
    }
  };




// ====== Class for ellips =====================================
class input_ell                 // ellips a,b   
  {
  private:
  real8         e2;                     // squared first  eccentricity (derived)
  real8         e2b;                    // squared second eccentricity (derived)
  // ______ Helpers should be private ______
  inline void set_ecc1st_sqr()          // first ecc.
    {e2=1.0-sqr(b/a);}//  faster than e2=(sqr(a)-sqr(b))/sqr(a);
  inline void set_ecc2nd_sqr()          // second ecc.
    {e2b=sqr(a/b)-1.0;}// faster than e2b=(sqr(a)-sqr(b))/sqr(b);

  public:
  real8         a;                      // semi major
  real8         b;                      // semi minor
  char          name[EIGHTY];
  inline void set_name(const char *s)    {strcpy(name,s);} 
  // ______ Default constructor ______
  input_ell()
    {a   = WGS84_A;
     b   = WGS84_B; 
     e2  = 0.006694379990141;
     e2b = 0.006739496742276;
     set_name("WGS84");
    }
  // ______ constructor ellips(a,b) ______
  input_ell(const real8 &semimajor, const real8 &semiminor)
    {a = semimajor;
     b = semiminor; 
     set_ecc1st_sqr();// compute e2 (not required for zero-doppler iter.)
     set_ecc2nd_sqr();// compute e2b (not required for zero-doppler iter.)
     //set_name("unknown");
     set_name("Default");
    }
  // ______ Copy constructor ______
  input_ell(const input_ell& ell)
    {a=ell.a; b=ell.b; e2=ell.e2; e2b=ell.e2b; strcpy(name,ell.name);}
  // ______ Destructor ______
  ~input_ell()
    {;}// nothing to destruct that isn't destructed automatically
  // ______ Public function in struct ______
  inline input_ell& operator = (const input_ell &ell)// assignment operator
    {
    if (this != &ell) 
      {a=ell.a; b=ell.b; e2=ell.e2; e2b=ell.e2b; strcpy(name,ell.name);}
    return *this;
    }
  inline void showdata() const
    {
    INFO << "ELLIPSOID: \tEllipsoid used (orbit, output): " << name << ".";
    INFO.print();
    INFO << "ELLIPSOID: a   = " << setw(15) << setprecision(13) << a;
    INFO.print();
    INFO << "ELLIPSOID: b   = " << setw(15) << setprecision(13) << b;
    INFO.print();
    INFO << "ELLIPSOID: e2  = " << e2;
    INFO.print();
    INFO << "ELLIPSOID: e2' = " << e2b;
    INFO.print();
    INFO.reset();
    }

  // ______ Convert xyz cartesian coordinates to ______
  // ______ Geodetic ellipsoid coordinates latlonh ______
  /****************************************************************
   *    xyz2ell                                                   *
   *                                                              *
   * Converts geocentric cartesian coordinates in the XXXX        *
   *  reference frame to geodetic coordinates.                    *
   *  method of bowring see globale en locale geodetische systemen*
   * input:                                                       *
   *  - ellipsinfo, xyz, (phi,lam,hei)                            *
   * output:                                                      *
   *  - void (updated lam<-pi,pi>, phi<-pi,pi>, hei)              *
   *                                                              *
   *    Bert Kampes, 05-Jan-1999                                  *
   ****************************************************************/
  inline void xyz2ell(const cn &xyz, real8 &phi, real8 &lambda, real8 &height) const
    {
    TRACE_FUNCTION("xyz2ell (BK 05-Jan-1999)");
    const real8 r    = sqrt(sqr(xyz.x)+sqr(xyz.y));
    const real8 nu   = atan2((xyz.z*a),(r*b));
    const real8 sin3 = pow(sin(nu),3);
    const real8 cos3 = pow(cos(nu),3);
    phi              = atan2((xyz.z+e2b*b*sin3),(r-e2*a*cos3));
    lambda           = atan2(xyz.y,xyz.x);
    const real8 N    = a / sqrt(1.0-e2*sqr(sin(phi)));
    height           = (r/cos(phi)) - N;
    } // END xyz2ell
  // --- Same but without height ---
  inline void xyz2ell(const cn &xyz, real8 &phi, real8 &lambda) const
    {
    TRACE_FUNCTION("xyz2ell (BK 05-Jan-1999)");
    const real8 r    = sqrt(sqr(xyz.x)+sqr(xyz.y));
    const real8 nu   = atan2((xyz.z*a),(r*b));
    const real8 sin3 = pow(sin(nu),3);
    const real8 cos3 = pow(cos(nu),3);
    phi              = atan2((xyz.z+e2b*b*sin3),(r-e2*a*cos3));
    lambda           = atan2(xyz.y,xyz.x);
    } // END xyz2ell


  /****************************************************************
   *    ell2xyz                                                   *
   *                                                              *
   * Converts wgs84 ellipsoid cn to geocentric cartesian coord.   *
   * input:                                                       *
   *  - phi,lam,hei (geodetic co-latitude, longitude, [rad] h [m] *
   * output:                                                      *
   *  - cn XYZ                                                    *
   *                                                              *
   *    Bert Kampes, 05-Jan-1999                                  *
   ****************************************************************/
  cn ell2xyz(const real8 &phi, const real8 &lambda, const real8 &height) const
    {
    const real8 N      = a / sqrt(1.0-e2*sqr(sin(phi)));
    const real8 Nph    = N + height;
    return cn(Nph * cos(phi) * cos(lambda),
              Nph * cos(phi) * sin(lambda),
             (Nph - e2*N)    * sin(phi));
    } // END ell2xyz


  // --- Test program for ellips class ----------------------------
  // --- This test is executed in inittest() --------------------------
  inline void test()
    {
    // constructors
    input_ell X;
    DEBUG << "ELLIPS: "  << X.name 
          << ": X.a="  << X.a  << "; X.b="   << X.b 
          << "; X.e2=" << X.e2 << "; X.e2b=" << X.e2b;
    DEBUG.print();
    }
  };







 /***************************************************************
 *                                                              *
 *        Class window + functions                              *
 *        image starts at 1, matrix starts at 0                 *
 *                                                              *
 *    Bert Kampes, 19-Dec-1999                                  *
 *    [MA] Window extend error check - 2008                     *
 ****************************************************************/
class window                            // window file={1,N,1,N}
  {
  public:
  uint          linelo,                         // min. line coord.
                linehi,                         // max. line coord.
                pixlo,                          // min. pix coord.
                pixhi;                          // max. pix coord.
  // ______ Constructors ______
  window()
    {linelo=0; linehi=0; pixlo=0; pixhi=0;}
  window(uint ll, uint lh, uint pl, uint ph)
    {
     TRACE_FUNCTION("window() (BK 19-Dec-1998)");

     linelo=ll; linehi=lh; pixlo=pl; pixhi=ph;

     // ______ Check window extend ______          [MA]
     if(  int32(linelo) > int32(linehi) )
       {
        ERROR << "Window: Impossible to continue... l0 coordinate [" << linelo << "] >  linehi [" << linehi << "].";
        ERROR.print();
        throw(usage_error) ;
       }
     else if( int32(pixlo) > int32(pixhi) )
       {
        ERROR << "Window: Impossible to continue... p0 coordinate [" << pixlo << "] > pixelhi [" << pixhi << "].";
        ERROR.print();
        throw(usage_error) ;
       }
     }

  // ______ Copy constructor ______
  window(const window& w)
    {linelo=w.linelo; linehi=w.linehi; pixlo=w.pixlo; pixhi=w.pixhi;}
  // ______ Destructor ______
  ~window()
    {;};// nothing to destruct that isn't destructed automatically
  // ______ Public function in struct ______
  inline void disp() const                      // show content
    {
    DEBUG << "window [l0:lN, p0:pN] = ["
         << linelo << ":" << linehi << ", " << pixlo  << ":" << pixhi << "]";
    DEBUG.print();
    }
  inline uint lines() const                     // return number of lines
    {return linehi-int32(linelo)+1;}
  inline uint pixels() const                    // return number of pixels
    {return pixhi-int32(pixlo)+1;}
  inline window& operator = (const window &X)// assignment operator
    {if (this != &X)
      {linelo=X.linelo;linehi=X.linehi;pixlo=X.pixlo;pixhi=X.pixhi;}
     return *this;}
  inline bool operator == (const window &X) const
    {return (linelo==X.linelo&&linehi==X.linehi &&
              pixlo==X.pixlo && pixhi==X.pixhi) ?   true : false;}
  inline bool operator != (const window &X) const
    {return (linelo==X.linelo&&linehi==X.linehi &&
              pixlo==X.pixlo && pixhi==X.pixhi) ?   false : true;}
  // --- Test program for coordinate class ----------------------------
  // --- This test is executed in inittest() --------------------------
  inline void test()
    {
    // constructors
    window X;
    DEBUG << "window X:      "; X.disp();
    window Y(11,21,103,114);
    DEBUG << "window Y(11,21,103,114): "; Y.disp();
    // functions
    X=Y;
    DEBUG << "X=Y:           "; X.disp();
    DEBUG << "X.lines():     " <<  X.lines(); DEBUG.print();
    DEBUG << "X.pixels():    " <<  X.pixels(); DEBUG.print();
    DEBUG << "X==Y:          " <<  (X==Y); DEBUG.print();
    DEBUG << "X!=Y:          " <<  (X!=Y); DEBUG.print();
    }
  };


// ====== Class slavewindow + functions ================================
// added by FvL
class slavewindow // window file={l00,p00,l0N,p0N,lN0,pN0,lNN,pNN}
  {
  public:
    real8 l00, p00, l0N, p0N, lN0, pN0, lNN, pNN;
  // ______ Constructors ______
  slavewindow()
    {l00=0; p00=0; l0N=0; p0N=0; lN0=0; pN0=0; lNN=0; pNN=0;}
  slavewindow(real8 ll00, real8 pp00, real8 ll0N, real8 pp0N, real8 llN0, real8 ppN0, real8 llNN, real8 ppNN)
    {l00=ll00; p00=pp00; l0N=ll0N; p0N=pp0N; lN0=llN0; pN0=ppN0; lNN=llNN; pNN=ppNN;}
  // ______ Copy constructor ______
  slavewindow(const slavewindow& w)
    {l00=w.l00; p00=w.p00; l0N=w.l0N; p0N=w.p0N; lN0=w.lN0; pN0=w.pN0; lNN=w.lNN; pNN=w.pNN;}
  // ______ Destructor ______
  ~slavewindow()
    {;};// nothing to destruct that isn't destructed automatically
};

#endif // CONSTANTS_H


