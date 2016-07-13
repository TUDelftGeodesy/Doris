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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/utilities.cc,v $  *
 * $Revision: 3.17 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of utility routines.                          *
 * - math utils                                                 * 
 * - solving small linear system                                *
 * - polynomial evaluation                                      *
 * - getorb call                                                *
 * - tiepoint computation                                       *
 ****************************************************************/


#include "matrixbk.hh"                  // my matrix class
#include "slcimage.hh"                  // my slc image class
#include "constants.hh"                 // global constants
#include "utilities.hh"                 // header file
#include "ioroutines.hh"                // error messages
#include "exceptions.hh"                 // my exceptions class

#include <strstream>                    // for memory stream
#include <iomanip>                      // setw
#include <cstring>                      // for strcat etc.
#include <cstdlib>                      // system
#include <cmath>                        // sqrt etc.
#include <cstdio>                       // some compilers, remove function
#include <ctime>                        // time functions
char *strptime(const char *s, const char  *format,  struct tm *tm);        



#ifndef NO_FASTTRIG
// ___
// ___ globals for lookuptable fast_sin() and fast_cos()
// ___ time for 10000000 std::sin(x):  2.73 sec.
// ___ time for LUT of lengths 128:65536 = 0.12sec
// ___ max_error is 1 bin (not 0.5 since we do not check for sign),
// ___ i.e., max slope =1--> maxerror=2pi/65536=0.0001
const uint32 LUT_LENGHT    = 65536;                                    // memory req: 65536*sizeof(Type=real8)= 524288 Bytes. check fast_sin definition before changing this value.
const uint16 LUT_ATAN_MAX  = 128;
const real8 LUT_DX_SIN     = real8(LUT_LENGHT/(2.0*PI));               // idx=int(x*LUT_STEP)
const real8 LUT_DX_ATAN    = real8(LUT_LENGHT)/real8(2*LUT_ATAN_MAX);  // idx=int(x*DX+0.5)
real8 LUT_SIN[LUT_LENGHT];                                             // lookup table for sin [0:2PI)
real8 LUT_ATAN[LUT_LENGHT];                                            // lookup table for atan [-ATAN_MAX:dx:ATAN_MAX-dx]
#endif




/****************************************************************
 *    getorb                                                    *
 *                                                              *
 * Calls getorb program to generate precise orbit file          *
 *  scratchorbit.                                               *
 * Datapoints with INTERVAL seconds in between.                 *
 * time1-BEFOREFIRST to timeN+BEFOREFIRST                       *
 * time1 in format as: 15-AUG-1996 09:49:16.913                 *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void getorb(
        const slcimage &image,
        const input_pr_orbits &inputorb,
        int16 ID)
  {
  TRACE_FUNCTION("getorb (BK 11-Dec-1998)")
  const int32 INTERVAL    = inputorb.timeinterval;
  const int32 BEFOREFIRST = inputorb.timebefore;
  char  orbdir[EIGHTY];
  char  dummyline[2*ONE27];
  char  startt[13];
  char  endt[13];
  char  *pstartt = startt;
  char  *pendt = endt;
  if (ID==MASTERID)
    strcpy(orbdir,inputorb.m_orbdir);
  else if (ID==SLAVEID)
    strcpy(orbdir,inputorb.s_orbdir);
  else
    {
    PRINT_ERROR("panic, impossible.");
    throw(unhandled_case_error);
    }

  // ______ Convert times to struct tm ______
  struct tm tijdstart;
  strptime(image.utc1,"%d-%b-%Y %T",&tijdstart);
  DEBUG << "value of image.utc1: " << image.utc1;
  DEBUG.print();

  // ______ Adjust interval ______
  int32 t1seconds = tijdstart.tm_sec + 60*tijdstart.tm_min + 3600*tijdstart.tm_hour;
  int32 t2seconds = t1seconds + BEFOREFIRST +
    int32(0.5+(image.originalwindow.linehi-image.originalwindow.linelo+1)/image.prf);
  t1seconds -= BEFOREFIRST;                             // adjust;
  if (t1seconds < 0)                                    // check 0:00'00" utc
    {
    WARNING.print("getorb: couldn't adjust time lower bound, using old value.");
    t1seconds += BEFOREFIRST;
    }
  tijdstart.tm_hour =  t1seconds / 3600;
  tijdstart.tm_min  = (t1seconds-3600*tijdstart.tm_hour) / 60;
  tijdstart.tm_sec  = (t1seconds-3600*tijdstart.tm_hour) % 60;
  DEBUG << "start time: tm struct: hour=" << tijdstart.tm_hour 
        << "; min=" << tijdstart.tm_min << "; sec=" << tijdstart.tm_sec;
  DEBUG.print();

  if (t2seconds / 3600 > 23)                            // hours + 1 day
    {
    WARNING.print("getorb: couldn't adjust time upper bound, using old value.");
    t2seconds -= BEFOREFIRST;
    }
  struct tm tijdend = tijdstart;                // copy year, month, day
  tijdend.tm_hour =  t2seconds / 3600;
  tijdend.tm_min  = (t2seconds-3600*tijdend.tm_hour) / 60;
  tijdend.tm_sec  = (t2seconds-3600*tijdend.tm_hour) % 60;
  DEBUG << "end time: tm struct: hour=" << tijdend.tm_hour 
        << "; min=" << tijdend.tm_min << "; sec=" << tijdend.tm_sec;
  DEBUG.print();

  // ______ Check timedifference ______
  //if (tijdend.tm_min - tijdstart.tm_min != 0 &&
  //    tijdend.tm_min - tijdstart.tm_min != 1 )
  if (tijdend.tm_min-tijdstart.tm_min >  10 ||
      tijdend.tm_min-tijdstart.tm_min < -50)
    WARNING.print("getorb: > 10 minutes of orbit ephemerides?");

  // ______convert time to strings[12]  YYMMDDHHMMSS(.S)______
  strftime(pstartt,13,"%y%m%d%H%M%S",&tijdstart);       // including \0
  strftime(pendt,13,"%y%m%d%H%M%S",&tijdend);           // including \0


  // ====== Get precise orbits ======
  // ______ Call getorb (unix) ______
  // ______Make string starttime YYMMDDHHMMSS(.S)______
  char strgetorb[2*ONE27];                                // unix system call
  ostrstream omemgetorb(strgetorb,2*ONE27);               // combine types in string
  omemgetorb.seekp(0);
  omemgetorb << "getorb " << startt << "," << endt << "," << INTERVAL << " "
             << orbdir << ">scratchorbit" << ends;
  // BK 4 july 2000: assume ODR files are already untarred
  //omemgetorb << "getorb " << startt << "," << endt << "," << INTERVAL << " ";
  //if (removeuntarred == false)        // .ODR file is in archive directory
  //  omemgetorb << orbdir << ">scratchorbit" << ends;
  //else                                // file is untarred to current dir
  //  omemgetorb << ". > scratchorbit" << ends;

 // 1.  The # of bits *captured* from a return code varies.
 //   Here are some examples:

 // some notes on status code:
 //       shell scripts           - capture 8-bits
 //       tcl/expect scripts      - capture 8-bits
 //       perl scripts            - capture 16-bits
 //       c binaries              - capture 16-bits
 //  see  http://www.etpenguin.com/cgi-bin/index2file.cgi?/pub/Reference/DOC_exit-codes.txt





  int32 status = 3*NaN;  // [MA] keep int32 due to 65280 value.
  for (register int32 i=0;i<3;i++)                      // try a few times.
    {
    DEBUG.print(strgetorb);
    status=(system(strgetorb));                         // run getorb  / return positive values -1 -->  65280
    DEBUG << "status of getorb: " << status;
    DEBUG.print();
    if (status == 0)
      {
      // ______Check really ok______
      DEBUG.print("Checking getorb result for errors.");
      //ifstream scratchorbit("scratchorbit",ios::in | ios::nocreate); // by getorb
      ifstream scratchorbit("scratchorbit",ios::in); // by getorb
      bk_assert(scratchorbit,"getorb: scratchorbit",__FILE__,__LINE__);

      int32 err=0;
      real8 dsec;
      bool getorbok=true;
      while(scratchorbit)
        {
        scratchorbit >> dsec >> err;
        DEBUG << "time: " << dsec << " err: " << err;
        DEBUG.print();
        if (err)
          {
          getorbok=false;
          WARNING.print("seems a problem with getorb, can you repair?");
          getanswer();
          break; // file?
          }
        scratchorbit.getline(dummyline,2*ONE27,'\n');
        }
      scratchorbit.close();
      if (getorbok) 
        {
        PROGRESS.print("getorb: program finished ok.");
        break;  // break try few times
        }
      }
    else if (status == -1 || status==-256 || status == 65280 )  // [MA] [-1] --> 65280
      {
      WARNING.print(strgetorb);
      WARNING.print("getorb: exit(-1): time within ODR, but outside precise part.");
      break;
      }
    else if (status ==  1 || status==256)
      {
      WARNING.print(strgetorb);
      WARNING.print("getorb: exit(1): error opening ODR file.\n\tPerhaps file is tarred?");
      getanswer();
      }
    else if (status ==  2 || status==512)
      {
      WARNING.print(strgetorb);
      PRINT_ERROR("code 905: getorb: exit(2): no ODR for requested epoch.")
      throw(some_error);
      }
    else if (status ==  3 || status==768)
      {
      WARNING.print(strgetorb);
      PRINT_ERROR("code 905: getorb: exit(3): too many records.")
      throw(some_error);
      }
    else 
      {
      WARNING.print(strgetorb);
      ERROR << "code 905: utilities.c: getorb: unknown exit code: " << status;
      PRINT_ERROR(ERROR.get_str())
      throw(some_error);
      }
    }

  // ______Tidy up______
  } // END getorb



/****************************************************************
 *    convertgetorbout                                          *
 *                                                              *
 * converts output of getorb (in "scratchorbit") to format:     *
 *  secsofday X Y Z                                             *
 * input:                                                       *
 *  - scratchorbit                                              *
 * output:                                                      *
 *  - scratchdatapoints, this file                              *
 *    might already exist due to readleader                     *
 *  - file is removed                                           *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void convertgetorbout(
        int16 FILEID,
        const char* file)                               // default=scratchorbit
  {
  TRACE_FUNCTION("convertgetorbout (BK 11-Dec-1998)")
  if (FILEID != MASTERID && FILEID!=SLAVEID) 
    {
    PRINT_ERROR("fileid not ok.")
    throw(some_error);
    }
  char                  dummyline[2*ONE27];

  // ______ Open files ______
  ifstream infile(file, ios::in);
  bk_assert(infile,file,__FILE__,__LINE__);

  ofstream outfile("scratchdatapoints", ios::out | ios::trunc); // do replace
  bk_assert(outfile,"convertgetorbout: scratchdatapoints",__FILE__,__LINE__);

  // ______Read line: sec85 err fie lambda h x y z______
  int32 err;
  real8 sec;
  real8 fie;
  real8 lam;
  real8 hei;
  char x[25];                                           // to be sure no round off
  char y[25];                                           // to be sure no round off
  char z[25];                                           // to be sure no round off

  // ______ Get number of points ______
  int32 numdatapoints = -1;
  int32 numerr        =  0;
  DEBUG.print("Counting number of lines in tmp file with orbit; checking err parameter");
  while(infile)
    {
    infile >> sec >> err >> fie >> lam >> hei >> x >> y >> z;
    infile.getline(dummyline,2*ONE27,'\n');
    DEBUG  << "scratchfile: sec=" << sec << " err=" << err 
          << " x=" << x << " y=" << y << " z=" << z;
    DEBUG.print();
    if (err!=0) numerr++;
    if (err==0) numdatapoints++;
    }
  infile.clear();// g++v3.2 has a problem if you don't to read it again
  //infile.close();// seems problem with g++v3.2 if stream is closed and opened?
  INFO << "Number of orbit ephemerides found: " << numdatapoints+1 << ends;
  INFO.print();
  INFO << "Number of erroneous orbit ephemerides found: " << numerr << ends;
  if (numerr==0) INFO.print();
  if (numerr>0)  WARNING.print(INFO.get_str());
  if ( numdatapoints == -1 )  // [MA] -1 --> no good pts]
    { 
     //WARNING.print(INFO.get_str()); 
     infile.close();
     if (remove("scratchdatapoints"))                                     // remove scratchdatapoints but keep scratchorbit
       WARNING.print("code 101: could not remove file.");
     ERROR << "No valid number of orbit ephemeris found! exiting..."; 
     PRINT_ERROR(ERROR.get_str())
     throw(some_error);
    }
  INFO.reset();// rewind

  // ___ Write results in parameter file via 2nd tmpfile ___
  outfile << "\n\n*******************************************************************"
          << "\n*_Start_";
          //<< "\n*_Start_precise_datapoints"
  if (FILEID==MASTERID)
    outfile << processcontrol[pr_m_porbits];
  if (FILEID==SLAVEID)
    outfile << processcontrol[pr_s_porbits];
  outfile << "\n*******************************************************************"
          << "\n\tt(s)\tX(m)\tY(m)\tZ(m)"
          << "\nNUMBER_OF_DATAPOINTS: \t\t\t"
          <<  numdatapoints
          << endl;

  // ______Set outfile in correct digits______
  int32 secint;                                         // integer part
  real8 secfrac;                                        // fractional part
  outfile.setf(ios::fixed);
  outfile.setf(ios::showpoint);
  //infile.open(file);
  infile.seekg(0, ios::beg);// rewind file to read again
  DEBUG.print("Rewinding file and reading orbit ephemerides again");
  for (register int32 i=0;i<numdatapoints+numerr;i++)// total number of lines
    {
    infile >> sec >> err >> fie >> lam >> hei >> x >> y >> z;
    infile.getline(dummyline,2*ONE27,'\n');               // go to next line
    DEBUG  << "scratchfile2: sec=" << sec << " err=" << err 
          << " x=" << x << " y=" << y << " z=" << z;
    DEBUG.print();

    // ______Convert time sec85 to sec of day______
    secint  = int32(sec)%86400;                         // 60*60*24 integer part
    secfrac = sec - int32(sec);                         // fractional part 
    if (err!=0) 
      {
      WARNING << "err == "<< err << " " << "err!=0: something went wrong in getorb? (skipping point)";
      WARNING.print();
      }
    if (err==0)
      outfile << secint+secfrac << "\t" << x << "\t" << y << "\t" << z << endl;
    }

  outfile << "\n*******************************************************************"
          //<< "\n* End_precise_datapoints:_NORMAL"
          << "\n* End_";
  if (FILEID==MASTERID)
    outfile << processcontrol[pr_m_porbits];
  if (FILEID==SLAVEID)
    outfile << processcontrol[pr_s_porbits];
  outfile << "_NORMAL"
          << "\n*******************************************************************\n";

  // ______Tidy up______
  infile.close();
  outfile.close();
  if (remove(file))                                     // remove file
    WARNING.print("code 101: could not remove file.");
  } // END convertgetorbout



/****************************************************************
 *    solve33                                                   *
 *                                                              *
 * Solves setof 3 equations by straightforward (no pivotting) LU*
 * y=Ax (unknown x)                                             *
 *                                                              *
 * input:                                                       *
 *  - matrix<real8> righthandside 3x1 (y)                       *
 *  - matrix<real8> partials 3x3 (A)                            *
 *                                                              *
 * output:                                                      *
 *  - matrix<real8> result 3x1 unknown                          *
 *                                                              *
 *    Bert Kampes, 22-Dec-1998                                  *
 ****************************************************************/
void solve33(
        matrix<real8>       &RESULT,
        const matrix<real8> &rhs,
        const matrix<real8> &A)
  {
  TRACE_FUNCTION("solve33 (BK 22-Dec-1998)")
#ifdef __DEBUG
  if (RESULT.lines() != 3 || RESULT.pixels() != 1)
    {
    PRINT_ERROR("solve33: input: size RES not 3x1.")
    throw(input_error);
    }
  if (A.lines() != 3 || A.pixels() != 3)
    {
    PRINT_ERROR("solve33: input: size of A not 33.")
    throw(input_error);
    }
  if (rhs.lines() != 3 || rhs.pixels() != 1)
    {
    PRINT_ERROR("solve33: input: size rhs not 3x1.")
    throw(input_error);
    }
#endif
   
  // ______  real8 L10, L20, L21: used lower matrix elements
  // ______  real8 U11, U12, U22: used upper matrix elements
  // ______  real8 b0,  b1,  b2:  used Ux=b
  const real8 L10 =  A(1,0)/A(0,0);
  const real8 L20 =  A(2,0)/A(0,0);
  const real8 U11 =  A(1,1)-L10*A(0,1);
  const real8 L21 = (A(2,1)-(A(0,1)*L20))/U11;
  const real8 U12 =  A(1,2)-L10*A(0,2);
  const real8 U22 =  A(2,2)-L20*A(0,2)-L21*U12;
   
  // ______ Solution: forward substitution ______
  const real8 b0  =  rhs(0,0);
  const real8 b1  =  rhs(1,0)-b0*L10;
  const real8 b2  =  rhs(2,0)-b0*L20-b1*L21;
   
  // ______ Solution: backwards substitution ______
  RESULT(2,0)     =  b2/U22;
  RESULT(1,0)     = (b1-U12*RESULT(2,0))/U11;
  RESULT(0,0)     = (b0-A(0,1)*RESULT(1,0)-A(0,2)*RESULT(2,0))/A(0,0);
  } // END solve33



/****************************************************************
 *    solve22                                                   *
 *                                                              *
 * Solves setof 2 equations by straightforward substitution     *
 * y=Ax (unknown x)                                             *
 *                                                              *
 * input:                                                       *
 *  - matrix<real8> righthandside 2x1 (y)                       *
 *  - matrix<real8> partials 2x2 (A)                            *
 * output:                                                      *
 *  - matrix<real8> result 2x1 unknown dx,dy,dz                 *
 *                                                              *
 * Future: LaTeX software doc.                                  *
 * Future: this has been tested                                 *
 *                                                              *
 *    Bert Kampes, 22-Dec-1998                                  *
 ****************************************************************/
matrix<real8> solve22(
        const matrix<real8> &y,
        const matrix<real8> &A)
  {
  TRACE_FUNCTION("solve22 (BK 22-Dec-1998)")
#ifdef __DEBUG
  // ______ Check input ______
  if (A.lines() != 2 || A.pixels() != 2)
    {
    PRINT_ERROR("error...")
    throw(input_error);
    }
  if (y.lines() != 2 || y.pixels() != 1)
    {
    PRINT_ERROR("error...")
    throw(input_error);
    }
#endif

  matrix<real8> Result(2,1);
  Result(1,0)= (y(0,0) - ((A(0,0)/A(1,0))*y(1,0))) / (A(0,1) - ((A(0,0)*A(1,1))/A(1,0)));
  Result(0,0)= (y(0,0) - A(0,1)*Result(1,0)) / A(0,0);
  return Result;
  } // END solve22



/****************************************************************
 *    nextpow2                                                  *
 *                                                              *
 * Returns 32 for 31.3 etc.                                     *
 *                                                              *
 *    Bert Kampes, 03-Feb-1999                                  *
 * AFter bug report from Bruno who discover bug.  b-1 instead   *
 * of b.  Did not                                               *
 * have effect on Doris since this function was effectively not *
 * used.                                                        *
 #%// BK 16-Apr-2002                                            *
 ****************************************************************/
uint nextpow2(
        real8 w)        
  {
  TRACE_FUNCTION("nextpow2 (BK 03-Feb-1999)")
  int32 b    = 0; 
  int32 *pnt = &b;
  real8 f    = frexp(w,pnt);
  if (f==0.5) 
    return uint(w);
  return (uint(pow(real8(2),b)));
  } // END nextpow2



/****************************************************************
 *    polyval                                                   *
 *                                                              *
 * evaluates polynomial at (x,y)                                *
 * input:                                                       *
 *  - x,y: point                                                *
 *  - polynomial coefficients                                   *
 * (a00 | 01 10 | 20 11 02 | 30 21 12 03 | ...)                 *
 *                                                              *
 * output:                                                      *
 *  - real8 functional value at (x,y).                          *
 *                                                              *
 *    Bert Kampes, 15-Mar-1999                                  *
 ****************************************************************/
real8 polyval(
        real8 x,
        real8 y,
        const matrix<real8> &coeff)
  {
  const int32 degreee = degree(coeff.size());
#ifdef __DEBUG
  TRACE_FUNCTION("polyval (BK 15-Mar-1999)")
  if (coeff.size() != coeff.lines())
    {
    PRINT_ERROR("polyval: require standing data vector.")
    throw(input_error);
    }
  if (degreee<0 || degreee > 1000)
    {
    PRINT_ERROR("polyval: degreee < 0")
    throw(input_error);
    }
#endif

// ______Check default arguments______
//  if (degreee==-1)                    // default applies
//    const int32 degreee = degree(coeff.size());

// ______Speed up______
// ______Evaluate polynomial______

  real8 sum = coeff(0,0);

  if (degreee == 1)
    {
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y );
    }

  else if (degreee == 2)
    {
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * sqr(x)
            + coeff(4,0) * x*y
            + coeff(5,0) * sqr(y) );
    }

  else if (degreee == 3)
    {
    register real8 xx = sqr(x);
    register real8 yy = sqr(y);
    sum += (  coeff(1,0) * x
          + coeff(2,0) * y
          + coeff(3,0) * xx
          + coeff(4,0) * x*y
          + coeff(5,0) * yy
          + coeff(6,0) * xx*x
          + coeff(7,0) * xx*y
          + coeff(8,0) * x*yy
          + coeff(9,0) * yy*y );
    }

  else if (degreee == 4)
    {
    register real8 xx  = sqr(x);
    register real8 xxx = xx*x;
    register real8 yy  = sqr(y);
    register real8 yyy = yy*y;
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * xx
            + coeff(4,0) * x*y
            + coeff(5,0) * yy
            + coeff(6,0) * xxx
            + coeff(7,0) * xx*y
            + coeff(8,0) * x*yy
            + coeff(9,0) * yyy
            + coeff(10,0)* xx*xx
            + coeff(11,0)* xxx*y
            + coeff(12,0)* xx*yy
            + coeff(13,0)* x*yyy
            + coeff(14,0)* yy*yy );
    }

  else if (degreee == 5)
    {
    register real8 xx   = sqr(x);
    register real8 xxx  = xx*x;
    register real8 xxxx = xxx*x;
    register real8 yy   = sqr(y);
    register real8 yyy  = yy*y;
    register real8 yyyy = yyy*y;
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * xx
            + coeff(4,0) * x*y
            + coeff(5,0) * yy
            + coeff(6,0) * xxx
            + coeff(7,0) * xx*y
            + coeff(8,0) * x*yy
            + coeff(9,0) * yyy
            + coeff(10,0)* xxxx
            + coeff(11,0)* xxx*y
            + coeff(12,0)* xx*yy
            + coeff(13,0)* x*yyy
            + coeff(14,0)* yyyy
            + coeff(15,0)* xxxx*x
            + coeff(16,0)* xxxx*y
            + coeff(17,0)* xxx*yy
            + coeff(18,0)* xx*yyy
            + coeff(19,0)* x*yyyy
            + coeff(20,0)* yyyy*y );
    }

  else if (degreee != 0)                // degreee > 5
    {
    sum = 0.0;
    register int32 coeffindex = 0;
    for (register int32 l=0; l<=degreee; l++)
      {
      for (register int32 k=0; k<=l ;k++)
        {
        sum += coeff(coeffindex,0) * pow(x,real8(l-k)) * pow(y,real8(k));
        coeffindex++;
        }
      }
    }

  return sum;  
  } // END polyval



/****************************************************************
 *    polyval                                                   *
 *                                                              *
 * evaluates polynomial at (x,y)                                *
 * input:                                                       *
 *  - x,y: point                                                *
 *  - polynomial coefficients                                   *
 * (a00 | 01 10 | 20 11 02 | 30 21 12 03 | ...)                 *
 *                                                              *
 * output:                                                      *
 *  - real8 functional value at (x,y).                          *
 *                                                              *
 *    Bert Kampes, 15-Mar-1999                                  *
 ****************************************************************/
real8 polyval(
        real8 x,
        real8 y,
        const matrix<real8> &coeff,
        int32 degreee)
  {
#ifdef __DEBUG
  TRACE_FUNCTION("polyval (BK 15-Mar-1999)")
  if (coeff.size() != coeff.lines())
    {
    PRINT_ERROR("polyval: require standing data vector.")
    throw(input_error);
    }
  if (degreee<0 || degreee > 1000)
    {
    PRINT_ERROR("polyval: degreee < 0")
    throw(input_error);
    }
#endif

// ______Check default arguments______
//  if (degreee==-1)                    // default applies
//    degreee = degree(coeff.size());

// ______Speed up______
// ______Evaluate polynomial______

  real8 sum = coeff(0,0);

  if (degreee == 1)
    {
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y );
    }

  else if (degreee == 2)
    {
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * sqr(x)
            + coeff(4,0) * x*y
            + coeff(5,0) * sqr(y) );
    }

  else if (degreee == 3)
    {
    register real8 xx = sqr(x);
    register real8 yy = sqr(y);
    sum += (  coeff(1,0) * x
          + coeff(2,0) * y
          + coeff(3,0) * xx
          + coeff(4,0) * x*y
          + coeff(5,0) * yy
          + coeff(6,0) * xx*x
          + coeff(7,0) * xx*y
          + coeff(8,0) * x*yy
          + coeff(9,0) * yy*y );
    }

  else if (degreee == 4)
    {
    register real8 xx  = sqr(x);
    register real8 xxx = xx*x;
    register real8 yy  = sqr(y);
    register real8 yyy = yy*y;
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * xx
            + coeff(4,0) * x*y
            + coeff(5,0) * yy
            + coeff(6,0) * xxx
            + coeff(7,0) * xx*y
            + coeff(8,0) * x*yy
            + coeff(9,0) * yyy
            + coeff(10,0)* xx*xx
            + coeff(11,0)* xxx*y
            + coeff(12,0)* xx*yy
            + coeff(13,0)* x*yyy
            + coeff(14,0)* yy*yy );
    }

  else if (degreee == 5)
    {
    register real8 xx   = sqr(x);
    register real8 xxx  = xx*x;
    register real8 xxxx = xxx*x;
    register real8 yy   = sqr(y);
    register real8 yyy  = yy*y;
    register real8 yyyy = yyy*y;
    sum += (  coeff(1,0) * x
            + coeff(2,0) * y
            + coeff(3,0) * xx
            + coeff(4,0) * x*y
            + coeff(5,0) * yy
            + coeff(6,0) * xxx
            + coeff(7,0) * xx*y
            + coeff(8,0) * x*yy
            + coeff(9,0) * yyy
            + coeff(10,0)* xxxx
            + coeff(11,0)* xxx*y
            + coeff(12,0)* xx*yy
            + coeff(13,0)* x*yyy
            + coeff(14,0)* yyyy
            + coeff(15,0)* xxxx*x
            + coeff(16,0)* xxxx*y
            + coeff(17,0)* xxx*yy
            + coeff(18,0)* xx*yyy
            + coeff(19,0)* x*yyyy
            + coeff(20,0)* yyyy*y );
    }

  else if (degreee != 0)                // degreee > 5
    {
    sum = 0.0;
    register int32 coeffindex = 0;
    for (register int32 l=0; l<=degreee; l++)
      {
      for (register int32 k=0; k<=l ;k++)
        {
        sum += coeff(coeffindex,0) * pow(x,real8(l-k)) * pow(y,real8(k));
        coeffindex++;
        }
      }
    }
  return sum;  
  } // END polyval


// MA explicit instantiation: w/o these two linesr, whole function needs to go to utilities.hh 
template matrix<real4> polyval<real4>( const matrix<real4>&, const matrix<real4>&, const matrix<real8>&, int32 );
template matrix<real8> polyval<real8>( const matrix<real4>&, const matrix<real4>&, const matrix<real8>&, int32 ); 

/****************************************************************
 *    polyval                                                   *
 *                                                              *
 * Evaluate 2d polynomial at regular grid.                              *
 * input:                                                       *
 *  - x,y: (standing) data axis                                 *
 *  - polynomial coefficients                                   *
 * (a00 | 01 10 | 20 11 02 | 30 21 12 03 | ...)                 *
 *                                                              *
 * output:                                                      *
 * ( - file ascii "TEST" with grid.)                            *
 *  - matrix<Type> with evaluated polynomial                   *
 *                                                              *
 * note optimixed for more y than x                             *
 *    Bert Kampes, 12-Mar-1999                                  *
 *    Mahmut Arikan, 04-06-2009 use template <>                 *
 ****************************************************************/
template <class Type> // Type real8 or real4
matrix<Type> polyval(
        const matrix<real4> &x,
        const matrix<real4> &y,
        const matrix<real8> &coeff,
        int32 degreee)                  // optional
  {
  TRACE_FUNCTION("polyval (BK 12-Mar-1999)")
#ifdef __DEBUG
  if (x.size() != x.lines())
    {
    PRINT_ERROR("polyval: require standing x vector.")
    throw(input_error);
    }
  if (y.size() != y.lines())
    {
    PRINT_ERROR("polyval: require standing y vector.")
    throw(input_error);
    }
  if (coeff.size() != coeff.lines())
    {
    PRINT_ERROR("polyval: require standing coeff. vector.")
    throw(input_error);
    }
  if (degreee<-1)
    {
    PRINT_ERROR("polyval: degreee < -1")
    throw(input_error);
    }
  if (x.size()>y.size())
    DEBUG.print("x larger than y, while optimized for y larger x");
#endif

// ______Check default arguments______
  if (degreee==-1)                              // default applies
    degreee = degree(coeff.size());

// ______Evaluate polynomial______
  matrix<Type> Result(x.size(),y.size());
  register int32 i,j;

  register Type c00;
  register Type c10;
  register Type c01;
  register Type c20;
  register Type c11;
  register Type c02;
  register Type c30;
  register Type c21;
  register Type c12;
  register Type c03;
  register Type c40;
  register Type c31;
  register Type c22;
  register Type c13;
  register Type c04;
  register Type c50;
  register Type c41;
  register Type c32;
  register Type c23;
  register Type c14;
  register Type c05;
  switch (degreee)
    {
    case 0:                                     // constant
      Result.setdata(Type(coeff(0,0)));
      break;

    case 1:
      c00 = coeff(0,0);
      c10 = coeff(1,0);
      c01 = coeff(2,0);
      for (j=0; j<Result.pixels(); j++)
        {
        Type c00pc01y1 = c00 + c01*y(j,0);
        for (i=0; i<Result.lines(); i++)
          {
          Result(i,j) =  c00pc01y1 + c10 * x(i,0);
          }
        }
      break;

    case 2:
        c00 = coeff(0,0);
        c10 = coeff(1,0);
        c01 = coeff(2,0);
        c20 = coeff(3,0);
        c11 = coeff(4,0);
        c02 = coeff(5,0);
        for (j=0; j<Result.pixels(); j++)
          {
          Type y1 = y(j,0);
          Type c00pc01y1 = c00 + c01 * y1;
          Type c02y2 = c02 * sqr(y1);
          Type c11y1 = c11 * y1;
          for (i=0; i<Result.lines(); i++)
            {
            Type x1 = x(i,0);
            Result(i,j) =  c00pc01y1 
                         + c10   * x1 
                         + c20   * sqr(x1) 
                         + c11y1 * x1 
                         + c02y2;
            }
          }
        break;

    case 3:
        c00 = coeff(0,0);
        c10 = coeff(1,0);
        c01 = coeff(2,0);
        c20 = coeff(3,0);
        c11 = coeff(4,0);
        c02 = coeff(5,0);
        c30 = coeff(6,0);
        c21 = coeff(7,0);
        c12 = coeff(8,0);
        c03 = coeff(9,0);
        for (j=0; j<Result.pixels(); j++)
          {
          Type y1 = y(j,0);
          Type y2 = sqr(y1);
          Type c00pc01y1 = c00 + c01 * y1;
          Type c02y2 = c02 * y2;
          Type c11y1 = c11 * y1;
          Type c21y1 = c21 * y1;
          Type c12y2 = c12 * y2;
          Type c03y3 = c03 * y1 * y2;
          for (i=0; i<Result.lines(); i++)
            {
            Type x1 = x(i,0);
            Type x2 = sqr(x1);
            Result(i,j) =  c00pc01y1 
                         + c10   * x1 
                         + c20   * x2 
                         + c11y1 * x1 
                         + c02y2
                         + c30   * x1 * x2
                         + c21y1 * x2
                         + c12y2 * x1
                         + c03y3;
            }
          }
        break;

    case 4:
        c00 = coeff(0,0);
        c10 = coeff(1,0);
        c01 = coeff(2,0);
        c20 = coeff(3,0);
        c11 = coeff(4,0);
        c02 = coeff(5,0);
        c30 = coeff(6,0);
        c21 = coeff(7,0);
        c12 = coeff(8,0);
        c03 = coeff(9,0);
        c40 = coeff(10,0);
        c31 = coeff(11,0);
        c22 = coeff(12,0);
        c13 = coeff(13,0);
        c04 = coeff(14,0);
        for (j=0; j<Result.pixels(); j++)
          {
          Type y1 = y(j,0);
          Type y2 = sqr(y1);
          Type c00pc01y1 = c00 + c01 * y1;
          Type c02y2 = c02 * y2;
          Type c11y1 = c11 * y1;
          Type c21y1 = c21 * y1;
          Type c12y2 = c12 * y2;
          Type c03y3 = c03 * y1 * y2;
          Type c31y1 = c31 * y1;
          Type c22y2 = c22 * y2;
          Type c13y3 = c13 * y2 * y1;
          Type c04y4 = c04 * y2 * y2;
          for (i=0; i<Result.lines(); i++)
            {
            Type x1 = x(i,0);
            Type x2 = sqr(x1);
            Result(i,j) =  c00pc01y1 
                         + c10   * x1 
                         + c20   * x2 
                         + c11y1 * x1 
                         + c02y2
                         + c30   * x1 * x2
                         + c21y1 * x2
                         + c12y2 * x1
                         + c03y3
                         + c40   * x2 * x2
                         + c31y1 * x2 * x1
                         + c22y2 * x2
                         + c13y3 * x1
                         + c04y4;
            }
          }
        break;

    case 5:
        c00 = coeff(0,0);
        c10 = coeff(1,0);
        c01 = coeff(2,0);
        c20 = coeff(3,0);
        c11 = coeff(4,0);
        c02 = coeff(5,0);
        c30 = coeff(6,0);
        c21 = coeff(7,0);
        c12 = coeff(8,0);
        c03 = coeff(9,0);
        c40 = coeff(10,0);
        c31 = coeff(11,0);
        c22 = coeff(12,0);
        c13 = coeff(13,0);
        c04 = coeff(14,0);
        c50 = coeff(15,0);
        c41 = coeff(16,0);
        c32 = coeff(17,0);
        c23 = coeff(18,0);
        c14 = coeff(19,0);
        c05 = coeff(20,0);
        for (j=0; j<Result.pixels(); j++)
          {
          Type y1 = y(j,0);
          Type y2 = sqr(y1);
          Type y3 = y2*y1;
          Type c00pc01y1 = c00 + c01 * y1;
          Type c02y2 = c02 * y2;
          Type c11y1 = c11 * y1;
          Type c21y1 = c21 * y1;
          Type c12y2 = c12 * y2;
          Type c03y3 = c03 * y3;
          Type c31y1 = c31 * y1;
          Type c22y2 = c22 * y2;
          Type c13y3 = c13 * y3;
          Type c04y4 = c04 * y2 * y2;
          Type c41y1 = c41 * y1;
          Type c32y2 = c32 * y2;
          Type c23y3 = c23 * y3;
          Type c14y4 = c14 * y2 * y2;
          Type c05y5 = c05 * y3 * y2;
          for (i=0; i<Result.lines(); i++)
            {
            Type x1 = x(i,0);
            Type x2 = sqr(x1);
            Type x3 = x1*x2;
            Result(i,j) =  c00pc01y1 
                         + c10   * x1 
                         + c20   * x2 
                         + c11y1 * x1 
                         + c02y2
                         + c30   * x3
                         + c21y1 * x2
                         + c12y2 * x1
                         + c03y3
                         + c40   * x2 * x2
                         + c31y1 * x3
                         + c22y2 * x2
                         + c13y3 * x1
                         + c04y4
                         + c50   * x3 * x2
                         + c41y1 * x2 * x2
                         + c32y2 * x3
                         + c23y3 * x2
                         + c14y4 * x1
                         + c05y5;
            }
          }
        break;

// ______ solve up to 5 efficiently, do rest in loop ______
      default:
        c00 = coeff(0,0);
        c10 = coeff(1,0);
        c01 = coeff(2,0);
        c20 = coeff(3,0);
        c11 = coeff(4,0);
        c02 = coeff(5,0);
        c30 = coeff(6,0);
        c21 = coeff(7,0);
        c12 = coeff(8,0);
        c03 = coeff(9,0);
        c40 = coeff(10,0);
        c31 = coeff(11,0);
        c22 = coeff(12,0);
        c13 = coeff(13,0);
        c04 = coeff(14,0);
        c50 = coeff(15,0);
        c41 = coeff(16,0);
        c32 = coeff(17,0);
        c23 = coeff(18,0);
        c14 = coeff(19,0);
        c05 = coeff(20,0);
        for (j=0; j<Result.pixels(); j++)
          {
          Type y1 = y(j,0);
          Type y2 = sqr(y1);
          Type y3 = y2*y1;
          Type c00pc01y1 = c00 + c01 * y1;
          Type c02y2 = c02 * y2;
          Type c11y1 = c11 * y1;
          Type c21y1 = c21 * y1;
          Type c12y2 = c12 * y2;
          Type c03y3 = c03 * y3;
          Type c31y1 = c31 * y1;
          Type c22y2 = c22 * y2;
          Type c13y3 = c13 * y3;
          Type c04y4 = c04 * y2 * y2;
          Type c41y1 = c41 * y1;
          Type c32y2 = c32 * y2;
          Type c23y3 = c23 * y3;
          Type c14y4 = c14 * y2 * y2;
          Type c05y5 = c05 * y3 * y2;
          for (i=0; i<Result.lines(); i++)
            {
            Type x1 = x(i,0);
            Type x2 = sqr(x1);
            Type x3 = x1*x2;
            Result(i,j) =  c00pc01y1 
                         + c10   * x1 
                         + c20   * x2 
                         + c11y1 * x1 
                         + c02y2
                         + c30   * x3
                         + c21y1 * x2
                         + c12y2 * x1
                         + c03y3
                         + c40   * x2 * x2
                         + c31y1 * x3
                         + c22y2 * x2
                         + c13y3 * x1
                         + c04y4
                         + c50   * x3 * x2
                         + c41y1 * x2 * x2
                         + c32y2 * x3
                         + c23y3 * x2
                         + c14y4 * x1
                         + c05y5;
            }
          }

// ______ do rest (>5) in loop ______
WARNING.print("polyval on grid for > degree 5 : this seems to be wrong ?? BK 9 feb.00 : Is it? just check results");
WARNING.print("polyval on grid for > degree 5 : this seems to be wrong ?? BK 9 feb.00 : Is it? just check results");
    const int32 STARTDEGREE = 6;
    const int32 STARTCOEFF  = Ncoeffs(STARTDEGREE-1);   // 5-> 21 6->28 7->36 etc.
    for (j=0; j<Result.pixels(); j++)
      {
      Type yy=y(j,0);
      for (i=0; i<Result.lines(); i++)
        {
        Type xx=x(i,0);        // ??? this seems to be wrong (BK 9-feb-00)
        Type sum = 0.;
        register int32 coeffindex = STARTCOEFF;
        for (register int32 l=STARTDEGREE; l<=degreee; l++)
          {
          for (register int32 k=0; k<=l ;k++)
            {
            sum += coeff(coeffindex,0) * pow(xx,Type(l-k)) * pow(yy,Type(k));
            coeffindex++;
            }
          }
        Result(i,j) += sum;  
        }
      }
    } // switch degreee

  return Result;
  } // END polyval


/****************************************************************
 *    polyval1d (real8, const matrix<real8> &)                  *
 *                                                              *
 * evaluates polynomial at (x)                                  *
 * input:                                                       *
 *  - x: point                                                  *
 *  - polynomial coefficients in standing vector (a0,a1,a2,..). *
 *  - degree of polynomial.                                     *
 *                                                              *
 * output:                                                      *
 *  - real8 functional value y = f(x).                          *
 *                                                              *
 *    Bert Kampes, 03-Jun-1999                                  *
 ****************************************************************/
real8 polyval1d(
        real8 x,
        const matrix<real8> &coeff)
  {
#ifdef __DEBUG
  TRACE_FUNCTION("polyval1d (BK 03-Jun-1999)")
  if (coeff.size() != coeff.lines())
    {
    PRINT_ERROR("polyval: require standing data vector.")
    throw(input_error);
    }
#endif
  real8 sum = 0.0;
  for (register int32 d=coeff.size()-1;d>=0;--d)
    {
    sum *= x;
    sum += coeff(d,0);
    }
  return sum;
  } // END polyval1d (real8, const matrix<real8> &)


/****************************************************************
 *    polyval1d (real8, const matrix<real8> &, uint)            *
 *                                                              *
 * evaluates derivative of polynomial at (x)                    *
 * input:                                                       *
 *  - x: point                                                  *
 *  - polynomial coefficients in standing vector (a0,a1,a2,..). *
 *  - derivative of polynomial                                  *
 *                                                              *
 * output:                                                      *
 *  - real8 functional value y = f(x).                          *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
real8 polyval1d(
	real8 x,
        const matrix<real8> &coeff,
	uint deriv)
  {
#ifdef __DEBUG
  TRACE_FUNCTION("polyval1d (HB 17-Mar-2011)")
  if (coeff.size() != coeff.lines())
    {
    PRINT_ERROR("polyval: require standing data vector.")
    throw(input_error);
    }
#endif
  real8 sum = 0.0;
  for (register int32 d=coeff.size()-1; d>=deriv ; d--)
    {
      int16 a = 1;
      for (register int32 p=d+1-deriv; p<=d; p++)
	a *= p;
      sum *= x;
      sum += a * coeff(d,0);
    }
  return sum;
  } // END polyval1d (real8, const matrix<real8> &, uint)


/****************************************************************
 *    normalize                                                 *
 *                                                              *
 * rescale data[min,max] to interval [-2,2]                     *
 *                                                              *
 * input:                                                       *
 *  - data                                                      *
 *  - min, max                                                  *
 * output:                                                      *
 *  - data matrix by reference                                  *
 *                                                              *
 *    Bert Kampes, 21-Jun-1999                                  *
 ****************************************************************/
void normalize(
        matrix<real4> &data,
        real8 min,
        real8 max
        )
  {
  TRACE_FUNCTION("normalize (BK 21-Jun-1999)")
  // ______ zero mean, rescale ______
  data -= real4(.5*(min+max));                          // [.5a-.5b, -.5a+.5b]
  data /= real4(.25*(max-min));                         // [-2 2]
  } // END normalize



/****************************************************************
 *    normalize                                                 *
 *                                                              *
 * rescale data to interval [-1,1]                              *
 * data -> X(data-min)/(max-min) ; data+SS                      *
 * X=1--1==EE-SS; SS=-1, EE=1                                   *
 *                                                              *
 * input:                                                       *
 *  - data                                                      *
 *  - min, max                                                  *
 * output:                                                      *
 *  - by reference                                              *
 *                                                              *
 *    Bert Kampes, 21-Jun-1999                                  *
 ****************************************************************/
void normalize(
        matrix<real8> &data,   // real8
        real8 min,
        real8 max
        )
  {
  TRACE_FUNCTION("ROUTINE: normalize (BK 21-Jun-1999)")
  // ______ zero mean, rescale ______
  data -= .5*(min+max);                         // [.5a-.5b, -.5a+.5b]
  data /= .25*(max-min);                                // [-2 2]
  } // END normalize



/****************************************************************
 *    BBparBperp                                                *
 * get Baseline, Bpar and Bperp based on 3 positions,           *
 * not efficient! See documentation for definitions             *
 *                                                              *
 #%// BK 22-Sep-2000                                            *
 ****************************************************************/
void BBparBperp(real8 &B, real8 &Bpar, real8 &Bperp,
           const cn Master, const cn Point, const cn Slave)
  {
  #ifdef __DEBUG
  TRACE_FUNCTION("BBparBperp (BK 22-Sep-2000)")
  #endif
  B = Master.dist(Slave);               // baseline. abs. value
  const real8 range1 = Master.dist(Point);
  const real8 range2 = Slave.dist(Point);
  Bpar = range1-range2;                 // parallel baseline, sign ok
  const cn r1 = Master.min(Point);
  const cn r2 = Slave.min(Point);
  const real8 theta = Master.angle(r1);
  Bperp = sqr(B)-sqr(Bpar);
  if (Bperp < 0.0) Bperp=0.0;
  else Bperp = (theta > Master.angle(r2)) ?     // perpendicular baseline, sign ok
     sqrt(Bperp) : -sqrt(Bperp);
  } // END BBparBperp



/****************************************************************
 *    BBhBv                                                     *
 * get Baseline, Bh and Bv based on 3 positions,                *
 * not efficient! See documentation for definitions             *
 * cannot get sign of Bh w/o Point P (?)                        *
 *                                                              *
 #%// BK 19-Oct-2000                                            *
 ****************************************************************/
void BBhBv(
        real8 &B,
        real8 &Bh,
        real8 &Bv,
        const cn Master,
        const cn Slave)
  {
  TRACE_FUNCTION("BBhBv (BK 19-Oct-2000)")
  WARNING.print("sign Bh not computed...");
  B                = Master.dist(Slave);// baseline. abs. value
  const real8 rho1 = Master.norm();
  const real8 rho2 = Slave.norm();
  Bv               = rho2-rho1;                 // Bv, sign ok
  Bh               = sqr(B)-sqr(Bv);
  Bh               = (Bh<0.0) ? 0.0 : sqrt(Bh);
  } // END BBhBv



/****************************************************************
 *    Btemp                                                     *
 * returns temporal baseline (days) based on 2 strings with UTC *
 * utc as: "30-AUG-1995 09:49:20.453"                           *
 *         "%d-%b-%Y %T"                                        *
 * (Btemp>0) if (Tslave later than Tmaster)                     *
 #%// BK 19-Oct-2000                                            *
 ****************************************************************/
int32 Btemp(
        const char *utc_master,
        const char *utc_slave)
  {
  TRACE_FUNCTION("Btemp (BK 19-Oct-2000)")
  struct tm tm_ref;     // reference time
  struct tm tm_master;
  struct tm tm_slave;
  char utc_ref[ONE27] = "01-JAN-1985 00:00:01.000";
  strptime(utc_ref,"%d-%b-%Y %T",&tm_ref);
  strptime(utc_master,"%d-%b-%Y %T",&tm_master);
  strptime(utc_slave, "%d-%b-%Y %T",&tm_slave);

  time_t time_ref      = mktime(&tm_ref);
  time_t time_master   = mktime(&tm_master);
  time_t time_slave    = mktime(&tm_slave);
  real8 master_since85 = difftime(time_master,time_ref);        // 1>2
  real8 slave_since85  = difftime(time_slave,time_ref);         // 1>2
  int32 Btemp = int32(floor(
    ((slave_since85-master_since85)/(60.0*60.0*24.0))+0.5));
  return Btemp;
  } // END Btemp



/****************************************************************
 *    BalphaBhBvBparBperpTheta                                  *
 * get Baseline, Bh and Bv based on 3 positions,                *
 * not efficient? See documentation for definitions             *
 * better use velocity vector, not P?                           *
 * input:                                                       *
 *  - coordinates x,y,z of Master, Slave en Point on ground     *
 * output:                                                      *
 *  - abs(B), alpha[-pi,pi], theta[approx23deg], Bh, Bv, Bpar   *
 *    Bperp                                                     *
 *                                                              *
 * theta  is the looking angle to point P.                      *
 * theta2 is the looking angle to point P for slave satellite.  *
 * cos theta  = (rho1^2+r1^2-p^2) / (2*rho1*r1)                 *
 * approx.                                                      *
 * cos theta2 = (rho2^2+r2^2-p^2) / (2*rho2*r2)                 *
 * if (theta1>theta2) then (Bperp>0)                            *
 * if (costheta1>costheta2) then (Bperp<0)                      *
 *                                                              *
 #%// BK 19-Oct-2000                                            *
 ****************************************************************/
void BalphaBhBvBparBperpTheta(
        real8 &B,
        real8 &alpha,
        real8 &Bh,
        real8 &Bv,
        real8 &Bpar,
        real8 &Bperp,
        real8 &theta,
        const cn     M,
        const cn     P,
        const cn     S)
  {
  TRACE_FUNCTION("BalphaBhBvBparBperpTheta (BK 19-Oct-2000)")
  const cn    r1        = P.min(M);
  const real8 P_2       = P.norm2();
  const real8 M_2       = M.norm2();
  const real8 r1_2      = r1.norm2();
  const real8 costheta1 = (M_2+r1_2-P_2)/(2*sqrt(M_2*r1_2));

  const real8 S_2       = S.norm2();
  const cn    r2        = P.min(S);
  const real8 r2_2      = r2.norm2();
  const real8 costheta2 = (S_2+r2_2-P_2)/(2*sqrt(S_2*r2_2));

  const cn    vecB      = S.min(M);             // (vector B points to S)
  const real8 B_2       = vecB.norm2();
  Bpar  = sqrt(r1_2)-sqrt(r2_2);                // parallel baseline, sign ok
  Bperp = (costheta1>costheta2) ?               // perpendicular bln, sign ok
    -sqrt(B_2-sqr(Bpar)) :
     sqrt(B_2-sqr(Bpar));
    
  theta = acos(costheta1);                      // <0,pi/2>, ok
  alpha = theta-atan2(Bpar,Bperp);              // <-pi,pi>
  B     = sqrt(B_2);                            // abs. value
  Bv    = B * sin(alpha);                       // Bv, sign ok
  Bh    = B * cos(alpha);                       // Bh, sign ok
  } // END BalphaBhBvBparBperpTheta















/****************************************************************
 *    shiftazispectrum                                          *
 * Shift spectrum of input matrix data either from fDC to zero, *
 * or from zero to fDC (Doppler centroid frequency).            *
 * Slcimage gives required polynomial fDC(column).              *
 * Abs(shift) gives the first range column for this polynomial. *
 *  First column is 1, so use this, not 0!                      *
 * (shift<0) indicates to shift towards zero, and               *
 * (shift>0) indicates to shift spectrum to fDC (back)          *
 * Spectrum is shifted by multiplication by a trend e^{iphase}  *
 * in the space domain.                                         *
 *                                                              *
 * If fDC is more or less equal for all columns, then an        *
 * approximation is used. If fDC is smaller then 2 percent of   *
 * PRF then no shifting is performed.                           *
 * memory intensive cause matrix is used.                       *
 *                                                              * 
 #%// BK 09-Nov-2000                                            *
 * this is not ok for ovs data *
 #%// BK 12-Aug-2005                                            *
 ****************************************************************/
void shiftazispectrum(
        matrix<complr4> &data,  // slcdata in space domain
        const slcimage  &slcinfo, // polynomial fDC
        const real4      shift) // abs(shift)==firstpix in slave system)
                                // if (shift>0) then shift from zero to fDC
  {
  // ______ Really debug this by defining: ______
  //#define REALLYDEBUG
  TRACE_FUNCTION("shiftazispectrum (BK 09-Nov-2000)")

  // ______ Evaluate fdc for all columns ______
  if (slcinfo.f_DC_a0==0 && slcinfo.f_DC_a1==0 && slcinfo.f_DC_a2==0) // no shift at all
    return;
  
  // ______ Some constants ______
  const int32 NLINES = data.lines();            // compute e^ix
  const int32 NCOLS  = data.pixels();           // width
  const int32 SIGN   = (shift<0) ? -1:1;        // e^{SIGN*i*x}
  // -1
  const real4 PIXLO  = abs(shift)-1;            // first column for fDC polynomial

  // ______ Compute fDC for all columns ______
  // ______ Create axis to evaluate fDC polynomial ______
  // ______ fDC(column) = fdc_a0 + fDC_a1*(col/RSR) + fDC_a2*(col/RSR)^2 ______
  // ______ offset slave,master defined as: cols=colm+offsetP ______
  matrix<real8> FDC(1,NCOLS);
  FDC = slcinfo.f_DC_a0;                                        // constant term
  if (slcinfo.f_DC_a1!=0 || slcinfo.f_DC_a2!=0)         // check to save time
    {
    matrix<real8> xaxis(1,NCOLS);                       // lying
    for (register int32 ii=0; ii<NCOLS; ++ii)
      xaxis(0,ii) = real8(ii+PIXLO);
    xaxis /= (slcinfo.rsr2x/2.0);
    FDC += (slcinfo.f_DC_a1 * xaxis);                   // linear term
    FDC += (slcinfo.f_DC_a2 * sqr(xaxis));              // cubic term
    #ifdef REALLYDEBUG
    cout << "REALLYDEBUG: matrix xaxis dumped.\n";
    cout << "REALLYDEBUG: SIGN   = " << SIGN << "\n";
    cout << "REALLYDEBUG: fda0   = " << slcinfo.f_DC_a0 << "\n";
    cout << "REALLYDEBUG: PIXLO  = " << PIXLO << "\n";
    cout << "REALLYDEBUG: shift  = " << shift << "\n";
    cout << "REALLYDEBUG: fda1   = " << slcinfo.f_DC_a1 << "\n";
    cout << "REALLYDEBUG: fda2   = " << slcinfo.f_DC_a2 << "\n";
    cout << "REALLYDEBUG: FDC(1) = " << FDC(0,0) << "\n";
    dumpasc("xaxis",xaxis);
    cout << "REALLYDEBUG: matrix FDC with doppler centroid dumped.\n";
    dumpasc("FDC",FDC);
    #endif
    }
  DEBUG << "fDC of first pixel: " << FDC(0,0);
  DEBUG.print();
  if (SIGN==1) DEBUG.print("Shifting from zero to fDC.");
  else         DEBUG.print("Shifting from fDC to zero.");


  // ====== Actually shift the azimuth spectrum ======
  matrix<real8> trend(NLINES,1);
  for (register int32 ii=0; ii<NLINES; ++ii)
    trend(ii,0) = (real8(2*ii)*PI) / slcinfo.prf;
  matrix<real8> P = trend * FDC;
  matrix<complr4> TREND =                               // forward or backward
    (SIGN==-1) ? mat2cr4(cos(P),-sin(P)) :
                 mat2cr4(cos(P), sin(P));


  // ______ check if this goes ok with matlab ______
  #ifdef REALLYDEBUG
  static int32 si_buffernum = 0;
  si_buffernum++;
  char OFILE[2*ONE27];
  ostrstream tmpomem(OFILE,2*ONE27);
  tmpomem.seekp(0);
  tmpomem << "datain." << si_buffernum << ends;
  cerr << "REALLYDEBUG: matrix data dumped to file " << OFILE << " (cr4).\n";
  cerr << "#lines, #pixs: " << data.lines() << ", " << data.pixels() << endl;
  ofstream of;
  of.open(OFILE); 
  of << data;
  of.close();
  tmpomem.seekp(0);
  tmpomem << "TREND." << si_buffernum << ends;
  cerr << "REALLYDEBUG: trend matrix dumped to file " << OFILE << " (cr4).\n";
  of.open(OFILE); 
  of << TREND;
  of.close();
  #endif

  data *= TREND;                                        // return by reference

  #ifdef REALLYDEBUG
  tmpomem.seekp(0);
  tmpomem << "dataout." << si_buffernum << ends;
  cerr << "REALLYDEBUG: matrix data dumped to file " << OFILE << " (cr4).\n";
  of.open(OFILE); 
  of << data;
  of.close();
  #endif
  } // END shiftazispectrum




/****************************************************************
 *    tiepoint                                                  *
 *                                                              *
 * For a given tiepoint, compute the pixel coordinates and      *
 * timing etc.  This can then be used in matlab to correct the  *
 * unwrapped phase by hand.                                     *
 * This routine is called only if tiepoint is given.            *
 *                                                              *
 * Input:                                                       *
 *  lat/lon/hei of a tiepoint                                   *
 * Output:                                                      *
 *  INFO the details                                            *
 *                                                              *
 #%// BK 14-May-2004                                            *
 ****************************************************************/
void tiepoint(
        const input_gen     &generalinput,
        const slcimage      &master,
        const slcimage      &slave,
              orbit         &masterorbit,
              orbit         &slaveorbit,
        const input_ell     &ellips)
  {
  TRACE_FUNCTION("tiepoint (Bert Kampes 14-May-2004)")
  if (abs(generalinput.tiepoint.x) < 0.001 &&
      abs(generalinput.tiepoint.y) < 0.001 &&
      abs(generalinput.tiepoint.z) < 0.001)
    {
    INFO.print("No tiepoint given.");
    return;
    }
  if (masterorbit.npoints()==0)
    {
    INFO.print("Exiting tiepoint, no orbit data master available.");
    return;
    }
  DEBUG.precision(30); // permanent; INFO << setp(13)<<q;// once
  DEBUG.width(14);
  INFO.precision(13); // permanent; INFO << setp(13)<<q;// once
  INFO.width(14);
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;
  const real8 m_minpi4cdivlam = (-4.0*PI*SOL)/master.wavelength;
  const real8 s_minpi4cdivlam = (-4.0*PI*SOL)/slave.wavelength;
  DEBUG << "Master wavelength = " << master.wavelength;
  DEBUG.print();
  DEBUG << "Slave wavelength  = " << slave.wavelength;
  DEBUG.print();


  ;// ______ Convert to X,Y,Z ______
  INFO.print("Computing radar coordinates and unwrapped phase of tie point");
  real8 tiepointphi    = deg2rad(generalinput.tiepoint.x);// lat in [dec.degrees]
  real8 tiepointlambda = deg2rad(generalinput.tiepoint.y);// lon in [dec.degrees]
  real8 tiepointheight = generalinput.tiepoint.z;//          height in [m]
  DEBUG << "TIEPOINT: lat/lon/hei [rad/rad/m]: " 
        << tiepointphi << " " << tiepointlambda << " " << tiepointheight;
  DEBUG.print();
  //cn tiepointpos;// compute phi/lambda/h-->x,y,z
  //ell2xyz(ellips,tiepointpos,tiepointphi,tiepointlambda,tiepointheight);
  cn tiepointpos = ellips.ell2xyz(tiepointphi,tiepointlambda,tiepointheight);
  INFO  << "TIEPOINT: x,y,z (WGS84):           "
        << tiepointpos.x << " " << tiepointpos.y << " " << tiepointpos.z;
  INFO.print();

  // ______ Convert to line/pix in master ______
  real8 tiepointline;// returned
  real8 tiepointpix;// returned
  int32 n_iter = xyz2lp(tiepointline, tiepointpix,
                        master, masterorbit, tiepointpos,
                        MAXITER, CRITERTIM);
  INFO  << "TIEPOINT line/pix in master:       " << tiepointline << " " << tiepointpix;
  INFO.print();

  // ______ Compute master azimuth/range time of this pixel______
  const real8 m_aztime = master.line2ta(tiepointline);
  const real8 m_trange = master.pix2tr(tiepointpix);
  DEBUG << "TIEPOINT master azimuth time:      " << m_aztime;
  DEBUG.print();
  DEBUG << "TIEPOINT master range time:        " << m_trange;
  DEBUG.print();

  // ______ Ellipsoid reference phase for this point master(l,p) ______
  // ______ Compute xyz for slave satellite from P ______
  cn P_ref;// point at H=0 (ellips)
  lp2xyz(tiepointline,tiepointpix,
         ellips,master,masterorbit,P_ref,
         MAXITER,CRITERPOS);
  INFO  << "TIEPOINT: same range, but on ellipsoid: x,y,z (WGS84): "
        << P_ref.x << " " << P_ref.y << " " << P_ref.z;
  INFO.print();

  // ______ Check if slave is there ______
  if (slaveorbit.npoints()==0)
    {
    INFO.print("tiepoint: no orbit data slave available cannot compute more.");
    DEBUG.reset();
    INFO.reset();
    }
  else
    {
    // ______ Compute xyz for slave satellite from P ______
    real8 s_aztime;                             // returned
    real8 s_trange;                             // returned
    xyz2t(s_aztime, s_trange,
          slave, slaveorbit, tiepointpos, MAXITER, CRITERTIM);
    DEBUG << "TIEPOINT slave azimuth time:       " << s_aztime;
    DEBUG.print();
    DEBUG << "TIEPOINT slave range time:         " << s_trange;
    DEBUG.print();
    INFO  << "TIEPOINT line/pix in slave:        " 
         << slave.ta2line(s_aztime) << " "
         << slave.tr2pix (s_trange);
    INFO.print();

    // ___ Phase of this point ___
    real8 phase = m_minpi4cdivlam*m_trange - s_minpi4cdivlam*s_trange;// [rad]
    INFO << "TIEPOINT unwrapped phase (incl. reference phase): "
         << phase << " [rad]";
    INFO.print();
  
    // _____ Times for slave image ______
    real8 s_tazi_ref;                               // returned
    real8 s_trange_ref;                             // returned
    xyz2t(s_tazi_ref, s_trange_ref,
          slave, slaveorbit, P_ref, MAXITER,CRITERTIM);

    // ___ report reference phase corrected unwrapped phase of tiepoint ___
    // ___ real8 m_trange_ref = m_trange;// by definition
    // ___ t_corrected = t_m - t_m_ref - (t_s - t_s_ref) = t_s_ref-t_s;
    DEBUG << "TIEPOINT slave range time to Pref: " << s_trange_ref;
    DEBUG.print();
    real8 phase_ref_corrected = s_minpi4cdivlam*(s_trange_ref-s_trange);// [rad]
    INFO << "TIEPOINT unwrapped phase (ellipsoid ref.phase corrected) [rad]: "
         << phase_ref_corrected;
    INFO.print();
    INFO.print("This phase can be used to correct the unwrapped interferogram");
    INFO.print("make sure that you compute the correct pixel position");
    INFO.print("the line/pix here are in original master coordinates.");
    INFO.print("And use Matlab to add this phase to the unwrapped interferogram.");
    INFO.print();


    // ______ Check if zero-Doppler iteration is correct by going back ______
    // ------ Go from P_ref to M and to S, then back from M to P_refcheckM
    // ------ and from S to P_refcheckS; Bert Kampes, 06-Apr-2005
    real8 t_az, t_rg;// tmp variables
    // --- from P_ref to orbits ---
    // ______ master: from P_ref to M ____
    n_iter     = xyz2t(t_az,t_rg,master,masterorbit,P_ref,MAXITER,CRITERTIM);// return t
    const cn M = masterorbit.getxyz(t_az);
    // ______ master: from M to P_ref_checkM ____
    cn P_ref_checkM;
    n_iter     = lp2xyz(master.ta2line(t_az), master.tr2pix(t_rg),
                    ellips, master, masterorbit, 
                    P_ref_checkM, MAXITER, CRITERPOS);
    // ______ slave: from P_ref to S ____
    n_iter     = xyz2t(t_az,t_rg,slave, slaveorbit, P_ref,MAXITER,CRITERTIM);// return t
    const cn S = slaveorbit.getxyz(t_az);
    // ______ slave: from S to P_ref_checkS ____
    cn P_ref_checkS;
    n_iter     = lp2xyz(slave.ta2line(t_az), slave.tr2pix(t_rg),
                    ellips, slave, slaveorbit, 
                    P_ref_checkS, MAXITER, CRITERPOS);
    // ______ finally back again from P_ref to slave and master ______
    cn M_check;
    n_iter     = xyz2orb(M_check,master,masterorbit,P_ref,MAXITER,CRITERTIM);
    cn S_check;
    n_iter     = xyz2orb(S_check,slave, slaveorbit, P_ref,MAXITER,CRITERTIM);
    // ______ Report: ______
    INFO.print("Checking zero-Doppler iteration by doing forward and inverse transform:");
    DEBUG << "M(x,y,z) = " << M.x     << "; " << M.y     << "; "<< M.z;
    DEBUG.print();
    DEBUG << "M'(x,y,z)= " << M_check.x     << "; " << M_check.y     << "; "<< M_check.z;
    DEBUG.print();
    DEBUG << "P(x,y,z) = " << P_ref.x << "; " << P_ref.y << "; "<< P_ref.z;
    DEBUG.print();
    DEBUG << "P'(x,y,z)= " << P_ref_checkM.x << "; " << P_ref_checkM.y << "; "<< P_ref_checkM.z;
    DEBUG.print();
    DEBUG << "P''(x,y,z)=" << P_ref_checkS.x << "; " << P_ref_checkS.y << "; "<< P_ref_checkS.z;
    DEBUG.print();
    DEBUG << "S(x,y,z) = " << S.x     << "; " << S.y     << "; "<< S.z;
    DEBUG.print();
    DEBUG << "S'(x,y,z)= " << S_check.x     << "; " << S_check.y     << "; "<< S_check.z;
    DEBUG.print();
    if (abs(M.dist(M_check))>0.001)// assume 1 mm or so?
      WARNING.print("check failed for M_check");
    if (abs(S.dist(S_check))>0.001)// assume 1 mm or so?
      WARNING.print("check failed for S_check");
    if (abs(P_ref.dist(P_ref_checkM))>0.001)// assume 1 mm or so?
      WARNING.print("check failed for P_ref_checkM");
    if (abs(P_ref.dist(P_ref_checkS))>0.001)// assume 1 mm or so?
      WARNING.print("check failed for P_ref_checkS");
    INFO.print("If no warnings seen then zero-Doppler is within 1 mm accurate.");
    }// slave orbit present

  // ______ Tidy up ______
  PROGRESS.print("Finished tiepoint.");
  DEBUG.reset();
  INFO.reset();
  } // END tiepoint


/****************************************************************
 * offsets2timing()                                             *
 *                                                              *
 *                                                              *
 * input:                                                       *
 *  - master info                                               *
 *  - simulated amplitude overall offsetL and offsetP           * 
 * output:                                                      *
 *  - master timing error in the azimuth and range direction    *
 *                                                              *
 *  : see also  slcimage.hh                                     *
 *                                                              *
 *    Mahmut Arikan, 13-Dec-2008                                *
 ****************************************************************/
void offsets2timing(
        const slcimage        &slcinfo, // normalization factors, ovs_rg/az
        const int32           &offsetL,
        const int32           &offsetP,
        real8                 &azTIME,
        real8                 &rgTIME)
  {
  TRACE_FUNCTION("offsets2timing (MA 13-DEC-2008)")
  
  // switch INFO-->DEBUG when safe.  
  INFO << "The given azimuth timing offset: " << offsetL << " lines.";
  INFO.print();
  INFO << "The given   range timing offset: " << offsetP << " pixels.";
  INFO.print();

  // if offsetL is NaN exit
  // TODO check if line pixel offsets are more than certian threshold ex: > 1000 in
  // lines > 200 in pixels
  // this should be a warning
  
  
  // ______ Initialize Variables ______
  const real8 &mPRF   =  slcinfo.prf   ;   // PRF  see slcimage.hh
  const real8 &mRSRx2 =  slcinfo.rsr2x ;   // 2xRSR
                                          // oversampling updates these values
                                          // variables see slcimage.cc
  DEBUG << "PRF: " << mPRF;
  DEBUG.print();
  DEBUG     << "RSR: " << mRSRx2/2.0 ;
  DEBUG.print();
  
  //slcinfo.showdata(); // debug
  
  // ______ Compute Time ______
  azTIME  = real8( offsetL/mPRF   );   // Azimuth timing, in seconds
  rgTIME  = real8( offsetP/mRSRx2 );   // Range timing,   in seconds ( two-way)


  INFO << "Original azimuth timing shift [sec]: " << slcinfo.az_timing_error << " sec.";
  INFO.print();
  INFO << "Original   range timing shift [sec]: " << slcinfo.r_timing_error <<  " sec.";
  INFO.print();

//  following already in mtiming, keep for future implementations
//  DEBUG << "Estimated master azimuth timing error [sec]: " << azTIME << " sec.";
//  DEBUG.print();
//  DEBUG << "Estimated master range timing error   [sec]: " << rgTIME << " sec.";
//  DEBUG.print();

  } // END offsets2timing

/****************************************************************
 *    Write Template KML                                        *
 *                                                              *
 *  Generate a KML file for displaying crop and interferograms  *
 * on Google Earth, based on 4 corner coordinates from orbits.  *
 *                                                              *
 * Batuhan Osmanoglu <batu@rsmas.miami.edu>                     *
 *                                                              *
 ****************************************************************/
void write_kml(
        const input_ell        &input_ellips,
        const window           &input_dbow,
        const slcimage         &master,
        orbit                  &masterorbit,
	      const char		       filename[4*ONE27])
  {
  TRACE_FUNCTION("write_kml (MA, FvL, BO 20100123)")
      real8 cornerPhi, cornerLam, cornerHei, One80overPi;
      real8 north=0;//initialize values...
      real8 south=0;
      real8 east=0;
      real8 west=0;
      double rotation=0;
      double rotationNumerator=0;
      double rotationDenominator=0;
      One80overPi=180/PI;
      char pngFilename[4*ONE27]; 
      strcpy(pngFilename, filename);
      strcat(pngFilename,".png");
	//corner 1 (l0,p0)
      lp2ell(input_dbow.linelo, input_dbow.pixlo, input_ellips, master, masterorbit, cornerPhi, cornerLam, cornerHei, 10, 1e-3);
      cornerPhi*=One80overPi;
      cornerLam*=One80overPi;
      DEBUG << "philamhei(l0,p0)= " << cornerPhi << ", " << cornerLam << ", " << cornerHei;
      DEBUG.print();
      north+=cornerPhi;			//Image should flip or flop automatically.   
      west+= cornerLam;	
      rotationNumerator=cornerPhi;
      rotationDenominator=cornerLam;
	//corner 2 (l0,pN)      
      lp2ell(input_dbow.linelo, input_dbow.pixhi, input_ellips, master, masterorbit, cornerPhi, cornerLam, cornerHei, 10, 1e-3);
      cornerPhi*=One80overPi;
      cornerLam*=One80overPi;
      DEBUG << "philamhei(l0,pN)= " << cornerPhi << ", " << cornerLam << ", " << cornerHei;
      DEBUG.print();
      north+=cornerPhi;north/=2;	//north=[ phi(l0,p0)+phi(l0,pN) ]/2
      east+=cornerLam;
      	//corner 3 (lN,p0)
      lp2ell(input_dbow.linehi, input_dbow.pixlo, input_ellips, master, masterorbit, cornerPhi, cornerLam, cornerHei, 10, 1e-3);
      cornerPhi*=One80overPi;
      cornerLam*=One80overPi;
      DEBUG << "philamhei(lN,p0)= " << cornerPhi << ", " << cornerLam << ", " << cornerHei;
      DEBUG.print();
      south+=cornerPhi;
      west+=cornerLam;west/=2;		//west= [ lam(l0,pN)+lam(lN,pN) ]/2
      rotationNumerator-=cornerPhi;
      rotationDenominator-=cornerLam;
      	//corner 4 (lN,pN)
      lp2ell(input_dbow.linehi, input_dbow.pixhi, input_ellips, master, masterorbit, cornerPhi, cornerLam, cornerHei, 10, 1e-3);
      cornerPhi*=One80overPi;
      cornerLam*=One80overPi;
      DEBUG << "philamhei(lN,pN)= " << cornerPhi << ", " << cornerLam << ", " << cornerHei;
      DEBUG.print();
      south+=cornerPhi;south/=2;	//south=[ phi(lN,p0)+phi(lN,pN) ]/2
      east+=cornerLam;east/=2;		//east= [ lam(l0,p0)+lam(lN,p0) ]/2
            
      DEBUG << "GoogleEarth North/South/East/West: " << north << " / " << south << " / " << east << " / " << west;
      DEBUG.print();
      
      rotation=-(90-atan2(rotationNumerator, rotationDenominator)*One80overPi);
      DEBUG << "GoogleEarth Rotation: " << rotation;
      DEBUG.print();
      // Write output KML. 
	//open output file
	ofstream kmlfile;
	openfstream(kmlfile, filename, true);
	bk_assert(kmlfile, filename, __FILE__, __LINE__);
	kmlfile << "<kml xmlns=\"http://earth.google.com/kml/2.1\">" << endl;
	kmlfile << "<Document>" << endl;
	kmlfile << "<GroundOverlay>" << endl;
        kmlfile << "        <visibility>1</visibility>" << endl;
        kmlfile << "        <LatLonBox>" << endl;
        kmlfile << "                <north>" << north << "</north>" << endl;
        kmlfile << "                <south>" << south << "</south>" << endl;
        kmlfile << "                <east>" << east << "</east>" << endl;
        kmlfile << "                <west>" << west << "</west>" << endl; 
        kmlfile << "                <rotation>" << rotation << "</rotation>" << endl;
        kmlfile << "        </LatLonBox>" << endl;
        kmlfile << "        <Icon>" << endl;
        kmlfile << "                <href>" << pngFilename << "</href>" << endl;
        kmlfile << "        </Icon>" << endl;
        kmlfile << "</GroundOverlay>" << endl;
	kmlfile << "</Document>" << endl;
	kmlfile << "</kml>" << endl;
	kmlfile.close();
      //End of -> Output Crop Corner Coordinates.
      //////////////////////////////////////////

  } //END write_kml

