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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixspecs.cc,v $
 * $Revision: 3.23 $
 * $Date: 2005/10/07 15:17:21 $
 * $Author: kampes $
 *
 * definition of routines for matrices (helpers)
 ****************************************************************/


#include "matrixbk.hh"
#include "utilities.hh"                 // fast_sin()
#include <strstream>                    // 
#include <algorithm>                    // max

#ifdef __DEBUGMAT1
#include <fstream>                      // type ofstream
#endif

//#if defined (__DEBUGMAT2) || defined (__DEBUGMAT1)
#include <cstdlib>                      // exit
//#endif



// ______ message objects, global, set in main ______
extern bk_messages matERROR;
extern bk_messages matDEBUG;



// ______ Prototype 1d fft in this file ______
// ______ Only change this routine for other library for fft _____
// ______ At the moment one can select VECLIB/FFTW/INTERNAL (slow?)
// ______ VECLIB/FFTW call optimized routines for 2dfft's
void four1(
        complr4 data[],
        int32 fftlength,
        int32 isign);

// ______ Prototype veclib functions ______
// ====== File matrixspecialization.c ======
#ifdef __USE_VECLIB_LIBRARY__
extern "C" { int32 sgemm (char* , char*, int32*, int32*, int32*,
              real4*, real4*, int32*,
              real4*, int32*, real4*,
              real4*, int32*, int32, int32); }
extern "C" { int32 dgemm (char* , char*, int32*, int32*, int32*,
              real8*, real8*, int32*,
              real8*, int32*, real8*,
              real8*, int32*, int32, int32); }
extern "C" { int32 cgemm (char* , char*, int32*, int32*, int32*,
              complr4*, complr4*, int32*,
              complr4*, int32*, complr4*,
              complr4*, int32*, int32, int32); }
extern "C" { int32 zgemm (char* , char*, int32*, int32*, int32*,
              complr8*, complr8*, int32*,
              complr8*, int32*, complr8*,
              complr8*, int32*, int32, int32); }
#endif

#ifdef __USE_VECLIB_LIBRARY__
extern "C" { int32 c1dfft (complr4*, int32*, real4*, int32*, int32*); }
extern "C" { int32 c2dfft (complr4*, int32*, int32*, int32*, int32*, int32*); }
//extern "C" { int32 s2dfft (real4*, real4*, int32*, int32*, int32*, int32*, int32*); }
#endif

// ______ USE FFTW IF AVAILABLE ______
// first include complex then this to have complex type known?
#ifdef __USE_FFTW_LIBRARY__
  #include <fftw3.h>                                  // fftw types and functions v3.0.1
#endif

#ifdef __USE_LAPACK_LIBRARY__
  // ___ FORTRAN version on SUN solaris (appends an underscore) ___
  #ifdef __DEBUGMAT2
    //matDEBUG.print("Using LAPACK library fortran version."); // MA this won't work here [DEL]
      #warning "Using LAPACK library fortran version."
  #endif
  #define  spotrf spotrf_
  #define  spotri spotri_
  #define  spotrs spotrs_
  #define  dpotrf dpotrf_
  #define  dpotri dpotri_
  #define  dpotrs dpotrs_
  #ifdef WIN32
  // ___ Windows version, Jia defined this ___
  // Bert Kampes, 24-Aug-2005
  #define  spotrf SPOTRF
  #define  spotri SPOTRI
  #define  spotrs SPOTRS
  #define  dpotrf DPOTRF
  #define  dpotri DPOTRI
  #define  dpotrs DPOTRS
  #endif
  // ___ C version on HP uses names "as is" ___
extern "C" { int32 spotrf (const char*, int32*, real4*, int32*, int32*, int32); } // [MA] fixed: deprecated conversion from string constant to char*
extern "C" { int32 spotri (const char*, int32*, real4*, int32*, int32*, int32); }
extern "C" { int32 spotrs (const char*, int32*, int32*, real4*, int32*,
                           real4*, int32*, int32*, int32); }
extern "C" { int32 dpotrf (const char*, int32*, real8*, int32*, int32*, int32); }
extern "C" { int32 dpotri (const char*, int32*, real8*, int32*, int32*, int32); }
extern "C" { int32 dpotrs (const char*, int32*, int32*, real8*, int32*,
                           real8*, int32*, int32*, int32); }
#endif




/****************************************************************
 * myispower2                                                   *
 *    Bert Kampes, 22-Mar-2000                                  *
 ****************************************************************/
inline bool myispower2(uint w)
{if (w==2    || w==4     || w==8    || w==16   ||
     w==32   || w==64    || w==128  || w==256  ||
     w==512  || w==1024  || w==2048 || w==4096 ||
     w==8192 || w==16384 || w==32768 )
   return true;
 return false;}



/****************************************************************
 *    matassert                                                 *
 * check if file is opened ok                                   *
 *    Bert Kampes, 16-Jun-1999                                  *
 ****************************************************************/
void matassert(
        const ifstream &str,
        const char* ifilename,
        const char* callingfile,
        int32 linenumber)
  {
  if (!str)
    {
    matERROR << "Cannot open file: \""  << ifilename
         << "\" (input). Source: \""          << callingfile
         << "\", line: "                      << linenumber;
    matERROR.print();
    }
  #ifdef __DEBUGMAT2
  matDEBUG << "OK input stream associated with file:  \""  
              << ifilename << "\"";
  matDEBUG.print();
  #endif
  } // END matassert
/****************************************************************
 *    matassert                                                 *
 * check if file is opened ok                                   *
 *    Bert Kampes, 16-Jun-1999                                  *
 ****************************************************************/
void matassert(
        const ofstream &str,
        const char* ofilename,
        const char* callingfile,
        int32 linenumber)
  {
  if (!str)
    {
    matERROR << "Cannot open file: \""  << ofilename
         << "\" (output). Source: \""         << callingfile
         << "\", line: "                      << linenumber
         << "; OVERWRIT OFF/non existing directory?";
    matERROR.print();
    }
  #ifdef __DEBUGMAT2
  matDEBUG << "OK output stream associated with file: \""  
              << ofilename << "\"";
  matDEBUG.print();
  #endif
  } // END matassert



/****************************************************************
 *    fileci2tomatcr4                                           *
 * reads file format complex short (or 2*short)                 *
 * fills matrix complex real4                                   *
 *  window specified in system of file. win(10:20) means        *
 *  start at 10th line after 1st Byte on file.                  *
 * currentwindow should be passed to account for any offset     *
 * between what is linenumber of first line in file and byte = 1*
 *    Bert Kampes, 14-Jan-1999                                  *
 ****************************************************************/
void fileci2tomatcr4(
        matrix<complr4> &Result,
        const char *file,
        uint filelines,
        window win,
        window winoffset)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("ROUTINE: fileci2tomatcr4.");
#endif

// ______First account for possible offsets of file______
  win.linelo = win.linelo - winoffset.linelo + 1;
  win.linehi = win.linehi - winoffset.linelo + 1;
  win.pixlo  = win.pixlo  - winoffset.pixlo  + 1;
  win.pixhi  = win.pixhi  - winoffset.pixlo  + 1;
          
  // ifstream ifile(file, ios::in | ios::ate | ios::binary);
  // g++ seems to haeve a prob. with this... BK 130700
  // ifstream ifile(file, ios::in | ios::app | ios::binary);
  //ifstream ifile(file, ios::in | ios::binary | ios::nocreate);
#ifdef __NO_IOS_BINARY__
  ifstream ifile(file, ios::in);
#else
  ifstream ifile(file, ios::in | ios::binary);
#endif
//  ifstream ifile;
//  openfstream(ifile,file);
  ifile.seekg(0,ios::end);                      // pointer to end...
  matassert(ifile,file,__FILE__,__LINE__);
// const uint sizefile = ifile.tellg();         // opened ate
  const streamoff &sizefile = ifile.tellg();            // opened ate, [MA]
 
  if (sizefile==0)
    matERROR.print("filesize==0...");


// ______Check input______
  const uint pixelsxbytes = sizefile/filelines;
  const uint SIZE = sizeof(compli16);           // size per element on disk (=4)
  const uint filepixels = pixelsxbytes/SIZE;
  const uint sizepixel  = pixelsxbytes/filepixels;
  const uint lines  = win.lines();
  const uint pixels = win.pixels();
  //const uint lines  = win.linehi - win.linelo + 1;
  //const uint pixels = win.pixhi - win.pixlo + 1;
  //const uint start  = ((win.linelo-1)*filepixels+win.pixlo-1)*sizepixel;
  const uint64 start  = (uint64)((win.linelo-1)*filepixels+win.pixlo-1)*sizepixel; //[MA], see matrixbk.hh for details on this line

#ifdef __DEBUGMAT2
  if (win.linelo<=0)
    matERROR.print("minimum line <= 0, (should be 1?).");
  if (win.pixlo<=0)
    matERROR.print("minimum pixel <= 0, (should be 1?).");
  if (win.linelo>win.linehi)
    matERROR.print("minimum line > max line.");
  if (win.pixlo>win.pixhi)
    matERROR.print("minimum pixel > max pixel.");
  if (win.linehi>filelines)
    matERROR.print("max. line is larger than on file.");
  if (win.pixhi>filepixels)
    matERROR.print("max. pixel is larger than on file.");
  if (SIZE != sizepixel)
    matERROR.print("Formatflag not consistent with file.");
#endif

#ifdef __DEBUGMAT2
  if (Result.lines() < lines || Result.pixels() < pixels)
    matERROR.print("code ?::fileci2tomatcr4: Result matrix not ok.");
  if (Result.lines() != lines || Result.pixels() != pixels)
    matDEBUG.print("debug info: fileci2tomatcr4: matrix not fully filled.");
#endif

  int16 el[2];                                          // pseudo complex type
// ______Convert io_____
  for (int32 lin=0; lin<lines; ++lin)
    {
    ifile.seekg(start+filepixels*lin*sizepixel,ios::beg);    // data at row: lin
    for (int32 pix=0; pix<pixels; ++pix)
      {
      ifile.read((char*)el,SIZE);                       // element compli16 = 4
      Result(lin,pix)=complr4(el[0],el[1]);
      //Result(lin,pix)=complr4(fileci16.real(),fileci16.imag());
      }
    }
  ifile.close();
  } // END fileci2tomatcr4



/****************************************************************
 *    magnitude                                                 *
 * Computes magnitude of complex matrix                         *
 *                                                              *
 *    Bert Kampes, 14-Jan-1999                                  *
 ****************************************************************/
matrix<real4> magnitude(
        const matrix<complr4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("magnitude.");
#endif
// _____Ensure stride one in memory______
  matrix<real4> Result(A.lines(),A.pixels());
  complr4 *pntA = A[0];
  real4 *pntR = Result[0];
  //for (register int32 i=0; i<Result.size(); i++)
  //  (*pntR++) = sqrt(sqr((*pntA).real()) + sqr((*pntA++).imag()));
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<Result.size(); i++)
    {
      *(pntR++) = sqrt(sqr((*pntA).real()) + sqr((*pntA).imag()));
        pntA++;                                         // [MA]
    }
  return Result;
  } // END magnitude



/****************************************************************
 *    intensity                                                 *
 * Computes intensity of complex matrix                         *
 *    Bert Kampes, 14-Jan-1999                                  *
 ****************************************************************/
matrix<real4> intensity(
        const matrix<complr4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("intensity.");
#endif
// _____Ensure stride one in memory______
  matrix<real4> Result(A.lines(),A.pixels());
  complr4 *pntA = A[0];
  real4 *pntR = Result[0];
  //for (register int32 i=0; i<Result.size(); i++)
  //  (*pntR++) = sqr((*pntA).real()) + sqr((*pntA++).imag());
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<Result.size(); i++)
    {
      *(pntR++) = sqr((*pntA).real()) + sqr((*pntA).imag());
       pntA++;
    }
  return Result;
  } // END intensity



#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * choles(A);   cholesky factorisation lapack                   *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void choles(
        matrix<real4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("choles: LAPACK");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrixes for cholesky.");
#endif

  int32 lda = A.pixels();                // leading dimension
  int32 n   = lda;                       // order of A
  int32 ierr;                           // 
  spotrf("U",&n,A[0],&lda,&ierr,1);     // Upper remains unchanged

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("choles: ok. lower contains factorisation.");
  else if (ierr <  0)
    matERROR.print("choles: the kth element is not ok.");
  else 
    matERROR.print("choles: not positive definite.");
#endif
  } // END choles


// ______ Else use simple intern algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * choles(A);   cholesky factorisation internal implementation  *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 * this one is a lot slower then veclib and there may be more   *
 * efficient implementations.                                   *
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void choles(matrix<real4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("choles: internal");
#endif
  const int32 N = A.lines();
  register real4 sum;
  for (register int32 i=0; i<N; ++i)
    {
    for (register int32 j=i; j<N; ++j)
      {
      sum = A[i][j];
      for (register int32 k=i-1; k>=0; --k)
        {
        sum -= A[i][k] * A[j][k];
        }
      if (i == j)
        {
        if (sum <= 0.) {matERROR.print("choles: internal: A not pos. def.");}
        A[i][i] = sqrt(sum);
        }
      else
        {
        A[j][i] = sum/A[i][i];
        }
      }
    }
  } // END choles internal, self bk
#endif


#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * choles(A);   cholesky factorisation lapack                   *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void choles(
        matrix<real8> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("choles: LAPACK");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrixes for cholesky.");
#endif

  int32 lda= A.pixels();                // leading dimension
  int32 n  = lda;                       // order of A
  int32 ierr;                           // 
  dpotrf("U",&n,A[0],&lda,&ierr,1);     // Upper remains unchanged

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("choles: ok. lower contains factorisation.");
  else if (ierr <  0)
    matERROR.print("choles: the kth element is not ok.");
  else 
    matERROR.print("choles: not positive definite.");
#endif
  } // END choles

// ______ Else use simple intern algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * choles(A);   cholesky factorisation internal implementation  *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 * this one is a lot slower then veclib and there may be more   *
 * efficient implementations.                                   *
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void choles(matrix<real8> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("choles: internal");
#endif
  const int32 N = A.lines();
  register real8 sum;
  for (register int32 i=0; i<N; ++i)
    {
    for (register int32 j=i; j<N; ++j)
      {
      sum = A[i][j];
      for (register int32 k=i-1; k>=0; --k)
        {
        sum -= A[i][k] * A[j][k];
        }
      if (i == j)
        {
        if (sum <= 0.) {matERROR.print("choles: internal: A not pos. def.");}
        A[i][i] = sqrt(sum);
        }
      else
        {
        A[j][i] = sum/A[i][i];
        }
      }
    }
  } // END choles internal, self bk
#endif


#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * solvechol(A,rhs); solution of AX=rhs                         *
 *  cholesky factorisation lapack                               *
 * A contains cholesky factorisation of A                       *
 * rhs contains estimated X on output                           *
 *                                                              *
 * known BUG: niet mogelijk meer dan 1 rechterlid simultaan op  *
 *  te lossen, waarschijnlijk door fortran <> c tegenstelling   *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void solvechol(
        const matrix<real4> &A,
        matrix<real4> &B)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("solvechol::LAPACK");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrices for cholesky.");
  if (A.lines() != B.lines())
    matERROR.print("solvechol: must same size a,b.");
  if (B.pixels() != 1)
    matERROR.print("solvechol: NOT possible simultaneous solution in this way.");
#endif

  int32 lda= A.pixels();                // leading dimension
  int32 ldb= B.lines();                 // leading dimension
  int32 nrhs= B.pixels();               // number of righthandsides
  int32 n  = lda;                       // order of A
  int32 ierr;                           // 
  spotrs("U",&n,&nrhs,A[0],&lda,B[0],&ldb,&ierr,1);     // factorisation stored in lower

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("solvechol: ok, rhs contains solution.");
  else if (ierr <  0)
    matERROR.print("solvechol: the kth element is not ok.");
  else 
    matERROR.print("solvechol: impossible ierr.");
#endif
  } // END solvechol

// ______ Else use simple (inaccurate, slow?) internal algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * solvechol(A,rhs); solution of AX=rhs                         *
 *  cholesky factorisation internal implemetnation              *
 * A contains cholesky factorisation of A                       *
 * rhs contains estimated X on output                           *
 * there may be more efficient implementations.                 *
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void solvechol(
        const matrix<real4> &A,
        matrix<real4> &B)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("solvechol::INTERNAL");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("solvechol: INTERNAL: symmetrical matrices for cholesky.");
  if (A.lines() != B.lines())
    matERROR.print("solvechol: must same size a,b.");
  if (B.pixels() != 1)
    matERROR.print("solvechol: NOT possible simultaneous solution in this way.");
#endif
  const int32 N = A.lines();
  register real4 sum;
  register int32 i,j;
// ______ Solve Ly=b, use B to store y ______
  for (i=0; i<N; ++i)
    {
    sum = B[i][0];
    for (j=i-1; j>=0; --j)
      {
      sum -= A[i][j]*B[j][0];
      }
    B[i][0] = sum/A[i][i];
    }
// ______ Solve Ux=y, use B to store unknowns ______
  for (i=N-1; i>=0; --i)
    {
    sum = B[i][0];
    for (j=i+1; j<N; ++j)
      {
      sum -= A[j][i]*B[j][0];
      }
    B[i][0] = sum/A[i][i];
    }
  } // END solvechol
#endif


#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * solvechol(A,rhs); solution of AX=rhs                         *
 *  cholesky factorisation lapack                               *
 * A contains cholesky factorisation of A                       *
 * rhs contains estimated X on output                           *
 *                                                              *
 * known BUG: niet mogelijk meer dan 1 rechterlid simultaan op  *
 *  te lossen, waarschijnlijk door fortran <> c tegenstelling   *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void solvechol(
        const matrix<real8> &A,
        matrix<real8> &B)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("solvechol: LAPACK: there was problem here, ldb,nrhs.");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrixes for cholesky.");
  if (A.lines() != B.lines())
    matERROR.print("solvechol: must same size A,B.");
  if (B.pixels() != 1)
    matERROR.print("solvechol: NOT possible simultaneous solution in this way.");
#endif

  int32 lda  = A.pixels();              // leading dimension
  int32 ldb  = B.lines();               // leading dimension
  int32 nrhs = B.pixels();              // number of righthandsides
  int32 n    = lda;                     // order of A
  int32 ierr;                           // 
  dpotrs("U",&n,&nrhs,A[0],&lda,B[0],&ldb,&ierr,1);     // factorisation stored in lower

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("solvechol: ok, rhs contains solution.");
  else if (ierr <  0)
    matERROR.print("solvechol: the kth element is not ok.");
  else 
    matERROR.print("solvechol: impossible ierr.");
#endif
  } // END solvechol

// ______ Else use simple (inaccurate, slow?) internal algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * solvechol(A,rhs); solution of AX=rhs                         *
 *  cholesky factorisation internal implemetnation              *
 * A contains cholesky factorisation of A                       *
 * rhs contains estimated X on output                           *
 * there may be more efficient implementations.                 *
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void solvechol(
        const matrix<real8> &A,
        matrix<real8> &B)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("solvechol::INTERNAL");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("solvechol: INTERNAL: symmetrical matrices for cholesky.");
  if (A.lines() != B.lines())
    matERROR.print("solvechol: must same size a,b.");
  if (B.pixels() != 1)
    matERROR.print("solvechol: NOT possible simultaneous solution in this way.");
#endif
  const int32 N = A.lines();
  register real8 sum;
  register int32 i,j;
// ______ Solve Ly=b, use B to store y ______
  for (i=0; i<N; ++i)
    {
    sum = B[i][0];
    for (j=i-1; j>=0; --j)
      {
      sum -= A[i][j]*B[j][0];
      }
    B[i][0] = sum/A[i][i];
    }
// ______ Solve Ux=y, use B to store unknowns ______
  for (i=N-1; i>=0; --i)
    {
    sum = B[i][0];
    for (j=i+1; j<N; ++j)
      {
      sum -= A[j][i]*B[j][0];
      }
    B[i][0] = sum/A[i][i];
    }
  } // END solvechol
#endif


#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * invertchol(A); inversion by cholesky factorisation lapack    *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void invertchol(
        matrix<real4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("invertchol: LAPACK");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrixes for cholesky.");
#endif

  int32 lda= A.pixels();                // leading dimension
  int32 n  = lda;                       // order of A
  int32 ierr;                           // 
  spotri("U",&n,A[0],&lda,&ierr,1);     // Upper remains unchanged

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("invertcholes: ok, lower contains inverse.");
  else if (ierr <  0)
    matERROR.print("invertcholes: the kth element is not ok.");
  else 
    matERROR.print("invertcholes: singular.");
#endif
  } // END invertchol


// ______ Else use simple (inaccurate? slow?) internal algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * invertchol(A); inversion by cholesky factorisation internal  *
 * lower triangle of A is changed on output                     *
 * upper remains un referenced                                  *
 * implementation with double loop probably can be improved easy*
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void invertchol(
        matrix<real4> &A)
  {
  const int32 N = A.lines();
  register real4 sum;
  register int32 i,j,k;
// ______ Compute inv(L) store in lower of A ______
  for (i=0; i<N; ++i)
    {
    A[i][i] = 1./A[i][i];
    for (j=i+1; j<N; ++j)
      {
      sum = 0.;
      for (k=i; k<j; ++k)
        {
        sum -= A[j][k]*A[k][i];
        }
      A[j][i]=sum/A[j][j];
      }
    }
// ______ Compute inv(A)=inv(LtL) store in lower of A ______
  for (i=0; i<N; ++i)
    {
    for (j=i; j<N; ++j)
      {
      sum = 0.;
      for (k=j; k<N; ++k)
        {
        sum += A[k][i]*A[k][j];                 // transpose
        }
      A[j][i] = sum;
      }
    }
  } // END invertchol BK
#endif


#ifdef __USE_LAPACK_LIBRARY__
/****************************************************************
 * invertchol(A); inversion by cholesky factorisation veclib    *
 * lower triangle of A is changed on output                     *
 * upper reamins un referenced                                  *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
void invertchol(
        matrix<real8> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("invertchol: LAPACK");
#endif
#ifdef __DEBUGMAT1
  if (A.lines() != A.pixels())
    matERROR.print("only symmetrical matrixes for cholesky.");
#endif

  int32 lda= A.pixels();                // leading dimension
  int32 n  = lda;                       // order of A
  int32 ierr;                           // 
  dpotri("U",&n,A[0],&lda,&ierr,1);     // Upper remains unchanged

#ifdef __DEBUGMAT1
  if      (ierr == 0)
    matDEBUG.print("invertcholes: ok, lower contains inverse.");
  else if (ierr <  0)
    matERROR.print("invertcholes: the kth element is not ok.");
  else 
    matERROR.print("invertcholes: singular.");
#endif
  } // END invertchol

// ______ Else use simple (inaccurate? slow?) internal algorithm ______
// #elif defined (INTERNAL)
#else
/****************************************************************
 * invertchol(A); inversion by cholesky factorisation internal  *
 * lower triangle of A is changed on output                     *
 * upper remains un referenced                                  *
 * implementation with double loop probably can be improved easy*
 *    Bert Kampes, 11-Oct-1999                                  *
 ****************************************************************/
void invertchol(
        matrix<real8> &A)
  {
  const int32 N = A.lines();
  register real8 sum;
  register int32 i,j,k;
// ______ Compute inv(L) store in lower of A ______
  for (i=0; i<N; ++i)
    {
    A[i][i] = 1./A[i][i];
    for (j=i+1; j<N; ++j)
      {
      sum = 0.;
      for (k=i; k<j; ++k)
        {
        sum -= A[j][k]*A[k][i];
        }
      A[j][i]=sum/A[j][j];
      }
    }
// ______ Compute inv(A)=inv(LtL) store in lower of A ______
  for (i=0; i<N; ++i)
    {
    for (j=i; j<N; ++j)
      {
      sum = 0.;
      for (k=j; k<N; ++k)
        {
        sum += A[k][i]*A[k][j];                 // transpose
        }
      A[j][i] = sum;
      }
    }
  } // END invertchol BK
#endif



/****************************************************************
 *  A=norm(B)                                                   *
 *    Bert Kampes, 08-Feb-1999                                  *
 ****************************************************************/
matrix<complr4> norm(
        const matrix<complr4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("norm");
#endif
// ______Ensure stride one memory______
  matrix<complr4> Result(A.lines(),A.pixels());
  complr4 *pntA = A[0];
  complr4 *pntR = Result[0];
  //for (register int32 i=0; i<Result.size(); i++)
  //  (*pntR++) = complr4(norm(*pntA++));
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<Result.size(); i++)
    {
      *(pntR++) = complr4(norm(*pntA));
        pntA++;                                         // [MA] see fast_dotmultconjphase 
    }
  return Result;
  } // END norm



/****************************************************************
 *    norm2                                                      *
 *    Bert Kampes, 08-Feb-1999                                  *
 ****************************************************************/
real4 norm2(
        const matrix<complr4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("norm2");
#endif
// ______Ensure stride one memory______
  real4 n=0.;
  complr4 *pntA = A[0];
  //for (register int32 i=0; i<A.size(); i++)
  //  n += norm(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
      n += norm(*pntA);
      pntA++;
    }
  return n;
  } // END norm2



/****************************************************************
 *    norm2                                                      *
 *    Bert Kampes, 08-Feb-1999                                  *
 ****************************************************************/
real4 norm2(
        const matrix<real4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("norm2");
#endif
// ______Ensure stride one memory______
  real4 n=0.;
  real4 *pntA = A[0];
  //for (register int32 i=0; i<A.size(); i++)
  //  n += sqr(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
      n += sqr(*pntA);
      pntA++;
    }
  return n;
  } // END norm2



/****************************************************************
 * A = abs(B)                                                   *
 * Bert Kampes, 16-Feb-1999                                     *
 ****************************************************************/
matrix<real4> abs(
        const matrix<real4> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("abs");                                          // TODO define <type>
#endif
  matrix<real4> Result(A.lines(),A.pixels());
// ______Ensure stride one memory______
  real4 *pntA   = A[0];
  real4 *pntR   = Result[0];
  //for (register int32 i=0; i<Result.size(); i++)
  // (*pntR++) = (abs((*pntA++)));
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<Result.size(); i++)
    {
    *(pntR++) = (abs((*pntA)));       
      pntA++;
    }
  return Result;
  } // END abs



/****************************************************************
 * A = abs(B)                                                   *
 * Bert Kampes, 02-Jun-1999                                     *
 ****************************************************************/
matrix<real8> abs(
        const matrix<real8> &A)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("abs");
#endif
  matrix<real8> Result(A.lines(),A.pixels());
// ______Ensure stride one memory______
  real8 *pntA   = A[0];
  real8 *pntR   = Result[0];
  //for (register int32 i=0; i<Result.size(); i++)
  //  (*pntR++) = (abs((*pntA++)));       
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<Result.size(); i++)
    {
    *(pntR++) = (abs((*pntA)));       
      pntA++;
    }
  return Result;
  } // END abs



/****************************************************************
 * A=real(B)                                                    *
 *    Bert Kampes, 02-Mar-1999                                  *
 ****************************************************************/
matrix<real4> real(
        const matrix<complr4> &A)
  {
  matrix<real4>  Result(A.lines(),A.pixels());
  real4         *pntR = Result[0];
  complr4       *pntA = A[0];
#ifdef __DEBUGMAT2
  matDEBUG.print("real. do not initialize memory (no alloc&and init).");
#endif
// ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = real(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = real(*pntA);
      pntA++;
    }
  return Result;
  } // END real



/****************************************************************
 * A=imag(B)                                                    *
 *    Bert Kampes, 02-Mar-1999                                  *
 ****************************************************************/
matrix<real4> imag(
        const matrix<complr4> &A)
  {
  matrix<real4> Result(A.lines(),A.pixels());
  real4         *pntR = Result[0];
  complr4       *pntA = A[0];
#ifdef __DEBUGMAT2
  matDEBUG.print("imag. do not initialize memory (no alloc).");
#endif
// ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = imag(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = imag(*pntA);
    pntA++;
    }
  return Result;
  } // END imag



/****************************************************************
 * A=angle(B)                                                   *
 *    Bert Kampes, 13-Apr-1999                                  *
 ****************************************************************/
matrix<real4> angle(
        const matrix<complr4> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("angle(complex) Bert Kampes, 13-Apr-1999.");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());// faster not to initialize
  real4         *pntR = Result[0];
  complr4       *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = arg(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = arg(*pntA);
      pntA++;
    }
  return Result;
  } // END angle


/****************************************************************
 * A=fast_angle(B)                                              *
 *    Bert Kampes, 06-Oct-2005                                  *
 ****************************************************************/
matrix<real4> fast_angle(
        const matrix<complr4> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("angle(complex) Bert Kampes, 06-Oct-2005.");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());// faster not to initialize
  real4         *pntR = Result[0];
  complr4       *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = fast_arg(*pntA++);// LUT
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = fast_arg(*pntA);// LUT
      pntA++;
    }
  return Result;
  } // END fast_angle



/****************************************************************
 * A=angle2cmplx(B)                                             *
 *    Bert Kampes, 06-Oct-1999                                  *
 ****************************************************************/
matrix<complr4> angle2cmplx(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("angle2cmplx() Bert Kampes, 06-Oct-1999");
  #endif
  matrix<complr4> Result(A.lines(),A.pixels());// faster not initialize matrix
  complr4       *pntR = Result[0];
  real4         *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = complr4(cos(*pntA),sin(*pntA++));
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = complr4(cos(*pntA),sin(*pntA));
     pntA++;
    }
  return Result;
  } // END angle2cmplx


/****************************************************************
 * A=fast_angle2cmplx(B)                                        *
 *    Bert Kampes, 06-Oct-2005                                  *
 ****************************************************************/
matrix<complr4> fast_angle2cmplx(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("fast_angle2cmplx() Bert Kampes, 06-Oct-2005");
  #endif
  matrix<complr4> Result(A.lines(),A.pixels()); // faster not initialize matrix
  complr4       *pntR = Result[0];
  real4         *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = complr4(fast_cos(*pntA),fast_sin(*pntA++));// LUT
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = complr4(fast_cos(*pntA),fast_sin(*pntA));// LUT
      pntA++;
    }
  return Result;
  } // END fast_angle2cmplx


/****************************************************************
 * A=fast_angle2cmplx(B)                                        *
 *    Bert Kampes, 06-Oct-2005                                  *
 *    Mahmut Arikan 2009
 ****************************************************************/
matrix<complr8> fast_angle2cmplx(
        const matrix<real8> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("fast_angle2cmplx() Bert Kampes, 06-Oct-2005");
  #endif
  matrix<complr8> Result(A.lines(),A.pixels()); // faster not initialize matrix
  complr8       *pntR = Result[0];
  real8         *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register int32 i=0; i<A.size(); i++)
  //  (*pntR++) = complr4(fast_cos(*pntA),fast_sin(*pntA++));// LUT
  // changed by FvL (for g++/gcc > 4.0):
  for (register int32 i=0; i<A.size(); i++)
    {
    *(pntR++) = complr8(fast_cos(*pntA),fast_sin(*pntA));// LUT
      pntA++;                                         // [MA] see fast_dotmultconjphase
    }
  return Result;
  } // END fast_angle2cmplx

/****************************************************************
 * dotmultconjphase(B,A)                                        *
 * subtract phase in real4 matrix A of complexr4 matrix B       *
 * by reference. (B .* conj(complr4(cos(phi,sin(phi))))         *
 *    Bert Kampes, 06-Oct-1999                                  *
 ****************************************************************/
void dotmultconjphase(
        matrix<complr4>     &cint,
        const matrix<real4> &refpha)
  {
  complr4       *pntC = cint[0];
  real4         *pntA = refpha[0];
#ifdef __DEBUGMAT2
  matDEBUG.print("dotmultconjphase.");
#endif
#if  defined (__DEBUGMAT1) || defined (__DEBUGMAT2)
  if (cint.lines() != refpha.lines())
    matERROR.print("dotmultconjphase: not same number of lines.");
  if (cint.pixels() != refpha.pixels())
    matERROR.print("dotmultconjphase: not same number of pixels.");
#endif
  // ______ Ensure stride one memory ______
  // --- looks a bit fishy: is pntA first increased from right to left,
  // --- or as i thought after sin is taken?
  // --- better would be to increase pntA in the for, together with i++
  // --- if we ever see a problem, it may be due to this?
  //for (register int32 i=0; i<cint.size(); i++)
  //  (*pntC++) *= complr4(cos(*pntA),-sin(*pntA++));
  //
  // Yes, indeed this gave a problem, for gcc 4.0 on a 64 bit
  // computer. Therefore, changed the code to this: (FvL)
  for (register int32 i=0; i<cint.size(); i++)
    {
      *(pntC++) *= complr4(cos(*pntA),-sin(*pntA));
        pntA++;
    }
  } // END dotmultconjphase


/****************************************************************
 * fast_dotmultconjphase(B,A)                                   *
 * subtract phase in real4 matrix A of complexr4 matrix B       *
 * by reference. (B .* conj(complr4(cos(phi,sin(phi))))         *
 *    Bert Kampes, 06-Oct-2005                                  *
 ****************************************************************/
void fast_dotmultconjphase(
        matrix<complr4>     &cint,
        const matrix<real4> &refpha)
  {
  complr4       *pntC = cint[0];
  real4         *pntA = refpha[0];
#ifdef __DEBUGMAT2
  matDEBUG.print("fast_dotmultconjphase.");
#endif
#if   defined (__DEBUGMAT1) || defined (__DEBUGMAT2)
  if (cint.lines() != refpha.lines())
    matERROR.print("dotmultconjphase: not same number of lines.");
  if (cint.pixels() != refpha.pixels())
    matERROR.print("dotmultconjphase: not same number of pixels.");
#endif
  // ______ Ensure stride one memory ______
  // --- looks a bit fishy: is pntA first increased from right to left,
  // --- or as i thought after sin is taken?
  // --- better would be to increase pntA in the for, together with i++
  // --- if we ever see a problem, it may be due to this?
  //for (register int32 i=0; i<cint.size(); i++)
  //  (*pntC++) *= complr4(fast_cos(*pntA),fast_min_sin(*pntA++));
  //
  // Yes, indeed this gave a problem, for gcc 4.0 on a 64 bit
  // computer. Therefore, changed the code to this: (FvL)
  for (register int32 i=0; i<cint.size(); i++)
    {
      *(pntC++) *= complr4(fast_cos(*pntA),fast_min_sin(*pntA));
       pntA++;                                                      // [MA] switched from *pntA++ --> pntA++ previously: 1. addr++ then 2. * is evaluated
    }
  } // END fast_dotmultconjphase


/****************************************************************  // see template in matrix.cc
 * fast_dotmultconjphase(B,A)                                   *
 * subtract phase in real4 matrix A of complexr4 matrix B       *
 * by reference. (B .* conj(complr4(cos(phi,sin(phi))))         *
 *    Bert Kampes, 06-Oct-2005                                  *
 *    Mahmut Arikan 2009 real8
 ****************************************************************/
void fast_dotmultconjphase(
        matrix<complr4>     &cint,
        const matrix<real8> &refpha)
  {
  complr4       *pntC = cint[0];
  real8         *pntA = refpha[0];
#ifdef __DEBUGMAT2
  matDEBUG.print("fast_dotmultconjphase.");
#endif
#if   defined (__DEBUGMAT1) || defined (__DEBUGMAT2)
  if (cint.lines() != refpha.lines())
    matERROR.print("dotmultconjphase: not same number of lines.");
  if (cint.pixels() != refpha.pixels())
    matERROR.print("dotmultconjphase: not same number of pixels.");
#endif
  // ______ Ensure stride one memory ______
  // --- looks a bit fishy: is pntA first increased from right to left,
  // --- or as i thought after sin is taken?
  // --- better would be to increase pntA in the for, together with i++
  // --- if we ever see a problem, it may be due to this?
  //for (register int32 i=0; i<cint.size(); i++)
  //  (*pntC++) *= complr4(fast_cos(*pntA),fast_min_sin(*pntA++));
  //
  // Yes, indeed this gave a problem, for gcc 4.0 on a 64 bit
  // computer. Therefore, changed the code to this: (FvL)
  for (register int32 i=0; i<cint.size(); i++)
    {
     *(pntC++) *= complr4(fast_cos(*pntA),fast_min_sin(*pntA));
       pntA++;                                                  // [MA] switched from *pntA++ --> pntA++ previously: 1. addr++ then 2. * is evaluated
    }
  } // END fast_dotmultconjphase

/****************************************************************
 * A=coherence(complexinterferogram,norms,estwindowsizel,p)     *
 * returned is matrix size=(orig-winsizel,orig.pixels)          *
 * This routine is time consuming so no subroutines are called. *
 * a fft version is tried but does not speed up                 *
 *    Bert Kampes, 19-Apr-1999                                  *
 ****************************************************************/
matrix<complr4> coherence(
        const matrix<complr4> &CINT, 
        const matrix<complr4> &NORMS,
        uint winL,
        uint winP)
  {
// #define PLACE
// #define FFTS
#define PLACE2                  // fastest


#ifdef __DEBUGMAT2
  matDEBUG.print("coherence");
  if (!winL>=winP)
    matDEBUG.print("coherence: estimator window size L<P not very efficiently programmed.");
#endif
#ifdef __DEBUGMAT1
  if (CINT.lines() != NORMS.lines() ||
      CINT.pixels() != NORMS.pixels() )
      matERROR.print("not same dimensions.");
#endif


  // ___ Allocate output matrix ___
  matrix<complr4> Result(CINT.lines()-winL+1,CINT.pixels());
  int32 leadingzeros  = (winP-1)/2;             // number of pixels=0 floor...
  int32 trailingzeros = (winP)/2;               // floor...

// --- METHOD 1: --------------------------------------------- //
#ifdef PLACE
  complr4 sum   = 0.;
  complr4 power = 0.;
  for (register int32 i=0; i<Result.lines(); i++)
    {
    for (register int32 j=leadingzeros; j<Result.pixels()-trailingzeros; j++)
      {
      for (register int32 k=i; k<i+winL; k++)
        {
        for (register int32 l=j-leadingzeros; l<j-leadingzeros+winP; l++)
          {
          sum   += CINT(k,l);
          power += NORMS(k,l);
          }
        }
      real4 p     = power.real() * power.imag();
      Result(i,j) = (p > 0.0) ?  sum/sqrt(p) : 0.0;
      sum         = 0.;
      power       = 0.;
      }
    }
#undef PLACE
#endif
// --- end METHOD 1: ----------------------------------------- //



// --- METHOD 2: --------------------------------------------- //
#ifdef PLACE2
  register int32   i,j,k,l;
  register complr4 sum;
  register complr4 power;
  for (j=leadingzeros; j<Result.pixels()-trailingzeros; j++)
    {
    sum   = complr4(0.);
    power = complr4(0.);
    // ______ Compute sum over first block ______
    for (k=0; k<winL; k++)
      {
      for (l=j-leadingzeros; l<j-leadingzeros+winP; l++)
        {
        sum   += CINT(k,l);
        power += NORMS(k,l);
        }
      }
    real4 p     = power.real() * power.imag();
    Result(0,j) = (p > 0.0) ?  sum/sqrt(p) : 0.0;

    // ______ Compute sum over rest of blocks ______
    for (i=0; i<Result.lines()-1; i++)
      {
      for (l=j-leadingzeros; l<j-leadingzeros+winP; l++)
        {
        sum   += (CINT(i+winL,l)  - CINT(i,l));
        power += (NORMS(i+winL,l) - NORMS(i,l));
        }
      real4 p = power.real() * power.imag();
      Result(i+1,j) = (p > 0.0) ?  sum/sqrt(p) : 0.0;
      }
    }
#undef PLACE2
#endif
// --- end METHOD 2: ----------------------------------------- //



// --- METHOD 3: --------------------------------------------- //
// ====== Compute coherence by fft's ======
#ifdef FFTS
  matDEBUG.print("THIS SEEMS TO GO WRONG SOMEWHERE, writes all zeros, computes ok?");
  bool  newline   = false;
  bool  nonewline = false;

  int32 SIZEBLOCKL = 256;                       // power2, must be larger than winL
  int32 SIZEBLOCKP = 1024;                      // power2, must be larger than winP
  int32 sizeL = SIZEBLOCKL;
  int32 sizeP = SIZEBLOCKP;
  if (CINT.lines() < SIZEBLOCKL)
    {
    SIZEBLOCKL = nextpow2(CINT.lines());
    sizeL = CINT.lines();
    nonewline = true;
    }
  if (CINT.pixels() < SIZEBLOCKP)
    {
    SIZEBLOCKP = nextpow2(CINT.pixels());
    sizeP = CINT.pixels();
    }

  const int32 twoL = 2*SIZEBLOCKL;
  const int32 twoP = 2*SIZEBLOCKP;
  const complr4 cr4one = complr4(1.,0);

  matrix<complr4> CINT2(twoL, twoP);
  matrix<complr4> NORMS2(twoL, twoP);
  matrix<complr4> BLOCK(twoL, twoP);
  register int32 i,j;
  for (i=0; i<winL; i++)
    for (j=0; j<winP; j++)
      BLOCK(i,j) = cr4one;
  fft2d(BLOCK);
  BLOCK = conj(BLOCK);                  // you should use theoretical expression.

  window win  = {0, sizeL-1, 0, sizeP-1};       // part from CINT,NORM
  window win2 = {0, sizeL-1, 0, sizeP-1};       // place where to put win
  register int32 i2 = 0;
  register int32 j2 = leadingzeros;

  for (;;) // ever
    {
#ifdef __DEBUGMAT2
  cout << "time for coherence init:\n";
  printcpu();
#endif

cout << "BERT: win:" 
     << win.linelo << " "
     << win.linehi << " "
     << win.pixlo << " "
     << win.pixhi << "\n";
cout << "BERT: win2:" 
     << win2.linelo << " "
     << win2.linehi << " "
     << win2.pixlo << " "
     << win2.pixhi << "\n";
cout << "BERT: i2,j2 sizeL,p:" 
     << i2 << " " 
     << j2 << " " 
     << sizeL << " " 
     << sizeP << "\n";


    CINT2.setdata(win2,CINT,win); 
    NORMS2.setdata(win2,NORMS,win); 
#ifdef __DEBUGMAT2
  cout << "time for setdata:\n";
  printcpu();
#endif
    fft2d(CINT2);
    fft2d(NORMS2);
#ifdef __DEBUGMAT2
  cout << "time for fft:\n";
  printcpu();
#endif

// ______ Convolution ______
    for (i=0; i<CINT2.lines(); i++)
      {
      for (j=0; j<CINT2.pixels(); j++)
        {
        CINT2(i,j)  *= BLOCK(i,j);
        NORMS2(i,j) *= BLOCK(i,j);
        }
      }

#ifdef __DEBUGMAT2
  cout << "time for convol.:\n";
  printcpu();
#endif
    ifft2d(CINT2);
    ifft2d(NORMS2);
#ifdef __DEBUGMAT2
  cout << "time for ifft.:\n";
  printcpu();
#endif

// ______ Coherence ______
    for (i=0; i<=sizeL-winL; i++)
      {
      for (j=0; j<=sizeP-winP; j++)
        {
        real4 p = NORMS2(i,j).real() * NORMS2(i,j).imag();
        Result(i2+i,j2+j) = (p>0.0) ? CINT2(i,j)/sqrt(p) : 0.0;
        }
      }

#ifdef __DEBUGMAT2
  cout << "time for coherence.:\n";
  printcpu();
#endif

// ______ Update window ______
    j2 += (SIZEBLOCKP - winP + 1);
    win.pixlo = win.pixhi + 1 - leadingzeros - trailingzeros;
    win.pixhi = win.pixlo + sizeP - 1;
    if (win.pixhi > CINT.pixels() - 1)
      {
cout << "BERT: win.pixhi > CINT.pixels() - 1\n";
      if (win.pixlo+leadingzeros > Result.pixels()-trailingzeros-1)
        newline = true;
      if (newline)
        {
cout << "BERT: newline==true\n";
        if (nonewline)
          break;                                        // break loop
        newline    = false;
        sizeP = SIZEBLOCKP;
        win.linelo = win.linehi + 1 - (winL-1)/2; 
        win.linehi = win.linelo + sizeL - 1;
        win.pixlo  = 0;
        win.pixhi  = win.pixlo + sizeP - 1;
        i2   += (SIZEBLOCKL - winL + 1);
        j2    = leadingzeros;
        if (win.linehi > CINT.lines() - 1)
          {
          nonewline  = true;                            // this is last block
          win.linehi = CINT.lines() - 1;                // no resizing for now
          sizeL      = win.linehi - win.linelo + 1;
          }
        }

      else      // last block of line
        {
cout << "BERT: newline==false\n";
        newline   = true;
        win.pixhi = CINT.pixels() - 1;                  // no resizing for now
        sizeP     = win.pixhi - win.pixlo + 1;
        }

      win2.linehi = sizeL - 1;
      win2.pixhi  = sizeP - 1;
      }
    }
#undef FFTS
#endif
// --- end METHOD 2: ----------------------------------------- //


  //cout << "Mean coherence: " << mean(Result); 
  return Result;
  } // END coherence



/****************************************************************
 * A=coherence(complexinterferogram,norms,estwindowsizel,p)     *
 * returned is matrix size=(orig-winsizel,orig.pixels)          *
 * real4 coherence                                              *
 * This routine is time consuming so no subroutines are called. *
 *    Bert Kampes, 19-Apr-1999                                  *
 ****************************************************************/
matrix<real4> coherence2(
        const matrix<complr4> &CINT, 
        const matrix<complr4> &NORMS,
        uint winL,
        uint winP)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("coherence2");
  if (!winL>=winP)
    matDEBUG.print("coherence: estimator window size L<P not very efficiently programmed.");
#endif
#ifdef __DEBUGMAT1
  if (CINT.lines() != NORMS.lines() ||
      CINT.pixels() != NORMS.pixels() )
      matERROR.print("coherence2::not same dimensions.");
#endif

  matrix<real4> Result(CINT.lines()-winL+1,CINT.pixels());
  int32 leadingzeros  = (winP-1)/2;             // number of pixels=0 floor...
  int32 trailingzeros = (winP)/2;               // floor...

  register int32   i,j,k,l;
  register complr4 sum;
  register complr4 power;
  for (j=leadingzeros; j<Result.pixels()-trailingzeros; j++)
    {
    sum   = complr4(0.);
    power = complr4(0.);

    // ______ Compute sum over first block ______
    for (k=0; k<winL; k++)
      {
      for (l=j-leadingzeros; l<j-leadingzeros+winP; l++)
        {
        sum   += CINT(k,l);
        power += NORMS(k,l);
        }
      }
    //Result(0,j) = norm(sum) / (power.real() * power.imag());
    real4 p = power.real() * power.imag();
    Result(0,j) = (p>0.0) ? sqrt(norm(sum)/p) : 0.0;

    // ______ Compute (relatively) sum over rest of blocks ______
    for (i=0; i<Result.lines()-1; i++)
      {
      for (l=j-leadingzeros; l<j-leadingzeros+winP; l++)
        {
        sum   += (CINT(i+winL,l)  - CINT(i,l));
        power += (NORMS(i+winL,l) - NORMS(i,l));
        }
      real4 p = power.real() * power.imag();
      Result(i+1,j) = (p>0.0) ? sqrt(norm(sum)/p) : 0.0;
      }
    }
  return Result;
  } // END coherence2




/****************************************************************
 * mysort2(A)                                                   *
 * sorts matrix first on column1, then on col2                  *
 * used to sort l,p,value matrices.                             *
 * calls to std c lib qsort, should better be c++ sort function?*
 * sorts in ascending order.                                    *
 * Bert Kampes, 10-Jan-2000                                     *
 ****************************************************************/
int32 mycomp2 (const void *x1, const void *x2)
  {
  //  float *pntf = (float*) x1;
  //  cout << " BB : " << *pntf << " " << *(pntf+1) << endl;
  //  cout << " BB : " << *(float*)x1 << " " << *((float*)x1+1) << endl;
  int ret = 0;
  if      (*(float*)x1 < *(float*)x2 ) ret = -1; // first field smaller
  else if (*(float*)x1 > *(float*)x2 ) ret =  1;
  else if (*((float*)x1+1) < *((float*)x2+1) ) ret = -1; // x1=x2
  else if (*((float*)x1+1) > *((float*)x2+1) ) ret =  1; // x1=x2
  // else ret=0; // same x1,x2 and x1+1,x2+1
  return ret;
  }
int32 mycomp231 (const void *x1, const void *x2)       // [MA] compare last 2 cols for getmodeoffset
  {
  //  float *pntf = (float*) x1;
  //  cout << " BB : " << *pntf << " " << *(pntf+1) << endl;
  //  cout << " BB : " << *(float*)x1 << " " << *((float*)x1+1) << endl;
  int ret = 0;
  if      (*((float*)x1+1) < *((float*)x2+1) ) ret = -1; // second field smaller
  else if (*((float*)x1+1) > *((float*)x2+1) ) ret =  1; // 
  else if (*((float*)x1+2) < *((float*)x2+2) ) ret = -1; // x2=x3 // third field smaller
  else if (*((float*)x1+2) > *((float*)x2+2) ) ret =  1; // x2=x3
  else if (*(float*)x1     < *(float*)x2 )     ret = -1; // first field smaller
  else if (*(float*)x1     > *(float*)x2 )     ret =  1;
  // else ret=0; // same x1,x2 and x1+1,x2+1 and x1+2,x2+2
  // else ret=0; // same x1,x2 and x1+1,x2+1
  return ret;
  }
//  /****************************************************************
void mysort2 (matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mysort2 (real4)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() < 2)
    matERROR.print("mysort sorts only min. 2 cols.");
  #endif
  qsort(&A[0][0],A.lines(),A.pixels()*sizeof(real4),mycomp2);
  } // END mysort2
//  /****************************************************************
void mysort2 (matrix<int32> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mysort2 (int)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() < 2)
    matERROR.print("mysort sorts only min. 2 cols.");
  #endif
  qsort(&A[0][0],A.lines(),A.pixels()*sizeof(int32),mycomp2);
  } // END mysort2
//  /****************************************************************
void mysort231 (matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mysort23 (real4) [MA]");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() < 3)
    matERROR.print("mysort sorts on  2nd, 3rd and 1st cols respectively.");
  #endif
  qsort(&A[0][0],A.lines(),A.pixels()*sizeof(real4),mycomp231);
  } // END mysort23


/****************************************************************
 * mysort2selcol(A,selcol)                                      *
 * first col=0, second col=1 and so on...                       *
 * sorts matrix on the selected col                             *
 * used to sort l,p,value matrices.                             *
 * Sorts using simple BubleSort, not very efficient: O(N^2)     *
 * sorts in ascending order.                                    *
 * Bert Kampes, 10-Jan-2000 (original mysort2)                  *
 * Batuhan Osmanoglu, Aug 2007                                  *
 ****************************************************************/
void mysort2selcol (matrix<real4> &A, int32 selcol)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mysort2selcol (real4)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() < 2)
    matERROR.print("mysort2selcol sorts only min. 2 cols.");
  #endif
  //qsort(&A[0][0],A.lines()+col*sizeof(real4),A.pixels()*sizeof(real4),mycomp2);
  for (int32 i=0; i<A.lines()-1; i++)
    {
    for (int32 r=0; r<A.lines()-1; r++)
      {
     // for (int32 c=0; c<width; c++)
     //   {
        real4 temp;
        if (A(r,selcol) > A(r+1,selcol))
          {
          //DEBUG << "Swapping Row " << r << " with " <<  r+1;
          //DEBUG.print();
          for (int32 c=0; c<A.pixels(); c++)
            {
            temp = A(r,c);
            A(r,c)= A(r+1,c);
            A(r+1,c) = temp;
            }
          }
      //  }
      }
    }
    DEBUG << " SORTED MATRIX BY COLUMN: " << selcol;
    DEBUG.print();
    for (int32 i=0; i< A.lines(); i++)
      {
        DEBUG << i << ", " << A(i,0) << ", " << A(i,1) << ", " << A(i,2);
        DEBUG.print();
      }
  } // END mysort2selcol
//  /****************************************************************
void mysort2selcol (matrix<int32> &A, int32 selcol)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mysort2 (int)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() < 2)
    matERROR.print("mysort sorts only min. 2 cols.");
  #endif
  //qsort(&A[0][0],A.lines(),A.pixels()*sizeof(int32),mycomp2);
  for (int32 i=0; i<A.lines()-1; i++)
    {
    for (int32 r=0; r<A.lines()-1; r++)
      {
     // for (int32 c=0; c<width; c++)
     //   {
        int32 temp;
        if (A(r,selcol) > A(r+1,selcol))
          {
          //DEBUG << "Swapping Row " << r << " with " <<  r+1;
          //DEBUG.print();
          for (int32 c=0; c<A.pixels(); c++)
            {
            temp = A(r,c);
            A(r,c)= A(r+1,c);
            A(r+1,c) = temp;
            }
          }
      //  }
      }
    }
    DEBUG << " SORTED MATRIX BY COLUMN: " << selcol;
    DEBUG.print();
    for (int32 i=0; i< A.lines(); i++)
      {
        DEBUG << i << ", " << A(i,0) << ", " << A(i,1) << ", " << A(i,2);
        DEBUG.print();
      }
  } // END mysort2selcol


//  /****************************************************************
//   * mysort(A)                                                    *
//   * sorts matrix on column col                                       *
//   * calls to std c lib qsort, should better be c++ sort function?*
//   * sorts in ascending order.                                    *
//   * Bert Kampes, 10-Jan-2000                                     *
//   ****************************************************************/
//  int32 mycomp1 (const void *x1, const void *x2)
//    {
//  //  float *pntf = (float*) x1;
//  //  cout << " BB : " << *pntf << " " << *(pntf+1) << endl;
//    cout << " BB : " << *(float*)x1 << " " << *((float*)x1+1) << endl;
//   
//    int ret = 0;
//    if      (*(float*)x1 < *(float*)x2 ) ret = -1; // first field smaller
//    else if (*(float*)x1 > *(float*)x2 ) ret =  1;
//    else if (*((float*)x1+columnnumber...) < *((float*)x2+columnnumber...) ) ret = -1; // x1=x2
//    else if (*((float*)x1+columnnumber...) > *((float*)x2+columnnumber...) ) ret =  1; // x1=x2
//    // else ret=0; // same x1,x2 and x1+1,x2+1
//    return ret;
//    }
//  void mysort1 (matrix<real4> &A, column)
//    {
//  #ifdef __DEBUGMAT2
//    matDEBUG.print("mysort2.");
//  #endif
//  #ifdef __DEBUGMAT1
//    if (A.pixels() < 2)
//      matERROR.print("mysort sorts only min. 2 cols.");
//  #endif
//    qsort(&A[0][0],A.lines(),A.pixels()*sizeof(real4),mycompare);
//    } // END mysort1





/****************************************************************
 * four1(complr4 *, length, isign)                              *
 * four1(&A[0][0], 128, 1)                                      *
 * either based on numerical recipes or veclib                  *
 * helper function for other fft routines, if no veclib         *
 * cooley-turkey, power 2, replaces input,                      *
 * isign=1: fft , isign=-1: ifft                                *
 * handling vectors should be simpler (lying, standing)         *
 * note that this is not a good implementation, only to get     *
 * doris software working without veclib.                       *
 *                                                              *
 * define SAMEASVECIB if you want the order of the coefficients *
 * of the fft the same as veclib. it seems this is not required *
 * for a good version of Doris, but in case of problems this    *
 * may be the solution.                                         *
 *                                                              *
 * VECLIB defines the FT same as matlab:                        * 
 *          N-1                                                 *
 *   X(k) = sum  x(n)*exp(-j*2*pi*k*n/N), 0 <= k <= N-1.        *
 *          n=0                                                 *
 *                                                              *
 * FFTW defines the same as Matlab, but inv. not normalized.    *
 * I don't know if the matrix must be allocated somehow, so for *
 * now we try only 1d ffts to build 2d too.                     *
 #%// BK 17-Sep-2003                                            *
 *                                                              *
 *    Bert Kampes, 12-Oct-1999                                  *
 ****************************************************************/
#ifndef __USE_VECLIB_LIBRARY__
#ifndef __USE_FFTW_LIBRARY__
#define SWAP(a,b) bert=(a);(a)=(b);(b)=bert
#endif
#endif
// --- common part to four1 for veclib,fftw, and internal ---
void four1(
        complr4 data[],
        int32 fftlength,
        int32 isign)
  {
  // ______ Repair original routine ______
  //#ifdef __DEBUGMAT2 // called very often...
    //matDEBUG.print("four1: fft internal");
  //#endif
  #ifdef __DEBUGMAT1
  if (!myispower2(fftlength))
    matERROR.print("four1: length fft must be power of 2");// not really...
  if (abs(isign) != 1)
    matERROR.print("four1: isign should be (-)1.");
  #endif

// ====== FFTW LIBRARY ROUTINE ======
#ifdef __USE_FFTW_LIBRARY__
  // if no plan yet, create it depending on direction.
  // use guru interface for created plan on different data.
  // first check if dimension of old plan is still valid, or that the size of data changed,
  // and we thus have to create a new plan?
  // FFTW_ESTIMATE does not overwrite data during planning.
  // the plan can be stored in matrixclass as I do now, or as a static here...
  // that would mean only one plan, and no changes to matrixbk.h, but...
// #define __NO_GURU_YET__
#ifdef __NO_GURU_YET__
  fftwf_plan fftwf_tmp_plan  = (isign==1) ?
    fftwf_plan_dft_1d(fftlength,
        reinterpret_cast<fftwf_complex*>(data), 
        reinterpret_cast<fftwf_complex*>(data), 
        FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT) :
    fftwf_plan_dft_1d(fftlength, r
        reinterpret_cast<fftwf_complex*>(data), 
        reinterpret_cast<fftwf_complex*>(data), 
        FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT) :
  fftwf_execute(fftwf_tmp_plan);
  if (isign==-1)
    for (int32 i=0; i<fftlength; ++i)
      #ifdef __GNUC__
      data[i] /= fftlength;
      #else
      data[i] /= complr4(fftlength);
      #endif
  fftwf_destroy_plan(fftwf_tmp_plan);// deallocs, cannot use destroy when static?


  // === use guru interface to re-use a static plan ===
#else
  if (isign == 1)
    {
    static int32 last_length_fwd = 0;
    static fftwf_plan fftwf_fwd_plan;// for fwd fft
    if (fftlength != last_length_fwd)
      {
      #ifdef __DEBUGMAT2
      matDEBUG.print("BK: FFTW: creating guru plan..."); 
      matDEBUG << "BK: FFTW: fftlength: "  << fftlength
               << "; last_length_fwd: "    << last_length_fwd;
      matDEBUG.print();
      #endif
      last_length_fwd = fftlength;
      fftwf_fwd_plan  = fftwf_plan_dft_1d(fftlength,
        reinterpret_cast<fftwf_complex*>(data), 
        reinterpret_cast<fftwf_complex*>(data), 
        FFTW_FORWARD,  
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      #ifdef __DEBUGMAT2
      fftwf_print_plan(fftwf_fwd_plan);// nerd-readable plan output
      matDEBUG.print(" ");
      #endif
      }
    // ___ Execute plan (guru interface) ___
    #ifdef __DEBUGMAT2
    matDEBUG << "BK: FFTW: executing plan..." << "last_length_fwd: " << last_length_fwd;
    matDEBUG.print();
    #endif
    fftwf_execute_dft(fftwf_fwd_plan, 
      reinterpret_cast<fftwf_complex*>(data),
      reinterpret_cast<fftwf_complex*>(data));
    }

  // ___ Inverse transform ___
  else
    {
    static int32 last_length_inv = 0;
    static fftwf_plan fftwf_inv_plan;// for inverse fft
    if (fftlength != last_length_inv)
      {
      #ifdef __DEBUGMAT2
      matDEBUG.print("BK: FFTW: creating guru plan..."); 
      matDEBUG << "BK: FFTW: fftlength: "  << fftlength
               << "; last_length_inv: "    << last_length_inv;
      matDEBUG.print();
      #endif
      last_length_inv = fftlength;
      fftwf_inv_plan  = fftwf_plan_dft_1d(fftlength,
        reinterpret_cast<fftwf_complex*>(data), 
        reinterpret_cast<fftwf_complex*>(data), 
        FFTW_BACKWARD,  
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
      #ifdef __DEBUGMAT2
      fftwf_print_plan(fftwf_inv_plan);// nerd-readable plan output
      matDEBUG.print(" ");
      #endif
      }
    // ___ Execute plan (guru interface) ___
    #ifdef __DEBUGMAT2
    matDEBUG << "BK: FFTW: executing plan..." << "last_length_inv: " << last_length_inv;
    matDEBUG.print();
    #endif
    fftwf_execute_dft(fftwf_inv_plan, 
      reinterpret_cast<fftwf_complex*>(data),
      reinterpret_cast<fftwf_complex*>(data));
    // ___ Inverse transform, normalization is not performed by fftw ___
    for (int32 i=0; i<fftlength; ++i)
      #ifdef __GNUC__
      data[i] /= fftlength;
      #else
      data[i] /= complr4(fftlength);
      #endif
    }
  }// END FOUR1
#endif // guru interface for an unknown number of transforms
#endif //fftw



  /**************************************************************
   * four1(complr4*, fftlength, iopt)                           *
   *    use veclib c1dfft                                       *
   *    data is complex vector (1,n)                            *
   *    isign=1: fft , isign=-1: ifft                           *
   **************************************************************/
// ______ INTERNAL ROUTINE BASED ON NRC ______
#ifndef __USE_FFTW_LIBRARY__
#ifndef __USE_VECLIB_LIBRARY__
  register int32 i,j,n,mmax,m,istep;
  register complr4 bert;                // for SWAP function
  register double theta;                // trigionometric expansions
  register complr8 w,wp;                // trigionometric expansions

  #define SAMEASVECLIB
  #ifdef SAMEASVECLIB
  // bit slower, but for sure OK, spectral filtering etc. order may be important
  // to get same output as veclib, not sure if it makes a difference in computations,
  // but order of coefficients is defined differently
  // see also swap for isign=1.
  // order veclib: 0, -1/Nd .. 1/Nd, and this algorithm gives 0, 1/Nd .. -1/Nd
  if (isign == -1)
    {for (i=1; i<fftlength/2; ++i)
      {SWAP(data[i],data[fftlength-i]);}}
  #endif

  n = fftlength << 1;
  j = 1;

  // ______ Check if input is standing or lying vector ______
  // index: j -> j/2 (floor) to correct for complex and first element
  // this should be repaired, made better later... within loop counter
  // define order differently
  for (i=1; i<n; i+=2)
    {
    if (j > i)
      {
      SWAP(data[j/2],data[i/2]);
      }
    m = n >> 1;
    while (m>=2 && j>m)
      {
      j  -= m;
      m >>= 1;
      }
    j += m;
    }
  mmax=2;
  while (n > mmax)
    {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wp    = complr8(-2.0*sqr(sin(0.5*theta)),sin(theta));
    w     = complr8(1.0,0.0);
    for (m=1; m<mmax; m+=2)
      {
      for (i=m; i<=n; i+=istep)
        {
        j = i+mmax;
        bert        = data[j/2];
  #ifdef __GNUC__
        bert       *= complr4(real4(w.real()),real4(w.imag())); // STUPID...
  #else                                         // assume aCC hp or so...
        bert       *= w; // not defined for g++ ???
  #endif
        data[j/2]  = data[i/2] - bert;
        data[i/2] += bert;
        }
      w += w * wp;
      }
    mmax = istep;
    }

  #ifdef SAMEASVECLIB
  // to get same output as veclib, not sure if it makes a difference in computations,
  // but order of coefficients is defined differently
  // this swap can be done better... above within loop somehow.
  // coefficients: 0 1 2 3 4 .. fftlength-1
  if (isign == 1)
    {for (i=1; i<fftlength/2; ++i)
      {SWAP(data[i],data[fftlength-i]);}}
  #endif


  // ______ Scale by fftlength if inverse transform ______
  if (isign == -1)
    for (i=0; i<fftlength; ++i)
#ifdef __GNUC__
      data[i] /= fftlength;
#else
      data[i] /= complr4(fftlength);
#endif
#undef SWAP
  }// END FOUR1
#endif //__USE_VECLIB_LIBRARY__
#endif //__USE_FFTW_LIBRARY__



// i prefer fftw for ffts
#ifndef __USE_FFTW_LIBRARY__
#ifdef __USE_VECLIB_LIBRARY__
  /**************************************************************
   * four1(complr4*, fftlength, iopt)                           *
   *    use veclib c1dfft                                       *
   *    data is complex vector (1,n)                            *
   *    isign=1: fft , isign=-1: ifft                           *
   *                                                            *
   *    If you want to include another library, intention is    *
   *    only to interchange this 1d fft routine.                *
   *                                                            *
   *    Bert Kampes, 22-Mar-2000                                *
   **************************************************************/
  int32 ierr = 0;
  static matrix<real4> work;                    // vector best lying
  if (int32(2.5*fftlength) != work.size())      // new size or first time??
    {
    work.resize(1,int32(2.5*fftlength));
    // ______ call initial work: iopt:=-3 ______
    int32 iopt = -3;                    // initialize work
    c1dfft(data,&fftlength,work[0],&iopt,&ierr);
    #ifdef __DEBUGMAT1
    switch (ierr)
      {
      case 0:  matDEBUG.print("c1dfft ok."); break;
      case -1: matERROR.print("length < 0"); break;
      case -2: matERROR.print("l not of required form"); break;
      case -3: matERROR.print("invalid value of iopt"); break;
      case -4: matERROR.print("insufficient dynamical memory available in work"); break;
      default: matERROR.print("four1: unrecognized ierr.");
      }
    #endif
    }
    // ______ Do actual 1d transform ____
  c1dfft(data,&fftlength,work[0],&isign,&ierr);
  if (ierr != 0) {cerr << "veclib: c1dfft: ierr = " << ierr << endl; exit(-1);}
  } // END four1
#endif // either internal or veclib complex 1d fft or fftw
#endif // either internal or veclib complex 1d fft or fftw





/****************************************************************
 * fft(A,dim)                                                   *
 *    forward 1dfft over dim of A is returned in A by reference *
 *    if dim=1 fft is over all columns of A, if 2 over rows.    *
 *    data is stored major row order in memory, so dim=2 is     *
 *    probably much faster.                                     *
 *    fftlength should be power of 2                            *
 *    Bert Kampes, 22-Mar-2000                                  *
 ****************************************************************/
void fft(matrix<complr4> &A, int32 dimension)
  {
  register int32 i;
  const int32 iopt = 1;                         // forward FFT
  switch (dimension)
    {
    case 1:
      {
#ifdef __DEBUGMAT2
  matDEBUG.print("1d fft over columns");
#endif
      const int32 fftlength = A.lines();
      for (i=0; i<A.pixels(); ++i)
        {
        // note that fftw can be used directly to do this w/o data copying...
        matrix<complr4> VECTOR = A.getcolumn(i);// tmp copy, not very efficient maybe
        four1(&VECTOR[0][0],fftlength,iopt);// but generic.
        A.setcolumn(i,VECTOR);
        }
      break;
      }
    case 2:
#ifdef __DEBUGMAT2
  matDEBUG.print("1d fft over rows");
#endif
      {
      const int32 fftlength = A.pixels();
      for (i=0; i<A.lines(); ++i)
        four1(&A[i][0],fftlength,iopt);
      break;
      }
    default:
      matERROR.print("fft: dimension != {1,2}");
    }
  } // END fft




/****************************************************************
 * ifft(A,dim)                                                  *
 *    inverse 1dfft over dim of A is returned in A by reference *
 *    if dim=1 ifft is over all columns of A, if 2 over rows.   *
 *    data is stored major row order in memory, so dim=2 is     *
 *    probably much faster.                                     *
 *    fftlength should be power of 2                            *
 *    Bert Kampes, 22-Mar-2000                                  *
 ****************************************************************/
void ifft(matrix<complr4> &A, int32 dimension)
  {
  register int32 i;
  const int32 iopt = -1;                        // inverse FFT (scaled)
  switch (dimension)
    {
    case 1:
      {
#ifdef __DEBUGMAT2
  matDEBUG.print("1d ifft over columns");
#endif
      const int32 fftlength = A.lines();
      for (i=0; i<A.pixels(); ++i)
        {
        matrix<complr4> VECTOR = A.getcolumn(i);
        four1(&VECTOR[0][0],fftlength,iopt);
        A.setcolumn(i,VECTOR);
        }
      break;
      }
    case 2:
#ifdef __DEBUGMAT2
  matDEBUG.print("1d ifft over rows");
#endif
      {
      const int32 fftlength = A.pixels();
      for (i=0; i<A.lines(); ++i)
        four1(&A[i][0],fftlength,iopt);
      break;
      }
    default:
      matERROR.print("fft: dimension != {1,2}");
    }
  } // END ifft



// common part to fft2d for veclib, fftw, internal.
/****************************************************************
 * fft2d(A);    2dfft veclib                                    *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
void fft2d(
        matrix<complr4> &A)
  {
#ifndef __USE_FFTW_LIBRARY__ // prefer FFTW
#ifdef __USE_VECLIB_LIBRARY__ // use VECLIB
  #ifdef __DEBUGMAT2
    matDEBUG.print("fft2d: veclib");
  #endif
  int32 l2 = A.lines();
  int32 l1 = A.pixels();
  int32 ldz= l1;                        // leading dimension
  int32 iopt = 1;                       // forward transform if >= 0
  int32 ierr;                           // 
  c2dfft(A[0],&l1,&l2,&ldz,&iopt,&ierr);
  #ifdef __DEBUGMAT1
  switch (ierr)
    {
    case 0:  matDEBUG.print("fft2d: ok"); break;
    case -1: matERROR.print("fft2d: l1 not of the required form"); break;
    case -2: matERROR.print("fft2d: l2 not of the required form"); break;
    case -3: matERROR.print("fft2d: ldz not of the required form"); break;
    case -4: matERROR.print("fft2d: probable error in ldz or dimension of z"); break;
    default: matERROR.print("fft2d: unrecognized error.");
    }
  #endif
  } // END fft2d
#endif
#endif



#ifdef __USE_FFTW_LIBRARY__ // prefer FFTW
  /**************************************************************
   * fft2d(A)                                                   *
   *    2dfft internal slow implementation,                     *
   *    four1 should be improved                                *
   *    optionally now fftw for 2d directly #%// BK 17-Sep-2003 *
   *    Bert Kampes, 12-Oct-1999                                *
   **************************************************************/
  // ______ fftw implementation 2d fft ______
  #ifdef __DEBUGMAT2
    matDEBUG.print("fft2d: fftw");
  #endif
  int32 nx = A.pixels();
  int32 ny = A.lines();
  fftwf_plan fftwf_tmp_plan;// plan created each time, stupid.  test.
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: creating tmp fwd 2d plan."); 
  #endif
  fftwf_tmp_plan  = fftwf_plan_dft_2d(
                          nx, ny,
                          reinterpret_cast<fftwf_complex*>(A[0]), 
                          reinterpret_cast<fftwf_complex*>(A[0]), 
                          FFTW_FORWARD,  
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  #ifdef __DEBUGMAT2
  fftwf_print_plan(fftwf_tmp_plan);// nerd-readable plan output
  matDEBUG.print(" ");
  #endif
  // ___ execute plan ___
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: executing tmp fwd 2d plan."); 
  #endif
  fftwf_execute(fftwf_tmp_plan);
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: destroying tmp fwd 2d plan."); 
  #endif
  fftwf_destroy_plan(fftwf_tmp_plan);// also associated data?
  } // END fft2d
#endif // select FFTW


// ====== Internal implementation using series of 1d fft-s ======
#ifndef __USE_FFTW_LIBRARY__
#ifndef __USE_VECLIB_LIBRARY__
  #ifdef __DEBUGMAT2
    matDEBUG.print("fft2d: internal (using multiple 1d ffts)");
  #endif
  fft(A,2);                             // forward transform over rows
  fft(A,1);                             // forward transform over columns
  } // END fft2d
#endif // select internal
#endif // select internal




// common part to ifft2d for veclib, fftw, internal.
/****************************************************************
 * ifft2d(A);   2dfft veclib                                    *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
void ifft2d(
        matrix<complr4> &A)
  {
#ifndef __USE_FFTW_LIBRARY__   // prefer fftw
#ifdef __USE_VECLIB_LIBRARY__  // but use veclib before internal implementation.
  #ifdef __DEBUGMAT2
    matDEBUG.print("ifft2d: veclib");
  #endif
  int32 l2 = A.lines();
  int32 l1 = A.pixels();
  int32 ldz= l1;                        // leading dimension
  int32 iopt = -1;                      // invers transform if < 0
  int32 ierr;                           // 
  c2dfft(A[0],&l1,&l2,&ldz,&iopt,&ierr);
  #ifdef __DEBUGMAT1
  switch (ierr)
    {
    case 0:  matDEBUG.print("ifft2d: ok"); break;
    case -1: matERROR.print("ifft2d: l1 not of the required form"); break;
    case -2: matERROR.print("ifft2d: l2 not of the required form"); break;
    case -3: matERROR.print("ifft2d: ldz not of the required form"); break;
    case -4: matERROR.print("ifft2d: probable error in ldz or dimension of z"); break;
    default: matERROR.print("ifft2d: unrecognized error.");
    }
  #endif
  } // END ifft2d
#endif
#endif



#ifdef __USE_FFTW_LIBRARY__   // prefer fftw
// ====== use internal slow implementation ======
/****************************************************************
 * ifft2d(A);   2dfft internal slow implementation,             *
 * should be improved                                           *
 *    Bert Kampes, 12-Oct-1999                                  *
 ****************************************************************/
  // ______ fftw implementation 2d fft ______
  #ifdef __DEBUGMAT2
      matDEBUG.print("ifft2d: fftw");
  #endif
  int32 nx = A.pixels();
  int32 ny = A.lines();
  fftwf_plan fftwf_tmp_plan;// plan created each time, stupid.  test.
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: creating tmp inv 2d plan."); 
  #endif
  fftwf_tmp_plan  = fftwf_plan_dft_2d(
                          nx, ny,
                          reinterpret_cast<fftwf_complex*>(A[0]), 
                          reinterpret_cast<fftwf_complex*>(A[0]), 
                          FFTW_BACKWARD,
                          FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  #ifdef __DEBUGMAT2
  fftwf_print_plan(fftwf_tmp_plan);// nerd-readable plan output
  matDEBUG.print(" ");
  #endif
  // ___ execute plan ___
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: executing tmp inv 2d plan."); 
  #endif
  fftwf_execute(fftwf_tmp_plan);
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: destroying tmp inv 2d plan."); 
  #endif
  fftwf_destroy_plan(fftwf_tmp_plan);// also associated data?
  // ___ normalize inverse transform, not done by fftw ___
  #ifdef __DEBUGMAT2
  matDEBUG.print("BK: FFTW: normalizing inverse transform."); 
  #endif
  A /= complr4(nx*ny);
  //for (int32 i=0; i<ny; ++i)
  //  for (int32 j=0; i<nx; ++j)
  //#ifdef __GNUC__
  //        A[i][j] /= nx*ny;
  //#else
  //        A[i][j] /= complr4(nx*ny);
  //#endif
  } // END ifft2d
#endif



// ====== Internal implementation using series of 1d fft-s ======
#ifndef __USE_FFTW_LIBRARY__
#ifndef __USE_VECLIB_LIBRARY__
  #ifdef __DEBUGMAT2
    matDEBUG.print("ifft2d: internal");
  #endif
  ifft(A,2);                            // inverse transform over rows
  ifft(A,1);                            // inverse transform over columns
  } // END ifft2d
#endif
#endif



/****************************************************************
 * B=oversample(A, factor);             harmonic interpolation  *
 *    Bert Kampes, 01-Feb-1999                                  *
 *    removed from code.                                        *
 * B=oversample(A, factorrow, factorcol);                       *
 *    2 factors possible, extrapolation at end.                 *
 *    no vectors possible.                                      *
 *    Bert Kampes, 28-Mar-2000                                  *
 ****************************************************************/
matrix<complr4> oversample(
        const matrix<complr4> &AA,// not by reference, changed by fft
        uint factorrow,
        uint factorcol)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("oversample");
#endif
  matrix<complr4> A = AA;// copy, AA is changed by in-place fft;
  const uint l      = A.lines();
  const uint p      = A.pixels();
  const uint halfl  = l/2;
  const uint halfp  = p/2;
  const uint L2     = factorrow*l;      // numrows of output matrix
  const uint P2     = factorcol*p;      // columns of output matrix

  #ifdef __DEBUGMAT1
  if (A.isvector())
    matERROR.print("OVERSAMPLE: only 2d matrices.");
  if (!myispower2(l) && factorrow != 1)
    matERROR.print("OVERSAMPLE: numlines != 2^n.");
  if (!myispower2(p) && factorcol != 1)
    matERROR.print("OVERSAMPLE: numcols != 2^n.");
  #endif

  #ifdef __GNUC__
  const real4 half = 0.5;
  #else
  const complr4 half = complr4(0.5);
  #endif
  matrix<complr4> Res(L2,P2);
  register int32 i,j;
  if (factorrow==1)
    {
    fft(A,2);                           // 1d fourier transform per row
    for (i=0; i<l; ++i)                 // divide by 2 'cause even fftlength
      A(i,halfp) *= half;               // complex *=?
    const window winA1(0, l-1, 0, halfp);       // zero padding windows
    const window winA2(0, l-1, halfp, p-1);
    const window winR2(0, l-1, P2-halfp, P2-1);
    Res.setdata(winA1,A,winA1); 
    Res.setdata(winR2,A,winA2); 
    ifft(Res,2);                        // inverse fft per row
    }
  else if (factorcol==1)
    {
    fft(A,1);                           // 1d fourier transform per column
    for (i=0; i<p; ++i)                 // divide by 2 'cause even fftlength
      A(halfl,i) *= half;
    const window winA1(0, halfl, 0, p-1);       // zero padding windows
    const window winA2(halfl, l-1, 0, p-1);
    const window winR2(L2-halfl, L2-1, 0, p-1);
    Res.setdata(winA1,A,winA1); 
    Res.setdata(winR2,A,winA2); 
    ifft(Res,1);                        // inverse fft per row
    }
  else
    {
    fft2d(A);                           // A=fft2d(A)
    for (i=0; i<l; ++i) A(i,halfp) *= half;
    for (i=0; i<p; ++i) A(halfl,i) *= half;
    const window winA1(0, halfl,   0, halfp);   // zero padding windows
    const window winA2(0, halfl,   halfp, p-1);
    const window winA3(halfl, l-1, 0, halfp);
    const window winA4(halfl, l-1, halfp, p-1);
    const window winR2(0, halfl,       P2-halfp, P2-1);
    const window winR3(L2-halfl, L2-1, 0, halfp);
    const window winR4(L2-halfl, L2-1, P2-halfp, P2-1);
    Res.setdata(winA1,A,winA1); 
    Res.setdata(winR2,A,winA2); 
    Res.setdata(winR3,A,winA3); 
    Res.setdata(winR4,A,winA4); 
    ifft2d(Res);
    }
  Res *= real4(factorrow*factorcol);
  return Res;
  } // END oversample



/****************************************************************
 * B=oversample(A, f1, f2);     calls complex oversampling      *
 *    Bert Kampes, 28-Mar-2000                                  *
 ****************************************************************/
matrix<real4> oversample(
        const matrix<real4> &A,
        uint factorrow,
        uint factorcol)
  {
  // not efficient, mat2cr4 allocates new, and oversample copies this.
  matrix<complr4> TMP = oversample(mat2cr4(A),factorrow,factorcol);
  return real(TMP);                     // imag == 0
  }



/****************************************************************
 * dotmult(complr4* pnt, const matrix<real4>&B, stride)         *
 * multiply memory pointed to by content of B, with mem. stride *
 *    Bert Kampes, 04-Apr-2000                                  *
 ****************************************************************/
void dotmult(
        complr4 *pntR,
        const matrix<real4> &B,
        int32 stride)
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("dotmult, no range checking");
#endif
  real4   *pntB   = B[0];
  for (register int32 i=0; i<B.size(); i++)
    {
#ifdef __GNUC__
    //(*pntR) *= (*pntB++);
    // changed by FvL (for g++/gcc > 4.0):
    {    
    (*pntR) *= (*pntB);
    *pntB++;
    }
#else
    //(*pntR) *= complr4(*pntB++);
    // changed by FvL (for g++/gcc > 4.0):
    {    
    (*pntR) *= complr4(*pntB);
    *pntB++;
    }
#endif
    pntR += stride;
    }
  } // dotmult



// ====== First declare (define?) specialized function ======
// ______ Before template function ______
#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C = A * B; real4                                             *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
matrix<real4> operator * (const matrix<real4>& A, const matrix<real4>& B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matrices * (veclib)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.lines())
    matERROR.print("matrix::operator *: multiplication not possible");
  if (!A.size())
    matERROR.print("matrix:: operator * with empty matrices.");
  #endif
  matrix<real4> Result(A.lines(),B.pixels());
  int32 n       = B.pixels();
  int32 m       = A.lines();
  int32 k       = A.pixels();
  real4 r4alpha = 1.;
  real4 r4beta  = 0.;
  sgemm("N","N",&n,&m,&k,&r4alpha,B[0],&n,A[0],&k,&r4beta,Result[0],&n,1,1);
  return Result; 
  } // END *
#endif


#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C = A * B; real8                                             *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
matrix<real8> operator * (const matrix<real8>& A, const matrix<real8>& B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matrices * (veclib)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.lines())
    matERROR.print("matrix::operator *: multiplication not possible");
  if (!A.size())
    matERROR.print("matrix:: operator * with empty matrices.");
  #endif
  matrix<real8> Result(A.lines(),B.pixels());
  int32 n       = B.pixels();
  int32 m       = A.lines();
  int32 k       = A.pixels();
  real8 r8alpha =1.;
  real8 r8beta = 0.;
  dgemm("N","N",&n,&m,&k,&r8alpha,B[0],&n,A[0],&k,&r8beta,Result[0],&n,1,1);
  return Result; 
  } // END *
#endif


#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C = A * B; complex real4                                             *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
matrix<complr4> operator * (const matrix<complr4>& A, const matrix<complr4>& B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matrices * (veclib)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.lines())
    matERROR.print("matrix::operator *: multiplication not possible");
  if (!A.size())
    matERROR.print("matrix:: operator * with empty matrices.");
  #endif
  matrix<complr4> Result(A.lines(),B.pixels());
  int32 n       = B.pixels();
  int32 m       = A.lines();
  int32 k       = A.pixels();
  complr4 c4alpha(1.,0.);
  complr4 c4beta(0.,0.);
  cgemm("N","N",&n,&m,&k,&c4alpha,B[0],&n,A[0],&k,&c4beta,Result[0],&n,1,1);
  return Result; 
  } // END *
#endif


#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C = A * B; complex real8                                             *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
matrix<complr8> operator * (const matrix<complr8>& A, const matrix<complr8>& B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matrices * (veclib)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.lines())
    matERROR.print("matrix::operator *: multiplication not possible");
  if (!A.size())
    matERROR.print("matrix:: operator * with empty matrices.");
  #endif
  matrix<complr8> Result(A.lines(),B.pixels());
  int32 n       = B.pixels();
  int32 m       = A.lines();
  int32 k       = A.pixels();
  complr8 c8alpha = 1.;
  complr8 c8beta  = 0.;
  zgemm("N","N",&n,&m,&k,&c8alpha,B[0],&n,A[0],&k,&c8beta,Result[0],&n,1,1);
  return Result; 
  } // END *
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<real4> matTxmat(const matrix<real4> &A, const matrix<real4> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("matTxmat, veclib::sgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines() != B.lines())
    matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
  #endif
  matrix<real4> Result(A.pixels(),B.pixels());
  int32 m = A.pixels();
  int32 n = B.pixels();
  int32 k = B.lines();
  real4 r4alpha =1.;
  real4 r4beta = 0.;
  sgemm("N","T",&n,&m,&k,&r4alpha,B[0],&n,A[0],&m,&r4beta,Result[0],&n,1,1);
  return Result;
  } // END matTxmat
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<real8> matTxmat(const matrix<real8> &A, const matrix<real8> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("matTxmat, veclib::dgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines() != B.lines())
    matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
  #endif
  matrix<real8> Result(A.pixels(),B.pixels());
  int32 m = A.pixels();
  int32 n = B.pixels();
  int32 k = B.lines();
  real8 r8alpha = 1.;
  real8 r8beta  = 0.;
  dgemm("N","T",&n,&m,&k,&r8alpha,B[0],&n,A[0],&m,&r8beta,Result[0],&n,1,1);
  return Result;
  } // END matTxmat
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<complr4> matTxmat(const matrix<complr4> &A, const matrix<complr4> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("matTxmat, veclib::cgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines() != B.lines())
    matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
  #endif
  matrix<complr4> Result(A.pixels(),B.pixels());
  int32 m = A.pixels();
  int32 n = B.pixels();
  int32 k = B.lines();
  complr4 c4alpha(1.,0.);
  complr4 c4beta(0.,0.);
  cgemm("N","T",&n,&m,&k,&c4alpha,B[0],&n,A[0],&m,&c4beta,Result[0],&n,1,1);
  return Result;
  } // END matTxmat
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<complr8> matTxmat(const matrix<complr8> &A, const matrix<complr8> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("matTxmat, veclib::zgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines() != B.lines())
    matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
  #endif
  matrix<complr8> Result(A.pixels(),B.pixels());
  int32 m = A.pixels();
  int32 n = B.pixels();
  int32 k = B.lines();
  complr8 c8alpha = 1.;
  complr8 c8beta  = 0.;
  zgemm("N","T",&n,&m,&k,&c8alpha,B[0],&n,A[0],&m,&c8beta,Result[0],&n,1,1);
  return Result;
  } // END matTxmat
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<real4> matxmatT(const matrix<real4> &A, const matrix<real4> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matxmatT, veclib::sgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.pixels())
    matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
  #endif
  matrix<real4> Result(A.lines(),B.lines());
  int32 m = A.lines();
  int32 n = B.lines();
  int32 k = B.pixels();
  real4 r4alpha =1.;
  real4 r4beta = 0.;
  sgemm("T","N",&n,&m,&k,&r4alpha,B[0],&k,A[0],&k,&r4beta,Result[0],&n,1,1);
  return Result;
  } // END matxmatT
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<real8> matxmatT(const matrix<real8> &A, const matrix<real8> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matxmatT, veclib::dgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.pixels())
    matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
  #endif
  matrix<real8> Result(A.lines(),B.lines());
  int32 m = A.lines();
  int32 n = B.lines();
  int32 k = B.pixels();
  real8 r8alpha =1.;
  real8 r8beta = 0.;
  dgemm("T","N",&n,&m,&k,&r8alpha,B[0],&k,A[0],&k,&r8beta,Result[0],&n,1,1);
  return Result;
  } // END matxmatT
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<complr4> matxmatT(const matrix<complr4> &A, const matrix<complr4> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matxmatT, veclib::cgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.pixels())
    matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
  #endif
  matrix<complr4> Result(A.lines(),B.lines());
  int32 m = A.lines();
  int32 n = B.lines();
  int32 k = B.pixels();
  complr4 c4alpha(1.,0.);
  complr4 c4beta(0.,0.);
  cgemm("T","N",&n,&m,&k,&c4alpha,B[0],&k,A[0],&k,&c4beta,Result[0],&n,1,1);
  return Result;
  } // END matxmatT
#endif



#ifdef __USE_VECLIB_LIBRARY__
/****************************************************************
 * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<complr8> matxmatT(const matrix<complr8> &A, const matrix<complr8> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matxmatT, veclib::zgemm");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.pixels())
    matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
  #endif
  matrix<complr8> Result(A.lines(),B.lines());
  int32 m = A.lines();
  int32 n = B.lines();
  int32 k = B.pixels();
  complr8 c8alpha = 1.;
  complr8 c8beta  = 0.;
  zgemm("T","N",&n,&m,&k,&c8alpha,B[0],&k,A[0],&k,&c8beta,Result[0],&n,1,1);
  return Result;
  } // END matxmatT
#endif






//#ifdef WIN32
// Jia defined a function which I think is identical.
// maybe it did not compile under windows.
// but, then it should be placed in matrixspecs.cc and
// a prototype below the class in matrixbk.hh
// Bert Kampes, 24-Aug-2005
//
/****************************************************************
 * C=diagxmat(vec,B) C=diag(vec) * B; (diag R4, mat CR4)        *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
matrix<complr4>  diagxmat    (const matrix<real4> &diag, const matrix<complr4> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("diagxmat: R4*CR4");
  #endif
  #ifdef __DEBUGMAT1
  if (min(diag.lines(),diag.pixels()) != 1)
    matERROR.print("diagxmat: sizes A,B: diag is vector.");
  if (diag.size() != B.lines())
    matERROR.print("diagxmat: sizes A,B: input is vector, matrix.");
  #endif

  matrix<complr4> Result=B;
  if (diag.lines() != 1)        // standing
    {
    for (register uint i=0; i<Result.lines(); i++)
      for (register uint j=0; j<Result.pixels(); j++)
        Result(i,j) *= diag(i,0);
    }
  else
    {
    for (register uint i=0; i<Result.lines(); i++)
      for (register uint j=0; j<Result.pixels(); j++)
        Result(i,j) *= diag(0,i);
    }
  return Result;
  } // END diagxmat



/****************************************************************
 * A=cos(B)                                                     *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<real4> cos(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("cos");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());
  real4 *pntR = Result[0];
  real4 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = cos(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = cos(*pntA);
    *pntA++;
    }
  return Result;
  } // END cos


#ifndef NO_FASTTRIG // [MA]
/****************************************************************
 * A=fast_cos(B)                                                *
 #%// BK 06-Oct-2005                                            *
 ****************************************************************/
matrix<real4> fast_cos(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("fast_cos");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());
  real4 *pntR = Result[0];
  real4 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = fast_cos(*pntA++);// LUT
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = fast_cos(*pntA);// LUT
    *pntA++;
    }
  return Result;
  } // END fast_cos
#endif

/****************************************************************
 * A=cos(B)                                                     *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<real8> cos(
        const matrix<real8> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("cos");
  #endif
  matrix<real8> Result(A.lines(),A.pixels());
  real8 *pntR = Result[0];
  real8 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = cos(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = cos(*pntA);
    *pntA++;
    }
  return Result;
  } // END cos



/****************************************************************
 * A=sin(B)                                                     *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<real4> sin(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("sin");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());
  real4 *pntR = Result[0];
  real4 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = sin(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = sin(*pntA);
    *pntA++;
    }
  return Result;
  } // END sin



#ifndef NO_FASTTRIG // [MA]
/****************************************************************
 * A=fast_sin(B)                                                *
 #%// BK 06-Oct-2005                                            *
 ****************************************************************/
matrix<real4> fast_sin(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("fast_sin");
  #endif
  matrix<real4> Result(A.lines(),A.pixels());
  real4 *pntR = Result[0];
  real4 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = fast_sin(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = fast_sin(*pntA);
    *pntA++;
    }
  return Result;
  } // END fast_sin
#endif


/****************************************************************
 * A=sin(B)                                                     *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<real8> sin(
        const matrix<real8> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("sin");
  #endif
  matrix<real8> Result(A.lines(),A.pixels());
  real8 *pntR = Result[0];
  real8 *pntA = A[0];
  // ______Ensure stride one memory______
  //for (register uint i=0; i<A.size(); i++)
  //  (*pntR++) = sin(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); i++)
    {
    (*pntR++) = sin(*pntA);
    *pntA++;
    }
  return Result;
  } // END sin



/****************************************************************
 * B=mat2cr4(A); conversion to complex<real4>                   *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
matrix<complr4> mat2cr4(
        const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mat2cr4. (real4)");
  #endif
  matrix<complr4> Result(A.lines(),A.pixels());
  real4   *pntA = A[0];
  complr4 *pntR = Result[0];
  //for (register uint i=0; i<Result.size(); i++)
  //(*pntR++) = complr4(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<Result.size(); i++)
    {
    (*pntR++) = complr4(*pntA);
    *pntA++;
    }
  return Result;
  } // END mat2cr4 (real4)



/****************************************************************
 * A=mat2cr4(B,C)                                               *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<complr4> mat2cr4(
        const matrix<real4>& A,
        const matrix<real4>& B)
  {
  #ifdef __DEBUGMAT1
  if (A.lines()!=B.lines() || A.lines()!=B.lines())
    matERROR.print("operator complr4, input matrices not same size");
  #endif
  #ifdef __DEBUGMAT2
  matDEBUG.print("operator complr4(r4,r4)");
  #endif
  matrix<complr4> Result(A.lines(),A.pixels());
  real4   *pntA = A[0];
  real4   *pntB = B[0];
  complr4 *pntR = Result[0];
  //for (register uint i=0; i<Result.size(); i++)
  //(*pntR++) = complr4(*pntA++,*pntB++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<Result.size(); i++)
    {
    (*pntR++) = complr4(*pntA,*pntB);
    *pntA++;
    *pntB++;
    }
  return Result;
  } // END mat2cr4 complr4 (real4,real4)



/****************************************************************
 * A=mat2cr4(r8B,r8C)                                           *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
matrix<complr4> mat2cr4(
        const matrix<real8>& A,
        const matrix<real8>& B)
  {
  #ifdef __DEBUGMAT1
  if (A.lines()!=B.lines() || A.lines()!=B.lines())
    matERROR.print("operator complr4, input matrices not same size");
  #endif
  #ifdef __DEBUGMAT2
  matDEBUG.print("operator complr4(r8,r8)");
  #endif
  matrix<complr4> Result(A.lines(),A.pixels());
  real8   *pntA = A[0];
  real8   *pntB = B[0];
  complr4 *pntR = Result[0];
  //for (register uint i=0; i<Result.size(); i++)
  //(*pntR++) = complr4(*pntA++,*pntB++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<Result.size(); i++)
    {
    (*pntR++) = complr4(*pntA,*pntB);
    *pntA++;
    *pntB++;
    }
  return Result;
  } // END mat2cr4(real8,real8)




/****************************************************************
 ****************************************************************
 ****************************************************************
 SPECIALIZATIONS MUST BE IN HERE DIRECTLY AFTER CLASS
 for g++-4.0 and up;
 Bert Kampes, 07-Oct-2005
 ****************************************************************
 ****************************************************************
 ****************************************************************/


#ifdef WIN32
  /********************************************************************
   * C *= A;      pointwise multiplication, CR4 R4              *
   #%// BK 26-Oct-2000                        *
   replace the origin pointwise multiplication, CR4 R4  Jia YouLiang  *
   ********************************************************************/
  matrix<complr4> timesCxR (matrix<complr4> &B, matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("*=");
  #endif
  #ifdef __DEBUGMAT1
  if (B.nrows != A.lines() || B.ncols != A.pixels())
    matERROR.print("matrix:: *= matrices must be same size.");
  #endif
  complr4 *pntmat = B[0];
  real4   *pntA   = A[0];
  //for (register uint i=0; i<B.size(); i++)
  //  (*pntmat++) *= (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<B.size(); i++)
    {
    (*pntmat++) *= (*pntA);
    *pntA++;
    }
  return B;
  } // END *= matrices CR4 with R4 by Jia YouLiang
#else


/****************************************************************
 For gcc4.0 it seemed template<> is required to specialize this.
 I guess this means it should be added everywhere, also for Veclib.
 // Bert Kampes, 06-Oct-2005
 ****************************************************************/
/****************************************************************
 * C *= A;      pointwise multiplication, CR4 R4                *
 #%// BK 26-Oct-2000
 ****************************************************************/
#if __GNUC__ < 4
//template<> // essential for specialization
//template<class>  // seems not allowed for g++2.95
//#if __GNUC_MINOR__ < 95
//#endif
//
#else 
#if __GNUC_MINOR__ > 0
template<> // essential for specialization
template<>  // seems required for specialization g++-4.0 upwards
#else
template<> // essential for specialization
template<class>  // seems required for specialization g++-4.0 upwards
#endif
#endif
// ****************************************************************
matrix<complr4>& matrix<complr4>::operator *= (const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("*=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: *= matrices must be same size.");
  #endif
  complr4 *pntmat = data[0];
  real4   *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) *= (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) *= (*pntA);
    *pntA++;
    }
  return *this;
  } // END *= matrices CR4 with R4


/****************************************************************
 * C /= A;      pointwise division, CR4 R4                      *
 #%// BK 27-Oct-2000                                            *
 ****************************************************************/
#if __GNUC__ < 4
//template<> // essential for specialization
//template<class>  // seems not allowed for g++2.95
//#if __GNUC_MINOR__ < 95
//#endif
//
#else   /* g++-4.0 upwards */
#if __GNUC_MINOR__ > 0
template<> // essential for specialization
template<>  // seems required for specialization g++-4.0 upwards
#else
template<> // essential for specialization
template<class>  // seems required for specialization g++-4.0 upwards
#endif
#endif
matrix<complr4>& matrix<complr4>::operator /= (const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("/=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: *= matrices must be same size.");
  #endif
  complr4 *pntmat = data[0];
  real4   *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) /= (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) /= (*pntA);
    *pntA++;
    }
  return *this;
  } // END /= matrices CR4 with R4


/****************************************************************
 * C += A;      pointwise addition, CR4 R4              *
 #%// BK 27-Oct-2000                                            *
 ****************************************************************/
#if __GNUC__ < 4
//template<> // essential for specialization
//template<class>  // seems not allowed for g++2.95
//#if __GNUC_MINOR__ < 95
//#endif
//
#else   /* g++-4.0 upwards */
#if __GNUC_MINOR__ > 0
template<> // essential for specialization
template<>  // seems required for specialization g++-4.0 upwards
#else
template<> // essential for specialization
template<class>  // seems required for specialization g++-4.0 upwards
#endif
#endif
matrix<complr4>& matrix<complr4>::operator += (const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("+=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: += matrices must be same size.");
  #endif
  complr4 *pntmat = data[0];
  real4   *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) += (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) += (*pntA);
    *pntA++;
    }
  return *this;
  } // END += matrices CR4 with R4


/****************************************************************
 * C -= A;      pointwise subtraction, CR4 R4           *
 #%// BK 27-Oct-2000                                            *
 ****************************************************************/
#if __GNUC__ < 4
//template<> // essential for specialization
//template<class>  // seems not allowed for g++2.95
//#if __GNUC_MINOR__ < 95
//#endif
//
#else   /* g++-4.0 upwards */
#if __GNUC_MINOR__ > 0
template<> // essential for specialization
template<>  // seems required for specialization g++-4.0 upwards
#else
template<> // essential for specialization
template<class>  // seems required for specialization g++-4.0 upwards
#endif
#endif
matrix<complr4>& matrix<complr4>::operator -= (const matrix<real4> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("-=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: -= matrices must be same size.");
  #endif
  complr4 *pntmat = data[0];
  real4   *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) -= (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) -= (*pntA);
    *pntA++;
    }
  return *this;
  } // END -= matrices CR4 with R4


/****************************************************************
 * dummy function: use the operators that are specialized, otherwise
 * they are not generated in the object code!
 * this is for g++-4.0
 #%// Bert Kampes, 07-Oct-2005
 ****************************************************************/
void dummy_calls_to_specs()
  {
  matrix<complr4> A(1,1);
  matrix<real4>   B(1,1);
  A*=B;
  A+=B;
  A-=B;
  A/=B;// impossible divide zero, but not executed
  } // dummy_calls_to_specs

#endif


