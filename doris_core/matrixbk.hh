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
/************************************************************************
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixbk.hh,v $
 * $Revision: 3.16 $                                                    *
 * $Date: 2005/10/18 13:46:51 $                                         *
 * $Author: kampes $                                                    *
 *                                                                      *
 *  template (base) class for matrices.                                 *
 *                                                                      *
 *  ====== Defines (compiler -D option, see Makefile) ======            *
 *  __USE_FFTW_LIBRARY__ do use the fftw lib (compile separately)       *
 *                          prefered over veclib fft.                   *
 *  __USE_VECLIB_LIBRARY__ use VECLIB (FFT & matrix multiplication)     *
 *                      else internal routines will be used             *
 *  __USE_LAPACK_LIBRARY__ use LAPACK (cholesky routines)               *
 *                      else internal routines will be used             *
 *  __DEBUGMAT1         check index in matrix, dimensions, etc.         *
 *  __DEBUGMAT2         give info, can savely be un-defined             *
 *                                                                      *
 * Better make a vector baseclass (stl?) and a 2dmatrix class,          *
 *  which are friends of eachother,                                     *
 *  then derive an image class and add functions there.                 *
 * Note that data has to be continuously in memory, because             *
 *  VECLIB/FFTW assumes this (and me as well...) so a vector of vectors *
 *  cannot be used.                                                     *
 ************************************************************************/

#ifndef MATRIXBK_H                              // guard
#define MATRIXBK_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                         // typedefs, window
#include <fstream>                              // ofstream type


//#ifdef __USE_FFTW_LIBRARY__
//  #include <fftw3.h>                                  // fftw types
//#endif

// ______ message objects, global, set in main ______
extern bk_messages matERROR;
extern bk_messages matDEBUG;

#ifdef __USE_LAPACK_LIBRARY__
  #define ilaver ilaver_
  extern "C" { int32 ilaver ( int32*, int32*, int32* ) ; } // [MA] get lapack version, used in processor.cc
#endif

/************************************************************************
friend functions are defined inside the class,
since if they are declared here, and defined in matrixbk.cc
then gcc gives errors like (for friend functions):
In file included from processor.c:23:
matrixbk.h:162: warning: friend declaration `void myswap(matrix<Type> &, matrix<Type> &)'
matrixbk.h:162: warning:   declares a non-template function
matrixbk.h:162: warning:   (if this is not what you intended, make sure
matrixbk.h:162: warning:   the function template has already been declared,
matrixbk.h:162: warning:   and add <> after the function name here)
matrixbk.h:162: warning:   -Wno-non-template-friend disables this warning.

Adding <> only makes the code less readable,
cause other compilers don't want it.
BK 17-Aug-2000
*************************************************************************/

/************************************************************************
Compiler g++-4.0 seems to require that for overloaded template
operators of the template matrix class the following syntax:
  template<> template<class> 
  matrix<complr4>& matrix<complr4>::operator += (const matrix<real4> &A) {code;}  
this function has to be in the same file as the class, below it.
unfortunately, g++2.95.2 does crash if this is done.  it does not
understand the "template<> template<class>" syntax, but more importantly,
it crashes if these functions are put below the class, saying they
are defined more than once!
There does not seem to be an easy solution to this problem.  What seems
to work is to keep them in matrixspecs.cc and to use a dummy call to
all specialized functions.  this guarantees that they are instantiated.
in short, this compiles with
gcc version 4.0.0 20050301 (prerelease) (Debian 4.0-0pre6ubuntu7)
#%// Bert Kampes, 07-Oct-2005
*************************************************************************/



// ====== Define template functions (no member no friend) ======
// ______ (matrix class is declared way below) ______
template <class Type> class matrix;

#ifdef __USE_VECLIB_LIBRARY__
// ______ See: matrixspecs.c ______
// ______ Compiler uses specialization if declared here ______
// ______ and not the template definition in matrixbk.cc ______
matrix<real4>    operator *  (const matrix<real4>& A,   const matrix<real4>& B);
matrix<real8>    operator *  (const matrix<real8>& A,   const matrix<real8>& B);
matrix<complr4>  operator *  (const matrix<complr4>& A, const matrix<complr4>& B);
matrix<complr8>  operator *  (const matrix<complr8>& A, const matrix<complr8>& B);
matrix<real4>    matTxmat    (const matrix<real4> &A,   const matrix<real4> &B);
matrix<real8>    matTxmat    (const matrix<real8> &A,   const matrix<real8> &B);
matrix<complr4>  matTxmat    (const matrix<complr4> &A, const matrix<complr4> &B);
matrix<complr8>  matTxmat    (const matrix<complr8> &A, const matrix<complr8> &B);
matrix<real4>    matxmatT    (const matrix<real4> &A,   const matrix<real4> &B);
matrix<real8>    matxmatT    (const matrix<real8> &A,   const matrix<real8> &B);
matrix<complr4>  matxmatT    (const matrix<complr4> &A, const matrix<complr4> &B);
matrix<complr8>  matxmatT    (const matrix<complr8> &A, const matrix<complr8> &B);
#endif


// ______ File: matrixspecs.c ______
//#if defined (__DEBUGMAT2) || defined (__DEBUGMAT1)    // extra checking
//void matDEBUG(char ch[ONE27]);
//#endif
//void matERROR(char ch[ONE27]);
void matassert(
        const ofstream &str,
        const char* ofilename,
        const char* callingfilename = "?",
        int32 linenumber = 0);
void matassert(
        const ifstream &str,
        const char* ifilename,
        const char* callingfilename = "?",
        int32 linenumber = 0);



// ====== Without VECLIB these are also implemented (slower) ======
// ______ NOW used is www.fftw.org routines ______
// ______ matrixspecs.c: VECLIB [cz]1dfft or internal routine ______
// ______ matrixspecs.c: FFTW  is used since 16-sep-2003 ______
  void            fft           (matrix<complr4> &A, int32 dimension);
  void            ifft          (matrix<complr4> &A, int32 dimension);
  void            fft2d         (matrix<complr4> &A);
  void            ifft2d        (matrix<complr4> &A);
  //void          fft2d         (matrix<real4> &A, matrix<real4> &B);
  //void          ifft2d        (matrix<real4> &A, matrix<real4> &B);

  matrix<complr4> oversample    (const matrix<complr4> &A,// new 8/05: by ref.
                                 uint frow, uint fcol);
  matrix<real4>   oversample    (const matrix<real4>   &A,
                                 uint frow, uint fcol);

// ______ LAPACK_LIBRARY ______
// ______ Without LAPACK, internal routines can be used (slower, inaccurate?) ______
  void            choles        (matrix<real4> &A);
  void            invertchol    (matrix<real4> &A);
  void            solvechol     (const matrix<real4> &A, matrix<real4> &B);
  void            choles        (matrix<real8> &A);
  void            invertchol    (matrix<real8> &A);
  void            solvechol     (const matrix<real8> &A, matrix<real8> &B);


  matrix<real4>   intensity     (const matrix<complr4> &A);
  matrix<real4>   magnitude     (const matrix<complr4> &A);

  matrix<real4>   real          (const matrix<complr4> &A);
  matrix<real4>   imag          (const matrix<complr4> &A);

  real4           norm2         (const matrix<complr4> &A);
  real4           norm2         (const matrix<real4> &A);
  matrix<complr4> norm          (const matrix<complr4> &A);

  matrix<real4>   abs           (const matrix<real4> &A);
  matrix<real8>   abs           (const matrix<real8> &A);

// ______ Read file complex<short> into matrix complex<real4> ______
  void fileci2tomatcr4(
        matrix<complr4>         &Result,
        const char              *file, 
        uint                     filelines,
        window                   win,
        window                   winoffset);

// ______ Casts ______
  // operator complr4 (matrix<real4>) seems better... (but how?)
  //  matrix<complr4> mat2cr4(
  //    const matrix<compli16>  &A);
  matrix<complr4> mat2cr4(
        const matrix<real4>     &A);
  matrix<complr4> mat2cr4(
        const matrix<real4>     &A,
        const matrix<real4>     &B);
  matrix<complr4> mat2cr4(
        const matrix<real8>     &A,
        const matrix<real8>     &B);

// ______ phase ______
  matrix<real4>   angle(
        const matrix<complr4>   &A);    // phase        
  matrix<complr4> angle2cmplx(
        const matrix<real4>     &A);    // phasor, a=1
  void dotmultconjphase(
        matrix<complr4>         &complexinterferogram,//        by ref.
        const matrix<real4>     &refphase);     // phasor, a=1
  // --- Using lookup table ---
  matrix<real4>   fast_angle(
        const matrix<complr4>   &A);    // phase        
  matrix<complr4> fast_angle2cmplx(
        const matrix<real4>     &A);    // phasor, a=1
  matrix<complr8> fast_angle2cmplx(     // MA
        const matrix<real8>     &A);    // phasor, a=1
  void fast_dotmultconjphase(
        matrix<complr4>         &complexinterferogram,//        by ref.
        const matrix<real4>     &refphase);     // phasor, a=1
  void fast_dotmultconjphase(
        matrix<complr4>         &complexinterferogram,//        by ref.
        const matrix<real8>     &refphase);     // phasor, a=1


// ______ complex coherence ______
  matrix<complr4> coherence(
        const matrix<complr4>   &complex_interferogram,
        const matrix<complr4>   &norm_image1_and_2,
        uint                     estimatorwinsizeL,
        uint                     estimatorwinsizeP);

// ______ complex coherence ______
  matrix<real4> coherence2(
        const matrix<complr4>   &complex_interferogram,
        const matrix<complr4>   &norm_image1_and_2,
        uint                     estimatorwinsizeL,
        uint                     estimatorwinsizeP);

// should be in matrix class?
// ______ sort rows of matrix on some column; uses qsort ______
//    void mysort1(matrix<real4> &A);   // sort matrix based on col. number;
// ______ (ascending) sort rows of matrix on first col, then second; uses qsort ______
  void mysort2(matrix<real4> &A);       // mysort(A);
  void mysort2(matrix<int32> &A);       // mysort(A);
  void mysort231(matrix<real4> &A);     // mysort(A);
  void mysort2selcol(matrix<real4> &A, int32 selcol);
  void mysort2selcol(matrix<int32> &A, int32 selcol);


// ______ multiply strike x memory with values in B ______
// ______ to multiply a column of A by vector B use:
// ______ complr4 pntA=&A[0][c], strike=numpixels(A), length(B)=numrows.
  void dotmult (complr4 *startaddress, const matrix<real4> &B, int32 strike);
  //void dotmultrow (matrix &, row, matrix<real4> &B);
  //void dotmultcol (matrix &, col, const matrix<real4> &B);
  //matrix<complr4> dotmult(const matrix<complr4> &B, const matrix<real4> &A);
  //matrix<complr4> dotmult(const matrix<real4> &A,   const matrix<complr4> &B);

// ______ Some trigonometric functions #%// BK 09-Nov-2000 ______
// add lookup table for fast trig! (someday..)
// Bert Kampes, 10-Apr-2005
  matrix<real4> cos (const matrix<real4> &A);
  matrix<real4> fast_cos (const matrix<real4> &A);// lookup table
  matrix<real8> cos (const matrix<real8> &A);
  matrix<real4> sin (const matrix<real4> &A);
  matrix<real4> fast_sin (const matrix<real4> &A);// lookup table
  matrix<real8> sin (const matrix<real8> &A);

// ______ Some casts for complex (Re,Im) ______ HOW?
//matrix<complr4> operator complr4  (const matrix<real4>& A);
//matrix<complr4> operator complr4  (const matrix<real8>& A,   const matrix<real8>& B);
//matrix<complr8> operator complr8  (const matrix<real8>& A,   const matrix<real8>& B);


// ______ See: matrixbk.cc for definition (bottom) ______
// ______ functions are no friends no members ______
template <class Type>
  matrix<Type>  operator * (const matrix<Type>& A, const matrix<Type>& B);
template <class Type>
  matrix<Type>  operator * (const matrix<Type>& A, Type scalar);
template <class Type>
  matrix<Type>  operator *  (Type  scalar, const matrix<Type> &A);
template <class Type>
  matrix<Type>  operator / (const matrix<Type>& A, Type B);
template <class Type>
  matrix<Type>  operator / (const matrix<Type>& A, const matrix<Type>& B);
template <class Type>
  matrix<Type>  operator - (const matrix<Type>& A, const matrix<Type>& B);
template <class Type>
  matrix<Type>  operator - (const matrix<Type>& A, Type scalar);
template <class Type>
  matrix<Type>  operator + (const matrix<Type>& A, const matrix<Type>& B);
template <class Type>
  matrix<Type>  operator + (const matrix<Type>& A, Type B);
template <class Type>
  Type          max(const matrix<Type> &A);
template <class Type>
  Type          max(const matrix<Type> &A, uint& line, uint& pixel);
template <class Type>
  Type          min(const matrix<Type> &A);
template <class Type>
  Type          min(const matrix<Type> &A, uint& line, uint& pixel);
template <class Type>
  matrix<Type>  matTxmat(const matrix<Type> &A, const matrix<Type> &B);
template <class Type>
  matrix<Type>  matxmatT(const matrix<Type> &A, const matrix<Type> &B);
template <class Type>
  void dumpasc(const char *file, const matrix<Type>& A);
template <class Type>
  matrix<Type>  dotmult     (const matrix<Type> &A, const matrix<Type> &B);
template <class Type>
  matrix<Type>  dotdiv      (const matrix<Type> &A, const matrix<Type> &B);
template <class Type>
  matrix<Type>  sqr         (const matrix<Type> &A);
template <class Type>
  matrix<Type>  conj        (const matrix<Type> &A);
template <class Type>
  matrix<Type>  diagxmat    (const matrix<Type> &diag, const matrix<Type> &B);
// ______ Specialisation R4*CR4 ______
// for windows it seems this is not supported.  Jia defined a new function
// I hope.
// Bert Kampes, 24-Aug-2005
#ifndef WIN32
matrix<complr4> diagxmat    (const matrix<real4> &diag, const matrix<complr4> &B);
#endif
template <class Type>
  matrix<Type>  multilook   (const matrix<Type> &A, uint factorL, uint factorP);
template <class Type, class Type2>
  void  converttype         (const matrix<Type> &A, matrix<Type2> &B);
template <class Type>
  matrix<real4> correlate   (const matrix<Type> &A, matrix<Type> Mask);
template <class Type>
  matrix<Type>  operator -  (const matrix<Type>& A);
template <class Type>
  real8         mean        (const matrix<Type> &A);
template <class Type>
  matrix<Type>  sum         (const matrix<Type> &A, int32 dim);
//template <class Type>
//  matrix<Type>  operator ^ (const matrix<Type>& A);
//template <class Type>
//  matrix<Type>  operator ^ (Type scalar);






// ====== ====== ====== ====== ====== ======
// ====== Start of class declaration: ======
// ====== ====== ====== ====== ====== ======
template <class Type>                           // class template
class matrix
  {
  private:
    // ______ Private data ______
    Type **data;                                // two dimensional array
    uint nrows;                                 // number of rows (lines)
    uint ncols;                                 // number of rows (pixels)
    uint nsize;                                 // lines*pixels

    // ______ Private functions ______
    void   allocate                                         (uint l, uint p);
    void   initialize                                       (uint l, uint p);

    #ifdef __DEBUGMAT1
    void   checkindex                                       (uint l, uint p) const;
    #endif

  public:
    // ______ Constructors ______
    matrix                      ();
    matrix                      (uint l, uint p);
    matrix                      (const matrix<Type>& A);
    matrix                      (window w, const matrix<Type>& A);

    // ______ Destructor ______
    ~matrix                     ();

    // ______ Data functions ______
    void          setdata       (Type w);
    void          setdata       (uint l, uint p, const matrix<Type>& A);
    void          setdata       (window winin, const matrix<Type> &A, window winA);
    void          setdata       (const matrix<Type> &A, window winA);
    matrix<Type>  getdata       (window w)                      const;
    matrix<Type>  getrow        (uint l)                        const;
    matrix<Type>  getcolumn     (uint p)                        const;
    // later... index matrix, mask, ... functies findlt findgt, etc.
    // later,,, void operator ()   (indexmatrix, mask, or whatever)
    // BK 23-Oct-2000
    //matrix<int32> find          (Type w)                              const;
    void          showdata      ()                              const;
    inline bool   isvector      ()                              const;
    inline uint   lines         ()                              const;
    inline uint   pixels        ()                              const;
    inline uint   size          ()                              const;
    void          resize        (uint l, uint p);
    void          clean         ();

    void          setrow        (uint l, const matrix<Type> &L);
    void          setrow        (uint l, Type scalar);
    void          setcolumn     (uint p, const matrix<Type> &C);
    void          setcolumn     (uint p, Type scalar);
    void          fliplr        ();
    void          flipud        ();

    inline Type*  operator []   (uint l)                        const;
    inline Type&  operator ()   (uint l, uint p)                const;
    matrix<Type>  operator ()   (const window &win)             const;
    matrix<Type>  operator ()   (const uint &l0, const uint &lN,
                                 const uint &p0, const uint &pN) const;

    matrix<Type>& operator =    (const matrix<Type>& A);
    matrix<Type>& operator =    (const Type scalar);
    matrix<Type>& operator -=   (Type scalar);
    matrix<Type>& operator -=   (const matrix<Type>& A);
    matrix<Type>& operator +=   (Type scalar);
    matrix<Type>& operator +=   (const matrix<Type>& A);
    matrix<Type>& operator *=   (Type scalar);
    matrix<Type>& operator *=   (const matrix<Type> &A);
    matrix<Type>& operator /=   (Type scalar);
    matrix<Type>& operator /=   (const matrix<Type> &A);
    bool          operator ==   (Type scalar)                   const;
    bool          operator ==   (const matrix<Type> &A)         const;
    bool          operator !=   (Type scalar)                   const;
    bool          operator !=   (const matrix<Type> &A)         const;
    void          conj          ();
    void          mypow         (Type s);               // better use operator ^

// for MS-Windows it seems this is not supported.  Jia defined a new function
// timesCxR instead of "*=" with complex and real input
// Bert Kampes, 24-Aug-2005
#ifdef WIN32
  // maybe this function needs to be moved out of class for compiler
  matrix<complr4> timesCxR (matrix<complr4> &B, matrix<real4> &A);
#else
    // ______ Functions to, e.g., multiply a CR4 by a R4 matrix ______
    // ______ Declare for all, but only define for a few in matrixspecs.cc ______
    // ______ for g++-4.0 this required some additions in matrixspecs.cc ______
    template <class TypeB>
      matrix<Type>& operator *=   (const matrix<TypeB> &A);
    template <class TypeB>
      matrix<Type>& operator /=   (const matrix<TypeB> &A);
    template <class TypeB>
      matrix<Type>& operator +=   (const matrix<TypeB> &A);
    template <class TypeB>
      matrix<Type>& operator -=   (const matrix<TypeB> &A);
#endif



// ================
//  Template friend functions have to be defined within class
//  (where they are declared) to satisfy with g++ preference (why??)
// ================


/****************************************************************
 * file << A;                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
friend ostream&      operator << (ostream& file, const matrix<Type>& A)
  {                                                              // uses copy constructor
  #ifdef __DEBUGMAT2
    matDEBUG.print("operator <<");
    matDEBUG <<  "outfile << [A] size: " << A.size() << " lines: " << A.lines() << " pixels: " << A.pixels()  << " address: " << &A;
    matDEBUG.print();
  #endif
  if ( A.nsize == 0)   // [MA] fix for empty matrices due to multilooking ... etc.
    {
    cout << "op <<   : Buffer was empty, nothing done.\n" ;  // TEMP
    DEBUG.print("input matrix buffer was empty, nothing is dumped to file.");
    return file; // do nothing
    }
  const uint sizeofType = sizeof(Type);
  for (register uint i=0;i<A.nrows;i++)
    for (register uint j=0;j<A.ncols;j++)
      file.write((char*)&A.data[i][j],sizeofType);
  return file;
  } // END <<



/****************************************************************
 * file >> A;                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
//friend istream&      operator >> (istream& file, const matrix<Type>& A)
// BK 25-Sep-2000
friend istream&      operator >> (istream& file, matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("operator >>");
  #endif
  const uint sizeofType = sizeof(Type);
  for (register uint i=0;i<A.nrows;i++)
    for (register uint j=0;j<A.ncols;j++)
      file.read((char*)&A.data[i][j],sizeofType);
  // too slow:  
  // file.read((char*)&A.data[0][0],A.nrows*A.ncols*sizeof(Type));
  return file;
  } // END >>


 
/****************************************************************
 * myswap(A,B)                                                  *
 *    interchange matrices of same size                         *
 *    somehow g++ does not like name swap, so myswap            *
 *    Bert Kampes, 14-Apr-1999                                  *
 ****************************************************************/
friend void          myswap      (matrix<Type> &A, matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("myswap.");
  #endif
  #ifdef __DEBUGMAT1
  if (A.nrows  != B.nrows ||
      A.ncols != B.ncols  )
    matERROR.print("swap: matrices must be same size (for now).");
  #endif

  Type **pntA = A.data;
  Type **pntB = B.data;
  A.data = pntB;
  B.data = pntA;
  } // END myswap



/****************************************************************
 * A = sqrt(B)                                                  *
 * Bert Kampes, 16-Feb-1999                                     *
 ****************************************************************/
friend matrix<Type>  sqrt        (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("sqrt");
  #endif
  matrix<Type> Result(A.lines(),A.pixels());
  // ______Ensure stride one memory______
  Type *pntA   = A.data[0];
  Type *pntR   = Result.data[0];
  //for (register uint i=0; i<Result.nsize; i++)
  //  (*pntR++) = (sqrt((*pntA++)));
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<Result.nsize; i++)
    {
    (*pntR++) = (sqrt((*pntA)));
    *pntA++;
    }
  return Result;
  } // END sqrt



/****************************************************************
 * A = transpose(B)                                             *
 * Bert Kampes, 16-Oct-2005                                     *
 ****************************************************************/
friend matrix<Type>  transpose        (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("transpose (Bert Kampes, 16-Oct-2005)");
  #endif
  matrix<Type> Result(A.pixels(),A.lines());
  // not very fast implementation; better use pointer
  for (register uint i=0; i<Result.lines(); i++)
    for (register uint j=0; j<Result.pixels(); j++)
      Result(i,j) = A(j,i);
  return Result;
  } // END transpose



/****************************************************************
 * readfile                                                     *
 *    Type must be same as on disk. Result is filled.           *
 *    Windows must be specified in same system.                 *
 *    Either use win[0:N] and offset[0:N] (where 0 is the first *
 *    line), or use win[1:N] and offset[1:N] (where 1 indicates *
 *    the first line.                                           *
 *    Since the information in the resultfiles is written       *
 *    starting at line 1 (not 0) a 'currentwindow' can be read  *
 *    by a call like:                                           *
 *    readfile(mat,name,lines,winstartingatline1,currentwindow);*
 *    Result is filled at R(0:L,0:P), L=win.lines.              *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
friend void readfile(matrix<Type> &Result, const char *file,
  uint filelines, window win, window winoffset)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("readfile (window).");
  #endif

// ______First account for possible offsets of file______
// BK 18/1/00: winoffset starts at line 1?
  win.linelo = win.linelo - winoffset.linelo + 1;
  win.linehi = win.linehi - winoffset.linelo + 1;
  win.pixlo  = win.pixlo  - winoffset.pixlo  + 1;
  win.pixhi  = win.pixhi  - winoffset.pixlo  + 1;
  //win.pixhi  -= (winoffset.pixlo - 1);

  // ______ Check input ______
  // ifstream ifile(file, ios::in | ios::ate | ios::binary);
  // g++ seems to have a prob. with this... BK 130700
  // ifstream ifile(file, ios::in | ios::app | ios::binary);
  //ifstream ifile(file, ios::in | ios::binary | ios::nocreate);
  // new compiler v3.2 handles nocreate automatically fine 
#ifdef __NO_IOS_BINARY__
  ifstream ifile(file, ios::in);
#else
  ifstream ifile(file, ios::in | ios::binary);
#endif
  matassert(ifile,file,__FILE__,__LINE__);
  ifile.seekg(0,ios::end);                              // pointer to end...
  //const uint sizefile     = ifile.tellg();              // opened ate
  const streamoff &sizefile     = ifile.tellg();              // opened ate, [MA] 64-bit pointer
  const uint pixelsxbytes = sizefile/filelines;
  const uint filepixels   = pixelsxbytes/sizeof(Type);  // Type on disk.
  const uint sizepixel    = pixelsxbytes/filepixels;

  #ifdef __DEBUGMAT2
  if (win.linelo<=0 || win.pixlo<=0)
    matERROR.print("minimum line(pixel) is 0, (should be 1?).");
  if (win.linelo>win.linehi || win.pixlo>win.pixhi)
    matERROR.print("minimum line(pixel) is larger than max.");
  if (win.linehi>filelines || win.pixhi>filepixels)
    matERROR.print("max. line (pixel) is larger then on file.");
  if (sizeof(Type)!=sizepixel)
    matERROR.print("Type on disk is different than type of matrix.");
  if (sizefile==0)
    matERROR.print("filesize==0...");
  #endif

  const uint lines  = win.lines();
  const uint pixels = win.pixels();
  //INFO << "lines " << lines <<"  "  << pixels;
  //INFO.print();
  //const uint start  = ((win.linelo-1)*filepixels+win.pixlo-1)*sizepixel; [MA]
  const uint64 start  = (uint64)((win.linelo-1)*filepixels+win.pixlo-1)*sizepixel; // both sides should have the same type to
                                                                                  //  detect/eliminate integer overflow [MA]

  #ifdef __DEBUGMAT1
  if (Result.lines() < lines || Result.pixels() < pixels)
    matERROR.print("readfile: matrix too small to contain window from file.");
  if (Result.lines() != lines || Result.pixels() != pixels)
    matDEBUG.print("debug info: readfile: matrix is not fully filled.");
  #endif
  for (register uint lin=0; lin<lines; ++lin)
    {
    // read data at row: lin
   // INFO << "lin  " << lin<<endl;
   // INFO << "seek " << start+filepixels*lin*sizepixel<<endl;
   // INFO << "read " << pixels*sizepixel<<endl;
   // INFO.print();
    
    ifile.seekg(start+filepixels*lin*sizepixel,ios::beg);
    ifile.read((char*)&Result.data[lin][0],pixels*sizepixel);
    }
  ifile.close();
  } // END readfile



/****************************************************************
 * writefile                                                    *
 * write matrix partially to file.                              *
 * win 0 = first line (matrix index)                            *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
friend void          writefile   (
  ofstream &file, const matrix<Type> &tobwritten, window win)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("writefile.");
  #endif
  #ifdef __DEBUGMAT1
  if (win.linelo>win.linehi || win.pixlo>win.pixhi)
    matERROR.print("minimum line(pixel) is larger than max.");
  if (win.linehi>=tobwritten.lines() || win.pixhi>=tobwritten.pixels())
    matERROR.print("window not ok with matrix.");
  #endif
  matassert(file,"writefile: file unknown",__FILE__,__LINE__);
  const uint pixels = win.pixels();
  const uint size   = pixels*sizeof(Type);
  for (uint line=win.linelo; line<=win.linehi; ++line)
    file.write((char*)&tobwritten.data[line][win.pixlo],size);
  } // END writefile



/****************************************************************
 * fftshift(A)                                                  *
 *    fftshift of vector A is returned in A by reference        *
 *    shift DC term to middle. p=ceil(m/2); A=A[p:m-1 0:p-1];   *
 *    Bert Kampes, 24-Mar-2000                                  *
 ****************************************************************/
friend void          fftshift    (matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("fftshift");
  #endif
  #ifdef __DEBUGMAT1
  if (!A.isvector())
    matERROR.print("fftshift: only vectors");
  #endif
  matrix<Type> Res(A.nrows,A.ncols);
  const int32 start  = int32(ceil(real8(A.nsize)/2));
  memcpy(Res.data[0],A.data[0]+start,sizeof(Type)*(A.nsize-start));
  memcpy(Res.data[0]+A.nsize-start,A.data[0],sizeof(Type)*start);
  myswap(A,Res);                // prevent copy
  } // END fftshift



/****************************************************************
 * ifftshift(A)                                                 *
 *    ifftshift of vector A is returned in A by reference       *
 *    undo effect of fftshift. ?p=floor(m/2); A=A[p:m-1 0:p-1]; *
 *    Bert Kampes, 24-Mar-2000                                  *
 ****************************************************************/
friend void          ifftshift   (matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("ifftshift");
  #endif
  #ifdef __DEBUGMAT1
  if (!A.isvector())
    matERROR.print("ifftshift: only vectors");
  #endif
  matrix<Type> Res(A.nrows,A.ncols);
  const int32 start  = int32(floor(real8(A.nsize)/2));
  memcpy(Res.data[0],A.data[0]+start,sizeof(Type)*(A.nsize-start));
  memcpy(Res.data[0]+A.nsize-start,A.data[0],sizeof(Type)*start);
  myswap(A,Res);
  } // END ifftshift



/****************************************************************
 * wshift(A,n)                                                  *
 *    circular shift of vector A by n pixels. positive n for    *
 *    right to left shift.                                      *
 *    implementation: WSHIFT(A,n) == WSHIFT(A,n-sizeA);         *
 *    A is changed itself!                                      *
 *    Bert Kampes, 02-Nov-2000                                  *
 ****************************************************************/
friend void          wshift       (matrix<Type> &A, int32 n)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("wshift");
  #endif
  #ifdef __DEBUGMAT1
  if (abs(n)>=A.nsize)
    matERROR.print("wshift: shift larger than matrix not implemented.");
  if (!A.isvector())
    matERROR.print("wshift: only vectors");
  #endif
  // positive only, use rem!  n = n%A.nsize;
  if (n==0) return;
  if (n<0)  n += A.nsize;
  matrix<Type> Res(A.nrows,A.ncols);
  // ______ n always >0 here ______
  memcpy(Res.data[0],A.data[0]+n,sizeof(Type)*(A.nsize-n));
  memcpy(Res.data[0]+A.nsize-n,A.data[0],sizeof(Type)*n);
  myswap(A,Res);                // prevent copy
  } // END wshift


  }; // END matrix class




// ______ Compilation with g++ seems impossible any other way? ______
#include "matrixbk.cc"



#endif // MATRIXBK_H guard

