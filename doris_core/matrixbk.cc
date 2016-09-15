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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixbk.cc,v $   *
 * $Revision: 3.14 $                                                    *
 * $Date: 2005/10/18 14:16:32 $                                         *
 * $Author: kampes $                                                    *
 *                                                                      *
 *  template (base) class for matrices.                                 *
 *                                                                      *
 *  ====== Defines (compiler -D option) ======                          *
 *  __USE_VECLIB_LIBRARY__      for using VECLIB library                *
 *  __USE_LAPACK_LIBRARY__      for using LAPACK library                *
 *  __USE_FFTW_LIBRARY__        for using FFTW library                  *
 *  __DEBUGMAT1         index checking, dimensions etc.                 *
 *  __DEBUGMAT2         info, can savely be un-defined                  *
 *                                                                      *
 * Note that data is continuously in memory, because                    *
 *  VECLIB/FFTW assumes this (and me as well...) so a vector of vectors *
 *  cannot be used.                                                     *
 ************************************************************************/

#include "constants.hh"                 // typedefs, window

#include <iostream>                     // cout etc.
#include <fstream>                      // ofstream type
#include <strstream>                    // memory stream
#include <iomanip>                      // setw etc.
#include <algorithm>                    // max
#include <cstring>                      // memset according to g++ (?)
#include <complex>
#ifdef __DEBUGMAT1                      // use index checking, alloc
  #include <new>                        // bad_alloc
#endif

// ______ Keep track of total allocated memory ______
// ______ This will work upto 2^31 bytes (2GB) ______
// ______ there is a problem here, but how to declare a global? ______
// ______ if called from different files, not correct bookkeeping ______
// ______ therefor we define it above main in processor.cc ______
#ifdef __DEBUGMAT2                              // use index checking, alloc
  extern uint totalallocated;                   // [B]
#endif

// ______ message objects, global, set in main ______
extern bk_messages matERROR;
extern bk_messages matDEBUG;




// ====== Start of class definition ======
// ====== Private functions ======
/****************************************************************
 * allocate                                                     *
 *    allocate 2d array in major row order, continuously in     *
 *    memory (required by VECLIB/FFTW). consider rewriting to vector    *
 *    of vectors of stl, but may not be cont. in mem.           *
 *    g++ w/o exception handling, ...                           *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::allocate(uint numlines, uint numpixels)      // allocator
  {
  #ifdef __DEBUGMAT1
  if (numlines==0 || numpixels==0)
    {
    matERROR << "Allocation impossible: size (l,p): "
         << numlines << ", " << numpixels;
    matERROR.print();
    }
  #endif
  nrows = numlines;
  ncols = numpixels;
  nsize = numlines*numpixels;
  // Bert Kampes, 07-Apr-2005: try/catch should work by now...
  try
    {
    data = new Type*[numlines];// get memory : make linear array of pointers
    }
  catch (bad_alloc)
    {
    matERROR << "code 502: first allocation failed, size: "
         << numlines << ", " << numpixels;
    matERROR.print();
    }
  try
    {
    data[0] = new Type[nsize];// get memory : get the first pointer for the memory block, later we'll fill in all the address
    }
  catch(bad_alloc)
    {
    matERROR << "code 502: second allocation failed: size: "
         << numlines << ", " << numpixels;
    matERROR.print();
    }
  for (register uint i=1; i<numlines; i++)
    data[i] = data[i-1]+numpixels;              // start at 0,0

  #ifdef __DEBUGMAT2
  uint allocated = sizeof(Type)*nsize;          // [B]
  totalallocated += allocated;                  // [B]
  matDEBUG << "allocated   matrix("
       << numlines << "," << numpixels << ") at: " << &data[0][0] << " ("
       << setw(10) << allocated << "B, total: "
       << setw(10) << totalallocated << "B)";
  matDEBUG.print();
  #endif
  } // END allocate



/****************************************************************
 * initialize = allocate + set 0                                *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::initialize(uint numlines, uint numpixels)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("initialization matrix.");
  #endif
  allocate(numlines,numpixels);
  clean();
  } // END initialize



#ifdef __DEBUGMAT1
/****************************************************************
 * checkindex (l,p)                                             *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::checkindex(uint line, uint pixel) const
  {
  if (int32(line) > int32(nrows)-1 || int32(pixel) > int32(ncols)-1)
    {
    matERROR << "Wrong index (l,p)=(" << line << "," << pixel
         << "); Matrix(" << nrows << "," << ncols
         << ") at " << &data[0][0];
    matERROR.print();
    }
  } // END checkindex
#endif



// ====== Public functions ======
// ====== Constructors ======
/****************************************************************  \\
 * matrix<real8> A;                                             *  \\
 * Bert Kampes, 11-Dec-1998                                     *  \\
 ****************************************************************/
template <class Type>
matrix<Type>::matrix()                          // constructor (0 arg)
  {
  nrows = 0;
  ncols = 0;
  nsize = 0;
  data  = 0;                                    // address of pointer array is set to null
  } // END constructor



/****************************************************************
 * matrix<real8> A(3,3);                                        *
 * Bert Kampes, 11-Dec-1998                                     *
 ****************************************************************/
template <class Type>
matrix<Type>::matrix(uint lines, uint pixels)
  {
  initialize(lines,pixels);                             // set to 0
  } // END constructor



/****************************************************************
 * matrix<real8> A=B;                                           *
 * copy constructor; avoids default for dynamical memory:       * 
 * bitwise copy.                                                *
 * Bert Kampes, 11-Dec-1998                                     *
 ****************************************************************/
template <class Type>
matrix<Type>::matrix(const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("copy constructor.");  // [MA] ex: copy constructor is called when operator << is used.
  #endif
  if (A.nsize)
    {
    allocate(A.nrows,A.ncols);
    memcpy(data[0],A.data[0],nsize*sizeof(Type));
    }
  else
  {                                      // [MA] allow for NULL matrix returns when dynamic memory  nsize =0 
    allocate(A.nrows,A.ncols);
    #ifdef __DEBUGMAT2  // [MA] delete this set in future
    matDEBUG << "new mtx nsize: " << nsize << " datasize: " << nsize*sizeof(Type) << " dataaddress: " << &data << " vs inputadd: " << &A  << " addressNULLmemblock: " << data[0];
    matDEBUG.print();
    #endif
    //memcpy(data[0],0,nsize*sizeof(Type)); // clean() // when nsize=0 and A.data[0] doesn't exist (data=0), such as a null mtx, thus memcpy crashes with "Caught SIGSEGV: Segmentation fault."
                                          // this is necessary otherwise it points to a memory block
   data[0]=0;
   }
  } // END constructor



/****************************************************************
 * matrix<int32> A(win,B)                                       *
 * Constructor (part)                                           *
 *    Bert Kampes, 01-Feb-1999                                  *
 *    Bert Kampes, 31-Mar-1999 memcpy ipv. for_j                *
 * window starts at 0 (matrix index)                            *
 #%// BK 25-Oct-2000                                            *
 ****************************************************************/
template <class Type>
matrix<Type>::matrix (window win, const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("constructor as part.");
  #endif
  // ______ Check arguments ______
  #ifdef __DEBUGMAT1
  if (win.linehi<win.linelo)
    matERROR.print("constructor (4uint,matrix): win.linehi<linelo ?");
  if (win.pixhi<win.pixlo)
    matERROR.print("constructor (4uint,matrix): win.pixhi<pixlo ?");
  A.checkindex(win.linehi,win.pixhi);
  #endif
  // ______ Allocate new matrix and fill ______
  const uint numlin  = win.lines();
  const uint numpix  = win.pixels();
  const uint sizelin = numpix*sizeof(Type);
  allocate(numlin,numpix);
  for(register uint i=0; i<numlin; i++)
    memcpy(data[i],A[win.linelo+i]+win.pixlo,sizelin);
  } // END constructor



// ======Destructor======
/****************************************************************
 * Destructor                                                   *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type>::~matrix()
  {
  if (data==0) return;
  #ifdef __DEBUGMAT2
  uint deallocated  = sizeof(Type)*nsize;       // [B]
  totalallocated   -= deallocated;              // [B]
  matDEBUG << "deallocated matrix("
       << nrows << "," << ncols << ") at: " << data[0] << " ("
       << setw(10) << deallocated    << "B; total: "
       << setw(10) << totalallocated << "B)";
  matDEBUG.print();
  #endif

  delete [] data[0];                    // deallocate
  delete [] data;                       // deallocate
  data=0;                               // set to null pointer
  } // END destructor



// ======Data functions======
/****************************************************************
 * A.setdata(w)                                                 *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setdata(Type w)      // sets matrix to constant
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("A.setdata(w)");
  #endif
  Type *pnt = data[0];
  for (register uint i=0;i<nsize;i++)
    (*pnt++) = Type(w);
  } // END setdata



/****************************************************************
 * setdata(i,j,A)                                               *
 *    put matrix A on l1,p1                                     *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setdata(uint l1, uint p1, const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT1
  checkindex(l1+A.nrows-1,p1+A.ncols-1);// check far most corner 
  #endif
  const uint sizelin = A.ncols*sizeof(Type);
  for (register uint i=0;i<A.nrows;i++)
    memcpy(data[i+l1]+p1,A[i],sizelin);
  } // END setdata



/****************************************************************
 * B.setdata(win, A, winA):                                     *
 *  set win of B to winA of A                                   *
 * if win==0 defaults to totalB, winA==0 defaults to totalA     *
 * first line matrix =0 (?)
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setdata(window winin, const matrix<Type> &A, window winA)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("setdata (win,A,win)");
  #endif

  // ______Check default request______
  if (winin.linehi == 0 && winin.pixhi == 0)
    {winin.linehi = nrows  -1;
     winin.pixhi  = ncols -1;}
  if  (winA.linehi == 0 &&  winA.pixhi == 0)
    {winA.linehi  = A.lines()  -1;
     winA.pixhi   = A.pixels() -1;}
#ifdef __DEBUGMAT1
  if (((winin.linehi - winin.linelo) != (winA.linehi - winA.linelo)) ||
       ((winin.pixhi  - winin.pixlo) != (winA.pixhi  - winA.pixlo))    )
    matERROR.print("code 901: wrong input.");
  if (winin.linehi<winin.linelo || winin.pixhi<winin.pixlo)
    matERROR.print("code 901: wrong input.1");
  if ((winin.linehi > nrows  -1) ||
      (winin.pixhi  > ncols -1)   )
    matERROR.print("code 901: wrong input.2");
  if ((winA.linehi > A.lines()  -1) ||
      (winA.pixhi  > A.pixels() -1)   )
    matERROR.print("code 901: wrong input.3");
#endif
  // ______ Fill data ______
  const uint sizelin = (winA.pixhi - winA.pixlo + 1)*sizeof(Type);
  for(register uint i=winin.linelo; i<=winin.linehi;i++)
    memcpy(data[i]+winin.pixlo,A[i-winin.linelo+winA.linelo]+winA.pixlo,sizelin);
  } // END setdata



/****************************************************************
 * B.setdata(A, winA):                                          *
 *    set total of B to winA of A                               *
 *    Bert Kampes, 17-Mar-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setdata(const matrix<Type> &A, window winA)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("setdata (A,win)");
  #endif
  #ifdef __DEBUGMAT1
  if (((nrows  -1) != (winA.linehi - winA.linelo)) ||
       ((ncols -1) != (winA.pixhi  - winA.pixlo))    )
    matERROR.print("code 901: wrong input.");
  if ((winA.linehi > A.lines()  -1) ||
      (winA.pixhi  > A.pixels() -1)   )
    matERROR.print("code 901: wrong input.3");
  #endif
  // ______ Fill data ______
  const uint sizelin = ncols*sizeof(Type);
  for(register uint i=0; i<nrows;i++)
    memcpy(data[i],A[i+winA.linelo]+winA.pixlo,sizelin);
  } // END setdata



/****************************************************************
 * A=B.getdata(win)                                             *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> matrix<Type>::getdata(window win) const
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("getdata.");
  #endif
  #ifdef __DEBUGMAT1
  if (win.linehi<win.linelo || win.pixhi<win.pixlo)
    matERROR.print(
    "code 501: matrix::getdata (win): arguments are wrong, l1<l2,p1<p2");
  checkindex(win.linehi,win.pixhi);
  #endif
  const uint numlin = win.lines();
  const uint numpix = win.pixels();
  matrix<Type> Result(numlin,numpix);                   // =data(;
  for(register uint i=0; i<numlin; i++)
    memcpy(Result[i],data[i+win.linelo]+win.pixlo,numpix*sizeof(Type));
  return Result;
  } // END getdata



/****************************************************************
 * A=B.getrow(row)                                              *
 * rows: 0 1 2 ..                                               *       
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> matrix<Type>::getrow(uint line) const
  {
  #ifdef __DEBUGMAT1
  checkindex(line,0);
  #endif
  matrix<Type> Result(1,ncols);
  memcpy(Result[0],data[line],ncols*sizeof(Type));
  return Result;
  } // END getrow



/****************************************************************
 * A=B.getcolumn(col)                                           *
 * cols: 0 1 2 ..                                               *       
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> matrix<Type>::getcolumn(uint pixel) const
  {
  #ifdef __DEBUGMAT1
    checkindex(0,pixel);
  #endif
  matrix<Type> Result(nrows,1);
  Type *pntA = data[0]+pixel;
  Type *pntR = Result[0];
  for (register uint i=0; i<nrows; ++i)
    {
    (*pntR++) = *pntA;
    pntA     += ncols;
    }
  return Result;
  } // END getcolumn



/****************************************************************
 * B.showdata()                                                 *
 *    Bert Kampes, 01-Feb-1999                                  *
 *    Mahmut Arikan, 22-May-2009 - Adjustment to see more digits*
 ****************************************************************/
template <class Type>
void matrix<Type>::showdata() const                     // show all data in matrix
  {
  #ifdef __DEBUGMAT1
    matDEBUG.print("showdata.");
  #endif
  #ifdef __DEBUGMAT2
  if (nrows>100 || ncols>15)
    {
    matDEBUG << "matrix ("
         << nrows << "," << ncols
         << "); only showing data (0:99,0:9).";
    matDEBUG.print();
    }
  #endif
  const uint L = (nrows<=100) ? nrows : 100;
  const uint P = (ncols<=10)  ? ncols : 10;
  matDEBUG.precision(11);  // [MA] 9 --> 11
  matDEBUG.width(12);      // 10 --> 12
  //matDEBUG.setf(ios::left, ios::adjustfield); // [MA] this failed by 4.01 version of messages.hh
  for (register uint i=0; i<L; i++)
    {
    for (register uint j=0; j<P; j++)
      {
      // matDEBUG << data[i][j] << " ";
      matDEBUG << left << data[i][j] << " ";
      }
    matDEBUG.print();// prevent too much in ostream
    }
  matDEBUG.reset();
  } // END showdata



/****************************************************************
 * bool v = A.isvector()                                        *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
inline
bool matrix<Type>::isvector() const
  {
  return (nrows==1 || ncols==1) ? true : false;
  } // END isvector



/****************************************************************
 * uint l = A.lines()                                           *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
inline
uint matrix<Type>::lines() const                // return number of lines
  {
  return nrows;
  } // END lines



/****************************************************************
 * uint p = A.pixels()                                          *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
inline
uint matrix<Type>::pixels() const               // return number of pixels
  {
  return ncols;
  } // END pixels



/****************************************************************
 * uint s = A.size()                                            *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
inline
uint matrix<Type>::size() const                         // return nsize
  {
  return nsize;
  } // END size



/****************************************************************
 * A.resize(l,p)                                                *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::resize(uint l1, uint p1)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("resize.");
  #endif
  if (l1 == nrows && p1 == ncols) return;
  else if (data!=0)                             // check for allocated memory
    {
    #ifdef __DEBUGMAT2
    uint deallocated  = sizeof(Type)*nsize;       // [B]
    totalallocated   -= deallocated;            // [B]
    matDEBUG << "deallocated matrix("
             << nrows << "," << ncols << ") at: " << data[0] << " ("
             << setw(10) << deallocated << "B; total: "
             << setw(10) << totalallocated << "B)";
    matDEBUG.print();
    #endif
    delete [] data[0];
    delete [] data;
    data=0;
    }
  initialize(l1,p1);                            // set to 0
  } // END resize



/****************************************************************
 * A.clean()                                                    *
 *    Bert Kampes, 01-Feb-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::clean()              // sets 2 zero
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mtx clean.");
  #endif
  memset(data[0],0,nsize*sizeof(Type));
  } // END clean



/****************************************************************
 * B.setrow(row, data)                                          *
 * rows: 0 1 2 .., should fit exactly                           *       
 * orientation of vector is disregarded.                        *
 *    Bert Kampes, 12-Oct-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setrow(uint line, const matrix<Type> &LINE)
  {
  #ifdef __DEBUGMAT1
  checkindex(line,0);
  if (!(LINE.nrows == 1 || LINE.ncols == 1))
    matERROR.print("setrow: only vector input.");
  if (LINE.nsize != ncols)
    matERROR.print("setrow: sizeofvector should be same as matrix.");
  #endif
  memcpy(data[line],LINE[0],ncols*sizeof(Type));
  } // END setrow



/****************************************************************
 * B.setrow(row, scalar)                                        *
 * rows: 0 1 2 .., should fit exactly                           *       
 * orientation of vector is disregarded.                        *
 *    Bert Kampes, 12-Oct-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setrow(uint line, Type scalar)
  {
  #ifdef __DEBUGMAT1
  checkindex(line,0);
  #endif
  for (register Type *pntB =&data[line][0];
                      pntB<=&data[line][ncols-1];
                      pntB++)
    *pntB = scalar;
  } // END setrow



/****************************************************************
 * B.setcolumn(col, COL)                                        *
 * cols: 0 1 2 ..                                               *       
 * orientation of vector is disregarded.                        *
 *    Bert Kampes, 12-Oct-1999                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::setcolumn(uint pixel, const matrix<Type> &COLUMN)
  {
  #ifdef __DEBUGMAT1
  checkindex(0,pixel);
  if (!(COLUMN.nrows == 1 || COLUMN.ncols == 1))
    matERROR.print("setcolumn: only vector input.");
  if (COLUMN.nsize != nrows)
    matERROR.print("setcolumn: sizeofvector should be same as matrix.");
  #endif
  Type *pntCOL = COLUMN[0];
  for (register Type *pntB =&data[0][pixel];
                      pntB<=&data[nrows-1][pixel];
                      pntB+=ncols)
    *pntB = *pntCOL++;
  } // END setcolumn



/****************************************************************
 * B.setcolumn(col, scalar)                                     *
 * cols: 0 1 2 ..                                               *       
 #%// BK 25-Sep-2000                                            *
 ****************************************************************/
template <class Type>
void matrix<Type>::setcolumn(uint pixel, Type scalar)
  {
  #ifdef __DEBUGMAT1
  checkindex(0,pixel);
  #endif
  for (register Type *pntB =&data[0][pixel];
                      pntB<=&data[nrows-1][pixel];
                      pntB+=ncols)
    *pntB = scalar;
  } // END setcolumn



/****************************************************************
 * B.fliplr()                                                   *
 * Mirror in center vertical (flip left right).                 *
 *    Bert Kampes, 23-Mar-2000                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::fliplr()
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("fliplr.");
  #endif
  if (nrows==1)
    {
    Type *pnt1 =  data[0];                              // first one
    Type *pnt2 =  data[0]+ncols-1;                      // last one
    Type tmp   = *pnt1;
    for (register int32 i=0; i<int32(ncols/2); ++i)     // floor
      {
      (*pnt1++) = *pnt2;
      (*pnt2--) =  tmp;
      tmp       = *pnt1;
      }
    }
  else
    {
    for (register int32 i=0; i<int32(ncols/2); ++i)     // floor
      {
      matrix<Type> tmp1 = getcolumn(i);
      matrix<Type> tmp2 = getcolumn(ncols-i-1);
      setcolumn(i,tmp2);
      setcolumn(ncols-i-1,tmp1);
      }
    }
  } // END fliplr



/****************************************************************
 * B.flipud()                                                   *
 * Mirror in center vertical (flip left right).                 *
 *    Actually move data around, not pointers, to be sure       *
 *    veclib works ok, data is cont. in memory.                 *
 *    Bert Kampes, 23-Mar-2000                                  *
 ****************************************************************/
template <class Type>
void matrix<Type>::flipud()
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("flipud.");
  #endif
  if (ncols==1)
    {
    Type *pnt1 =  data[0];                              // first one
    Type *pnt2 =  data[0]+nrows-1;                      // last one
    Type tmp   = *pnt1;
    for (register int32 i=0; i<int32(nrows/2); ++i)     // floor
      {
      (*pnt1++) = *pnt2;
      (*pnt2--) =  tmp;
      tmp       = *pnt1;
      }
    }
  else
    {
    for (register int32 i=0; i<int32(ncols/2); ++i)     // floor
      {
      matrix<Type> tmp1 = getrow(i);
      matrix<Type> tmp2 = getrow(nrows-i-1);
      setrow(i,tmp2);
      setrow(nrows-i-1,tmp1);
      }
    }
  } // END flipud







/****************************************************************
 ****************************************************************
 * OVERLOADED OPS                                               *
 ****************************************************************
 ****************************************************************/




/****************************************************************
 * a = A[5][2];                                                 *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
inline
Type* matrix<Type>::operator [] (uint line) const
  {
  #ifdef __DEBUGMAT1
  checkindex(line,0);
  #endif
  return data[line];
  } // END []



/****************************************************************
 * a = A(i,j);                                                  *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
inline
Type& matrix<Type>::operator () (uint line, uint pixel) const
  {
  #ifdef __DEBUGMAT1
  checkindex(line,pixel);
  #endif
  return data[line][pixel];
  } // END ()



/****************************************************************
 * matrix<T> B = A(window);                                     *
 * may not be to efficient cause twice allocation?              *
 * Bert Kampes, 31-Mar-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> matrix<Type>::operator () (const window &win) const
  {
  matrix<Type> Res(win,*this);
  return Res;
  } // END (win)



/****************************************************************
 * matrix<T> B = A(uint,uint,uint,uint);                        *
 * Bert Kampes, 31-Mar-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> matrix<Type>::operator () (
        const uint &l0, const uint &lN,
        const uint &p0, const uint &pN) const
  {
  const window win(l0,lN,p0,pN);
  matrix<Type> Res(win,*this);
  return Res;
  } // END (4 uint)



/****************************************************************
 *  =                                                           *
 * Bert Kampes, 14-Jan-1999                                     *
 * Mahmut Arikan, 19-May-2009 Null matrix  assignment           *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator = (const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("operator =");
  #endif
  //cerr << "operator =\n";
  if (this != &A)                               // prevent copy to itself
    {
    if (A.nsize)                                // if allocated : disable this to return null matrix 
      {
      if (data != 0)                            // if allocated
        {
        #ifdef __DEBUGMAT2
        uint deallocated  = sizeof(Type)*nsize; // [B]
        totalallocated   -= deallocated;        // [B]
        matDEBUG << "deallocated matrix("
             << nrows << "," << ncols << ") at: " << data[0] << " ("
             << setw(10) << deallocated << "B; total: "
             << setw(10) << totalallocated << "B)";
        matDEBUG.print();
        #endif
        delete [] data[0];
        delete [] data;
        data = 0;
        }
      allocate(A.nrows,A.ncols);
      memcpy(data[0],A.data[0],nsize*sizeof(Type));
      }
    else
      {                                          // [MA] allow for NULL matrix returns when dynamic memory  nsize =0 
        cerr << "op =    : NULL matrix is assigned.";
        allocate(A.nrows,A.ncols);
        data[0]=0;  // point to no memory block
       }
    }
  return *this;
  } // END =


/****************************************************************
 *  =                                                           *
 #%// BK 09-Nov-2000                                            *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator = (const Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("operator = (scalar)");
  #endif
  setdata(scalar);
  return *this;
  } // END = (scalar)




/****************************************************************
 * C *= 5.0;                                                    *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator *= (Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("*=");
  #endif
  #ifdef __DEBUGMAT1
  if (!nsize)
    matERROR.print("matrix:: *= with empty matrix.");
  #endif
  Type *pntmat = data[0];
  for (register uint i=0; i<nsize; ++i)
    (*pntmat++) *= scalar;      
  return *this;
  } // END *= scalar




/****************************************************************
 * C *= A;      pointwise multiplication, a,c same size         *
 * Bert Kampes, 06-Oct-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator *= (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("*= pointwise");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: *= matrices must be same size.");
  #endif
  Type *pntmat = data[0];
  Type *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) *= (*pntA++); 
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) *= (*pntA);
    *pntA++;
    }
  return *this;
  } // END *= matrices




/****************************************************************
 * C /= 5.0;                                                    *
 * Bert Kampes, 12-Oct-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator /= (Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("/= scalar");
  #endif
  #ifdef __DEBUGMAT1
  if (!nsize)
    matERROR.print("matrix:: /= with empty matrix.");
  #endif
    Type *pntmat = data[0];
    for (register uint i=0; i<nsize; i++)
      (*pntmat++) /= scalar;    
    return *this;
  } // END /= scalar




/****************************************************************
 * C /= A;      pointwise division, a,c same size               *
 * Bert Kampes, 06-Oct-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator /= (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("/=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.lines() || ncols != A.pixels())
    matERROR.print("matrix:: /= matrices must be same size.");
  #endif
  Type *pntmat = data[0];
  Type *pntA   = A[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) /= (*pntA++); 
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) /= (*pntA);
    *pntA++;
    }
  return *this;
  } // END /= matrices



/****************************************************************
 * C -= A;                                                      *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator -= (const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("-= mat");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.nrows || ncols != A.ncols)
    matERROR.print("error dimensions.");
  if (!A.nsize)
    matERROR.print("matrix:: -= with empty matrices.");
  #endif
  Type *pntmat = data[0];
  Type *pntA   = A.data[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) -= (*pntA++); 
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) -= (*pntA);     
    *pntA++;
    }
  return *this;
  } // END -=



/****************************************************************
 * C -= 5.0;                                                    *
 * Bert Kampes, 26-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator -= (Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("-= scalar");
  #endif
  #ifdef __DEBUGMAT1
  if (!nsize)
    matERROR.print("matrix:: -= with empty matrix.");
  #endif
  Type *pntmat = data[0];
  for (register uint i=0; i<nsize; i++)
    (*pntmat++) -= scalar;      
  return *this;
  } // END -= scalar




/****************************************************************
 * C += A;                                                      *
 * Bert Kampes, 14-Jan-1999                                     *
 #%// Bert Kampes, 10-Apr-2005  (why const?)
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator += (const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("+=");
  #endif
  #ifdef __DEBUGMAT1
  if (nrows != A.nrows || ncols != A.ncols)
    matERROR.print("error dimensions.");
  if (!A.nsize)
    matERROR.print("matrix:: += with empty matrices.");
  #endif
  Type *pntmat = data[0];
  Type *pntA   = A.data[0];
  //for (register uint i=0; i<nsize; i++)
  //  (*pntmat++) += (*pntA++); 
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<nsize; i++)
    {
    (*pntmat++) += (*pntA);
    *pntA++;
    }
  return *this;
  } // END +=




/****************************************************************
 * C += 5.0;                                                    *
 * Bert Kampes, 26-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>& matrix<Type>::operator += (Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("+= scalar");
  #endif
  #ifdef __DEBUGMAT1
  if (!nsize)
    matERROR.print("matrix:: += with empty matrix.");
  #endif
  Type *pntmat = data[0];
  for (register uint i=0; i<nsize; i++)
    (*pntmat++) += scalar;      
  return *this;
  } // END += scalar



/****************************************************************
 * A.conj();complex conjugated                                  *
 * Bert Kampes, 18-Oct-1999                                     *
 ****************************************************************/
template <class Type>
void matrix<Type>::conj()
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("conj"); 
  #endif
  Type *pntmat = data[0];
  for (register uint i=0; i<nsize; ++i)
    {
    (*pntmat) = Type((*pntmat).real(), -(*pntmat).imag());
    pntmat++;
    }
  } // END conj()



/****************************************************************
 * if (A==scalar) all elements equal scalar                     *
 * Bert Kampes, 08-Oct-1999                                     *
 ****************************************************************/
template <class Type>
bool matrix<Type>::operator == (Type scalar) const
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("== (scalar)");
  #endif
  if (nsize == 0)
    return false;
  bool same = true;
  Type *pnt = data[0];
  for (register int32 i=0; i<nsize; ++i)
    {
    if ((*pnt) != scalar)
      {
      same = false;
      break;
      }
    }
  return same;
  } // END ==



/****************************************************************
 * if (A==B)                                                    *
 * Bert Kampes, 08-Oct-1999                                     *
 ****************************************************************/
template <class Type>
bool matrix<Type>::operator == (const matrix<Type> &A) const
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("==");
#endif
  if ((A.lines() != nrows) || (A.pixels() != ncols))
    return false;
  bool same  = true;
  Type *pnt  = data[0];
  Type *pntA = A[0];
  for (register uint i=0; i<nsize; ++i)
    {
    if ((*pnt) != (*pntA))
      {
      same = false;
      break;
      }
    }
  return same;
  } // END ==



/****************************************************************
 * if (A!=scalar)                                               *
 * Bert Kampes, 08-Oct-1999                                     *
 ****************************************************************/
template <class Type>
bool matrix<Type>::operator != (Type scalar) const
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("!= (scalar)");
#endif
  if (*this == scalar)
    return false;
  else
    return true;
  } // END !=



/****************************************************************
 * if (A!=B)                                                    *
 * Bert Kampes, 08-Oct-1999                                     *
 ****************************************************************/
template <class Type>
bool matrix<Type>::operator != (const matrix<Type> &A) const
  {
#ifdef __DEBUGMAT2
  matDEBUG.print("!=");
#endif
  if (*this == A)
    return false;
  else
    return true;
  } // END !=




// ++++++++++++++++++++++++++
// template functions, stupid place, but else not found?
// specializations in matrixspecs.cc are NOT used before these?
// BK 16-Aug-2000
// ++++++++++++++++++++++++++
/****************************************************************
 * C = A * B;                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator * (const matrix<Type>& A, const matrix<Type>& B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("matrices * (no veclib)");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.lines())
    matERROR.print("matrix::operator *: multiplication not possible");
  if (!A.size())
    matERROR.print("matrix:: operator * with empty matrices.");
  #endif
  // ______ Straightforward, no veclib _____
  // ______ it is probably worth to make this faster, e.g.,
  // ______ using pointers (using A[i][j] is same speed)
  // ______ and not initializing Result
  matrix<Type>  Result(A.lines(),B.pixels());
  register Type sum = Type(0.0);
  // --- use straightforward notation for slow matrix access --------
  #define NO_POINTERS // this is likely a bit slower
  #ifdef NO_POINTERS // straightforward notation
  for (register uint i=0; i<Result.lines(); i++) 
    {
    for (register uint j=0; j<Result.pixels(); j++) 
      {
      for (register uint k=0; k<A.pixels(); k++) 
        {
        sum += A(i,k) * B(k,j); 
        }
      Result(i,j) = sum; 
      sum         = Type(0.0);// complex requires this
      }
    }
  // --- use pointers for faster matrix access --------------------
  // Bert Kampes, 13-Oct-2005: tested OK.
  #else  // use pointers for faster matrix access
  Type *pntR = Result[0];// use to point from first to last element 
  for (register uint i=0; i<Result.lines(); i++) 
    {
    for (register uint j=0; j<Result.pixels(); j++) 
      {
      Type *pntA = A[i];// point to first element this row of A
      Type *pntB = B[0]+j;// point to first element in column of B
      for (register uint k=0; k<A.pixels(); k++) 
        {
        //sum  += (*pntA++) * (*pntB); // pointer over row A
        // changed by FvL (for g++/gcc > 4.0):
        sum  += (*pntA) * (*pntB); // pointer over row A
        *pntA++;
        pntB += B.pixels();// point to next element in this column of B
        }
      (*pntR++) = sum;// matrix is linear in memory, major-row
      sum       = Type(0.0);// complex requires this
      }
    }
  #endif
  return Result; 
  } // END *



/****************************************************************
 * C = 5.0 * B;                                                 *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator * (const matrix<Type>& A, Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("scalar *");
  #endif
  // ______Perform multiplication______
  matrix<Type>  Result=A;
  return Result *= scalar;              // checks in *=
  } // END * scalar



/****************************************************************
 * C = B * 5.0;                                                 *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator *  (Type  scalar, const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("* scalar.");
  #endif
  return A*scalar;                      // checks in *=
  } // END scalar *



/****************************************************************
 * C = A / 5;                                                   *
 * Bert Kampes, 04-Apr-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator / (const matrix<Type>& A, Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("/");
  #endif
  matrix<Type> Result = A;
  return Result      /= scalar;         // checks are performed here.
  } // END /



/****************************************************************
 * C = A / B;                                                   *
 * Bert Kampes, 04-Apr-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator / (const matrix<Type>& A, const matrix<Type>& B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("/");
  #endif
  matrix<Type> Result = A;
  return Result      /= B;      // checks are performed here.
  } // END /



/****************************************************************
 * C = A - B;                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator - (const matrix<Type>& A, const matrix<Type>& B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("-");
  #endif
  matrix<Type> Result = A;
  return Result      -= B;              // checks are performed here.
  } // END - (binary)



/****************************************************************
 * C = A - 5;                                                   *
 * Bert Kampes, 04-Apr-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator - (const matrix<Type>& A, Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("-");
  #endif
  matrix<Type> Result = A;
  return Result      -= scalar;         // checks are performed here.
  } // END - (binary)



/****************************************************************
 * C = A + B;                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator + (const matrix<Type>& A, const matrix<Type>& B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("+");
  #endif
  matrix<Type> Result = B;
  return Result      += A;              // checks are in +=
  } // END + (binary)



/****************************************************************
 * C = A + 5;                                                   *
 * Bert Kampes, 04-Apr-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type> operator + (const matrix<Type>& A, Type scalar)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("+");
  #endif
  matrix<Type> Result = A;
  return Result      += scalar;          // checks are performed here.
  } // END + (binary)



/****************************************************************
 * a = max(A)                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
Type max(const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("max");
  #endif
  Type m=A(0,0);
  for (register uint i=0; i<A.lines(); ++i)
    for (register uint j=0; j<A.pixels(); ++j)
      if (A(i,j)>m) m=A(i,j);
  return m;     
  } // END max



/****************************************************************
 * a = max(A,linemax,pixelmax)                                  *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
Type max(const matrix<Type> &A, uint& line, uint& pixel)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("max");
  #endif
  Type m=A(0,0);
  for (register uint i=0; i<A.lines(); ++i)
    for (register uint j=0; j<A.pixels(); ++j)
      if (A(i,j)>=m)
        {
        m     = A(i,j);
        line  = i;
        pixel = j;
        }
  return m;     
  } // END max



/****************************************************************
 * a = min(A)                                                   *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
Type min(const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("min");
  #endif
  Type m=A(0,0);
  for (register uint i=0; i<A.lines(); ++i)
    for (register uint j=0; j<A.pixels(); ++j)
      if (A(i,j)<m) m=A(i,j);
  return m;     
  } // END min



/****************************************************************
 * a = min(A,linemax,pixelmax)                                  *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
Type min(const matrix<Type> &A, uint& line, uint& pixel)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("min");
  #endif
  Type m=A(0,0);
  for (register int32 i=0; i<A.lines(); i++)
    for (register int32 j=0; j<A.pixels(); j++)
      if (A(i,j)<=m) 
        {
        m     = A(i,j);
        line  = i;
        pixel = j;
        }
  return m;     
  } // END min



/****************************************************************
 * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> matTxmat(const matrix<Type> &A, const matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matTxmat: no veclib");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines() != B.lines())
    matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
  #endif
  matrix<Type> Result(A.pixels(),B.pixels());
  register Type sum = Type(0.0);
  // --- use straightforward notation for slow matrix access --------
  #define NO_POINTERS // this is likely a bit slower
  #ifdef NO_POINTERS // straightforward notation
  for (register uint i=0; i<Result.lines(); i++)
    {
    for (register uint j=0; j<Result.pixels(); j++)
      {
      for (register uint k=0; k<A.lines(); k++)
        {
        sum += A(k,i) * B(k,j);
        }
      Result(i,j) = sum;
      sum         = Type(0.0);
      }
    }
  // --- use pointers for faster matrix access --------------------
  #else
  Type *pntR = Result[0];// use to point from first to last element 
  for (register uint i=0; i<Result.lines(); i++) 
    {
    for (register uint j=0; j<Result.pixels(); j++) 
      {
      Type *pntA = A[0]+i;// point to first element this row of A
      Type *pntB = B[0]+j;// point to first element in column of B
      for (register uint k=0; k<A.lines(); k++) 
        {
        sum  += (*pntA) * (*pntB); // pointer over row A
        pntA += A.pixels();
        pntB += B.pixels();// point to next element in this column of B
        }
      (*pntR++) = sum;// matrix is linear in memory, major-row
      sum       = Type(0.0);// complex requires this
      }
    }
  #endif
  return Result;
  } // END matTxmat



/****************************************************************
 * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> matxmatT(const matrix<Type> &A, const matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
    matDEBUG.print("matxmatT: no veclib");
  #endif
  #ifdef __DEBUGMAT1
  if (A.pixels() != B.pixels())
    matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
  #endif
  register Type sum = Type(0.0);
  matrix<Type> Result(A.lines(),B.lines());
  // --- use straightforward notation for slow matrix access --------
  #define NO_POINTERS // this is likely a bit slower
  #ifdef NO_POINTERS // straightforward notation
  for (register uint i=0; i<Result.lines(); i++)
    {
    for (register uint j=0; j<Result.pixels(); j++)
      {
      for (register uint k=0; k<A.pixels(); k++)
        {
        sum += A(i,k) * B(j,k);
        }
      Result(i,j) = sum;
      sum         = Type(0.0);
      }
    }
  // --- use pointers for faster matrix access --------------------
  #else
  Type *pntR = Result[0];// use to point from first to last element 
  for (register uint i=0; i<Result.lines(); i++) 
    {
    for (register uint j=0; j<Result.pixels(); j++) 
      {
      Type *pntA = A[i];// point to first element this row of A
      Type *pntB = B[j];// point to first element in column of B
      for (register uint k=0; k<A.pixels(); k++) 
        {
        //sum  += (*pntA++) * (*pntB++); // pointer over row A
        // changed by FvL (for g++/gcc > 4.0):
        sum  += (*pntA) * (*pntB); // pointer over row A
        *pntA++;
        *pntB++;
        }
      (*pntR++) = sum;// matrix is linear in memory, major-row
      sum       = Type(0.0);// complex requires this
      }
    }
  #endif
  return Result;
  } // END matxmatT



/****************************************************************
 * dumpasc(file,A);                                             *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
void dumpasc(const char *file, const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("dumpasc to file");
  #endif
  ofstream fo(file,ios::out | ios::trunc);
  matassert(fo,file,__FILE__,__LINE__);
  fo.precision(3);
  fo.width(11);
  fo.setf(ios::fixed);
  for (register int32 i=0; i<A.lines(); ++i)
    {
    for (register int32 j=0; j<A.pixels(); ++j)
      {
      fo << A(i,j) << " ";
      }
    fo << endl;
    }
  fo.close();
  } // END dumpasc




/****************************************************************
 * C = dotmult(A,B) = A .* B                                    *
 * Bert Kampes, 26-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  dotmult     (const matrix<Type> &A, const matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("dotmult");
  #endif
  matrix<Type> Result = A;
  Result             *= B;                      // checks are here
  return Result;
  } // END dotmult



/****************************************************************
 * C = dotdiv(A,B) = A/B                                        *
 * Bert Kampes, 26-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  dotdiv      (const matrix<Type> &A, const matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("dotdiv");
  #endif
  matrix<Type> Result = A;
  Result             /= B;                      // checks are here
  return Result;
  } // END dotdiv



/****************************************************************
 * A = sqr(B)                                                   *
 * Bert Kampes, 16-Feb-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  sqr         (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("sqr");
  #endif
  matrix<Type> Result = dotmult(A,A);
  return Result;
  } // END sqr



/****************************************************************
 *    B = conj(A)                                               *
 *    Bert Kampes, 02-Mar-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type> conj       (const matrix<Type> &A)
  {
  matrix<Type> Result = A;
  Result.conj();
  return Result;
  } // END conj



/****************************************************************
 * C=diagxmat(vec,B) C=diag(vec) * B;                           *
 *    Bert Kampes, 22-Feb-1999                                  *
 ****************************************************************/
template <class Type>
matrix<Type>  diagxmat    (const matrix<Type> &diag, const matrix<Type> &B)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("diagxmat");
  #endif
  #ifdef __DEBUGMAT1
  if (min(diag.lines(),diag.pixels()) != 1)
    matERROR.print("diagxmat: sizes A,B: diag is vector.");
  if (diag.size() != B.lines())
    matERROR.print("diagxmat: sizes A,B: input is vector, matrix.");
  #endif
 
  matrix<Type> Result=B;
  if (diag.lines() != 1)        // standing
    {
    for (register int32 i=0; i<int32(Result.lines()); i++)
      for (register int32 j=0; j<int32(Result.pixels()); j++)
        Result(i,j) *= diag(i,0);
    }
  else
    {
    for (register int32 i=0; i<int32(Result.lines()); i++)
      for (register int32 j=0; j<int32(Result.pixels()); j++)
        Result(i,j) *= diag(0,i);
    }
  return Result;
  } // END diagxmat



/****************************************************************
 * multilook A with factors l,p                                 *
 * lines(A) have to be multiple of factorL.                     *
 * size of output is (Al/fl,Ap/fp)                              *
 *  multilooked is averaged by factorLP.                        *
 * Bert Kampes, 19-Apr-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  multilook   (const matrix<Type> &A, uint factorL, uint factorP)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("multilook.");
    if (A.lines()%factorL)                                 // [MA] we can handle this 
      //matERROR.print("lines A must be multiple of factorL.");
      matDEBUG.print("For this buffer, lines A is not multiple of factorL.");
  #endif

  if (factorL==1 && factorP==1)
    {
    //matrix<Type> R=A;
    //return R;
    return A;        // [MA]: fastest solution
    }
 
  #ifdef __DEBUGMAT2
    matDEBUG <<  "multilook input [A] size: " << A.size() << " lines: " << A.lines() << " pixels: " << A.pixels()  << " address: " << &A ; // A.[0] can't be printed if 0
    matDEBUG.print();
  #endif

  if ( A.lines()/factorL == 0 || A.pixels()/factorP == 0 ) // [MA] fix for extra buffer when lines < mlfactor or ...
    {
    DEBUG.print("Multilooking was not necessary for this buffer: buffer.lines() < mlL or buffer.pixels < mlP");
    matrix<Type> R; //=A; // see initialize()
    //R.resize(1,1); // fill with 0 
    #ifdef __DEBUGMAT2
      matDEBUG <<  "multilook return [R] size: " << R.size() << " lines: " << R.lines() << " pixels: " << R.pixels() << " address: " << &R << endl;
    matDEBUG.print();
    #endif
    return R; // NULL 
    }
 
  Type sum;
  Type factorLP = Type(factorL * factorP);
  //cerr << "multilook: "<< A.lines()/factorL << " " << A.pixels()/factorP << endl;
  matrix<Type> Result(A.lines()/factorL,A.pixels()/factorP);
  for (register uint i=0; i<Result.lines(); i++)
    {
    for (register uint j=0; j<Result.pixels(); j++)
      {
      sum = Type(0.0);
      for (register uint k=i*factorL; k<(i+1)*factorL; k++)
        {
        for (register uint l=j*factorP; l<(j+1)*factorP; l++)
          {
          sum += A(k,l);
          }
        }
      Result(i,j) = sum/factorLP;
      }
    }
  //cerr << "multilook: l,p " << Result.lines() << " " << Result.pixels() << " size " << Result.size() << endl;
  return Result;
  } // END multilook



/****************************************************************
 * Correlate A with maskB, return C(sizeA)                      *
 *  egde is set to zero, not pad with zeros                     *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<real4> correlate   (const matrix<Type> &A, matrix<Type> Mask)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("Correlate.");
  matDEBUG.print("not yet correct for complex?: returns real4");
  if (Mask.lines()<2 || Mask.pixels()<2)
    matERROR.print("very small mask.");
  #endif
  #ifdef __DEBUGMAT1
  if (A.lines()<Mask.lines() || A.pixels()<Mask.pixels())
    matERROR.print("matrix input smaller than mask.");
  #endif

  real8 varM = 0.;                                // variance of Mask
  Mask      -= mean(Mask);
  Type *pntMsk = Mask[0];
  //for (register uint ii=0; ii<Mask.size(); ii++)
  //  varM += sqr(*pntMsk++);                         // 1/N later
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint ii=0; ii<Mask.size(); ii++)
    {
    varM += sqr(*pntMsk);                         // 1/N later
    *pntMsk++;
    }

// ______Compute correlation at these points______
  uint beginl = (Mask.lines()-1)  / 2 ;             // floor
  uint beginp = (Mask.pixels()-1) / 2 ;             // floor
  matrix<real4> Result(A.lines(),A.pixels());       // init to 0

// ______First window of A, updated at end of loop______
  window winA   (0, Mask.lines()-1, 0, Mask.pixels()-1);
  window windef (0,0,0,0);// defaults to total Am

// ______Correlate part of Result______
  matrix<Type> Am(Mask.lines(),Mask.pixels());
  for (register uint i=beginl; i<A.lines()-beginl; i++)
    {
    for (register uint j=beginp; j<A.pixels()-beginp; j++)
      {
      Am.setdata(windef,A,winA);                // Am no allocs.
      Am -= mean(Am);                           // center around mean
      real8 covAM  = 0.;                        // covariance A,Mask
      real8 varA   = 0.;                        // variance of A(part)
      Type  *pntM  = Mask[0];
      Type  *pntAm = Am[0];
      for (register uint l=0; l<Mask.size(); l++)
        {
        //covAM += ((*pntM++) * (*pntAm));        // wait for move pnt
        //varA  += sqr(*pntAm++);                 // pnt ++
        // changed by FvL (for g++/gcc > 4.0):
        covAM += ((*pntM) * (*pntAm));        // wait for move pnt
        varA  += sqr(*pntAm);                 // pnt ++
        *pntM++;
        *pntAm++;
        }
      // Result(i,j) = covAM / sqrt(varM*varA);
      Result(i,j) = real4(covAM / sqrt(varM*varA)); // [BO]
      winA.pixlo++;
      winA.pixhi++;
      }
    winA.linelo++;
    winA.linehi++;
    winA.pixlo = 0;
    winA.pixhi = winA.pixlo + Mask.pixels() - 1;
    }
  return Result;
  } // END correlate



/****************************************************************
 * C = -A;                                                      *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  operator - (const matrix<Type>& A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("operator -");
  #endif
  #ifdef __DEBUGMAT1
  if (!A.size())
    matERROR.print("matrix:: unary minus with empty matrix");
  #endif
  matrix<Type> Result(A.lines(),A.pixels());
  Type *pntA = A[0];
  Type *pntR = Result[0];
  //for (register uint i=0; i<Result.size(); ++i)
  //  (*pntR++) = -(*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<Result.size(); ++i)
    {
    (*pntR++) = -(*pntA);
    *pntA++;
    }
  return Result;
  } // END - (unary)



/****************************************************************
 * a = mean(A)                                                  *
 * Bert Kampes, 14-Jan-1999                                     *
 ****************************************************************/
template <class Type>
real8         mean        (const matrix<Type> &A)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("mean");
  #endif
  real8 sum=0.;
  // ______Ensure stride one memory______
  Type *pntA   = A[0];
  //for (register uint i=0; i<A.size(); ++i)
  //  sum += (*pntA++);
  // changed by FvL (for g++/gcc > 4.0):
  for (register uint i=0; i<A.size(); ++i)
    {
    sum += (*pntA);
    *pntA++;
    }
  return sum/real8(A.size());
  } // END mean



/****************************************************************
 * matrix<Type> S = sum(A,dim)                                  *
 * return sum over dim of A. dim=1 by default.                  *
 * always returns a matrix, not very handy...                   *
 * A=[1 2 3]  then: sum(A,1) = [5 7 9]; sum(A,2)=[6]            *
 *   [4 5 6]                                     [15]           *
 * Bert Kampes, 28-Mar-2000                                     *
 ****************************************************************/
template <class Type>
matrix<Type>  sum         (const matrix<Type> &A, int32 dim)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("sum");
  #endif
  Type sum = Type(0);
  matrix<Type> Res;                             // may be 1x1 ...
  if (A.isvector())
    {
    Res.resize(1,1);
    Type *pntA   = A[0];
    //for (uint i=0; i<A.size(); i++)
    //  sum += (*pntA++);
    // changed by FvL (for g++/gcc > 4.0):
    for (uint i=0; i<A.size(); i++)
      {
      sum += (*pntA);
      *pntA++;
      }
    Res(0,0) = sum;
    }
  else // no vector
    {
    switch (dim)
      {
      // ______ sum over rows ______
      case 1:
        {
        Res.resize(1,A.pixels());
        for (uint i=0; i<A.pixels(); ++i)
          {
          for (uint j=0; j<A.lines(); ++j)
            {
            sum += A(j,i);
            }
          Res(0,i) = sum;
          sum = Type(0);
          }
        }
        break;
      // ______ sum over columns, may be done by pointers for speed ______
      case 2:
        {
        Res.resize(A.lines(),1);
        for (uint i=0; i<A.lines(); ++i)
          {
          for (uint j=0; j<A.pixels(); ++j)
            {
            sum += A(i,j);
            }
          Res(i,0) = sum;
          sum = Type(0);
          }
        }
        break;
      default: matERROR.print("sum: dim!={1,2}");
      }
    } // ifelse vector
  return Res;
  } // END sum



/****************************************************************
 * A.mypow(scalar)                                              *
 * NOmatrix<Type> S = pow(A,scalar)                             *
 * return pointwize power.                                      *
 * always returns a matrix, not very handy...                   *
 #%// BK 26-Oct-2000
 ****************************************************************/
template <class Type>
void matrix<Type>::mypow(Type s)
  {
  for (register uint i=0; i<nrows; ++i)
    for (register uint j=0; j<ncols; ++j)
      data[i][j] = pow(data[i][j],s);
  } // END sum


/****************************************************************
 * convert matrix A to type of matrix B                         *
 * Mahmut Arikan, 07-Jun-2009                                   *
 * TODO convert_type whole stuff should goto operator = one day.*
 ****************************************************************/
template <class Type, class Type2>
void  convert_type   (const matrix<Type> &A, const matrix<Type2> &B)
  {
  TRACE_FUNCTION("convert_type(matrix<Type>A-->matrix<Type2>B) (MA 07-Jun-2009)")
//  #ifdef __DEBUGMAT2
//    matDEBUG.print("convert_type(A-->B).");
//  #endif

  Type             *pntA = A[0];
  Type2            *pntB = B[0];

  if ( A.lines()!=B.lines() || A.pixels()!=B.pixels() )
    {
    DEBUG.print("convert_type aborted since the number of lines or/and pixels of the matrices are not equal.");
    return;
    }
  else if (  A.lines()==0 || A.pixels() == 0 )
    {
    DEBUG.print("convert_type aborted since the number of lines or/and pixels of the input matrix is 0.");
    return;
    }

  //if ( sizeof(Type) == sizeof(Type2) )
  if ( getformat(*pntA) == getformat(*pntB) )
    {
    cerr << "==| convert_type input and output types are the same types" << "==| "  << getformat(*pntA) << " to "  << getformat(*pntB) << endl;
    DEBUG.print("convert_type was not necessary since both input and output types are the same.");
    memcpy(pntB,pntA,A.size()*sizeof(Type));
    return;
    }

  for (register uint32 i=0; i<A.size(); i++)
    {
      *pntB = Type2( *pntA ); // less ambiguous
       pntB++; 
       pntA++;
    }

  } // END convert_type


