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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixbk_test.cc,v $      *
 * $Revision: 3.11 $                                            *
 * $Date: 2005/10/18 13:46:51 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * Test program for matrix class.                               *
 * Compile with Doris Makefile.                                 *
 #%// BK 17-Aug-2000                                            *
 ****************************************************************/


#include "constants.hh"                 // typedefs
#include "matrixbk.hh"                  // matrix functions
bk_messages TRACE;
bk_messages DEBUG;
bk_messages INFO;
bk_messages PROGRESS;
bk_messages WARNING;
bk_messages ERROR;
bk_messages matDEBUG;
bk_messages matERROR;
// ______ used in matrix_bk.cc to keep track of total allocated ______
uint totalallocated=0;  // [B] matrixdebuging





/****************************************************************
 *    main                                                      *
 #%// BK 17-Aug-2000                                            *
 ****************************************************************/
int main(
        int argc,
        char* argv[])
  {
  const int L = 4;      // numlines A
  const int P = 2;      // numcols  A
  int i,j;

  matrix<complr4> A(L,P);
  matrix<complr4> B(P,L+1);
  for (i=0; i<A.lines(); ++i)
    for (j=0; j<A.pixels(); ++j)
      A(i,j) = complr4(float(i)/2.3+j,float(j)/2+i);
  for (i=0; i<B.lines(); ++i)
    for (j=0; j<B.pixels(); ++j)
      B(i,j) = complr4(float(j+i)/2.3+i,float(i)/(2.2+i-j*j));

  // ______ Show matrices for testing ______
  cout << "matrix A(" << A.lines() << "," << A.pixels() << "):\n";
  A.showdata();
  cout << "matrix B(" << B.lines() << "," << B.pixels() << "):\n";
  B.showdata();
  cout << "\n\n---------------------\n\n";

  // ______ Test - ______
  B = complr4(2.) * -B;
  cout << "B = 2*-B(" << B.lines() << "," << B.pixels() << "):\n";
  B.showdata();
  cout << "\n\n---------------------\n\n";

  // ______ Test * ______
  matrix<complr4> C = A*B;
  cout << "matrix C(" << C.lines() << "," << C.pixels() << ") = A*B\n";
  C.showdata();
  cout << "\n\n---------------------\n\n";

  // ______ Test matTxmat ______
  matrix<complr4> C1 = matTxmat(transpose(A),B);
  cout << "matrix C(" << C1.lines() << "," << C1.pixels() << ") = A*B (using matTxmat)\n";
  C1.showdata();
  cout << "\n\n---------------------\n\n";

  // ______ Test matTxmat ______
  matrix<complr4> C2 = matxmatT(A,transpose(B));
  cout << "matrix C(" << C2.lines() << "," << C2.pixels() << ") = A*B (using matxmatT)\n";
  C2.showdata();
  cout << "\n\n---------------------\n\n";
  // --- check result ---
  cout << "checksum1 = " << max(magnitude(C-C1)) << endl;
  cout << "checksum2 = " << max(magnitude(C-C2)) << endl;



  // --- Conjugated ---
  matrix<complr4> Bconj = conj(B);
  cout << "conjugated:\n";
  Bconj.showdata();
  cout << "\n\n---------------------\n\n";

  // ______ Test 1d fft ______
  fft(A,1);
  cout << "1d fft over columns of matrix A(" << A.lines() << "," << A.pixels() << "):\n";
  A.showdata();
  cout << "\n\n---------------------\n\n";

//    // ______ Test operator = ______
//    for (i=0; i<10; ++i)
//      {
//      cout << i << ": A=B\n";
//      A=B;
//      cout << i << ": B=B\n";
//      B=B;
//      cout << i << ": C=A\n";
//      matrix<complr4> C=A;
//      cout << "\n";
//      }

  // ______ Test find function ______
  // matrix<int32> indexNaN = A.find(2);

  // ______ Test mypow function ______
  matrix<float> R(L,P);
  for (i=0; i<R.lines(); ++i)
    for (j=0; j<R.pixels(); ++j)
      R(i,j) = float(i/2.3+j*i*j);
  
  cout << "\n\n-RRRRRRRRRRRRRRRRRRR-\n\n";
  R.showdata();
  R.mypow(1.5);
  R.showdata();
// ?? but does B*= A^C work?
// ?? but does B*= A.mypow(C) work?

  // ______ Test CR4*R4 function ______
  cout << "\n\n-RRRRRRRRRRRRRRRRRRR-\n\n";
  matrix<complr4> Q(L,P);
  for (i=0; i<Q.lines(); ++i)
    for (j=0; j<Q.pixels(); ++j)
      Q(i,j) = complr4(float(i)/2.3+j,float(j)/2+i);
  cout << "\n\n- Q.showdata()\n\n";
  Q.showdata();
  cout << "\n\n- R.showdata()\n\n";
  R.showdata();
  Q*=R;
  cout << "\n\n- (Q*=R).showdata()\n\n";
  Q.showdata();

  // ______ Test conjugated function ______
  cout << "\n\n- CONJUGATED: Q.conj()\n\n";
  Q.conj();
  Q.showdata();
  cout << "\n\n- Q2 = CONJ(Q)\n\n";
  matrix<complr4> Q2 = conj(Q);
  Q2.showdata();
return 0;

  // ______ Test smooth function ______
  // R.smooth(1);
  // ______ Test pow ^ operator ______
  // R = R^2.0;

  // ______ Test wshift function ______
  cout << "\n\n- WSHIFT -\n\n";
  matrix<real4> AA(1,7);
  for (i=0; i<AA.lines(); ++i)
    for (j=0; j<AA.pixels(); ++j)
      AA(i,j) = i+j;
  AA.showdata();
  cout << "wshift AA, -2\n";
  wshift(AA,-2);
  AA.showdata();

  // ______ Test diagxmat R4*CR4 ______
  matrix<complr4> QQ(7,7);
  for (i=0; i<QQ.lines(); ++i)
    for (j=0; j<QQ.pixels(); ++j)
      QQ(i,j) = complr4(float(i)/2.3+j,float(j)/2+i);
  QQ.showdata();
  matrix<complr4> BB = diagxmat(AA,QQ);
  cout << "BB=diagxmat(AA,QQ)\n";
  BB.showdata();

  // ______ Test cos/sin complr4 operator ______
  cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n";
  matrix<real8> r8A(3,5);
  r8A = 5.3;
  cout << "r8A=5.3:\n";
  r8A.showdata();

  for (i=0; i<r8A.lines(); ++i)
    for (j=0; j<r8A.pixels(); ++j)
      r8A(i,j) = (3.14/180.0)*((0.5+i-j)*100);  // rad
  cout << "r8A=radians:\n";
  r8A.showdata();
  matrix<real8> cosr8A = cos(r8A);
  cout << "\ncosr8A:\n";
  cosr8A.showdata();
  //
  matrix<real8> r8B = r8A+0.2;
  matrix<real8> sinr8B = sin(r8B);
  cout << "\nsinr8B:\n";
  sinr8B.showdata();
  //
  matrix<complr4> cr4AB = mat2cr4(r8A,r8B);
  cout << "\ncomplr4(r8A,r8B):\n";
  cr4AB.showdata();

  // ______ Test mysort rows on first col ______
  cout << "\n\nmysort(real4) matrix ascending on rows\n";
  matrix<real4> SS(5,3);
  SS(0,0)=2.1;  SS(0,1)=10.1;  SS(0,2)=500.1;
  SS(1,0)=1.1;  SS(1,1)=20.1;  SS(1,2)=400.1;
  SS(2,0)=3.1;  SS(2,1)=30.1;  SS(2,2)=300.1;
  SS(3,0)=5.1;  SS(3,1)=40.1;  SS(3,2)=200.1;
  SS(4,0)=1.1;  SS(4,1)=50.1;  SS(4,2)=100.1;
  cout << "original data (SS)\n";
  SS.showdata();
  cout << "mysort2(SS)\n";
  mysort2(SS);
  SS.showdata();

  matrix<int32> SSi(5,3);
  SSi(0,0)=2;  SSi(0,1)=10;  SSi(0,2)=500;
  SSi(1,0)=1;  SSi(1,1)=20;  SSi(1,2)=400;
  SSi(2,0)=3;  SSi(2,1)=30;  SSi(2,2)=300;
  SSi(3,0)=5;  SSi(3,1)=40;  SSi(3,2)=200;
  SSi(4,0)=1;  SSi(4,1)=50;  SSi(4,2)=100;
  cout << "original data (SSi)\n";
  SSi.showdata();
  cout << "mysort2(SSi)\n";
  mysort2(SSi);
  SSi.showdata();

  // ______ Test oversample ______
  cout << "\n\nmysort(real4) matrix ascending on rows\n";
  matrix<real4> OVS(4,4);
  OVS(0,0)=2.1;  OVS(0,1)=10.1;  OVS(0,2)=500.1;
  OVS(1,0)=1.1;  OVS(1,1)=20.1;  OVS(1,2)=400.1;
  OVS(2,0)=3.1;  OVS(2,1)=30.1;  OVS(2,2)=300.1;
  OVS(3,0)=5.1;  OVS(3,1)=40.1;  OVS(3,2)=200.1;
  cout << "original data (OVS)\n";
  OVS.showdata();
  cout << "mysort2(SS)\n";
  matrix<real4> OVS_ = oversample(OVS,2,2);
  OVS_.showdata();

  // ______ Test fliplr/flipud ______
  cout << "\n\nOVS.showdata()\n";
  OVS.showdata();
  cout << "OVS.flipud()\n";
  OVS.flipud();
  OVS.showdata();
  cout << "OVS.fliplr() (ie both up and lr)\n";
  OVS.fliplr();
  OVS.showdata();

  cout << "\n\nNormal termination.\n\n";
  return 0;
  } // END main

