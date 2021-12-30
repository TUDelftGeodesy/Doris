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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/referencephase.cc,v $ *
 * $Revision: 3.22 $                                            *
 * $Date: 2006/05/18 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * -computation flat earth correction.                          *
 * -computation radarcoding dem + interpolation to (1,1) grid   *
 ****************************************************************/

#include "matrixbk.hh"
#include "orbitbk.hh"
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class
#include "constants.hh"
#include "referencephase.hh"            // proto types
#include "ioroutines.hh"                // error etc.
#include "utilities.hh"                 // isodd; ones(), etc.
#include "coregistration.hh"            // distribute points
#include "exceptions.hh"                 // my exceptions class

#include <iomanip>                      // setw only for test..
#include <cstdlib>                      // system
#include <algorithm>                    // max
#ifdef WIN32
  // Jia defined this.
  // Bert Kampes, 24-Aug-2005
  #include "winsock2.h"
#else
  #include <netinet/in.h>                 // ntohl byteorder x86-HP unix
#endif

// Using A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator
// from  Jonathan Richard Shewchuk
// Some definition for triangulate call
#define VOID int
#define REAL double
#define ANSI_DECLARATORS
extern "C" {
#include "triangle.h"
}


/****************************************************************
 *    flatearth                                                 *
 *                                                              *
 * Compute polynomial model for 'flat earth' correction.        *
 *  fie(l,p) = sumj=0:d sumk=0:d Ajk l^j p^k (NOT bert 8sept99) *
 * precise orbits are used to compute delta range for Npoints   *
 * after which the polynomial model is fitted (LS).             *
 *                                                              *
 * input:                                                       *
 *  - inputoptions                                              *
 *  - info structs                                              *
 *  - platform data points                                      *
 * output:                                                      *
 *  - void (result to file "scratchresflat")                    *
 *  - coefficients normalized wrt. original window of master    *
 *                                                              *
 *    Bert Kampes, 09-Mar-1999                                  *
 *    Bert Kampes, 26-Oct-1999 normalization of coeff.,         *
 *    dump to logfile: var(unknowns) == diag(inv(AtA))          *
 ****************************************************************/
void flatearth(
        const input_comprefpha &comprefphainput,
        const input_ell        &ellips,
        const slcimage         &master,
        const slcimage         &slave,
        const productinfo      &interferogram,
        orbit                  &masterorbit,
        orbit                  &slaveorbit)
  {
  TRACE_FUNCTION("flatearth (BK 26-Oct-1999)")
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;
  char                  dummyline[2*ONE27];

  INFO << "FLATEARTH: MAXITER: "   << MAXITER   << "; "
                  << "CRITERPOS: " << CRITERPOS << " m; "
                  << "CRITERTIM: " << CRITERTIM << " s";
  INFO.print();

  // ______ Normalization factors for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;
  INFO << "flatearth: polynomial normalized by factors: "
       << minL << " " << maxL << " " << minP << " " << maxP
       << " to [-2,2]";
  INFO.print();

  // ______Handling of input______
  const real8 m_minpi4cdivlam = (-4.0*PI*SOL)/master.wavelength;
  const real8 s_minpi4cdivlam = (-4.0*PI*SOL)/slave.wavelength;
  DEBUG << "master wavelength = " << master.wavelength;
  DEBUG.print();
  DEBUG << "slave  wavelength = " << slave.wavelength;
  DEBUG.print();
  const int32 DEGREE    = comprefphainput.degree;
  const int32 Nunk      = Ncoeffs(DEGREE);              // Number of unknowns
  bool pointsrandom = true;
  if (specified(comprefphainput.ifpositions))
    pointsrandom = false;                       // only use those points



  // ______ Distribute points wel distributed over win ______
  // ______ or read from ascii file ______
  // ______(i,0): line, (i,1): pixel, (i,2) flagfromdisk______
  //matrix<uint> Position;
  // [FvL] for correct folding of points outside overlap window when inserted by file
  matrix<int> Position;
  const uint  Npoints = comprefphainput.Npoints;
  register int32 i,j,k,index;

  if (pointsrandom)                             // no filename specified
    {
    Position = distributepoints(Npoints,interferogram.win);
    }
  else // read from file
    {
    Position.resize(Npoints,3);
    //ifstream ifpos(comprefphainput.ifpositions, ios::in);
    ifstream ifpos;
    openfstream(ifpos,comprefphainput.ifpositions);
    bk_assert(ifpos,comprefphainput.ifpositions,__FILE__,__LINE__);
    uint ll,pp;
    for (i=0; i<Npoints; ++i)
      {
      ifpos >> ll >> pp;
      //Position(i,0) = uint(ll);
      //Position(i,1) = uint(pp);
      //Position(i,2) = uint(1);                // flag from file
      // [FvL]
      Position(i,0) = int(ll);
      Position(i,1) = int(pp);
      Position(i,2) = int(1);                // flag from file
      ifpos.getline(dummyline,2*ONE27,'\n');       // goto next line.
      }
    ifpos.close();

    // ______ Check last point ivm. EOL after last position in file ______
    if (Position(Npoints-1,0) == Position(Npoints-2,0) &&
        Position(Npoints-1,1) == Position(Npoints-2,1))
      {
      Position(Npoints-1,0) = uint(.5*(minL + maxL) + 27);      // random
      Position(Npoints-1,1) = uint(.5*(minP + maxP) + 37);      // random
      WARNING << "refpha: there should be no EOL after last point in file: "
           << comprefphainput.ifpositions;
      WARNING.print();
      }

    // ______ Check if points are in overlap ______
    // ______ no check for uniqueness of points ______
    }

  matrix<real8>         y(Npoints,1);                   // observation
  matrix<real8>         y_h2ph(Npoints,1);              // observation, h2ph factors, added by FvL
  matrix<real8>         A(Npoints,Nunk);                // designmatrix

  // ______Check redundancy______
  if (Npoints < Nunk)
    {
    PRINT_ERROR("flatearth: Number of points is smaller than parameters solved for.");
    throw(input_error);
    }



  // ======Compute delta r for all points======
  for (i=0; i<Npoints; ++i)
    {
    const real8 line  = Position(i,0);
    const real8 pixel = Position(i,1);
   
    // ______ Compute azimuth/range time of this pixel______
    //const real8 m_trange = pix2tr(pixel,master.t_range1,master.rsr2x);
    const real8 m_trange = master.pix2tr(pixel);
    const real8 m_tazi = master.line2ta(line); // added by FvL

    // ______ Compute xyz of this point P from position in image ______
    cn P;                                       // point, returned by lp2xyz
    lp2xyz(line,pixel,ellips,master,masterorbit,
           P,MAXITER,CRITERPOS);
 
    // ______ Compute xyz for slave satelite from P ______
    real8 s_tazi;                               // returned
    real8 s_trange;                             // returned
    xyz2t(s_tazi,s_trange,slave,
          slaveorbit,
          P,MAXITER,CRITERTIM);
   

// ______Compute delta range ~= phase______
// ______ real8 dr = dist(m_possat,pospoint) - dist(s_possat,pospoint);
// ______ real8 phase = -pi4*(dr/LAMBDA);
// ______  dr    == M-S         want if no flatearth M-S - flatearth = M-S-(M-S)=0
// ______  phase == -4pi*dr/lambda == 4pi*(S-M)/lambda
// BK: 24-9: actually defined as: phi = +pi4/lambda * (r1-r2) ???
    // real8 phase = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
    //y(i,0) = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
    //y(i,0) = pi4divlam*(s_possat.dist(pospoint)-m_possat.dist(pospoint));
    //y(i,0) = minpi4cdivlam * (m_trange - s_trange);
    y(i,0) = m_minpi4cdivlam*m_trange - s_minpi4cdivlam*s_trange;
    DEBUG << "l="   << line     << " p="  << pixel 
          << " t1=" << m_trange << " t2=" << s_trange
          << " fe=" << y(i,0) << " [rad]";
    DEBUG.print();

    // ____________________________________________________________________________________
    // _____________ Vector with h2ph factors for random number of points by FvL __________
    //_____________________________________________________________________________________
    
    cn Psat_master = masterorbit.getxyz(m_tazi);
    cn Psat_slave = slaveorbit.getxyz(s_tazi);
    real8 B    = Psat_master.dist(Psat_slave);      // abs. value
    // const real8 Bpar = P.dist(M) - P.dist(S);    // sign ok
    real8 Bpar = SOL*(m_trange-s_trange);   // sign ok

    // ______ if (MP>SP) then S is to the right of slant line, then B perp is positive.
    cn r1 = Psat_master.min(P);
    cn r2 = Psat_slave.min(P);
    // real8 theta = Psat_master.angle(r1);  // look angle
    real8 theta       = P.angle(r1);            // incidence angle
    real8 theta_slave = P.angle(r2);            // incidence angle slave
    real8 Bperp = (theta > theta_slave ) ?     // sign ok
                   sqrt(sqr(B)-sqr(Bpar)) :
                  -sqrt(sqr(B)-sqr(Bpar))  ;

    y_h2ph(i,0) = Bperp/(m_trange*SOL*sin(theta));
    
    // ____________________________________________________________________________________
    // _____________ End added part by FvL ________________________________________________
    //_____________________________________________________________________________________

    // ______Set up system of equations______
    // ______Order unknowns: A00 A10 A01 A20 A11 A02 A30 A21 A12 A03 for degree=3______
    // ______  normalize data [-2,2] ______
    real8 posL          = normalize(line,minL,maxL);
    real8 posP          = normalize(pixel,minP,maxP);
    index = 0;
    for (j=0; j<=DEGREE; j++)
      {
      for (k=0; k<=j; k++)
        {
        A(i,index) = pow(posL,real8(j-k)) * pow(posP,real8(k));
        index++;
        }
      }
    }


  // ======Compute polynomial for these phases (LS)======
  // ______Compute Normalmatrix, rghthandside______
  matrix<real8> N   = matTxmat(A,A);
  matrix<real8> rhs = matTxmat(A,y);
  matrix<real8> rhs_h2ph = matTxmat(A,y_h2ph); // Added by FvL, same A matrix can be used


  // ______Compute solution______
  matrix<real8> Qx_hat = N;
  choles(Qx_hat);               // Cholesky factorisation normalmatrix
  solvechol(Qx_hat,rhs);        // Estimate of unknowns in rhs
  solvechol(Qx_hat,rhs_h2ph);   // Estimate of unknowns in rhs_h2ph, added by FvL
  invertchol(Qx_hat);           // Covariance matrix


  // ______Test inverse______
  for (i=0; i<Qx_hat.lines(); i++)
    for (j=0; j<i; j++)
      Qx_hat(j,i) = Qx_hat(i,j);// repair Qx_hat
  const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
  INFO << "flatearth: max(abs(N*inv(N)-I)) = " << maxdev;
  INFO.print();
  if (maxdev > .01) 
    {
    ERROR << "Deviation too large. Decrease degree or number of points?";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  else if (maxdev > .001) 
    {
    WARNING << "Deviation quite large. Decrease degree or number of points?";
    WARNING.print();
    }
  else
    {
    INFO.print("Deviation is OK.");
    }

  // ______Some other stuff, scale is ok______
  //  matrix<real8> Qy_hat        = A * (matxmatT(Qx_hat,A));
  matrix<real8> y_hat           = A * rhs;
  matrix<real8> y_hat_h2ph      = A * rhs_h2ph; // added by FvL
  matrix<real8> e_hat           = y - y_hat;
  matrix<real8> e_hat_h2ph      = y_h2ph - y_hat_h2ph; // added by FvL


  // ______Overall model test (variance factor)______
  // ... ?





  // ______ Wrap offset ______ 
  // BK 30/9/99 do not do this, later absolute ref. phase is used.
  // in s2h rodriguez for example.
  // it does not change anything for compinterfero etc.
  //  rhs(0,0) = remainder(rhs(0,0),2*PI);



  // ______Write results to file______
  ofstream scratchlogfile("scratchlogflat", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"flatearth: scratchlogflat",__FILE__,__LINE__);

  scratchlogfile << "\n\n*******************************************************************"
                 << "\n* FLATEARTH: "
                 //<< "\n*_Start_" << processcontrol[pr_i_comprefpha]
                  << "\n*******************************************************************"
                 << "\nDegree_flat:\t" << DEGREE
                 << "\nEstimated coefficients:\n"
                 << "\nx_hat \tstd:\n";
  for (i=0; i<Nunk; i++)
    scratchlogfile
                 << setiosflags(ios::fixed)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::right)
                 << setw(8) << setprecision(4) 
                 <<  rhs(i,0) << " \t" << sqrt(Qx_hat(i,i)) << endl;

  // ___________________ added by FvL _________________________________________________________
  
  scratchlogfile << "\n" << "\nDegree_h2ph:\t" << DEGREE
     << "\nEstimated coefficients:\n"
     << "\nx_hat \tstd:\n";
  for (i=0; i<Nunk; i++)
    scratchlogfile
     << setiosflags(ios::fixed)
     << setiosflags(ios::showpoint)
     << setiosflags(ios::right)
     << setw(8) << setprecision(4) 
     <<  rhs_h2ph(i,0) << " \t" << sqrt(Qx_hat(i,i)) << endl;
  // ___________________ end added by FvL _________________________________________________________

  scratchlogfile << "\nCovariance matrix estimated parameters:"
                 << "\n---------------------------------------\n";
  for (i=0; i<Nunk; i++)
    {
    for (j=0; j<Nunk; j++)
      {
      scratchlogfile 
                 << setiosflags(ios::fixed)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::right)
                 << setw(8) << setprecision(4) 
                 << Qx_hat(i,j) << " ";
      }
    scratchlogfile << endl;
    }

  scratchlogfile << "\nMaximum deviation N*inv(N):"
                 << setiosflags(ios::scientific)
                 << maxdev
                 << "\nSome more info for each observation:"
                 << "\nline \tpixel \tobs \t\tobs_hat \t\t err_hat\n";
  for (i=0; i<Npoints; i++)
    scratchlogfile <<  Position(i,0) << "\t" <<  Position(i,1) << "\t"
                   << y(i,0) << "\t" << y_hat(i,0) << "\t"
                   << setiosflags(ios::fixed)
                   << setiosflags(ios::showpoint)
                   << setiosflags(ios::right)
                   << setw(8) << setprecision(4) 
                   << e_hat(i,0) << "\n";
  scratchlogfile << "\nMaximum absolute error: \t\t\t"
                 <<  max(abs(e_hat))
                 << "\n*******************************************************************\n";
  //  scratchlogfile.close(); // commented out for coordinates dumped below

  ofstream scratchresfile("scratchresflat", ios::out | ios::trunc);
  bk_assert(scratchresfile,"flatearth: scratchresflat",__FILE__,__LINE__);

  scratchresfile.setf(ios::scientific, ios::floatfield);
  scratchresfile.setf(ios::right, ios::adjustfield);
  scratchresfile.precision(8);
  scratchresfile.width(18);

  scratchresfile << "\n\n*******************************************************************"
                 //<< "\n*_Start_flat_earth"
                 << "\n*_Start_" << processcontrol[pr_i_comprefpha]
                 << "\n*******************************************************************"
                 << "\nDegree_flat:\t" << DEGREE
                 << "\nEstimated_coefficients_flatearth:\n";
  int32 coeffL = 0;
  int32 coeffP = 0;
  for (i=0; i<Nunk; i++)
    {
    if (rhs(i,0) < 0.)
      scratchresfile <<         rhs(i,0);
    else
      scratchresfile << " " <<  rhs(i,0);

    // ______ Add coefficient number behind value ______
    scratchresfile << " \t" <<  coeffL << " " << coeffP << "\n";
    coeffL--;
    coeffP++;
    if (coeffL == -1)
      {
      coeffL = coeffP;
      coeffP = 0;
      }
    }
  
  //_________ added by FvL _______________________________________  
  scratchresfile << "\n" << "\nDegree_h2ph:\t" << DEGREE
     << "\nEstimated_coefficients_h2ph:\n";
  coeffL = 0;
  coeffP = 0;
  for (i=0; i<Nunk; i++)
    {
    if (rhs_h2ph(i,0) < 0.)
      scratchresfile <<         rhs_h2ph(i,0);
    else
      scratchresfile << " " <<  rhs_h2ph(i,0);

  // ______ Add coefficient number behind value ______
    scratchresfile << " \t" <<  coeffL << " " << coeffP << "\n";
    coeffL--;
    coeffP++;
    if (coeffL == -1)
      {
      coeffL = coeffP;
      coeffP = 0;
      }
    }
  //_________ end added by FvL _______________________________________  

  scratchresfile << "*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_comprefpha] << "_NORMAL"
                 << "\n*******************************************************************\n";


  // ====== Compute coordinates of corners of interferogram here ======
  // ______ (though better place this somewhere else ....
  real8 phi,lambda,height;                              // rad, returned
  lp2ell(interferogram.win.linelo,interferogram.win.pixlo,
         ellips,master,masterorbit,
         phi,lambda,height,MAXITER,CRITERPOS);
  INFO << "Coordinates of corner interferogram: "
       << interferogram.win.linelo << ", " << interferogram.win.pixlo
       << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
  INFO.print();
  
  scratchlogfile << "\n\n********************************************"
                 <<   "\n* [Lat_Long] coordinates of crop [START]   *"
                 <<   "\n********************************************\n";

  scratchlogfile << "\nCoords_of_ifg_corner [l,p] : [phi,lam]: "
                 << interferogram.win.linelo << " , " << interferogram.win.pixlo
                 << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

  lp2ell(interferogram.win.linehi,interferogram.win.pixlo,
         ellips,master,masterorbit,
         phi,lambda,height,MAXITER,CRITERPOS);
  INFO << "\nCoordinates of corner interferogram: "
       << interferogram.win.linehi << ", " << interferogram.win.pixlo
       << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
  INFO.print();

  scratchlogfile << "\nCoords_of_ifg_corner [l,p] : [phi,lam]: "
                 << interferogram.win.linehi << " , " << interferogram.win.pixlo
                 << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

  lp2ell(interferogram.win.linelo,interferogram.win.pixhi,
         ellips,master,masterorbit,
         phi,lambda,height,MAXITER,CRITERPOS);
  INFO << "\nCoordinates of corner interferogram: "
       << interferogram.win.linelo << ", " << interferogram.win.pixhi
       << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
  INFO.print();

  scratchlogfile << "\nCoordinates of ifgs  [low,hi] corner [l,p] : [phi,lam]: "
                 << interferogram.win.linelo << " , " << interferogram.win.pixhi
                 << " = " << rad2deg(phi) << " , " << rad2deg(lambda);

  lp2ell(interferogram.win.linehi,interferogram.win.pixhi,
         ellips,master,masterorbit,
         phi,lambda,height,MAXITER,CRITERPOS);
  INFO << "\nCoordinates of corner interferogram: "
       << interferogram.win.linehi << ", " << interferogram.win.pixhi
       << " = " << rad2deg(phi) << ", " << rad2deg(lambda);
  INFO.print();

  scratchlogfile <<  "\nCoords_of_ifg_corner [l,p] : [phi,lam]: "
                 << interferogram.win.linehi << " , " << interferogram.win.pixhi
                 << " = " << rad2deg(phi) << " , " << rad2deg(lambda);
  
  scratchlogfile << "\n\n******************************************"
                 <<   "\n* [Lat_Long] coordinates of crop [STOP]  *"
                 <<   "\n******************************************\n";
  // ______Tidy up______
  scratchresfile.close(); // close res file
  scratchlogfile.close(); // close log.out file
  } // END flatearth



/****************************************************************
 *    demassist                                                 * 
 *                                                              *
 * Coregistration based on DEM (SRTM)                          *
 * DEM on equiangular grid (lat/lon) assumed                    *
 * DEM seems stored from North to South                         *
 *                                                              *
 * Freek van Leijen, Liu Guang, Mahmut Arikan, 21-Sep-2007      *
 ****************************************************************/
void demassist(
        const input_gen        &generalinput,
        const input_ell        &ellips,
        const input_demassist  &demassistinput,
        const slcimage         &master,
        const slcimage         &slave,
        orbit                  &masterorbit,
        orbit                  &slaveorbit)
  {
  TRACE_FUNCTION("demassist (FvL 21-SEP-2007)")

  const string STEP="DAC: ";
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;

  const real8 lat0file    = demassistinput.demlatleftupper;// first pix on disk w02090
  const real8 lon0file    = demassistinput.demlonleftupper;// first pix on disk
  const real8 DEMdeltalat = demassistinput.demdeltalat;    // in radians
  const real8 DEMdeltalon = demassistinput.demdeltalon;    // in radians
  const int32 numberoflonpixels = demassistinput.demcols;  // NCOLS on file
  const int32 numberoflatpixels = demassistinput.demrows;  // NROWS on file
  const real8 NODATA      =  demassistinput.demnodata;     // (BK 4 may 2001)
  const bool outputdemi   =  specified(demassistinput.fodemi);// if spec. then output
  const bool outputrefdemhei   =  specified(demassistinput.forefdemhei);

  ////const real8 m_min4picdivlam = (-4.0*PI*SOL)/master.wavelength;
  ////const real8 s_min4picdivlam = (-4.0*PI*SOL)/slave.wavelength;

  const real8 latNfile = lat0file-DEMdeltalat*(numberoflatpixels-1); // upper=max. lat value
  const real8 lonNfile = lon0file+DEMdeltalon*(numberoflonpixels-1); // left=min. lon value

  // ______ Extra info ______
  INFO << "DEM input: w/e/s/n:          \t"
       << rad2deg(lon0file) << "/" << rad2deg(lonNfile) << "/"
       << rad2deg(latNfile) << "/" << rad2deg(lat0file);
  INFO.print();

  // ______ Get corners of master (approx) to select DEM ______
  // ______ in radians (if height were zero)______
  real8 extralat = (1.5*DEMdeltalat + deg2rad(4.0/25.0));
  real8 extralong = (1.5*DEMdeltalon + deg2rad(4.0/25.0));

  real8 phimin;
  real8 phimax;
  real8 lambdamin;
  real8 lambdamax;
  int32 indexphi0DEM;
  int32 indexphiNDEM;
  int32 indexlambda0DEM;
  int32 indexlambdaNDEM;
  getcorners(master.currentwindow.linelo,master.currentwindow.linehi,
             master.currentwindow.pixlo,master.currentwindow.pixhi,
             extralat,extralong,lat0file,lon0file,
             DEMdeltalat,DEMdeltalon,numberoflatpixels,numberoflonpixels,
             ellips,master,masterorbit,phimin,phimax,lambdamin,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

  // ______ Extra info ______
  INFO << "DEM input required: w/e/s/n: \t"
       << rad2deg(lambdamin) << "/" << rad2deg(lambdamax) << "/"
       << rad2deg(phimin)    << "/" << rad2deg(phimax);
  INFO.print();
  INFO << "For window (l0,lN,p0,pN):    \t"
       << master.currentwindow.linelo << " "
       << master.currentwindow.linehi << " "
       << master.currentwindow.pixlo << " "
       << master.currentwindow.pixhi;
  INFO.print();


  // ______ Check corners of DEM ______
  // check if DEM is appropriate for master crop
  // DEM should at least partially cover master crop
  // note: phi is [90:-90]
  if (phimax <= latNfile)// DEM is more north than master
    {
    ERROR << "master crop outside DEM: most South latitude: " << rad2deg(latNfile)
         << " [deg]; master crop requires: " << rad2deg(phimax) 
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    //throw(some_error);
    }
  // DEM is more south than master crop
  if (phimin >= lat0file)// largest latitude at first line of file
    {
    ERROR << "master crop outside DEM: most North latitude: " << rad2deg(lat0file)
         << " [deg]; master crop requires: " << rad2deg(phimax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    //throw(some_error);
    }
  if (lambdamax <= lon0file)
    {
    ERROR << "master crop outside DEM: most West longitude: " << rad2deg(lon0file)
         << " [deg]; master crop window requires: " << rad2deg(lambdamax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    //throw(some_error);
    }
  if (lambdamin >= lonNfile)
    {
    ERROR << "master crop outside DEM: most East longitude: " << rad2deg(lonNfile)
         << " [deg]; master crop window requires: " << rad2deg(lambdamin)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    //throw(some_error);
    }


  //===================================================================
  //============ First loop: radarcode DEM ============================
  //============ (DEM geometry)            ============================
  //===================================================================

  int32 numvalid        = 0;// number of good values, not NODATA in buffer
  int32 numNODATA       = 0;// number of NODATA values in buffer
  real8 meancroppedDEM  = 0.0;// to detect byte order problems, formats
  real8 min_input_dem   =  100000.0;// stats
  real8 max_input_dem   = -100000.0;// stats

  // ______ Compute buffer size radarcoding DEM______
  const real8 BUFFERMEMSIZE = generalinput.memory;// Bytes
  int32 NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
  int32 NrowsDEM = indexphiNDEM-indexphi0DEM+1;
  const uint32 NcolsDEMlog = NcolsDEM;       // since NcolsDEM updated after getcorners() call.
  const uint32 NrowsDEMlog = NrowsDEM;
  const real8 Nrows_possible_DEM = BUFFERMEMSIZE / (5*8*NcolsDEM);
  int32 bufferlines = int32(ceil(Nrows_possible_DEM));   // [MA] checked ok. Sinces  SLC is not multilooked, see comprefdem for solution
  if (bufferlines>NrowsDEM) bufferlines=NrowsDEM;
  int32 numfullbuffers = NrowsDEM / bufferlines;
  int32 restlines      = NrowsDEM % bufferlines;
  int32 extrabuffer = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "DEM output total pixels: " << NcolsDEM;
  INFO.print();
  INFO << "DEM output total lines : " << NrowsDEM;
  INFO.print();
  INFO << "Radar coding of DEM in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();
  


  // ______ Open (temporary) output files ______
  // DEM heights 
  ofstream demofile;
  openfstream(demofile,demassistinput.fodem,generalinput.overwrit); // dem_crop radarcoded
  bk_assert(demofile,demassistinput.fodem,__FILE__,__LINE__);

  // master line coordinates of DEM
  ofstream masterdemlineoutfile("dac_m_demline.temp", ios::out | ios::trunc); 
  bk_assert(masterdemlineoutfile,"dac_m_demline.temp",__FILE__,__LINE__);
  
  // master pixel coordinates of DEM
  ofstream masterdempixeloutfile("dac_m_dempixel.temp", ios::out | ios::trunc); 
  bk_assert(masterdempixeloutfile,"dac_m_dempixel.temp",__FILE__,__LINE__);

  // delta line coordinates of DEM  ( slave-master )
  ofstream deltademlineoutfile("dac_delta_demline.temp", ios::out | ios::trunc); 
  bk_assert(deltademlineoutfile,"dac_delta_demline.temp",__FILE__,__LINE__);
  
  // delta pixel coordinates of DEM
  ofstream deltadempixeloutfile("dac_delta_dempixel.temp", ios::out | ios::trunc); 
  bk_assert(deltadempixeloutfile,"dac_delta_dempixel.temp",__FILE__,__LINE__);



  // ______ DEM loop per buffer ______
  register int32 j,i;// DEM index grid counter, register j first to ensure allocation
  for (register int32 buffer=0; buffer<numfullbuffers+extrabuffer; ++buffer)
    {

     // Determine indices for buffer
    const int32 indexphi0BUFFER = indexphi0DEM+buffer*bufferlines;
    const int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
    const int32 indexphiNBUFFER = indexphi0BUFFER+(blines-1);
    matrix<real4> DEM(blines,NcolsDEM);

    // ______ Extra info ______
    PROGRESS << STEP << "Buffer# [l0:lN, p0:pN]: " << buffer+1 << " ["
   << indexphi0BUFFER  << ": " << indexphiNBUFFER  << ", "
   << indexlambda0DEM << ": " << indexlambdaNDEM << "]";
    PROGRESS.print();

    // ______ lat/lon for first pixel in matrix read from file ______
    // ______ upper is max. latitude, left is min. longitude ______
    const real8 upperleftphi    = lat0file-indexphi0BUFFER*DEMdeltalat;
    const real8 upperleftlambda = lon0file+indexlambda0DEM*DEMdeltalon;

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0BUFFER,indexphiNBUFFER,
                        indexlambda0DEM,indexlambdaNDEM);

    // ______ Read in grdfile of DEM in matrix R4 (raw data, no header) _______
    // ______ added formats (BK 4-May-2001) ______
    PROGRESS << STEP << "Reading crop of DEM for buffer: " << buffer+1;
    PROGRESS.print();
    DEBUG.print("Reading input DEM into real4 matrix (buffer).");
    switch (demassistinput.iformatflag)
      {
      // ______ Read as short BE, then convert to host order ______
      case FORMATI2_BIGENDIAN:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj)));// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATI2:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMi2(iii,jjj);// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATR4:
        readfile(DEM,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        INFO.print("Read crop of input DEM: format: REAL4.");
        break;
      case FORMATR8:
        {
        matrix<real8> DEMr8(blines,NcolsDEM);
        readfile(DEMr8,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMr8(iii,jjj);// cast to real4
        DEMr8.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: REAL8.");
        break;
        }
      default:
        PRINT_ERROR("totally impossible, checked input.")
        //throw(unhandled_case_error);
      }


    // ----- Loop over DEM for stats ------------------------
    real8 min_dem_buffer =  100000.0;
    real8 max_dem_buffer = -100000.0;
    for (i=0; i<DEM.lines(); ++i)
      {
      // ----- Loop over oversampled matrix in x ------
      for (j=0; j<DEM.pixels(); ++j)
        {
        if(DEM(i,j)!=NODATA)
          {
          numvalid++;
          meancroppedDEM += DEM(i,j);// divide by numvalid later
          if (DEM(i,j)<min_dem_buffer) min_dem_buffer=DEM(i,j);//buffer
          if (DEM(i,j)>max_dem_buffer) max_dem_buffer=DEM(i,j);// stats
          }
        else
          {
          numNODATA++;
          }
        }//loop dem for stats
      }//loop dem for stats
    min_input_dem = min(min_input_dem,min_dem_buffer);//global stats
    max_input_dem = max(max_input_dem,max_dem_buffer);//global stats


    // ====== Radarcoding DEM ==============================
    // ______ DEM contains values from leftupper with ______
    // ______ spacing (DEMdeltalat,DEMdeltalon) ______
    // ______ Transform DEM to l,p,refphase ______
    PROGRESS.print("Converting DEM to radar system for this buffer.");
    const int32 NpointsDEM = DEM.size();
    const int32 NpixelsDEM = DEM.pixels();
    // ______ Extra info ______
    INFO << "Number of points in DEM: "
         << NpointsDEM;
    INFO.print();

    matrix<real8> masterDEMline(DEM.lines(),DEM.pixels());
    matrix<real8> masterDEMpixel(DEM.lines(),DEM.pixels());
    matrix<real8> deltaDEMline(DEM.lines(),DEM.pixels());
    matrix<real8> deltaDEMpixel(DEM.lines(),DEM.pixels());

    // --- Loop DEM ---
    real8 phi,lambda,height,m_l,m_p,s_l,s_p;


    phi = upperleftphi;
    for (i=0; i<DEM.lines(); ++i)
      {
      if ((i%100)==0)
        {
        // ______ Extra info ______
        PROGRESS << STEP << "Radarcoding buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " DEM line: " << i << " ("
             << floor(.5+(100.*real8(i)/real8(DEM.lines())))
             << "%)";
        PROGRESS.print();
        }

      lambda = upperleftlambda;
      for (j=0; j<DEM.pixels(); ++j)
        {
    height = DEM(i,j);
    ell2lp(m_l,m_p,ellips,master,masterorbit,phi,lambda,height,MAXITER,CRITERTIM);
    ell2lp(s_l,s_p,ellips,slave,slaveorbit,phi,lambda,height,MAXITER,CRITERTIM);
    masterDEMline(i,j) = m_l;
    masterDEMpixel(i,j) = m_p;
          deltaDEMline(i,j) = s_l-m_l;
          deltaDEMpixel(i,j) = s_p-m_p;

          lambda += DEMdeltalon;
        } // loop DEM pixels

      // ______ update latitude of next line ______
      phi -= DEMdeltalat;           // upper left is max. value
      } // loop DEM lines


    // Write results to output files 
    PROGRESS << STEP << "Writing radar coded DEM to file, buffer: " << buffer+1 << " of " << numfullbuffers + extrabuffer ;
    PROGRESS.print();

    demofile << DEM;                         
    masterdemlineoutfile << masterDEMline;
    masterdempixeloutfile << masterDEMpixel;
    deltademlineoutfile << deltaDEMline;
    deltadempixeloutfile << deltaDEMpixel;

    masterDEMline.resize(1,1); //deallocate
    masterDEMpixel.resize(1,1); //deallocate
    deltaDEMline.resize(1,1); //deallocate
    deltaDEMpixel.resize(1,1); //deallocate
    DEM.resize(1,1); //deallocate
    } // buffer loop

  demofile.close();
  masterdemlineoutfile.close();
  masterdempixeloutfile.close();
  deltademlineoutfile.close();
  deltadempixeloutfile.close();


  //===================================================================
  //============ End first loop: radarcode DEM ========================
  //============ (DEM geometry)            ============================
  //===================================================================


  //===================================================================
  //============ Second loop: interpolation               =============
  //============ (radar geometry)                         =============
  //===================================================================
  
  INFO << STEP << "Start interpolation...";
  INFO.print();

  // ______ Line/pixel of first point in original master coordinates ______
  const int32 Nlinesml    = master.currentwindow.lines();
  const int32 Npixelsml   = master.currentwindow.pixels();
  
  const real8 veryfirstline = real8(master.currentwindow.linelo);
  const real8 verylastline = real8(master.currentwindow.linehi);
  const real8 firstpixel = real8(master.currentwindow.pixlo);
  const real8 lastpixel = real8(master.currentwindow.pixhi);


  //Determine range-azimuth spacing ratio, needed for proper triangulation
  cn P1, P2 , P3, P4;
  lp2xyz(veryfirstline,firstpixel,ellips,master,masterorbit,
           P1,MAXITER,CRITERPOS);
  lp2xyz(veryfirstline,lastpixel,ellips,master,masterorbit,
           P2,MAXITER,CRITERPOS);
  lp2xyz(verylastline,firstpixel,ellips,master,masterorbit,
           P3,MAXITER,CRITERPOS);
  lp2xyz(verylastline,lastpixel,ellips,master,masterorbit,
           P4,MAXITER,CRITERPOS);

  const real8 r_spacing  = ( (P1.min(P2)).norm() + (P3.min(P4)).norm() ) / 2 /(lastpixel - firstpixel) ;
  const real8 az_spacing = ( (P1.min(P3)).norm() + (P2.min(P4)).norm() ) /2 /(verylastline - veryfirstline); 
  const real8 r_az_ratio = r_spacing/az_spacing;

  INFO << "Master azimuth spacing: " << az_spacing;
  INFO.print();
  INFO << "Master range spacing: " << r_spacing;
  INFO.print();
  INFO << "Range-azimuth spacing ratio: " << r_az_ratio;
  INFO.print();

  // ______ Compute buffer size interpolation______
  const real8 Nlinesml_possible = BUFFERMEMSIZE / (6*8*Npixelsml);
  bufferlines = int32(ceil(Nlinesml_possible));
  if (bufferlines > Nlinesml) bufferlines=Nlinesml;
  numfullbuffers = Nlinesml / bufferlines;
  restlines      = Nlinesml % bufferlines; // the number of lines in extra buffer
  extrabuffer = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "Interpolation in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();


  // ______ Open output files ______
  ofstream deltalineofile("dac_delta_line.raw", ios::out | ios::trunc); 
  bk_assert(deltalineofile,"dac_delta_line.raw",__FILE__,__LINE__);
  
  ofstream deltapixelofile("dac_delta_pixel.raw", ios::out | ios::trunc); 
  bk_assert(deltapixelofile,"dac_delta_pixel.raw",__FILE__,__LINE__);
  
  // if request for height in radar coordinates l,p
  ofstream refdemheiofile;
  if (outputrefdemhei==true)
    {
    openfstream(refdemheiofile,demassistinput.forefdemhei,generalinput.overwrit);
    bk_assert(refdemheiofile,demassistinput.forefdemhei,__FILE__,__LINE__);
    }
  
  // ______ interpolation loop per buffer ______
  for (register int32 buffer = 0; buffer < numfullbuffers + extrabuffer; ++buffer)
    {

    // Determine indices for buffer
    const int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
    const real8 firstline_buffer = veryfirstline+buffer*bufferlines;
    const real8 lastline_buffer = firstline_buffer+blines-1;

    // ______ Extra info ______
    PROGRESS << STEP << "Interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " [l0:lN, p0:pN]: " << " ["
   << firstline_buffer << ": " << lastline_buffer << ", "
   << firstpixel << ": " << lastpixel << "]";
    PROGRESS.print();

    // Get corners of buffer
    real8 phimin_az;
    real8 phimax_az;
    real8 lambdamin_az;
    real8 lambdamax_az;
    getcorners(firstline_buffer,lastline_buffer,
             firstpixel,lastpixel,
             extralat,extralong,phimax,lambdamin,
             DEMdeltalat,DEMdeltalon,NrowsDEM,NcolsDEM,
             ellips,master,masterorbit,phimin_az,phimax_az,lambdamin_az,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0DEM,indexphiNDEM,
                        indexlambda0DEM,indexlambdaNDEM);
    const int32 NrowsDEM_buffer = indexphiNDEM-indexphi0DEM+1;
    const int32 NcolsDEM_buffer = indexlambdaNDEM-indexlambda0DEM+1;

    PROGRESS << STEP << "Reading input for interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer;
    PROGRESS.print();

    // read x,y
    matrix<real8> DEMline_buffer(NrowsDEM_buffer,NcolsDEM_buffer);
    matrix<real8> DEMpixel_buffer(NrowsDEM_buffer,NcolsDEM_buffer);

    readfile(DEMline_buffer,"dac_m_demline.temp",NrowsDEM,winfromfile,zerooffset);
    readfile(DEMpixel_buffer,"dac_m_dempixel.temp",NrowsDEM,winfromfile,zerooffset);
    
    // read z (multiple, number can easily be increased, e.g. simulated intensity)
    int32 Nz = 2; //number of z
    matrix<real8> input_buffer(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);
    matrix<real8> temp_input_buffer(NrowsDEM_buffer,NcolsDEM_buffer);
    if (outputrefdemhei==true)
      {
      Nz += 1;
      input_buffer.resize(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);
      }

    readfile(temp_input_buffer,"dac_delta_demline.temp",NrowsDEM,winfromfile,zerooffset);
    input_buffer.setdata(0, 0, temp_input_buffer);
    readfile(temp_input_buffer,"dac_delta_dempixel.temp",NrowsDEM,winfromfile,zerooffset);
    input_buffer.setdata(NrowsDEM_buffer, 0, temp_input_buffer);
    Nz = 2;
    if (outputrefdemhei==true)
      {
        Nz += 1;
        /// i would like to use real4, test later on
        matrix<real4> dem_input(NrowsDEM_buffer,NcolsDEM_buffer);
        readfile(dem_input,demassistinput.fodem,NrowsDEM,winfromfile,zerooffset);
        for (register int32 i =0 ; i < NrowsDEM_buffer ; i ++)
          for(register int32 j = 0; j < NcolsDEM_buffer; j++)
            temp_input_buffer(i,j) = real8(dem_input(i,j));
        input_buffer.setdata(NrowsDEM_buffer * (Nz-1), 0, temp_input_buffer);
      }
    
    // initialize output array
    Nz = 2;  
    matrix<real8> output_buffer(blines * Nz, Npixelsml);
    if (outputrefdemhei==true)
      {
        Nz += 1;
        output_buffer.resize(blines * Nz, Npixelsml);
      }

    // interpolation
    griddatalinear(DEMline_buffer,DEMpixel_buffer,input_buffer,
                   firstline_buffer,lastline_buffer,firstpixel,lastpixel,
                   1,1,r_az_ratio,0,NODATA,output_buffer);

    deltalineofile << output_buffer(window(0, blines - 1, 0, Npixelsml -1 ));
    deltapixelofile << output_buffer(window(blines , 2 * blines - 1, 0, Npixelsml -1 ));
    Nz = 2;
    if (outputrefdemhei==true)
      {
        Nz += 1;
        refdemheiofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
      }
    
    DEMline_buffer.resize(1,1);
    DEMpixel_buffer.resize(1,1);
    input_buffer.resize(1,1);
    temp_input_buffer.resize(1,1);
    output_buffer.resize(1,1);

  } // end loop azimuth direction

  INFO << "Closing output files";
  INFO.print();

  deltalineofile.close();
  deltapixelofile.close();
  if (outputrefdemhei==true) // Radarcoded DEM
    refdemheiofile.close();

  //===================================================================
  //============ End second loop: interpolation           =============
  //============ (radar geometry)                         =============
  //===================================================================

  
  //===================================================================
  //============ Determine inverse transformation         =============
  //============ (slave corners only, needed for overlap) =============
  //===================================================================
  
  real8 line, pixel;
  real8 deltaline_slave00,deltapixel_slave00,
    deltaline_slave0N,deltapixel_slave0N,
    deltaline_slaveN0,deltapixel_slaveN0,
    deltaline_slaveNN,deltapixel_slaveNN;
  real8 phimin_az,phimax_az,lambdamin_az,lambdamax_az;

  
  for (register int16 corner = 0 ; corner < 4 ; corner ++)
    {
      
      PROGRESS << "Radarcoding slave corner: " << corner+1;
      PROGRESS.print();

      switch (corner)
        {
        case 0:
          {
            line=slave.currentwindow.linelo;
            pixel=slave.currentwindow.pixlo;
            break;
          }
        case 1:
          {
            line=slave.currentwindow.linelo;
            pixel=slave.currentwindow.pixhi;
            break;
          }
        case 2:
          {
            line=slave.currentwindow.linehi;
            pixel=slave.currentwindow.pixlo;
            break;
          }
        case 3:
          {
            line=slave.currentwindow.linehi;
            pixel=slave.currentwindow.pixhi;
            break;
          }
        default:
          PRINT_ERROR("totally impossible, checked input.");
          
        }

      //use getcorners with line,line,pixel,pixel for single point
      getcorners(line,line,pixel,pixel,
                 extralat,extralong,lat0file,lon0file,
                 DEMdeltalat,DEMdeltalon,numberoflatpixels,numberoflonpixels,
                 ellips,slave,slaveorbit,
                 phimin_az,phimax_az,lambdamin_az,lambdamax,
                 indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);
      
      
      NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
      NrowsDEM = indexphiNDEM-indexphi0DEM+1;
      const real8 upperleftphi    = lat0file-indexphi0DEM*DEMdeltalat;
      const real8 upperleftlambda = lon0file+indexlambda0DEM*DEMdeltalon;
      
      window zerooffset  (0,0,0,0);
      window winfromfile (indexphi0DEM,indexphiNDEM,
                          indexlambda0DEM,indexlambdaNDEM);

      
      // ______ Read in DEM in matrix R4 (raw data, no header) _______

      matrix<real4> DEM(NrowsDEM,NcolsDEM);
     
      switch (demassistinput.iformatflag)
        {
          // ______ Read as short BE, then convert to host order ______
        case FORMATI2_BIGENDIAN:
          {
            matrix<int16> DEMi2(NrowsDEM,NcolsDEM);
            readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
            for (int32 iii=0; iii<DEM.lines(); ++iii)
              for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
                DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj)));// cast to real4
            DEMi2.resize(1,1);// dealloc...
            INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
            break;
          }
          
        case FORMATI2:
          {
            matrix<int16> DEMi2(NrowsDEM,NcolsDEM);
            readfile(DEMi2,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
            for (int32 iii=0; iii<DEM.lines(); ++iii)
              for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
                DEM(iii,jjj) = DEMi2(iii,jjj);// cast to real4
            DEMi2.resize(1,1);// dealloc...
            INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
            break;
          }
          
        case FORMATR4:
          readfile(DEM,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
          INFO.print("Read crop of input DEM: format: REAL4.");
          break;
        case FORMATR8:
          {
            matrix<real8> DEMr8(NrowsDEM,NcolsDEM);
            readfile(DEMr8,demassistinput.firefdem,numberoflatpixels,winfromfile,zerooffset);
            for (int32 iii=0; iii<DEM.lines(); ++iii)
              for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
                DEM(iii,jjj) = DEMr8(iii,jjj);// cast to real4
            DEMr8.resize(1,1);// dealloc...
            INFO.print("Read crop of input DEM: format: REAL8.");
            break;
          }
        default:
          PRINT_ERROR("totally impossible, checked input.");
          //throw(unhandled_case_error);
        }


      // radarcode dem
      matrix<real8> slaveDEMline(DEM.lines(),DEM.pixels());
      matrix<real8> slaveDEMpixel(DEM.lines(),DEM.pixels());
      matrix<real8> deltaDEMline(DEM.lines(),DEM.pixels());
      matrix<real8> deltaDEMpixel(DEM.lines(),DEM.pixels());
      
      // --- Loop DEM ---
      real8 phi,lambda,height,m_l,m_p,s_l,s_p;
      
      
      phi = upperleftphi;
      for (i=0; i<DEM.lines(); ++i)
        {
          
          lambda = upperleftlambda;
          for (j=0; j<DEM.pixels(); ++j)
            {
              height = DEM(i,j);
              ell2lp(m_l,m_p,ellips,master,masterorbit,phi,lambda,height,MAXITER,CRITERTIM);
              ell2lp(s_l,s_p,ellips,slave,slaveorbit,phi,lambda,height,MAXITER,CRITERTIM);
              slaveDEMline(i,j) = s_l;
              slaveDEMpixel(i,j) = s_p;
              deltaDEMline(i,j) = m_l-s_l;
              deltaDEMpixel(i,j) = m_p-s_p;
              
              lambda += DEMdeltalon;
            } // loop DEM pixels
          
          // ______ update latitude of next line ______
          phi -= DEMdeltalat;           // upper left is max. value
        } // loop DEM lines


      
      // interpolate to slave corner
      matrix<real8> input_buffer(DEM.lines()*2,DEM.pixels());
      input_buffer.setdata(0, 0, deltaDEMline);
      input_buffer.setdata(DEM.lines(), 0, deltaDEMpixel);

      matrix<real8> output_buffer(2,1);

      griddatalinear(slaveDEMline,slaveDEMpixel,input_buffer,
                     line,line,pixel,pixel,
                     1,1,r_az_ratio,0,NODATA,output_buffer);
      
      switch (corner)
        {
        case 0:
          {
            deltaline_slave00 = output_buffer(0,0);
            deltapixel_slave00 = output_buffer(1,0);
            INFO << "Deltaline_slave00: " << deltaline_slave00;
            INFO.print();
            INFO << "Deltapixel_slave00: " << deltapixel_slave00;
            INFO.print();
            break;
          }
        case 1:
          {
            deltaline_slave0N = output_buffer(0,0);
            deltapixel_slave0N = output_buffer(1,0);
            INFO << "Deltaline_slave0N: " << deltaline_slave0N;
            INFO.print();
            INFO << "Deltapixel_slave0N: " << deltapixel_slave0N;
            INFO.print();
            break;
          }
        case 2:
          {
            deltaline_slaveN0 = output_buffer(0,0);
            deltapixel_slaveN0 = output_buffer(1,0);
            INFO << "Deltaline_slaveN0: " << deltaline_slaveN0;
            INFO.print();
            INFO << "Deltapixel_slaveN0: " << deltapixel_slaveN0;
            INFO.print();
            break;
          }
        case 3:
          {
            deltaline_slaveNN = output_buffer(0,0);
            deltapixel_slaveNN = output_buffer(1,0);
            INFO << "Deltaline_slaveNN: " << deltaline_slaveNN;
            INFO.print();
            INFO << "Deltapixel_slaveNN: " << deltapixel_slaveNN;
            INFO.print();
            break;
          }
        default:
          PRINT_ERROR("totally impossible, checked input.");
          
        }

    }
  
  //===================================================================
  //============ End determine inverse transformation     =============
  //============ (slave corners only, needed for overlap) =============
  //===================================================================


  // ====== Write output information ======
  char croppeddemi[4*ONE27];
  strcpy(croppeddemi,"NO output requested");
  if (outputdemi) strcpy(croppeddemi,demassistinput.fodemi);
  INFO << "Min. value of input DEM covering master: " << min_input_dem; 
  INFO.print();
  INFO << "Max. value of input DEM covering master: " << max_input_dem; 
  INFO.print();

  ofstream scratchlogfile("scratchlogdemassist", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"demassist: scratchlogdemassist",__FILE__,__LINE__);
  scratchlogfile
    << "\n*******************************************************************"
    << "\n* " << processcontrol[pr_i_demassist]
    << "\n*******************************************************************"
    << "\n1) DEM source file:                   \t" <<  demassistinput.firefdem
    << "\nFormat:                               \t";
    switch (demassistinput.iformatflag)
      {
      case FORMATI2:
        {
        scratchlogfile << "SHORT SIGNED INTEGER (HOST ENDIANNESS)";
        break;
        }
      case FORMATI2_BIGENDIAN:
        {
        scratchlogfile << "SHORT SIGNED INTEGER, BIG ENDIAN";
        break;
        }
      case FORMATR4:
        {
        scratchlogfile << "REAL4 SIGNED FLOAT";
        break;
        }
      case FORMATR8:
        {
        scratchlogfile << "REAL8 SIGNED DOUBLE";
        break;
        }
      default:
        {
        scratchlogfile << "UNKNOWN? IMPOSSIBLE...";
        break;
        }
      }
  scratchlogfile
    << "\nByte order:                           \t" <<  "check it yourself..."
    << "\nNumber of lines:                      \t" <<  numberoflatpixels
    << "\nNumber of pixels:                     \t" <<  numberoflonpixels
    << "\nResolution latitude:                  \t" << rad2deg(DEMdeltalat) << " [deg]"
    << "\nResolution longitude:                 \t" << rad2deg(DEMdeltalon) << " [deg]"
    << "\nMost West point in input DEM:         \t" << rad2deg(lon0file)
    << "\nMost East point in input DEM:         \t" << rad2deg(lonNfile)
    << "\nMost South point in input DEM:        \t" << rad2deg(latNfile)
    << "\nMost North point in input DEM:        \t" << rad2deg(lat0file)
    << "\nMin. value of input DEM covering master: " << min_input_dem
    << "\nMax. value of input DEM covering master: " << max_input_dem
    << "\n2) Output file cropped DEM:           \t" <<  demassistinput.fodem //[FVL
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines               :        \t" << NrowsDEMlog 
    << "\nNumber of pixels              :        \t" << NcolsDEMlog
    << "\nDEM extend w/e/s/n            :      \t" << rad2deg(lambdamin) << "/" 
    << rad2deg(lambdamax) << "/" << rad2deg(phimin)    << "/" << rad2deg(phimax)
//    << "\nMean value:                           \t" <<  meancroppedDEM
    << "\n3) Output file interpolated crop DEM: \t" << croppeddemi
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):        \t" << Nlinesml
    << "\nNumber of pixels (multilooked):       \t" << Npixelsml
    << "\nDeltaline_slave00_dem:                    \t" << deltaline_slave00
    << "\nDeltapixel_slave00_dem:                   \t" << deltapixel_slave00
    << "\nDeltaline_slave0N_dem:                    \t" << deltaline_slave0N
    << "\nDeltapixel_slave0N_dem:                   \t" << deltapixel_slave0N
    << "\nDeltaline_slaveN0_dem:                    \t" << deltaline_slaveN0
    << "\nDeltapixel_slaveN0_dem:                   \t" << deltapixel_slaveN0
    << "\nDeltaline_slaveNN_dem:                    \t" << deltaline_slaveNN
    << "\nDeltapixel_slaveNN_dem:                   \t" << deltapixel_slaveNN
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();


  ofstream scratchresfile("scratchresdemassist", ios::out | ios::trunc);
  bk_assert(scratchresfile,"demassist: scratchresdemassist",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_demassist]
    << "\n*******************************************************************";
  scratchresfile
    << "\nDEM source file:                      \t" << demassistinput.firefdem
    << "\nMin. of input DEM:                    \t" << min_input_dem
    << "\nMax. of input DEM:                    \t" << max_input_dem
    << "\nFirst_line (w.r.t. original_master):   \t"
    <<  master.currentwindow.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  master.currentwindow.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  master.currentwindow.pixlo
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  master.currentwindow.pixhi
    << "\nNumber of lines:        \t" << Nlinesml
    << "\nNumber of pixels:       \t" << Npixelsml
    << "\nDeltaline_slave00_dem:      \t" << deltaline_slave00
    << "\nDeltapixel_slave00_dem:     \t" << deltapixel_slave00
    << "\nDeltaline_slave0N_dem:      \t" << deltaline_slave0N
    << "\nDeltapixel_slave0N_dem:     \t" << deltapixel_slave0N
    << "\nDeltaline_slaveN0_dem:      \t" << deltaline_slaveN0
    << "\nDeltapixel_slaveN0_dem:     \t" << deltapixel_slaveN0
    << "\nDeltaline_slaveNN_dem:      \t" << deltaline_slaveNN
    << "\nDeltapixel_slaveNN_dem:     \t" << deltapixel_slaveNN
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_demassist] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();


 } // END demassist

 
/****************************************************************
 *    radarcodedem (a.k.a comprefdem)                           *
 *                                                              *
 * Compute reference phase based on DEM (SRTM)                  *
 * DEM on equiangular grid (lat/lon) assumed                    *
 * DEM seems stored from North to South                         *
 *                                                              *
 * Freek van Leijen, 26-Sep-2007                                *
 ****************************************************************/
void radarcodedem(
        const input_gen        &generalinput,
        const input_ell        &ellips,
        const input_comprefdem &refdeminput,
        const slcimage         &master,
        const slcimage         &slave,
        const productinfo      &interferogram,
        orbit                  &masterorbit,
        orbit                  &slaveorbit)
  {
  TRACE_FUNCTION("radarcodedem (FvL 26-Sep-2007)")

  const string STEP="CRD: ";
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;

  const real8 lat0file    = refdeminput.demlatleftupper;      // first pix on disk w02090
  const real8 lon0file    = refdeminput.demlonleftupper;      // first pix on disk
  const real8 DEMdeltalat = refdeminput.demdeltalat;          // in radians
  const real8 DEMdeltalon = refdeminput.demdeltalon;          // in radians
  const int32 numberoflonpixels = refdeminput.demcols;        // NCOLS on file
  const int32 numberoflatpixels = refdeminput.demrows;        // NROWS on file
  const real8 NODATA      =  refdeminput.demnodata;           // (BK 4 may 2001)
  bool onlyrefphasetopo   = !refdeminput.includerefpha;       // true: phase DEM w.r.t. ellipsoid
  const bool outputdemi   =  specified(refdeminput.fodemi);   // if spec. then output
  const bool outputh2ph   =  specified(refdeminput.foh2ph);   // if spec. then output, added by FvL
  const bool outputrefdemhei   =  specified(refdeminput.forefdemhei);
  
  // _____ start added by MA _____
  bool mlookedIFG         =  false;                           // true: ifg is multilooked
 
  int32 mlL               = interferogram.multilookL;         // initialize multilookfactor
  int32 mlP               = interferogram.multilookP;
  const int32 &ifgmlL     = interferogram.multilookL;         // multilookfactor of interferogram
  const int32 &ifgmlP     = interferogram.multilookP;         // multilookfactor of interferogram
  if ( ifgmlL != 1 || ifgmlP != 1 )                           // [MA] additional entry for Coherence comptation using refdem.
    {                                                         //  always do computation without multilooking
      mlL        = 1;                                         // set multilookfactor for interpolation
      mlP        = 1;                                         // set multilookfactor for interpolation
      mlookedIFG = true;                                      // dealing with mlooked ifg.
    }
  // _____ end added by MA _____

  
  
  const real8 m_min4picdivlam = (-4.0*PI*SOL)/master.wavelength;
  const real8 s_min4picdivlam = (-4.0*PI*SOL)/slave.wavelength;
  DEBUG << "master wavelength = " << master.wavelength;
  DEBUG.print();
  DEBUG << "slave  wavelength = " << slave.wavelength;
  DEBUG.print();

  const real8 latNfile = lat0file-DEMdeltalat*(numberoflatpixels-1); // upper=max. lat value
  const real8 lonNfile = lon0file+DEMdeltalon*(numberoflonpixels-1); // left=min. lon value

  // ______ Extra info ______
  INFO << "DEM input: w/e/s/n:          \t"
       << rad2deg(lon0file) << "/" << rad2deg(lonNfile) << "/"
       << rad2deg(latNfile) << "/" << rad2deg(lat0file);
  INFO.print();

  // ______ Get corners of interferogram (approx) to select DEM ______
  // ______ in radians (if height were zero)______
  real8 extralat = (1.5*DEMdeltalat + deg2rad(4.0/25.0));
  real8 extralong = (1.5*DEMdeltalon + deg2rad(4.0/25.0));

  real8 phimin;
  real8 phimax;
  real8 lambdamin;
  real8 lambdamax;
  int32 indexphi0DEM;
  int32 indexphiNDEM;
  int32 indexlambda0DEM;
  int32 indexlambdaNDEM;
  const uint &ifglinelo = interferogram.win.linelo;   // [MA] win no-mlooked master coords 
  const uint &ifglinehi = interferogram.win.linehi;
  const uint &ifgpixlo  = interferogram.win.pixlo;
  const uint &ifgpixhi  = interferogram.win.pixhi;

  getcorners(ifglinelo,ifglinehi,
             ifgpixlo,ifgpixhi,
             extralat,extralong,lat0file,lon0file,
             DEMdeltalat,DEMdeltalon,numberoflatpixels,numberoflonpixels,
             ellips,master,masterorbit,phimin,phimax,lambdamin,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

  // ______ Extra info ______
  INFO << "DEM input required: w/e/s/n: \t"
       << rad2deg(lambdamin) << "/" << rad2deg(lambdamax) << "/"
       << rad2deg(phimin)    << "/" << rad2deg(phimax);
  INFO.print();
  INFO << "For window (l0,lN,p0,pN):    \t"
       << ifglinelo << " "
       << ifglinehi << " "
       << ifgpixlo << " "
       << ifgpixhi;
  INFO.print();


  // ______ Check corners of DEM ______
  // check if DEM is appropriate for interferogram
  // DEM should at least partially cover IFG
  // note: phi is [90:-90]
  if (phimax <= latNfile)// DEM is more north than IFG
    {
    ERROR << "IFG outside DEM: most South latitude: " << rad2deg(latNfile)
         << " [deg]; IFG requires: " << rad2deg(phimax) 
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  // DEM is more south than IFG
  if (phimin >= lat0file)// largest latitude at first line of file
    {
    ERROR << "IFG outside DEM: most North latitude: " << rad2deg(lat0file)
         << " [deg]; IFG requires: " << rad2deg(phimax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  if (lambdamax <= lon0file)
    {
    ERROR << "IFG outside DEM: most West longitude: " << rad2deg(lon0file)
         << " [deg]; IFG window requires: " << rad2deg(lambdamax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  if (lambdamin >= lonNfile)
    {
    ERROR << "IFG outside DEM: most East longitude: " << rad2deg(lonNfile)
         << " [deg]; IFG window requires: " << rad2deg(lambdamin)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }


  //===================================================================
  //============ First loop: radarcode DEM ============================
  //============ (DEM geometry)            ============================
  //===================================================================

  int32 numvalid        = 0;// number of good values, not NODATA in buffer
  int32 numNODATA       = 0;// number of NODATA values in buffer
  real8 meancroppedDEM  = 0.0;// to detect byte order problems, formats
  real8 min_input_dem   =  100000.0;// stats
  real8 max_input_dem   = -100000.0;// stats

  // ______ Compute buffer size radarcoding DEM______
  const real8 BUFFERMEMSIZE = generalinput.memory;// Bytes
  const int32 NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
  const int32 NrowsDEM = indexphiNDEM-indexphi0DEM+1;
  const real8 Nrows_possible_DEM = BUFFERMEMSIZE / (5*8*NcolsDEM);
  int32 bufferlines = int32(ceil(Nrows_possible_DEM));
  INFO << "Possible max. buffer lines: " << bufferlines << " for " << BUFFERMEMSIZE << " memory size.";
  INFO.print();
  if (bufferlines>NrowsDEM) bufferlines=NrowsDEM;
  int32 numfullbuffers = NrowsDEM / bufferlines;
  int32 restlines      = NrowsDEM % bufferlines;
  int32 extrabuffer = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "DEM output total pixels: " << NcolsDEM;
  INFO.print();
  INFO << "DEM output total lines : " << NrowsDEM;
  INFO.print();
  INFO << "Radar coding of DEM in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();



  // ______ Open (temporary) output files ______
  // DEM heights 
  INFO<<refdeminput.fodem << endl;
  INFO.print();
  ofstream demofile;
  openfstream(demofile,refdeminput.fodem,generalinput.overwrit);
  bk_assert(demofile,refdeminput.fodem,__FILE__,__LINE__);

  // master line coordinates of DEM
  ofstream masterdemlineoutfile("crd_m_demline.temp", ios::out | ios::trunc); 
  bk_assert(masterdemlineoutfile,"crd_m_demline.temp",__FILE__,__LINE__);
  
  // master pixel coordinates of DEM
  ofstream masterdempixeloutfile("crd_m_dempixel.temp", ios::out | ios::trunc); 
  bk_assert(masterdempixeloutfile,"crd_m_dempixel.temp",__FILE__,__LINE__);

  // ref phase in DEM geometry
  ofstream demrefphaseoutfile("crd_dem_refphase.temp", ios::out | ios::trunc); 
  bk_assert(demrefphaseoutfile,"crd_dem_refphase.temp",__FILE__,__LINE__);

  // h2ph factor in DEM geometry
  ofstream demh2phoutfile; 
  if (outputh2ph==true)
    {
      openfstream(demh2phoutfile,"crd_dem_h2ph.temp",generalinput.overwrit);
      bk_assert(demh2phoutfile,"crd_dem_h2ph.temp",__FILE__,__LINE__);
    }

  // ______ DEM loop per buffer ______
  register int32 j,i;// DEM index grid counter, register j first to ensure allocation
  for (register int32 buffer=0; buffer<numfullbuffers+extrabuffer; ++buffer)
    {

     // Determine indices for buffer
    const int32 indexphi0BUFFER = indexphi0DEM+buffer*bufferlines;
    const int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
    const int32 indexphiNBUFFER = indexphi0BUFFER+(blines-1);
    matrix<real4> DEM(blines,NcolsDEM);

    // ______ Extra info ______
    PROGRESS << STEP << "Buffer# [l0:lN, p0:pN]: " << buffer+1 << " ["
   << indexphi0BUFFER  << ": " << indexphiNBUFFER  << ", "
   << indexlambda0DEM << ": " << indexlambdaNDEM << "]";
    PROGRESS.print();

    // ______ lat/lon for first pixel in matrix read from file ______
    // ______ upper is max. latitude, left is min. longitude ______
    const real8 upperleftphi    = lat0file-indexphi0BUFFER*DEMdeltalat;
    const real8 upperleftlambda = lon0file+indexlambda0DEM*DEMdeltalon;

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0BUFFER,indexphiNBUFFER,
                        indexlambda0DEM,indexlambdaNDEM);

    // ______ Read in grdfile of DEM in matrix R4 (raw data, no header) _______
    // ______ added formats (BK 4-May-2001) ______
    PROGRESS << STEP << "Reading crop of DEM for buffer: " << buffer+1;
    PROGRESS.print();
    DEBUG.print("Reading input DEM into real4 matrix (buffer).");
     INFO<< "file info: name: "     << refdeminput.firefdem <<",  nof flat pixels, " << numberoflatpixels << endl;
     INFO<< "file info, Format : "  << refdeminput.iformatflag<<endl;
     INFO.print();
    switch (refdeminput.iformatflag)
      {
      // ______ Read as short BE, then convert to host order ______
      case FORMATI2_BIGENDIAN:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj)));// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATI2:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMi2(iii,jjj);// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATR4:
        readfile(DEM,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        INFO.print("Read crop of input DEM: format: REAL4.");
        break;
      case FORMATR8:
        {
        matrix<real8> DEMr8(blines,NcolsDEM);
        readfile(DEMr8,refdeminput.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMr8(iii,jjj);// cast to real4
        DEMr8.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: REAL8.");
        break;
        }
      default:
        PRINT_ERROR("totally impossible, checked input.")
        //throw(unhandled_case_error);
      }


    // ----- Loop over DEM for stats ------------------------
    real8 min_dem_buffer =  100000.0;
    real8 max_dem_buffer = -100000.0;
    for (i=0; i<DEM.lines(); ++i)
      {
      // ----- Loop over oversampled matrix in x ------
      for (j=0; j<DEM.pixels(); ++j)
        {
        if(DEM(i,j)!=NODATA)
          {
          numvalid++;
          meancroppedDEM += DEM(i,j);// divide by numvalid later
          if (DEM(i,j)<min_dem_buffer) min_dem_buffer=DEM(i,j);//buffer
          if (DEM(i,j)>max_dem_buffer) max_dem_buffer=DEM(i,j);// stats
          }
        else
          {
          numNODATA++;
          }
        }//loop dem for stats
      }//loop dem for stats
    min_input_dem = min(min_input_dem,min_dem_buffer);//global stats
    max_input_dem = max(max_input_dem,max_dem_buffer);//global stats


    // ====== Radarcoding DEM ==============================
    // ______ DEM contains values from leftupper with ______
    // ______ spacing (DEMdeltalat,DEMdeltalon) ______
    // ______ Transform DEM to l,p,refphase ______
    PROGRESS.print("Converting DEM to radar system for this buffer.");
    const int32 NpointsDEM = DEM.size();
    const int32 NpixelsDEM = DEM.pixels();
    // ______ Extra info ______
    INFO << "Number of points in DEM: "
         << NpointsDEM;
    INFO.print();

    matrix<real8> masterDEMline(DEM.lines(),DEM.pixels());
    matrix<real8> masterDEMpixel(DEM.lines(),DEM.pixels());
    matrix<real8> ref_phase_array(DEM.lines(),DEM.pixels());
    matrix<real8> h2ph_array(DEM.lines(),DEM.pixels());

    // --- Loop DEM ---
    cn P;
    real8 phi,lambda,height,l,p,ref_phase;


    phi = upperleftphi;
    for (i=0; i<DEM.lines(); ++i)
      {
      if ((i%100)==0)
        {
        // ______ Extra info ______
        PROGRESS << STEP << "Radarcoding buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " DEM line: " << i << " ("
             << floor(.5+(100.*real8(i)/real8(DEM.lines())))
             << "%)";
        PROGRESS.print();
        }

      lambda = upperleftlambda;
      for (j=0; j<DEM.pixels(); ++j)
        {
    height = DEM(i,j);
    ell2lp(l,p,ellips,master,masterorbit,phi,lambda,height,MAXITER,CRITERTIM);
    masterDEMline(i,j) = l;
    masterDEMpixel(i,j) = p;
    P = ellips.ell2xyz(phi,lambda,height);     // returns P(x,y,z)

    real8 t_range_master;
    real8 t_azi_master;
    xyz2t(t_azi_master,t_range_master,
              master,masterorbit,
              P,MAXITER,CRITERTIM);
    real8 t_azi_slave;
    real8 t_range_slave;
    xyz2t(t_azi_slave,t_range_slave,                // t returned
        slave,slaveorbit,
        P,MAXITER,CRITERTIM);
    
  if (outputh2ph==true)
{
    // compute h2ph factor
    cn Psat_master = masterorbit.getxyz(t_azi_master);
    cn Psat_slave = slaveorbit.getxyz(t_azi_slave);
    real8 B = Psat_master.dist(Psat_slave);
    real8 Bpar = SOL*(t_range_master-t_range_slave);
    cn r1 = Psat_master.min(P);
    cn r2 = Psat_slave.min(P);
    // real8 theta = Psat_master.angle(r1);  // look angle
    real8 theta       = P.angle(r1);            // incidence angle
    real8 theta_slave = P.angle(r2);            // incidence angle slave
    real8 Bperp = (theta > theta_slave ) ? // sign ok
                   sqrt(sqr(B)-sqr(Bpar)) :
                  -sqrt(sqr(B)-sqr(Bpar))  ;

    h2ph_array(i,j) = Bperp/(t_range_master*SOL*sin(theta));
}


    if (onlyrefphasetopo)                 // do not include flat earth phase
      {
        lp2xyz(l,p,ellips,                      // h==0
         master, masterorbit,
         P,MAXITER,CRITERPOS);                        // P returned

        real8 t_range_flatearth,t_azi_dummy;

        xyz2t(t_azi_dummy,t_range_flatearth,
        slave,slaveorbit,
        P,MAXITER,CRITERTIM);                 // P on h=0
        ref_phase = s_min4picdivlam*t_range_flatearth-
    s_min4picdivlam*t_range_slave;
      }
    else // include flatearth, ref.pha = phi_topo+phi_flatearth
      {
        ref_phase = 
    m_min4picdivlam*master.pix2tr(p)-
    s_min4picdivlam*t_range_slave;
      }
          
        ref_phase_array(i,j) = ref_phase;

        lambda += DEMdeltalon;
        } // loop DEM pixels

      // ______ update latitude of next line ______
      phi -= DEMdeltalat;           // upper left is max. value
      } // loop DEM lines


    // Write results to output files 
    PROGRESS << STEP << "Writing radar coded DEM to file, buffer: " << buffer+1 << " of " << numfullbuffers + extrabuffer ;
    PROGRESS.print();

    demofile << DEM;
    masterdemlineoutfile << masterDEMline;
    masterdempixeloutfile << masterDEMpixel;
    demrefphaseoutfile << ref_phase_array;
if (outputh2ph==true)
  demh2phoutfile << h2ph_array;

    masterDEMline.resize(1,1); //deallocate
    masterDEMpixel.resize(1,1); //deallocate
    DEM.resize(1,1); //deallocate
    ref_phase_array(1,1); //deallocate
    h2ph_array(1,1); //deallocate
    } // buffer loop

  demofile.close();
  masterdemlineoutfile.close();
  masterdempixeloutfile.close();
  demrefphaseoutfile.close();
if (outputh2ph==true)
  demh2phoutfile.close();


  //===================================================================
  //============ End first loop: radarcode DEM ========================
  //============ (DEM geometry)            ============================
  //===================================================================


  //===================================================================
  //============ Second loop: interpolation               =============
  //============ (radar geometry)                         =============
  //===================================================================
  
  INFO << STEP << "Start interpolation...";
  INFO.print();

  // ______ Line/pixel of first point in original master coordinates ______
  // ______ maybe this should be changed to be x+(ml/2) ?? but depends on
  // ______ definition of range_to_first_bin is to center or start..
  // Bert Kampes, 08-Apr-2005: chose center by adding ml/2
  const int32   Nlinesml      = interferogram.win.lines()  / mlL;    // ifg lines when mlL = 1 (no multilooking)
  const int32   Npixelsml     = interferogram.win.pixels() / mlP;
  const int32   ifgNlinesml   = interferogram.win.lines()  / ifgmlL; // for the result file when mlL != 1
  const int32   ifgNpixelsml  = interferogram.win.pixels() / ifgmlP;
  const real8 offset = 0;

//cerr << "xNFO:   linesnoml:    " << Nlinesml <<  " pixnoml:    " <<  Npixelsml << endl;
//cerr << "xNFO:   ifglinesml: " << ifgNlinesml <<  " ifgpixml: " <<  ifgNpixelsml << endl;
  
  const real8 veryfirstline = real8(ifglinelo) + (real8(mlL)-1.0)/2.0;
  const real8 verylastline  = veryfirstline + real8((Nlinesml-1)*mlL);
  const real8 firstpixel    = real8(ifgpixlo) + (real8(mlP)-1.0)/2.0;
  const real8 lastpixel     = firstpixel + real8((Npixelsml-1)*mlP);

//cerr << "xNFO:   vl0:vlN " << veryfirstline << ":"  << verylastline << " p1:pN " <<  firstpixel << ":" << lastpixel << endl;

  //Determine range-azimuth spacing ratio, needed for proper triangulation
  cn P1, P2 , P3, P4;
  lp2xyz(veryfirstline,firstpixel,ellips,master,masterorbit,
           P1,MAXITER,CRITERPOS);
  lp2xyz(veryfirstline,lastpixel,ellips,master,masterorbit,
           P2,MAXITER,CRITERPOS);
  lp2xyz(verylastline,firstpixel,ellips,master,masterorbit,
           P3,MAXITER,CRITERPOS);
  lp2xyz(verylastline,lastpixel,ellips,master,masterorbit,
           P4,MAXITER,CRITERPOS);

  const real8 r_spacing  = ( (P1.min(P2)).norm() + (P3.min(P4)).norm() ) / 2 /(lastpixel - firstpixel) ;
  const real8 az_spacing = ( (P1.min(P3)).norm() + (P2.min(P4)).norm() ) /2 /(verylastline - veryfirstline); 
  const real8 r_az_ratio = r_spacing/az_spacing;

  INFO << "Interferogram azimuth spacing: " << az_spacing;
  INFO.print();
  INFO << "Interferogram range spacing: " << r_spacing;
  INFO.print();
  INFO << "Range-azimuth spacing ratio: " << r_az_ratio;
  INFO.print();

  // ______ Compute buffer size interpolation______
  const real8 Nlinesml_possible = BUFFERMEMSIZE / (6*8*Npixelsml);
  bufferlines = int32(ceil(Nlinesml_possible));  // initialized

  if ( mlookedIFG == true)                                // if ifg is multilooked by a factor
    {
      bufferlines = int32( floor(Nlinesml_possible/ifgmlL) * ifgmlL ); // [HB] Hermann provided the fix:
                                                                       // Successive bufferlines must have a size which is multiple of multilooking factor
                                                                       // unless data can fit completely to an initial single buffer. 
                                                                       // Extra buffer will scale correctly and rounding due to multilooking 
                                                                       // will yield correct number lines for the output file.
                                                                       // Ex: floor(2097/25)*2 + floor(1578/25)    =  229 lines (wrong)  (bufferlines/mlL)*Nbuffers+extrabufferlines/mlL 
    }                                                                  // Ex: floor(2075/25)*25*2 + floor(1622/25) =  230 lines (correct)
  else // no-multilooking
    {
      bufferlines = int32( floor(Nlinesml_possible) );    // [MA] instead of ceil, prefered floor to use less mem
    }
  if (bufferlines > Nlinesml) bufferlines=Nlinesml;  // if bufferlines > Nlines then shrink bufferlines to Nlines, no extra buffer requested.
  INFO << "Possible max. buffer lines: " << bufferlines << " for " << BUFFERMEMSIZE << " memory size.";
  INFO.print();
  numfullbuffers = Nlinesml / bufferlines;
  restlines      = Nlinesml % bufferlines;  // the number of lines in extra buffer
  extrabuffer    = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "Interpolation in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();

INFO << "OutputFile forefdem: "<< refdeminput.forefdem;
INFO.print();
  // ______ Open output files ______
  ofstream refdemofile;       // refdem phase
  openfstream(refdemofile,refdeminput.forefdem,generalinput.overwrit);
  bk_assert(refdemofile,refdeminput.forefdem,__FILE__,__LINE__);
  
  // _____ start added by MA _____
  ofstream refdemofilenoML;                                 // [MA] refdem phase no-multilooked
  { // local scope practice
    string fname = string(refdeminput.forefdem) + ".noML";  // new name as m_s_refdemphase.raw.noML
    if ( mlookedIFG == true)                                // if ifg is multilooked by a factor
      {
        openfstream(refdemofilenoML,fname.c_str(),generalinput.overwrit);
        bk_assert(refdemofilenoML,fname.c_str(),__FILE__,__LINE__);
      }
    else // no-multilooking
      {
        if(!remove(fname.c_str()))                           // when success report removed. 
          {
          WARNING << "Removed existing " << fname << "file.";
          WARNING.print();
          }
      }
  }
  // _____ end added by MA _____
  
  // if request for height in radar coordinates l,p
  ofstream refdemheiofile;// Radarcoded DEM (Z.Perski)
  if (outputrefdemhei==true)
    {
      INFO << "OutputFile forefdemhei: "<< refdeminput.forefdemhei;
      INFO.print();

    openfstream(refdemheiofile,refdeminput.forefdemhei,generalinput.overwrit);
    bk_assert(refdemheiofile,refdeminput.forefdemhei,__FILE__,__LINE__);
    }
  
  // if request for h2ph in radar coordinates l,p
  ofstream h2phofile;
  if (outputh2ph==true)
    {
    openfstream(h2phofile,refdeminput.foh2ph,generalinput.overwrit);
    bk_assert(h2phofile,refdeminput.foh2ph,__FILE__,__LINE__);
    }

  // ______ interpolation loop per buffer ______
  for (register int32 buffer = 0; buffer < numfullbuffers + extrabuffer; ++buffer)
    {

    // Determine indices for buffer
    const int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
    const real8 firstline_buffer = veryfirstline+buffer*bufferlines*mlL;
    const real8 lastline_buffer = firstline_buffer+(blines-1)*mlL;

    // ______ Extra info ______
    PROGRESS << STEP << "Interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " [l0:lN, p0:pN]: " << " ["
   << firstline_buffer << ": " << lastline_buffer << ", "
   << firstpixel << ": " << lastpixel << "]";
    PROGRESS.print();

    // Get corners of buffer
    real8 phimin_az;
    real8 phimax_az;
    real8 lambdamin_az;
    real8 lambdamax_az;
    getcorners(firstline_buffer+offset,lastline_buffer+offset,
             firstpixel+offset,lastpixel+offset,
             extralat,extralong,phimax,lambdamin,
             DEMdeltalat,DEMdeltalon,NrowsDEM,NcolsDEM,
             ellips,master,masterorbit,phimin_az,phimax_az,lambdamin_az,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0DEM,indexphiNDEM,
                        indexlambda0DEM,indexlambdaNDEM);
    const int32 NrowsDEM_buffer = indexphiNDEM-indexphi0DEM+1;
    const int32 NcolsDEM_buffer = indexlambdaNDEM-indexlambda0DEM+1;

    PROGRESS << STEP << "Reading input for interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer;
    PROGRESS.print();

    // read x,y
    matrix<real8> DEMline_buffer(NrowsDEM_buffer,NcolsDEM_buffer);
    matrix<real8> DEMpixel_buffer(NrowsDEM_buffer,NcolsDEM_buffer);

    readfile(DEMline_buffer,"crd_m_demline.temp",NrowsDEM,winfromfile,zerooffset);
    readfile(DEMpixel_buffer,"crd_m_dempixel.temp",NrowsDEM,winfromfile,zerooffset);
    
    // read z (multiple, number can easily be increased, e.g. simulated intensity)
    int32 Nz = 1; //number of z
    matrix<real8> input_buffer(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);
    matrix<real8> temp_input_buffer(NrowsDEM_buffer,NcolsDEM_buffer);
    if (outputrefdemhei==true)
      {
      Nz += 1;
      input_buffer.resize(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);
      }
    if (outputh2ph==true)
      {
      Nz += 1;
      input_buffer.resize(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);
      }

    readfile(temp_input_buffer,"crd_dem_refphase.temp",NrowsDEM,winfromfile,zerooffset);
    input_buffer.setdata(0, 0, temp_input_buffer);
    Nz = 1;
    if (outputrefdemhei==true)
      {
        Nz += 1;
        /// i would like to use real4, test later on
        matrix<real4> dem_input(NrowsDEM_buffer,NcolsDEM_buffer);
        readfile(dem_input,refdeminput.fodem,NrowsDEM,winfromfile,zerooffset);
        for (register int32 i =0 ; i < NrowsDEM_buffer ; i ++)
          for(register int32 j = 0; j < NcolsDEM_buffer; j++)
            temp_input_buffer(i,j) = real8(dem_input(i,j));
        input_buffer.setdata(NrowsDEM_buffer * (Nz-1), 0, temp_input_buffer);
      }
    if (outputh2ph==true)
      {
        Nz += 1;
        readfile(temp_input_buffer,"crd_dem_h2ph.temp",NrowsDEM,winfromfile,zerooffset);
        input_buffer.setdata(NrowsDEM_buffer * (Nz-1) , 0, temp_input_buffer);
      }
    
    // initialize output array
    Nz = 1;  
    matrix<real8> output_buffer(blines * Nz, Npixelsml);

    if (outputrefdemhei==true)
      {
        Nz += 1;
        output_buffer.resize(blines * Nz, Npixelsml);
      }
    if (outputh2ph==true)
      {
        Nz += 1;
        output_buffer.resize(blines * Nz, Npixelsml);
      }

     
    // interpolation
    griddatalinear(DEMline_buffer,DEMpixel_buffer,input_buffer,
                   firstline_buffer,lastline_buffer,firstpixel,lastpixel,
                   mlL,mlP,r_az_ratio,offset,NODATA,output_buffer);


    //MA multilooking will start here.
    //MA cast all output files to type real4 
    matrix<real4> output_layer(blines,Npixelsml); 

    Nz = 1;          // matrix 3rd dimension counter, 1 --> PHASE

    for (register int32 i =0 ; i < blines ; i++)
      for(register int32 j = 0; j < Npixelsml; j++)
        output_layer(i,j) = real4(output_buffer(i,j));   // real8 --> real4 
   //     convert_type(output_buffer,output_layer); // TODO MA, should replace above 3 lines but this one doesn't work properly yet

    cerr << "refphase: blines: " << blines << " Npixelsml " << Npixelsml << endl;
    INFO <<  "size buffer: lines:  " <<  output_layer.lines()<< ",  " << output_layer.pixels() <<endl;
    INFO.print();
  // _____ start added by MA _____
    if ( mlookedIFG == true)                             // [MA] if ifg is multilooked by a factor
      {
        refdemofile     << multilook(output_layer, ifgmlL, ifgmlP);    // multilook to output
        refdemofilenoML << output_layer;                               // default output: PHASE
      }
    else
      {
        refdemofile << output_layer;
        //refdemofile << output_buffer(window(0, blines - 1, 0, Npixelsml -1 ));
      }
  // _____ end added by MA _____

    if (outputrefdemhei==true)
      {
        Nz += 1;
        for (register int32 i = blines * (Nz-1) ; i < blines * Nz  ; i++)
          for(register int32 j = 0; j < Npixelsml; j++)
            output_layer(i-(blines * (Nz-1)),j) = real4(output_buffer(i,j)); // real8 --> real4
        //refdemheiofile << output_layer;            // output reference dem heights
        (mlookedIFG == true) ? refdemheiofile << multilook(output_layer, ifgmlL, ifgmlP)  // [MA]
                             : refdemheiofile << output_layer ;
        //refdemheiofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
      }
    if (outputh2ph==true)
      {
        Nz += 1;
        for (register int32 i = blines * (Nz-1)  ; i < blines * Nz  ; i++)
          for(register int32 j = 0; j < Npixelsml; j++)
            output_layer(i-(blines * (Nz-1)),j) = real4(output_buffer(i,j)); // real8 --> real4
        //h2phofile  << output_layer;            // output h2ph matrix
        (mlookedIFG == true) ? h2phofile  << multilook(output_layer, ifgmlL, ifgmlP)    // [MA]
                             : h2phofile  << output_layer;
        //h2phofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
      }
    
    DEMline_buffer.resize(1,1);           // deallocate
    DEMpixel_buffer.resize(1,1);
    input_buffer.resize(1,1);
    temp_input_buffer.resize(1,1);
    output_buffer.resize(1,1);

  } // end loop azimuth direction

  INFO << "Closing output files";
  INFO.print();

  refdemofile.close();          // has the same multilook as interferrogram
  if (mlookedIFG==true)
    refdemofilenoML.close();    // [MA] if interferogram is mlooked then this is 
                                // generated as non-multilooked for coherence estimation.
  if (outputrefdemhei==true)    // For Zbigniew Perski
    refdemheiofile.close();
  if (outputh2ph==true)
    h2phofile.close();

  //===================================================================
  //============ End second loop: interpolation           =============
  //============ (radar geometry)                         =============
  //===================================================================

  // === Clean up temporary outfiles === [MA]

  if (remove("crd_m_demline.temp"))                                        // remove files
    WARNING.print("code 101: could not remove file: crd_m_demline.temp.");
  if (remove("crd_m_dempixel.temp"))                                       
    WARNING.print("code 101: could not remove file: crd_m_dempixel.temp.");
  if (remove("crd_dem_refphase.temp"))                                        
    WARNING.print("code 101: could not remove file: crd_dem_refphase.temp.");
  if ( outputh2ph==true && remove("crd_dem_h2ph.temp"))   // output available                                    
    WARNING.print("code 101: could not remove file: crd_dem_h2ph.temp.");

  // === Clean up Done. ===

  // ====== Write output information ======
  char croppeddemi[4*ONE27];
  strcpy(croppeddemi,"NO output requested");
  if (outputdemi) strcpy(croppeddemi,refdeminput.fodemi);
  INFO << "Min. value of input DEM covering interferogram: " << min_input_dem; 
  INFO.print();
  INFO << "Max. value of input DEM covering interferogram: " << max_input_dem; 
  INFO.print();

  ofstream scratchlogfile("scratchlogcomprefdem", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"comprefdem: scratchlogcomprefdem",__FILE__,__LINE__);
  scratchlogfile
    << "\n*******************************************************************"
    << "\n* " << processcontrol[pr_i_comprefdem]
    << "\n*******************************************************************"
    << "\n1) DEM source file:                   \t" <<  refdeminput.firefdem
    << "\nFormat:                               \t";
    switch (refdeminput.iformatflag)
      {
      case FORMATI2:
        {
        scratchlogfile << "SHORT SIGNED INTEGER (HOST ENDIANNESS)";
        break;
        }
      case FORMATI2_BIGENDIAN:
        {
        scratchlogfile << "SHORT SIGNED INTEGER, BIG ENDIAN";
        break;
        }
      case FORMATR4:
        {
        scratchlogfile << "REAL4 SIGNED FLOAT";
        break;
        }
      case FORMATR8:
        {
        scratchlogfile << "REAL8 SIGNED DOUBLE";
        break;
        }
      default:
        {
        scratchlogfile << "UNKNOWN? IMPOSSIBLE...";
        break;
        }
      }
  scratchlogfile
    << "\nByte order:                           \t" <<  "check it yourself..."
    << "\nNumber of lines:                      \t" <<  numberoflatpixels
    << "\nNumber of pixels:                     \t" <<  numberoflonpixels
    << "\nResolution latitude:                  \t" << rad2deg(DEMdeltalat) << " [deg]"
    << "\nResolution longitude:                 \t" << rad2deg(DEMdeltalon) << " [deg]"
    << "\nMost West point in input DEM:         \t" << rad2deg(lon0file)
    << "\nMost East point in input DEM:         \t" << rad2deg(lonNfile)
    << "\nMost South point in input DEM:        \t" << rad2deg(latNfile)
    << "\nMost North point in input DEM:        \t" << rad2deg(lat0file)
    << "\nMin. value of input DEM covering interferogram: " << min_input_dem
    << "\nMax. value of input DEM covering interferogram: " << max_input_dem
    << "\n2) Output file cropped DEM:           \t" <<  refdeminput.fodem //[FVL
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):         \t" << NcolsDEM
    << "\nNumber of pixels (multilooked):        \t" << NrowsDEM
    << "\nDEM extend w/e/s/n            :      \t" << rad2deg(lambdamin) << "/" 
    << rad2deg(lambdamax) << "/" << rad2deg(phimin)    << "/" << rad2deg(phimax)
//    << "\nMean value:                           \t" <<  meancroppedDEM
    << "\n3) Output file interpolated crop DEM: \t" <<  croppeddemi
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\n4) Output file synthetic phase:       \t" <<  refdeminput.forefdem
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):        \t" <<  ifgNlinesml
    << "\nNumber of pixels (multilooked):       \t" <<  ifgNpixelsml

// this is not correct, only stats per buffer...
//    << "\n\n----- Other STATS -----"
//    << "\nTotal points in cropped DEM:          \t" << numpoints
//    << "\nNumber of valid points in DEM:        \t" << numvalid
//    << " (" << 100*numvalid/numpoints << "%)"
//    << "\nNumber of NODATA points in DEM:       \t" << numNODATA
//    << " (" << 100*numNODATA/numpoints << "%)"
//    << "\nMean height in meters at valid points:\t" << meancroppedDEM
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();


  ofstream scratchresfile("scratchrescomprefdem", ios::out | ios::trunc);
  bk_assert(scratchresfile,"comprefdem: scratchrescomprefdem",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_comprefdem]
    << "\n*******************************************************************";
  if (onlyrefphasetopo==true) scratchresfile                 // [don]
    << "\nInclude_flatearth:                 \tNo";
  else scratchresfile
    << "\nInclude_flatearth:                 \tYes";
  scratchresfile
    << "\nDEM source file:                      \t" << refdeminput.firefdem
    << "\nMin. of input DEM:                    \t" << min_input_dem
    << "\nMax. of input DEM:                    \t" << max_input_dem
    << "\nData_output_file:                     \t" <<  refdeminput.forefdem
    << "\nData_output_format:                   \t" << "real4"
    << "\nFirst_line (w.r.t. original_master):  \t"
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  interferogram.win.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  interferogram.win.pixhi
    << "\nMultilookfactor_azimuth_direction:    \t" <<  ifgmlL
    << "\nMultilookfactor_range_direction:      \t" <<  ifgmlP
    << "\nNumber of lines (multilooked):        \t" <<  ifgNlinesml
    << "\nNumber of pixels (multilooked):       \t" <<  ifgNpixelsml
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_comprefdem] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();


  } // END radarcodedem


/****************************************************************
 *    getcorners                                                *
 *                                                              *
 * Get corners of window (approx) to select DEM in radians (if  *
 * height were zero)                                            *
 *                                                              *
 * Implementation:                                              *
 * 1) calculate phi, lambda of corners                          *
 * 2) select the extreme values                                 *
 * 3) add extra overlap                                         *
 * 4) determine the indices in the file                         *
 *                                                              *
 *    Freek van Leijen, 07-AUG-2006                             *
 *                                                              *
 ****************************************************************/
void getcorners(
                const real8            &l0,                   // master.currentwindow.linelo
                const real8            &lN,                   // master.currentwindow.
                const real8            &p0,                   // master.currentwindow.
                const real8            &pN,                   // master.currentwindow.
                const real8            &extralat,             // 
                const real8            &extralong,            //
                const real8            &lat0,                 // lat0file
                const real8            &long0,                // lon0file
                const real8            &DEMdeltalat,          //
                const real8            &DEMdeltalong,         //
                const int32            &Nlatpixels,           // numberoflatpixels
                const int32            &Nlongpixels,          // numberoflonpixels
                const input_ell        &ellips,
                const slcimage         &master,
                orbit                  &masterorbit,
                real8            &phimin,
                real8            &phimax,
                real8            &lambdamin,
                real8            &lambdamax,
                int32            &indexphi0DEM,              // returned
                int32            &indexphiNDEM,              // returned
                int32            &indexlambda0DEM,           // returned
                int32            &indexlambdaNDEM)           // returned
  {
  TRACE_FUNCTION("getcorners (FvL 07-AUG-2006)")

  DEBUG << "(getcorners) l0 :" << l0;
  DEBUG.print();
  DEBUG << "(getcorners) lN :" << lN;
  DEBUG.print();
  DEBUG << "(getcorners) p0 :" << p0;
  DEBUG.print();
  DEBUG << "(getcorners) pN :" << pN;
  DEBUG.print();
  DEBUG << "(getcorners) extralat :" << extralat;
  DEBUG.print();
  DEBUG << "(getcorners) extralong :" << extralong;
  DEBUG.print();
  DEBUG << "(getcorners) lat0 :" << lat0;
  DEBUG.print();
  DEBUG << "(getcorners) long0 :" << long0;
  DEBUG.print();
  DEBUG << "(getcorners) DEMdeltalat :" << DEMdeltalat;
  DEBUG.print();
  DEBUG << "(getcorners) DEMdeltalong :" << DEMdeltalong;
  DEBUG.print();
  DEBUG << "(getcorners) Total Nlatpixels :" << Nlatpixels;
  DEBUG.print();
  DEBUG << "(getcorners) Total Nlongpixels :" << Nlongpixels;
  DEBUG.print();
  
  real8 phi;
  real8 lambda;
  real8 height;  
  lp2ell(l0,p0,
         ellips, master, masterorbit,
         phi, lambda, height);          // returned
  real8 phil0p0    = phi;
  real8 lambdal0p0 = lambda;

  lp2ell(lN,p0,
         ellips, master, masterorbit,
         phi, lambda, height);          // returned
  real8 philNp0    = phi;
  real8 lambdalNp0 = lambda;

  lp2ell(lN,pN,
         ellips, master, masterorbit,
         phi, lambda, height);          // returned
  real8 philNpN    = phi;
  real8 lambdalNpN = lambda;

  lp2ell(l0,pN,
         ellips, master, masterorbit,
         phi, lambda, height);          // returned
  real8 phil0pN    = phi;
  real8 lambdal0pN = lambda;

  // ______ Select DEM values based on rectangle outside l,p border ______
  phimin    = min(min(min(phil0p0,philNp0),philNpN),phil0pN);
  phimax    = max(max(max(phil0p0,philNp0),philNpN),phil0pN);
  lambdamin = min(min(min(lambdal0p0,lambdalNp0),lambdalNpN),lambdal0pN);
  lambdamax = max(max(max(lambdal0p0,lambdalNp0),lambdalNpN),lambdal0pN);

  // ______ a little bit extra at edges to be sure ______ 
  // TODO this should change based on ascending or descending [MA] 
  phimin    -= extralat;
  phimax    += extralat;
  lambdamax += extralong;
  lambdamin -= extralong;

  DEBUG << "(getcorners) phimin :" << phimin;
  DEBUG.print();
  DEBUG << "(getcorners) phimax :" << phimax;
  DEBUG.print();
  DEBUG << "(getcorners) lambdamin :" << lambdamin;
  DEBUG.print();
  DEBUG << "(getcorners) lambdamax :" << lambdamax;
  DEBUG.print();

  // ______ Get indices of DEM needed ______
  // ______ Index boundary: [0:numberofx-1] ______
 
  indexphi0DEM = int32(floor((lat0-phimax)/DEMdeltalat));
  if (indexphi0DEM < 0) 
    {
    WARNING << "(getcorners) indexphi0DEM: " << indexphi0DEM; WARNING.print();
    indexphi0DEM=0;   // default start at first
    WARNING.print("DEM does not cover entire interferogram.");
    WARNING.print("input DEM should be extended to the North.");
    }
  indexphiNDEM = int32(ceil((lat0-phimin)/DEMdeltalat));
  if (indexphiNDEM > Nlatpixels-1)
    {
      WARNING << "(getcorners) indexphiNDEM: " << indexphiNDEM; WARNING.print();
    indexphiNDEM=Nlatpixels-1;
    WARNING.print("DEM does not cover entire interferogram.");
    WARNING.print("input DEM should be extended to the South.");
    }
  indexlambda0DEM = int32(floor((lambdamin-long0)/DEMdeltalong));
  if (indexlambda0DEM < 0)
    {
      WARNING << "(getcorners) indexlambda0DEM: " << indexlambda0DEM; WARNING.print();
    indexlambda0DEM=0;    // default start at first
    WARNING.print("DEM does not cover entire interferogram.");
    WARNING.print("input DEM should be extended to the West.");
    }
  indexlambdaNDEM = int32(ceil((lambdamax-long0)/DEMdeltalong));
  if (indexlambdaNDEM > Nlongpixels-1)
    {
      WARNING << "(getcorners) indexlambdaNDEM: " << indexlambdaNDEM; WARNING.print();
    indexlambdaNDEM=Nlongpixels-1;
    WARNING.print("DEM does not cover entire interferogram.");
    WARNING.print("input DEM should be extended to the East.");
    }
  
    DEBUG << "(getcorners) indexphi0DEM :" << indexphi0DEM;
    DEBUG.print();
    DEBUG << "(getcorners) indexphiNDEM :" << indexphiNDEM;
    DEBUG.print();
    DEBUG << "(getcorners) indexlambda0DEM :" << indexlambda0DEM;
    DEBUG.print();
    DEBUG << "(getcorners) indexlambdaNDEM :" << indexlambdaNDEM;
    DEBUG.print();


  } // END getcorners


/****************************************************************
 *    griddatalinear (naming after Matlab function)             *
 *                                                              *
 *    Implementation after GMT function triangulate.c           *
 ****************************************************************/
  void griddatalinear(
          const matrix<real8>  &x_in,
          const matrix<real8>  &y_in,
          const matrix<real8>  &z_in,
          const real8          &x_min,
                      const real8          &x_max,
                      const real8          &y_min,
                      const real8          &y_max,
          const int32          &x_inc,
          const int32          &y_inc,
          const real8          &r_az_ratio,
          const real8          &offset,
          const real8          &NODATA,
          matrix<real8>        &grd
          )
{
  TRACE_FUNCTION("griddatalinear (LG&FvL 13-AUG-2006)")
    
  INFO << "griddataLinear interpolation.";
  INFO.print();

  int32 i, j, k, ij, p, i_min, i_max, j_min, j_max;
  int32 n, nx, ny, zLoops, zLoop, zBlockSize,indexFirstPoint,zInterpolateBlockSize;
  real8 vx[4], vy[4], xkj, xlj, ykj, ylj, zj, zk, zl, zlj, zkj, xp, yp;
  real8 f, *a = NULL , *b = NULL , *c = NULL ;    // linear interpolation parameters
  struct triangulateio In, Out, vorOut;
  
  
  // Initialize variables
  
  n   = zBlockSize = x_in.size();    // block size of x and y coordination 
  
  /* How many groups of z value should be interpolated */
  if (( z_in.size() % zBlockSize ) != 0)
    {
      INFO << "The input of the DEM buffer and z is not the same...";
      INFO.print();
      return;
    }
  else
    zLoops = z_in.size()/x_in.size();
  
  a   = new real8[zLoops];
  b   = new real8[zLoops];
  c   = new real8[zLoops];
  if  (a  ==  NULL  ||  b ==  NULL  ||  c ==  NULL)
    {
      ERROR << "Memory ERROR in source file: " << __FILE__
      << " at line: " << __LINE__;
      PRINT_ERROR(ERROR.get_str());
      throw(memory_error);
    }
  
  nx  = grd.lines()/zLoops;
  ny  = grd.pixels();
  zInterpolateBlockSize = grd.size()/zLoops;

  /* Set everything to 0 and NULL */
  memset ((void *)&In,   0, sizeof (struct triangulateio));
  memset ((void *)&Out, 0, sizeof (struct triangulateio));
  memset ((void *)&vorOut, 0, sizeof (struct triangulateio));
  
  /* Allocate memory for input points */
  In.numberofpoints = n ;
  In.pointlist = new real8 [2 * n];
  
  /* Copy x,y points to In structure array */
  
  for (i = j = 0; i < n; i++) 
    {
      In.pointlist[j++] = *(x_in[0] + i);
      In.pointlist[j++] = (*(y_in[0] + i)) * r_az_ratio ; 
      // to eliminate the effect of difference in range and azimuth spacing;
    }
  
  /* Call Jonathan Shewchuk's triangulate algorithm.  This is 64-bit safe since
   * all the structures use 4-byte ints (longs are used internally). */
  
  triangulate ("zIQB", &In, &Out, &vorOut);
  //triangulate ("zIBV", &In, &Out, &vorOut); // [MA] for verbosing
  
  int32 *link = Out.trianglelist; /* List of node numbers to return via link */
  int32  np   = Out.numberoftriangles;

  for (k = ij = 0; k < np; k++) 
    {
      DEBUG << "k of np, ij: " << k << " of " << np << ", :" << ij;
      DEBUG.print();
      //Store the Index of the first Point of this triangle.
      indexFirstPoint = ij;
      
      vx[0] = vx[3] = *(x_in[0] + link[ij]);  vy[0] = vy[3] = *(y_in[0] + link[ij]); ij++;
      vx[1] = *(x_in[0] + link[ij]);  vy[1] = *(y_in[0]+link[ij]); ij++;
      vx[2] = *(x_in[0] + link[ij]);  vy[2] = *(y_in[0]+link[ij]); ij++;
      
      if ( vx[0] == NODATA || vx[1] == NODATA || vx[2] == NODATA ) continue;
      if ( vy[0] == NODATA || vy[1] == NODATA || vy[2] == NODATA ) continue;
      
      /* Compute grid indices the current triangle may cover.*/
      xp = min (min (vx[0], vx[1]), vx[2]); i_min = x_to_i (xp, x_min, x_inc, offset, nx);
      //INFO << "xp: " << xp;
      //INFO.print();
      xp = max (max (vx[0], vx[1]), vx[2]); i_max = x_to_i (xp, x_min, x_inc, offset, nx);
      //INFO << "xp: " << xp;
      //INFO.print();
      yp = min (min (vy[0], vy[1]), vy[2]); j_min = y_to_j (yp, y_min, y_inc, offset, ny);
      //INFO << "yp: " << yp;
      //INFO.print();
      yp = max (max (vy[0], vy[1]), vy[2]); j_max = y_to_j (yp, y_min, y_inc, offset, ny);
      //INFO << "yp: " << yp;
      //INFO.print();
      /* Adjustments for triangles outside -R region. */
      /* Triangle to the left or right. */
      if ((i_max < 0) || (i_min >= nx)) continue;
      /* Triangle Above or below */
      if ((j_max < 0) || (j_min >= ny)) continue;
      /* Triangle covers boundary, left or right. */
      if (i_min < 0) i_min = 0;       if (i_max >= nx) i_max = nx - 1;
      /* Triangle covers boundary, top or bottom. */
      if (j_min < 0) j_min = 0;       if (j_max >= ny) j_max = ny - 1;
      // for (kk = 0; kk<npar;kk++) {  //do for each parameter LIUG
      // read zj, zk, zl (instead of above) LIUG
      /* Find equation for the plane as z = ax + by + c */
      xkj = vx[1] - vx[0]; ykj = vy[1] - vy[0]; 
      xlj = vx[2] - vx[0]; ylj = vy[2] - vy[0]; 
      
      f = 1.0 / (xkj * ylj - ykj * xlj);
      
      for( zLoop = 0 ; zLoop < zLoops; zLoop++ )
  {
    zj = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint]);
    zk = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint + 1]);
    zl = *(z_in[0] + zLoop * zBlockSize + link[indexFirstPoint + 2]);
    zkj = zk - zj;  zlj = zl - zj;
    a[zLoop] = -f * (ykj * zlj - zkj * ylj);
    b[zLoop] = -f * (zkj * xlj - xkj * zlj);
    c[zLoop] = -a[zLoop] * vx[1] - b[zLoop] * vy[1] + zk;
  }
      
      for (i = i_min; i <= i_max; i++) 
        {
    xp = i_to_x (i, x_min, x_max, x_inc, offset, nx);
    p = i * ny + j_min;
    for (j = j_min; j <= j_max; j++, p++) 
            {
              yp = j_to_y (j, y_min, y_max, y_inc, offset, ny);
              if (!pointintriangle(vx, vy, xp, yp)) continue; /* Outside */

        for( zLoop = 0 ; zLoop < zLoops; zLoop++ )
    *(grd[0] + zLoop * zInterpolateBlockSize + p) =  a[zLoop] * xp + b[zLoop] * yp + c[zLoop] ;
            } //LIUG
        }
    }
  
  if(a) delete[]  a;
  if(b) delete[]  b;
  if(c) delete[]  c;
  if(In.pointlist)  delete[]  In.pointlist; //  only this use delete
  if(Out.pointlist) free(Out.pointlist);
  if(Out.trianglelist) free(Out.trianglelist);
  
  
} //END griddatalinear


/****************************************************************
 *    pointintriangle                                           *
 *                                                              *
 *    Liu, Guang and Freek van Leijen, 16-AUG-2006              *
 *                                                              *
 *                                                              *
 ****************************************************************/
int pointintriangle(double *xt,double * yt,double x ,double y)
{
  
  int iRet0 = ((xt[2] - xt[0]) * (y - yt[0])) >  ((x - xt[0]) * ( yt[2] - yt[0])) ? 1:-1; 
  int iRet1 = ((xt[0] - xt[1]) * (y - yt[1])) >  ((x - xt[1]) * ( yt[0] - yt[1])) ? 1:-1;
  int iRet2 = ((xt[1] - xt[2]) * (y - yt[2])) >  ((x - xt[2]) * ( yt[1] - yt[2])) ? 1:-1;
  
  if ((iRet0 >0 && iRet1 > 0 && iRet2 > 0 ) || (iRet0 <0 && iRet1 < 0 && iRet2 < 0 ))
    return 1;
  else 
    return 0;
  
} //END pointintriangle

