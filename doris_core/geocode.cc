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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/geocode.cc,v $    *
 * $Revision: 3.16 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of computation of endproducts (DEM, defo.map) *
 * -slant range 2 height (schwabisch)                           *
 * -slant range 2 height (rodriguez, exact, some problems)      *
 * -slant range 2 height (ambiguity)                            *
 * -geocode heightmatrix to phi,lambda                          *
 ****************************************************************/


#include "constants.hh"                 // global constants
#include "matrixbk.hh"                  // my matrix class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class
#include "geocode.hh"                   // header file
#include "utilities.hh"                 // utils
#include "ioroutines.hh"                // error messages
#include "coregistration.hh"            // distributepoints
#include "exceptions.hh"                 // my exceptions class

#include <strstream>                    // for memory stream
#include <cctype>                       // isspace
#include <algorithm>                    // max



/****************************************************************
 *    slant2h eight (schwabisch)                                *
 *                                                              *
 * compute height in radar coded system (master):               *
 *                                                              *
 * Input:                                                       *
 *  -                                                           *
 * Output:                                                      *
 *  -                                                           *
 *                                                              *
 * See thesis swabisch for method                               *
 * 1.  compute reference phase for h=0,2000,4000 in Npoints     *
 * 2.  solve system: h(phi) = a_0 + a_1*phi + a_2*phi*phi       *
 *      (2nd degree 1D polynomial) for all Npoints              *
 * 3.  compute a_i (l,p) = DEGREE2D 2D polynomial               *
 * 4.0 set offset to the one of first pixel , add this to all   *
 *     this step is skipped, phase is w.r.t. h=0, ref. is subtracted *
 * 4.1 evaluate polynomial of 3. for all points (l,p) of        *
 *      (multilooked) unwrapped interferogram                   *
 * solution to system for betas seems not very stable?          *
 *                                                              *
 *    Bert Kampes, 02-Jun-1999                                  *
 ****************************************************************/
void slant2hschwabisch(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage      &master,
        const slcimage      &slave,
        const productinfo   &unwrappedinterf,
        orbit               &masterorbit,
        orbit               &slaveorbit)
  {
  TRACE_FUNCTION("slant2hschwabisch (BK 02-Jun-1999)")

  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;
  const real8 m_minpi4cdivlambda = (-4.*PI*SOL)/master.wavelength;
  const real8 s_minpi4cdivlambda = (-4.*PI*SOL)/slave.wavelength;

  const int32 Npoints   = slant2hinput.Npoints;         // where ref.phase is evaluated
  const int32 DEGREE1D  = slant2hinput.degree1d;        // only possible now.
  const int32 DEGREE2D  = slant2hinput.degree2d;
  const int32 Nheights  = slant2hinput.Nheights;

  const int32 MAXHEIGHT = 5000;                         // max hei for ref.phase
  const int32 TEN       = 10;                           // used in pointer
  if (DEGREE1D - 1 > TEN)
    {
    PRINT_ERROR("panic, programmers problem: increase TEN.")
    throw(some_error);
    }

  const int32 HEIGHTSTEP = MAXHEIGHT / (Nheights - 1);  // heights to eval ref.phase

  // ______ Matrices for storing phase for all ref. ellipsoids ______
  // ______ PHASE(i,0)  phase for height 0
  // ______ PHASE(i,1)  phase for height Heigthsep * 1
  // ______ PHASE(i,Nh) phase for height 4000
  matrix<real8>         PHASE(Npoints,Nheights);        // pseudo-observation


  // ______ Distribute points in original master system (not multilooked) ______
  // ______ (i,0): line, (i,1): pixel, (i,2) flagfromdisk (not used here) ______
  //matrix<uint> Position = distributepoints(Npoints,unwrappedinterf.win);
  // [FvL] for correct folding of points outside overlap window when inserted by file
  matrix<int> Position = distributepoints(Npoints,unwrappedinterf.win);


// ====== STEP 1 ======
// ====== Compute reference phase in N points for height (numheight) ======
  PROGRESS.print("S2H: schwabisch: STEP1: compute reference phase for Nheights.");
  register int32 numheight;
  register int32 i,j,k,l,index;
  register int32 alfa;
  for (numheight=0; numheight<Nheights; numheight++)
    {
    cn pospoint;
    const int32 HEIGHT = numheight * HEIGHTSTEP;
    const input_ell ELLIPS(ellips.a+HEIGHT, ellips.b+HEIGHT);
    
    // ______ Compute delta r for all points ______
    for (i=0;i<Npoints;i++)
      {
      const real8 line  = Position(i,0);
      const real8 pixel = Position(i,1);
      const real8 m_trange = master.pix2tr(pixel);
  
      // ______ Compute xyz of point P on ELLIPS for this line,pixel ______
      lp2xyz(line,pixel,ELLIPS,master,masterorbit,
             pospoint,MAXITER,CRITERPOS);
  
      // ______ Compute xyz of slave satelite in orbit_slave from P ______
      real8 s_tazi;                             // returned
      real8 s_trange;                           // returned, unused?
      xyz2t(s_tazi,s_trange,slave, slaveorbit,
            pospoint,MAXITER,CRITERTIM);
  
      // ______ Compute delta range ~= phase, store in matrix ______
      // ______ real8 dr = dist(m_possat,pospoint) - dist(s_possat,pospoint);
      // ______ real8 phase = -pi4*(dr/LAMBDA);
      // ______  dr    == M-S   cause if no flatearth M-S - flatearth = M-S-(M-S)=0
      // ______  phase == -4pi*dr/lambda == 4pi*(S-M)/lambda
      PHASE(i,numheight) = m_minpi4cdivlambda*m_trange-
                           s_minpi4cdivlambda*s_trange;
      }
    }


  // ______ Subtract ref. phase at h=0 for all point ______
  // ______ this is the same as adding reference phase for all in uint ______
  // ______ BK tested 4/oct/99 ______
  // ______ This possibly needs to be CHANGED later for other situations ______
  for (i=0; i<Npoints; ++i)
    {
    real8 offset = PHASE(i,0);
    PHASE(i,0) = 0.;
    for (j=1; j<numheight; ++j)
      {
      PHASE(i,j) -= offset;
      }
    }



// ====== STEP 2 ======
  PROGRESS.print("S2H: schwabisch: STEP2: estimate coefficients 1d polynomial.");
  // ====== Compute alpha coefficients of polynomials for these points ======
  // ______ h = sum (alpha_i phi^i); ______
  // ______ bk 26/10/99 rescale phi -> phi[0,1] ______
  matrix<real8>         DESIGN(Nheights,DEGREE1D+1);    // design matrix
  matrix<real8>         ALPHAS(Npoints,DEGREE1D+1);     // pseudo-observation
  matrix<real8>         HEI(Nheights,1);
  for (i=0; i<Nheights; i++)
    HEI(i,0) = i*HEIGHTSTEP;                            // 0, .., 5000

  // ______ normalize data to [0,1] ______
  const real8 minphi = min(PHASE);
  const real8 maxphi = max(PHASE);
  normalize(PHASE,minphi,maxphi);                       // regrid data

  for (i=0; i<Npoints; i++)     // solve system for all points
    {

    // ______ Set up design matrix ______
    for (j=0; j<Nheights; j++)
      for (k=0; k<=DEGREE1D; k++)
        DESIGN(j,k) = pow(PHASE(i,j),real8(k));                 // PHASE is normalized


    // ______ Solve by cholesky (even the exactly determined case) ______
    matrix<real8> N   = matTxmat(DESIGN,DESIGN);
    matrix<real8> rhs = matTxmat(DESIGN,HEI);
    matrix<real8> Qx_hat = N;
    choles(Qx_hat);             // Cholesky factorisation normalmatrix
    solvechol(Qx_hat,rhs);      // Estimate of unknowns (alphas) in rhs


    // ______ Test inverse ______
    invertchol(Qx_hat);         // Covariance matrix (lower)
    for (uint k=0; k<Qx_hat.lines(); k++)
      for (uint j=0; j<k; j++)
        Qx_hat(j,k) = Qx_hat(k,j);              // was only stored in lower
    const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
    INFO << "s2h schwaebisch: max(abs(N*inv(N)-I)) = " << maxdev;
    INFO.print();
    if (maxdev > .01)
      WARNING.print("wrong solution for 1d polynomial? (decrease d1d or nhei)");
    // ______ Scale back unknowns: alpha_i <= alpha_i * (scale)^i ______
    // ______ Store solution in ALPHAS ______
    for (alfa=0; alfa<=DEGREE1D; alfa++)
      ALPHAS(i,alfa) = rhs(alfa,0);
    } // loop over all points



// ====== STEP 3 ======
  PROGRESS.print("S2H: schwabisch: STEP3: estimate coefficients for 2d polynomial.");
  // ====== Compute alpha_i coefficients of polynomials as function of (l,p) ======
  // ______ alpha_i = sum(k,l) beta_kl l^k p^l; ______
  // ______ Solve simultaneous for all betas ______
  // ______ this does not seem to be possibly with my routine, so do per alfa_i ______
  const int32 Nunk     = Ncoeffs(DEGREE2D);              // Number of unknowns

  // ______ Check redundancy is done before? ______
  if (Npoints < Nunk)
    {
    PRINT_ERROR("slant2hschwabisch: N_observations<N_unknowns (increase S2H_NPOINTS or decrease S2H_DEGREE2D.")
    throw(input_error);
    }

  matrix<real8>         A(Npoints,Nunk);                // designmatrix


  // ______ Set up system of equations ______
  // ______ Order unknowns: B00 B10 B01 B20 B11 B02 B30 B21 B12 B03 for degree=3 ______
  const real8 minL = min(Position.getcolumn(0)); 
  const real8 maxL = max(Position.getcolumn(0)); 
  const real8 minP = min(Position.getcolumn(1)); 
  const real8 maxP = max(Position.getcolumn(1)); 
  for (i=0; i<Npoints; i++)
    {
    // ______ normalize coordinates ______
    const real8 posL = normalize(real8(Position(i,0)),minL,maxL);
    const real8 posP = normalize(real8(Position(i,1)),minP,maxP);
   
    index = 0;
    for (j=0; j<=DEGREE2D; j++)
      {
      for (k=0; k<=j; k++)
        {
        A(i,index) = pow(posL,real8(j-k)) * pow(posP,real8(k));
        index++;
        }
      }
    }

  // ______ Solve 2d polynomial system for alfas at these points _____
  matrix<real8> N   = matTxmat(A,A);
  matrix<real8> rhs = matTxmat(A,ALPHAS);
  matrix<real8> Qx_hat = N;
  choles(Qx_hat);                       // Cholesky factorisation normalmatrix
/////  solvechol(Qx_hat,rhs);        // Estimate of unknowns (betas) in rhs, NOT OK!


  // ______ Solve the normal equations for all alpha_i ______
  // ______ Simultaneous solution doesn't work somehow ______
  for (uint i=0; i<rhs.pixels(); ++i)
    {
    matrix<real8> rhs_alphai = rhs.getcolumn(i);
    solvechol(Qx_hat,rhs_alphai);       // Solution in rhs_alphai
    rhs.setcolumn(i,rhs_alphai);        // place solution back
    }


  // ______ Test solution by inverse ______
  invertchol(Qx_hat);                   // Covariance matrix (lower)
  for (uint i=0; i<Qx_hat.lines(); i++)
    for (uint j=0; j<i; j++)
      Qx_hat(j,i) = Qx_hat(i,j);        // was only stored in lower
  const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
  INFO << "s2h schwaebisch: max(abs(N*inv(N)-I)) = " << maxdev;
  INFO.print();
  if (maxdev > 0.01)
    {
    WARNING << "slant2h: possibly wrong solution. deviation from unity AtA*inv(AtA) = "
         << maxdev << " > 0.01";
    WARNING.print();
    }



// ====== STEP 4 ======
  PROGRESS.print("S2H: schwabisch: STEP4: compute height for all pixels.");
  // ====== Evaluate for all points interferogram h=f(l,p,phase) ======
  // ______ recon with multilook, degree1D, degree2D free
  // ______ Multilook factors ______
  const real8 multiL = unwrappedinterf.multilookL;
  const real8 multiP = unwrappedinterf.multilookP;


  // ______ Number of lines/pixels of multilooked unwrapped interferogram ______
  const int32 mllines = int32(floor(real8(
    unwrappedinterf.win.linehi-unwrappedinterf.win.linelo+1) / multiL));
  const int32 mlpixels = int32(floor(real8(
    unwrappedinterf.win.pixhi-unwrappedinterf.win.pixlo+1)   / multiP));


  // ______ Line/pixel of first point in original master coordinates ______
  const real8 veryfirstline = real8(unwrappedinterf.win.linelo) +
                                (real8(multiL) - 1.) / 2.;
  const real8 firstpixel    = real8(unwrappedinterf.win.pixlo)  +
                                (real8(multiP) - 1.) / 2.;


  // ______ Constant axis of pixel coordinates ______
  matrix<real4> p_axis(mlpixels,1);
  for (i=0; i<mlpixels; i++)
  //    p_axis(i,0) = real8((firstpixel + i*multiP) - minP) / (maxP - minP);
    p_axis(i,0) = firstpixel + i*multiP;
  normalize(p_axis,minP,maxP);

  const int32 NUMMAT = 1 + (1 + DEGREE1D);              // number of heavy matrices
  //int32 bufferlines = generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)));
  int32 bufferlines = int32(ceil( real8( generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)) )) )); // [MA]
  if (bufferlines > mllines)                            // whole image fits in BUFFER
    bufferlines = mllines;

  const int32 FULLBUFFERS = mllines / bufferlines;
  const int32 RESTLINES   = mllines % bufferlines;
  const int32 EXTRABUFFER = RESTLINES ? 1 : 0;
  
  // ______ Window to be read into BUFFER from file in multilooked system ______
  const uint dummy = 999999;                    // large to force error if not ok
  window bufferwin(1, bufferlines, 1, mlpixels);        // initial
  window offsetbuffer(1,dummy,1,dummy); // dummy not used in readfile, no offset

  // ______ Open output file ______
  ofstream ofile;
  openfstream(ofile,slant2hinput.fohei,generalinput.overwrit);
  bk_assert(ofile,slant2hinput.fohei,__FILE__,__LINE__); 


  // ====== Process BUFFERS ======
  for (register int32 buffer=1; buffer<=FULLBUFFERS+EXTRABUFFER; buffer++)
    {

    // ______ Give progress ______
    PROGRESS << "SLANT2H: " 
         << 100*(buffer-1)/(FULLBUFFERS+EXTRABUFFER) << "%";
    PROGRESS.print();

    // ______ In original master coordinate system ______
    const real8 firstline = veryfirstline + (buffer-1) * bufferlines * multiL;

    // ______ Set indices for loading / check last buffer ______
    bufferwin.linelo = 1 + (buffer-1) * bufferlines;    // Update window to be read from file
    if (buffer == FULLBUFFERS+1)
      {
      bufferlines = RESTLINES;
      //BUFFER.resize(bufferlines,mlpixels);
      }
    bufferwin.linehi = bufferwin.linelo + bufferlines - 1; // window 2b read from file

    // ______ Read in phase buffer of unwrapped interferogram ______
    matrix<real4> BUFFER = unwrappedinterf.readphase(bufferwin);


    // ______ Evaluate polynomial coefficients for these points ______
    // ______ Compute first line of current buffer in master coordinates ______
    matrix<real4> l_axis(bufferlines,1);
    for (k=0; k<bufferlines; k++)
      l_axis(k,0) = firstline + k*multiL;
    normalize(l_axis,minL,maxL);

    

    // ______ Lookup table because not known in advance what DEGREE1D is ______
    // BK: usage:
    //  pntALPHA[0]->showdata();
    //  (*pntALPHA[0][0]).showdata();
    // containing a grid of 1D coefficients
    matrix<real4> *pntALPHA[TEN];
    for (k=0; k<=DEGREE1D; k++)
      {
      matrix<real8> beta(Ncoeffs(DEGREE2D),1);
      for (l=0; l<Ncoeffs(DEGREE2D); l++)
        {
        beta(l,0) = rhs(l,k);           // solution stored in rhs
        }
      pntALPHA[k] = new matrix<real4> (l_axis.lines(),p_axis.lines());
      (*pntALPHA[k]) = polyval<real4>(l_axis, p_axis, beta, DEGREE2D); // [MA]
      }

    // ______ Evaluate h=f(l,p,phi) for all points in grid in BUFFER ______
    matrix<real8> coeff_thispoint(DEGREE1D+1,1);
    for (uint line=0; line<BUFFER.lines(); line++)
      {
      for (uint pixel=0; pixel<BUFFER.pixels(); pixel++)
        {
        // ______ Check if unwrapped ok, else compute h ______
        if (BUFFER(line,pixel) != NaN)          // else leave NaN
          {
          for (k=0; k<DEGREE1D+1; k++)
            {
            coeff_thispoint(k,0) = (*pntALPHA[k])[line][pixel];
            }
          BUFFER(line,pixel) = polyval1d(normalize(
                    BUFFER(line,pixel),minphi,maxphi),  // regrid phi [0,1]
                    coeff_thispoint);
          }
        }
      }

    // ______ Write computed heights to file ______
    ofile << BUFFER;
    // ______ new Matrix should be deleted ______
    // ______ if errors occur sigsegv maybe because of this ______
    DEBUG.print("deleting new matrix, memory errors could be caused by this");
    for (k=0; k<=DEGREE1D; k++)
      delete pntALPHA[k];// correct?
    }// loop over BUFFERS
  ofile.close();



// ====== Write to result files ======
  ofstream scratchlogfile("scratchlogslant2h", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"slant2h: scratchlogslant2h",__FILE__,__LINE__);
  scratchlogfile << "\n\n*******************************************************************"
                 << "\n* " << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\tschwabisch"
                 << "\nData_output_file: \t\t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nEllipsoid (name,a,b): \t\t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b 
                 << endl;
  scratchlogfile.close();


  ofstream scratchresfile("scratchresslant2h", ios::out | ios::trunc);
  bk_assert(scratchresfile,"slant2h: scratchresslant2h",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\tschwabisch"
                 << "\nData_output_file:                     \t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nFirst_line (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.linelo
                 << "\nLast_line (w.r.t. original_master):   \t"
                 <<  unwrappedinterf.win.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  unwrappedinterf.win.pixlo
                 << "\nLast_pixel (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.pixhi
                 << "\nMultilookfactor_azimuth_direction:    \t"
                 <<  unwrappedinterf.multilookL
                 << "\nMultilookfactor_range_direction:      \t"
                 <<  unwrappedinterf.multilookP

                 << "\nEllipsoid (name,a,b):                 \t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b
                 << "\n*******************************************************************"
                 //<< "\n* End_slant2h:_NORMAL"
                 << "\n* End_" << processcontrol[pr_i_slant2h] << "_NORMAL"
                 << "\n*******************************************************************\n";

  scratchresfile.close();
  PROGRESS.print("finished slant2hschwabisch.");
  } // END slant2h



/****************************************************************
 *    slant2h ambiguity method                                  *
 *                                                              *
 * compute height in radar coded system (master):               *
 * use Bhor Bver                                                *
 * (constant per line for parallel orbits, otherwize nearly)    *
 * and transformation model to get slave position               *
 *                                                              *
 * Input:                                                       *
 *  -                                                           *
 * Output:                                                      *
 *  -                                                           *
 *                                                              *
 *    Bert Kampes, 07-Jul-1999                                  *
 ****************************************************************/
void slant2hambiguity(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage      &master,
        const slcimage      &slave,
        const productinfo   &unwrappedinterf,
        orbit               &masterorbit,
        orbit               &slaveorbit,
        const BASELINE      &baseline)
  {
  TRACE_FUNCTION("slant2hambiguity (BK 07-Jul-1999)")

  const int32 MAXITER   = 10;                   // iterations for lp2xyz
  const real8 CRITERPOS = 1e-6;                 // stop criterium for lp2xyz
  const real8 CRITERTIM = 1e-10;                // stop criterium for lp2xyz

  const int32 MAXITERHERE = 4;                  // iterations for h
  const real8 CRITERHERE  = 0.05;               // m (delta h in iterations)


  // ______ Multilook factors ______
  const real8 multiL = unwrappedinterf.multilookL;
  const real8 multiP = unwrappedinterf.multilookP;

  // ______ Number of lines/pixels of multilooked unwrapped interferogram ______
  const int32 mllines = int32(floor(real8(
    unwrappedinterf.win.linehi-unwrappedinterf.win.linelo+1) / multiL));
  const int32 mlpixels = int32(floor(real8(
    unwrappedinterf.win.pixhi-unwrappedinterf.win.pixlo+1)   / multiP));

  // ______ Line/pixel of first point in original master coordinates ______
  const real8 veryfirstline = real8(unwrappedinterf.win.linelo) +
                                (real8(multiL) - 1.) / 2.;
  const real8 firstpixel    = real8(unwrappedinterf.win.pixlo)  +
                                (real8(multiP) - 1.) / 2.;


  // ====== Compute number of buffers required ======
  // BK 8/10/99: why use buffer? just read in line by line is as efficient...
  const int32 NUMMAT = 3;       // number of large matrices BUFFER PHI LAMBDA
  //int32 bufferlines = generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)));
  int32 bufferlines = int32(ceil( real8( generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)) )) ));
  if (bufferlines > mllines)                            // whole image fits in BUFFER
    bufferlines = mllines;

  const int32 FULLBUFFERS = mllines / bufferlines;
  const int32 RESTLINES   = mllines % bufferlines;
  const int32 EXTRABUFFER = RESTLINES ? 1 : 0;
  

// ______ Window to be read into BUFFER from file in multilooked system ______
//  matrix<real4> BUFFER(bufferlines,mlpixels);         // unwrapped phase and height
//  matrix<real4> PHI(BUFFER.lines(),BUFFER.pixels());  // geocoded
//  matrix<real4> LAMBDA(BUFFER.lines(),BUFFER.pixels());       // geocoded

  window bufferwin(1, bufferlines, 1, mlpixels);        // initial


  // ______ Open output file ______
  ofstream fohei;
  openfstream(fohei,slant2hinput.fohei,generalinput.overwrit);
  bk_assert(fohei,slant2hinput.fohei,__FILE__,__LINE__); 

  ofstream fophi;
  openfstream(fophi,slant2hinput.fophi,generalinput.overwrit);
  bk_assert(fophi,slant2hinput.fophi,__FILE__,__LINE__); 

  ofstream folambda;
  openfstream(folambda,slant2hinput.folam,generalinput.overwrit);
  bk_assert(folambda,slant2hinput.folam,__FILE__,__LINE__); 

  // ______ Local variables ______
  cn M;                                 // coordinates of master in orbit
  cn Mdot;                              // velocity of master in orbit
  cn S;                                 // coordinates of slave in orbit
  cn P;                                 // coordinates of point on earth
  real8 sintheta  = 0.0;                // sine looking angle
  real8 costheta  = 0.0;                // cosine looking angle
  real8 inc_angle = 0.0;                // incidence angle by KKM
  matrix<real8> observations(3,1);      // setofeq
  matrix<real8> partials(3,3);          // partial derivatives
  matrix<real8> solution(3,1);          // solution of setofeq.

  
  // ====== Process BUFFERS ======
  for (register int32 buffer=1; buffer<=FULLBUFFERS+EXTRABUFFER; ++buffer)
    {
    // ______ Give progress ______
    PROGRESS << "SLANT2H: " 
         << 100*(buffer-1)/(FULLBUFFERS+EXTRABUFFER) << "%";
    PROGRESS.print();

    // ====== First read in data ======
    // ______ firstline of this buffer in master coordinate system ______
    real8 firstline = veryfirstline + real8((buffer-1) * bufferlines * multiL); 

    // ______ Set indices to be read from file / check if last buffer ______
    bufferwin.linelo = 1 + (buffer-1) * bufferlines;
    if (buffer == FULLBUFFERS+1)
      {
      bufferlines = RESTLINES;
      }
    bufferwin.linehi = bufferwin.linelo + bufferlines - 1;

    // ______ Read in buffer of unwrapped interferogram ______
    matrix<real4> BUFFER = unwrappedinterf.readphase(bufferwin);
    matrix<real4> PHI(BUFFER.lines(),BUFFER.pixels());          // geocoded
    matrix<real4> LAMBDA(BUFFER.lines(),BUFFER.pixels());       // geocoded


    // ====== Actually compute h for all points ======
    // ====== (better use inverse function in baseline, ie., schwabisch method.)
    input_ell ELLIPS = ellips;                  // to put P at height above ellips
    register real8 line = firstline - multiL;   // in master coordinate system
    for (uint i=0; i<BUFFER.lines(); i++)
      {
      line                += multiL;            // in master coordinate system
      real8 currentheight  = 0.0;               // this iteration
      real8 lastheight     = 0.0;               // last iteration
      register real8 pixel = firstpixel - multiP;// in master coordinate system

      // ______ Evaluate position M,S,P(ell(h)) for this pixel(l,p) ______
      // ______ only dependent on line ______
      const real8 m_tazi = master.line2ta(line);
      M    = masterorbit.getxyz(m_tazi);
      Mdot = masterorbit.getxyzdot(m_tazi);
      const real8 norm2M = M.norm2();                   // (sqr) constant per line

      // ______ compute baseline for point on ellips (h) ______
      // ______ Compute P(x,y,z), also iteratively ______
      //      real8 middlepixel = pixel+(.5*real8(BUFFER.pixels())*multiP);
      real8 tempphase = 0.0;    // to find out baseline parameters
      real8 temppixel = 0.0;  // to find out baseline parameters
      real8 Bpar      = 0.0;
      real8 Bperp     = 0.0;

      bool lineokunwrapped = false;
      for (uint t1=0; t1<BUFFER.pixels(); ++t1)
        if (BUFFER(i,t1) != NaN) lineokunwrapped = true;

      if (lineokunwrapped)
        {
        // ______ Get phase for a pixel to compute h to obtain baseline parameters ______
        int32 middle = BUFFER.pixels()/2;       // floor
        for (int32 j=0; j<middle+1; ++j)
          {
          if (BUFFER(i,middle-j) != NaN)
            {
            tempphase = BUFFER(i,middle-j);
            temppixel = firstpixel + (middle-j)*multiP;
            break; // for loop
            }
          else if (BUFFER(i,middle+j) != NaN)
            {
            tempphase = BUFFER(i,middle+j);
            temppixel = firstpixel + (middle+j)*multiP;
            break; // for loop
            }
          }

//for loop.... to find out h of point for correct baseline computation
//might be done by transformation model as well.(better?)
// ______ Compute h iteratively to get ok baseline parameters ______
        int32 j;
        for (j=0; j<=MAXITERHERE; ++j)
          {

bool DO_NEW_METHOD=false;// no time to test this, methdo seems wrong??
if (DO_NEW_METHOD==false)// old method, not using BASELINE class:
{
          real8 s_tazi;                                 // returned
          real8 s_trange;                               // returned, unused?
          ELLIPS.a = ellips.a + currentheight;          //next height
          ELLIPS.b = ellips.b + currentheight;          //next height
          lp2xyz(line,temppixel,
                 ELLIPS,master,masterorbit,  // intersect with ellips+hei
                 P,MAXITER,CRITERPOS);                  // P returned

          // ______ Compute S(x,y,z) ______ 
          xyz2t(s_tazi,s_trange,slave,slaveorbit,
                P,MAXITER,CRITERTIM);
          S = slaveorbit.getxyz(s_tazi);

          // ====== The baseline parameters, derived from the positions (x,y,z) ======
          // ====== Compute Bhor Bver (assumed contant per line) ======
          // ______ theta is angle (M,M-P) ______
          const real8 B = M.dist(S);                    // abs. value
          Bpar = P.dist(M) - P.dist(S); // sign ok


          // ______ Bperp>0 if (MP>SP) then S is to the right of slant line
          costheta = ((-(P.norm2()) + norm2M +          // cosine law
              sqr(master.pix2range(temppixel))) /
             (2*sqrt(norm2M)*master.pix2range(temppixel)));
          sintheta = sqrt(1-sqr(costheta));             // cos^2 + sin^2 = 1
 
          const cn r2 = S.min(P);                       // vector; to find sign
          const real8 costheta2 = M.in(r2) / (M.norm()*r2.norm());
          Bperp = (costheta < costheta2) ?              // sign ok
             sqrt(sqr(B)-sqr(Bpar)) :
            -sqrt(sqr(B)-sqr(Bpar));
}
else
{
// --- New method with BASELINE class and incidence angle ---
Bpar     = baseline.get_bpar(line,temppixel,currentheight);
Bperp    = baseline.get_bperp(line,temppixel,currentheight);
costheta = cos(baseline.get_theta_inc(line,temppixel,currentheight));
sintheta = sqrt(1.0-sqr(costheta));
inc_angle = baseline.get_theta_inc(line,temppixel,currentheight);//KKM
}

          // --- Update height ---
          lastheight         = currentheight;
          const real8 m_tr   = master.pix2tr(temppixel);
          const real8 tempr1 = SOL*m_tr;
          currentheight      = (-master.wavelength*tempr1*sin(inc_angle)*tempphase)/
                               (4.*PI*Bperp);//KKM added this


          // ______ Check convergence ______
          if (abs(currentheight-lastheight) < CRITERHERE)
            break; // iterate to get h
          }                             // loop iterations (j)

        // ______ Check number of iterations ______
        if (j >= MAXITERHERE)
          {
          WARNING << "slant2hambiguity: maxiter reached. "
               << "MAXITER: " << MAXITERHERE
               << "CRITER: "  << CRITERHERE << "m "
               << "last correction: " << currentheight-lastheight;
          WARNING.print();
          }
        } // check lineokunwrapped (all NaNs)

      // ______ The baseline parameters that are used foreach pixel ______
      // this is not correct, Bhor not same for each pixel? 
      // BK 09-Aug-2000
      const real8 Bhor  = Bperp*costheta + Bpar*sintheta;
      const real8 Bver  = Bperp*sintheta - Bpar*costheta;
      currentheight     = 0.0;                  // this iteration
      lastheight        = 0.0;                  // last iteration


// ====== Start loop over pixels ======
      for (uint j=0; j<BUFFER.pixels(); j++)
        {
        pixel += multiP;                        // in master coordinate system

        // ______ Check if conversion is necessary _______
        if (BUFFER(i,j) == NaN)                 // leave NaN in buffer
          {
          PHI(i,j)    = NaN;                    // not geocoded
          LAMBDA(i,j) = NaN;                    // not geocoded
          }

        else
          {
          // ______ Compute some constants (per pixel) ______
          const real8 m_trange       = master.pix2tr(pixel);
          const real8 normr1         = SOL*m_trange;
          const real8 partnumerator  = norm2M + sqr(normr1);
          const real8 denominator    = 2.0*sqrt(norm2M)*normr1;
          const real8 m_lamr1phidiv4pi = master.wavelength*normr1*BUFFER(i,j) / (4.0*PI);

          // ______ Iteratively solve for height ______
          register int32 iteration;
          for (iteration=0; iteration<=MAXITERHERE; ++iteration)
            {
            // ______ Evaluate position P(for l,p) at h ______
            // use function in orbitbk.cc cause eq1 there private f
            // BK 09-Aug-2000
            ELLIPS.a = ellips.a + currentheight;//next height
            ELLIPS.b = ellips.b + currentheight;//next height
            lp2xyz(line,pixel,
                   ELLIPS,master,masterorbit,
                   P,MAXITER,CRITERPOS);                        // P returned

            // ______ Compute h(theta,B,phi) ______
            // ______ dh=-lambdaover4pi*R1*(sintheta/(Bhor*costheta+Bver*sintheta))*(phi1-phi2)
            costheta      = (-(P.norm2()) + partnumerator) / denominator;// cosine law
            sintheta      = sqrt(1-sqr(costheta));              // cos^2 + sin^2 = 1
            lastheight    = currentheight;
            inc_angle     = baseline.get_theta_inc(line,pixel,currentheight);//KKM
            currentheight = (-m_lamr1phidiv4pi * sin(inc_angle)) /
                            (Bhor*costheta + Bver*sintheta);//KKM

            // ______ Check convergence ______
            if (abs(currentheight-lastheight) < CRITERHERE)
              break;
            }                           // loop iterations

          // ______ Check number of iterations ______
          if (iteration >= MAXITERHERE)
            {
            WARNING << "slant2hambiguity: maxiter reached. "
                 << "MAXITER: " << MAXITERHERE
                 << "CRITER: "  << CRITERHERE << "m "
                 << "last correction: " << currentheight-lastheight;
            WARNING.print();
            }


          // ______ Put computed height in buffer ______
          BUFFER(i,j) = currentheight;

          // ______ Geocode P(x,y,z) --> lat/lon/hei (h already known) ______
          real8 tmp_phi, tmp_lambda;
          ellips.xyz2ell(P, tmp_phi, tmp_lambda);// WGS84 coordinates
          PHI(i,j)    = rad2deg(tmp_phi);
          LAMBDA(i,j) = rad2deg(tmp_lambda);
          }                             // unwrapped ok
        }                               // loop over pixels
      }                                 // loop over lines
    fohei << BUFFER;                    // write output
    fophi << PHI;                       // write output
    folambda << LAMBDA;                 // write output
    }                                   // loop over buffers


// ====== Write to result files ======
  ofstream scratchlogfile("scratchlogslant2h", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"slant2h: scratchlogslant2h",__FILE__,__LINE__);
  scratchlogfile << "\n\n*******************************************************************"
                 << "\n* " << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\tambiguity"
                 << "\nData_output_file: \t\t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nData_output_file_phi: \t\t"
                 <<  slant2hinput.fophi
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nData_output_file_lam: \t\t"
                 <<  slant2hinput.folam
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nEllipsoid (name,a,b): \t\t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b 
                 << endl;
  scratchlogfile.close();


  ofstream scratchresfile("scratchresslant2h", ios::out | ios::trunc);
  bk_assert(scratchresfile,"slant2h: scratchresslant2h",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod:                               \t"
                 << "ambiguity"
                 << "\nData_output_file:                     \t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nData_output_file_phi:                 \t"
                 <<  slant2hinput.fophi
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nData_output_file_lam:                 \t"
                 <<  slant2hinput.folam
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nFirst_line (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.linelo
                 << "\nLast_line (w.r.t. original_master):   \t"
                 <<  unwrappedinterf.win.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  unwrappedinterf.win.pixlo
                 << "\nLast_pixel (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.pixhi
                 << "\nMultilookfactor_azimuth_direction:    \t"
                 <<  unwrappedinterf.multilookL
                 << "\nMultilookfactor_range_direction:      \t"
                 <<  unwrappedinterf.multilookP

                 << "\nEllipsoid (name,a,b):                 \t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b
                 << "\n*******************************************************************"
                 //<< "\n* End_slant2h:_NORMAL"
                 << "\n* End_" << processcontrol[pr_i_slant2h] << "_NORMAL"
                 << "\n*******************************************************************\n";

// ====== Tidy up ======
  scratchresfile.close();
  PROGRESS.print("finished slant2hambiguity.");
  } // END slant2hambiguity



/****************************************************************
 *    slant2h rodriguez92 method                                *
 *                                                              *
 * compute height in radar coded system (master):               *
 * use own derivations, check carefully                         *
 *                                                              *
 * Input:                                                       *
 *  -                                                           *
 * Output:                                                      *
 *  -                                                           *
 *                                                              *
 *    Bert Kampes, 30-Sep-1999                                  *
 ****************************************************************/
void slant2hrodriguez(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage      &master,
        const slcimage      &slave,
        const productinfo   &unwrappedinterf,
        const matrix<real8> &coeff_flatearth,
        orbit               &masterorbit,
        orbit               &slaveorbit,
        const BASELINE      &baseline)
  {
  TRACE_FUNCTION("slant2hrodriguez (BK 30-Sep-1999)")
  WARNING.print("this method is based on wrong approximations. (?)");

  const int32 MAXITER   = 10;                   // iterations for lp2xyz
  const real8 CRITERPOS = 1e-6;                 // stop criterium for lp2xyz
  const real8 CRITERTIM = 1e-10;                // stop criterium for lp2xyz

  const int32 degreecfe = degree(coeff_flatearth.size());

  // ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;


  // ______ Multilook factors ______
  const real8 multiL = unwrappedinterf.multilookL;
  const real8 multiP = unwrappedinterf.multilookP;
  const real8 m_pi4divlambda = (4.*PI)/master.wavelength;

  // ______ Number of lines/pixels of multilooked unwrapped interferogram ______
  const int32 mllines = int32(floor(real8(
    unwrappedinterf.win.linehi-unwrappedinterf.win.linelo+1) / multiL));
  const int32 mlpixels = int32(floor(real8(
    unwrappedinterf.win.pixhi-unwrappedinterf.win.pixlo+1)   / multiP));

  // ______ Line/pixel of first point in original master coordinates ______
  const real8 veryfirstline = real8(unwrappedinterf.win.linelo) +
                                (real8(multiL) - 1.) / 2.;
  const real8 firstpixel    = real8(unwrappedinterf.win.pixlo)  +
                                (real8(multiP) - 1.) / 2.;


  // ====== Compute number of buffers required ======
  const int32 NUMMAT = 1;       // number of large matrices to determine size of matrices
  //int32 bufferlines = generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)));
  int32 bufferlines = int32(ceil( real8( generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)) )) ));
  if (bufferlines > mllines)                            // whole image fits in BUFFER
    bufferlines = mllines;

  const int32 FULLBUFFERS = mllines / bufferlines;
  const int32 RESTLINES   = mllines % bufferlines;
  const int32 EXTRABUFFER = RESTLINES ? 1 : 0;
  

  // ______ Window to be read into BUFFER from file in multilooked system ______
  const uint dummy = 999999;                    // large to force error if not ok
  window bufferwin(1, bufferlines, 1, mlpixels);        // initial
  window offsetbuffer(1,dummy,1,dummy); // dummy not used in readfile, no offset


  // ______ Open output file ______
  ofstream fohei;
  openfstream(fohei,slant2hinput.fohei,generalinput.overwrit);
  bk_assert(fohei,slant2hinput.fohei,__FILE__,__LINE__); 



  
// ====== Process BUFFERS ======
  for (register int32 buffer=1; buffer<=FULLBUFFERS+EXTRABUFFER; ++buffer)
    {
    // ______ Give progress ______
    PROGRESS << "SLANT2H: " 
         << 100*(buffer-1)/(FULLBUFFERS+EXTRABUFFER) << "%";
    PROGRESS.print();

    // ====== First read in data ======
    // ______ firstline of this buffer in master coordinate system ______
    real8 firstline = veryfirstline + real8((buffer-1) * bufferlines * multiL); 

    // ______ Set indices to be read from file / check if last buffer ______
    bufferwin.linelo = 1 + (buffer-1) * bufferlines;
    if (buffer == FULLBUFFERS+1)
      {
      bufferlines = RESTLINES;
      //BUFFER.resize(bufferlines,mlpixels);
      }
    bufferwin.linehi = bufferwin.linelo + bufferlines - 1;

    // ______ Read in buffer of unwrapped interferogram ______
    matrix<real4> BUFFER = unwrappedinterf.readphase(bufferwin);


    // ====== Actually compute h for all points ======
    //register int32 i,j;
    register real8 line = firstline - multiL;           // in master coordinate system
    for (uint i=0; i<BUFFER.lines(); i++)
      {
      line += multiL;                                   // in master coordinate system
      register real8 pixel = firstpixel - multiP;       // in master coordinate system

      // ______ Evaluate position M,S,P(ell(h)) for this pixel(l,p) ______
      // ______ only dependent on line ______
      const real8 m_tazi = master.line2ta(line);
      cn M = masterorbit.getxyz(m_tazi);
      cn P;                                     // coordinates of point on earth

      // ______ compute baseline for point on ellips (h) ______
      // ______ Compute P(x,y,z) only to get baseline and alpha for this line ______
      const real8 middlepixel = pixel+(.5*BUFFER.pixels()*multiP);
      lp2xyz(line,middlepixel,
             ellips,master,masterorbit,
             P,MAXITER,CRITERPOS);              // P returned

      // ====== Compute which azimuth time this is for slave ======
      // ______ and compute S(x,y,z) for this time ______
      // ______ Compute S(x,y,z) could be by transf. model coregistration ______ 
      real8 s_tazi;                             // returned
      real8 s_trange;                           // returned, unused?
      xyz2t(s_tazi,s_trange,slave, slaveorbit,
            P,MAXITER,CRITERTIM);
      cn S = slaveorbit.getxyz(s_tazi);


// !! this position is not correct, depends on P(h) if orbits are not parallel.
// you should compute S by transformation model?
// error is probably very small and since H is not compute exact I will leave it for now.
// M -> P(h0) -> S -> B -> h1 = H - r costheta;
//   -> P(h1) -> S -> B -> h2
//   -> P(h2) -> S -> B -> h3 etc.


      // ====== The baseline parameters, derived from the positions (x,y,z) ======
      // ====== Compute B and alpha (contant per line) ======
      // ______ theta is angle (M,M-P) ______
      const real8 B        = M.dist(S);                         // abs. value
      const real8 Bpartemp = P.dist(M) - P.dist(S);             // sign ok

      // ______ if (MP>SP) then S is to the right of slant line, then B perp is positive.
      const real8 rho1sqr = M.norm2();                          // (sqr) constant per line
      const real8 costhetatemp = ((-(P.norm2()) + rho1sqr +             // cosine law
                   sqr(master.pix2range(middlepixel))) /
                   (2*sqrt(rho1sqr)*master.pix2range(middlepixel)));
      //sintheta = sqrt(1-sqr(costhetatemp));           // cos^2 + sin^2 = 1

      const cn r2 = S.min(P);                                   // vector; to find out sign
      const real8 Bperp = (acos(costhetatemp) > M.angle(r2)) ?  // sign ok
         sqrt(sqr(B)-sqr(Bpartemp)) :
        -sqrt(sqr(B)-sqr(Bpartemp));

      const real8 alpha  = acos(costhetatemp) - atan2(Bpartemp,Bperp);

      // ______ not used here ______
      //  const real8 Bhor  = Bperp*costhetatemp + Bpartemp*sintheta;
      //  const real8 Bver  = Bperp*sintheta - Bpartemp*costhetatemp;
      //  const real8 Bhor  = B * cos(alfa);
      //  const real8 Bver  = B * sin(alfa);


      // ______ Find out height of satellite ______
      // ______ assume radius to P is equal to this Radius. ______
      const real8 rho1  = sqrt(rho1sqr);
      real8 satphi,satlambda,satheight;
      ellips.xyz2ell(M, satphi,satlambda,satheight);
      const real8 Radius = rho1 - satheight;


      // ====== Compute per pixel: phi->Bpar, r, ->theta,p,mu,H,h ======
      for (uint j=0; j<BUFFER.pixels(); j++)
        {
        pixel += multiP;                        // in master coordinate system

        // ______ Check if conversion is necessary _______
        if (BUFFER(i,j) != NaN)                 // leave NaN in buffer
          {

          // ______ Compute some constants (per pixel) ______
          const real8 r1    = master.pix2range(pixel);
          const real8 r1sqr = sqr(r1);
          // ______ Compute reference phase to obtain Bpar ______
          // const real8 phiref= polyval(line,pixel,coeff_flatearth,degreecfe);
          const real8 phiref= polyval(normalize(line,minL,maxL),
                                      normalize(pixel,minP,maxP),
                                      coeff_flatearth,degreecfe);
          //const real8 Bpar  = -(BUFFER(i,j)+phiref)/pi4divlambda;
          // master [or slave] wavelength???
          const real8 Bpar  = -(BUFFER(i,j)+phiref)/m_pi4divlambda;

          const real8 sinthetaminalpha =
            (sqr(r1)+sqr(B)-sqr(r1-Bpar))/(2*r1*B);

// ______ two possible solutions, find out later how to do this more efficient ______
// ______ probably if Bperp is <0 then choose pi-theta
// ______ for now use theta is approx 20 degrees., (theta=1)=57 degrees
          real8 theta = asin(sinthetaminalpha)+alpha;
          if (theta<0.0 || theta>1.0)
            theta = PI-asin(sinthetaminalpha)+alpha;

          const real8 costheta = cos(theta);
          const real8 psqr  = rho1sqr+r1sqr - 2*rho1*r1*costheta;       // cosine law
          const real8 p     = sqrt(psqr);
          const real8 cosmu = ((-r1sqr+psqr+rho1sqr) / (2*p*rho1));     // cosine law
          const real8 H = rho1 - Radius*cosmu;


          // ______ Put computed height in buffer ______
          BUFFER(i,j) = H - r1*costheta;
          }                             // unwrapped ok
        }                               // loop over pixels
      }                                 // loop over lines
    fohei << BUFFER;                    // write output
    }                                   // loop over buffers


// ====== Write to result files ======
  ofstream scratchlogfile("scratchlogslant2h", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"slant2h: scratchlogslant2h",__FILE__,__LINE__);
  scratchlogfile << "\n\n*******************************************************************"
                 << "\n* " << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\trodriguez"
                 << "\nData_output_file: \t\t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nEllipsoid (name,a,b): \t\t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b 
                 << endl;
  scratchlogfile.close();


  ofstream scratchresfile("scratchresslant2h", ios::out | ios::trunc);
  bk_assert(scratchresfile,"slant2h: scratchresslant2h",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_slant2h]
                 << "\n*******************************************************************"
                 << "\nMethod:                               \t"
                 << "rodriguez"
                 << "\nData_output_file:                     \t"
                 <<  slant2hinput.fohei
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nFirst_line (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.linelo
                 << "\nLast_line (w.r.t. original_master):   \t"
                 <<  unwrappedinterf.win.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  unwrappedinterf.win.pixlo
                 << "\nLast_pixel (w.r.t. original_master):  \t"
                 <<  unwrappedinterf.win.pixhi
                 << "\nMultilookfactor_azimuth_direction:    \t"
                 <<  unwrappedinterf.multilookL
                 << "\nMultilookfactor_range_direction:      \t"
                 <<  unwrappedinterf.multilookP

                 << "\nEllipsoid (name,a,b):                 \t"
                 <<  ellips.name << " " 
                 <<  ellips.a    << " " 
                 <<  ellips.b
                 << "\n*******************************************************************"
                 //<< "\n* End_slant2h:_NORMAL"
                 << "\n* End_" << processcontrol[pr_i_slant2h] << "_NORMAL"
                 << "\n*******************************************************************\n";


// ====== Tidy up ======
  scratchresfile.close();
  PROGRESS.print("finished slant2hrodriguez.");
  } // END slant2rodriguez



/****************************************************************
 *    geocode                                                   *
 *                                                              *
 * compute phi, lambda, (height is input)                       *
 *                                                              *
 * Input:                                                       *
 *  - slant2h done, h filename in interferogram struct          *
 * Output:                                                      *
 *  - geocoded image                                            *
 *                                                              *
 * See thesis swabisch for method 3eq.                          *
 *  ellips is added known height,                               *
 *  point(x,y,z) is evaluated at that ellips.                   *
 * This is converted to ell. coord. (bowrings method)           *
 *                                                              *
 *    Bert Kampes, 02-Jun-1999                                  *
 ****************************************************************/
void geocode(
        const input_gen     &generalinput,
        const input_geocode &geocodeinput,
        const input_ell     &ellips,
        const slcimage      &master,
        const productinfo   &heightinradarsystem,
        orbit               &masterorbit)
  {
  TRACE_FUNCTION("geocode (BK 02-Jun-1999)")
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  //const real8 CRITERTIM = 1e-10;


  // ====== Evaluate for all points interferogram pos=f(l,p,height) ======
  // ______ recon with multilook, buffers
  // ______ Multilook factors ______
  const real8 multiL = heightinradarsystem.multilookL;
  const real8 multiP = heightinradarsystem.multilookP;

  // ______ Number of lines/pixels of multilooked unwrapped interferogram ______
  const int32 mllines = int32(floor(real8(
    heightinradarsystem.win.linehi-heightinradarsystem.win.linelo+1) / multiL));
  const int32 mlpixels = int32(floor(real8(
    heightinradarsystem.win.pixhi-heightinradarsystem.win.pixlo+1)   / multiP));

  // ______ Line/pixel of first point in original master coordinates ______
  const real8 veryfirstline = real8(heightinradarsystem.win.linelo) +
                                (multiL - 1.) / 2.;
  const real8 firstpixel    = real8(heightinradarsystem.win.pixlo)  +
                                (multiP - 1.) / 2.;


// ====== Compute number of buffers required ======
  const int32 NUMMAT = 3;       // number of large matrices to determine size of matrices
  //int32 bufferlines = generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)));
  int32 bufferlines = int32(ceil( real8( generalinput.memory / (NUMMAT * (mlpixels * sizeof(real4)) )) ));
  if (bufferlines > mllines)                            // whole image fits in BUFFER
    bufferlines = mllines;

  const int32 FULLBUFFERS = mllines / bufferlines;
  const int32 RESTLINES   = mllines % bufferlines;
  const int32 EXTRABUFFER = RESTLINES ? 1 : 0;
  

  // ______ Window to be read into BUFFER from file in multilooked system ______
  matrix<real4> PHI(bufferlines,mlpixels);              // 
  matrix<real4> LAMBDA(bufferlines,mlpixels);           // 
  matrix<real4> HEIGHT(bufferlines,mlpixels);           // also output
  const uint dummy = 999999;                    // large to force error if not ok
  window bufferwin(1, bufferlines, 1, mlpixels);        // initial
  window offsetbuffer(1,dummy,1,dummy); // dummy not used in readfile, no offset


  // ______ Open output files ______
  ofstream fophi;
  openfstream(fophi,geocodeinput.fophi,generalinput.overwrit);
  bk_assert(fophi,geocodeinput.fophi,__FILE__,__LINE__); 

  ofstream folam;
  openfstream(folam,geocodeinput.folam,generalinput.overwrit);
  bk_assert(folam,geocodeinput.folam,__FILE__,__LINE__); 

  
// ====== Process BUFFERS ======
  for (register int32 buffer=1; buffer<=FULLBUFFERS+EXTRABUFFER; buffer++)
    {
    // ______ Give progress ______
    PROGRESS << "GEOCODE: " 
         << 100*(buffer-1)/(FULLBUFFERS+EXTRABUFFER) << "%";
    PROGRESS.print();

    // ______ firstline of this buffer in master coordinate system ______
    real8 firstline = veryfirstline + (buffer-1) * bufferlines * multiL; 

    // ______ Set indices to be read from file / check if last buffer ______
    bufferwin.linelo = 1 + (buffer-1) * bufferlines;
    if (buffer == FULLBUFFERS+1)
      {
      bufferlines = RESTLINES;
      PHI.resize(bufferlines,mlpixels);
      LAMBDA.resize(bufferlines,mlpixels);
      HEIGHT.resize(bufferlines,mlpixels);
      }
    bufferwin.linehi = bufferwin.linelo + bufferlines - 1;
  
    // ______ Read in buffer of s2h ______
    switch (heightinradarsystem.formatflag)
      {
      case FORMATR4:
        readfile (HEIGHT, heightinradarsystem.file, mllines, bufferwin,
        offsetbuffer);
        break;
  
      default:
        PRINT_ERROR("geocode format flag on file heights (s2h output) only real4 possible.");
        throw(unhandled_case_error);
      } // switch formatflag
  
// ====== Compute xyz for all points on their height ======
    register real8 line = firstline - multiL;           // in master coordinate system
    cn pospoint;                                        // returned by lp2xyz
    input_ell ELLIPS = ellips;                          // to correct height of ellips
    real8 r;                                            // for conversion xyz2philambda
    real8 nu;                                           // for conversion xyz2philambda
    real8 sin3;                                         // for conversion xyz2philambda
    real8 cos3;                                         // for conversion xyz2philambda
    for (uint i=0; i<HEIGHT.lines(); i++)
      {
      line += multiL;                                   // in master coordinate system
      register real8 pixel = firstpixel - multiP;       // in master coordinate system
  
      for (uint j=0; j<HEIGHT.pixels(); j++)
        {
        pixel += multiP;                                // in master coordinate system
  
        // ______ Check if conversion is necessary _______
        if (HEIGHT(i,j) == NaN)
          {
          PHI(i,j)    = NaN;
          LAMBDA(i,j) = NaN;
          }

        else
          {
          // ______ Compute point on this ellips (Bowring) ______
          // ______ Adjust height of ellips to intersect zero-Doppler _______
          ELLIPS.a    = ellips.a + HEIGHT(i,j);
          ELLIPS.b    = ellips.b + HEIGHT(i,j);
          lp2xyz(line,pixel,ELLIPS,master,masterorbit,
                 pospoint,MAXITER,CRITERPOS);
          // _____ then convert cartesian xyz to geodetic WGS84 ______
          real8 tmp_phi, tmp_lambda;
          ellips.xyz2ell(pospoint, tmp_phi, tmp_lambda);// WGS84 coordinates
          PHI(i,j)    = rad2deg(tmp_phi);
          LAMBDA(i,j) = rad2deg(tmp_lambda);
          } //else
        } // loop j
      } // loop i
    fophi << PHI;
    folam << LAMBDA;
    } // loop buffers
  

  fophi.close();
  folam.close();
  
// ====== Write to result files ======
  ofstream scratchlogfile("scratchloggeocode", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"geocode: scratchloggeocode",__FILE__,__LINE__);
  scratchlogfile << "\n\n*******************************************************************"
                 << "\n* " << processcontrol[pr_i_geocoding]
                 << "\n*******************************************************************"
                 << "\nData_output_file_hei (slant2h): "
                 <<  heightinradarsystem.file
                 << "\nData_output_file_phi: \t\t"
                 <<  geocodeinput.fophi
                 << "\nData_output_file_lambda: \t"
                 <<  geocodeinput.folam
                 << endl;
  scratchlogfile.close();


  ofstream scratchresfile("scratchresgeocode", ios::out | ios::trunc);
  bk_assert(scratchresfile,"geocode: scratchresgeocode",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_geocoding]
                 << "\n*******************************************************************"
                 << "\nData_output_file_hei (slant2h): "
                 <<  heightinradarsystem.file
                 << "\nData_output_file_phi:                 \t"
                 <<  geocodeinput.fophi
                 << "\nData_output_file_lamda:               \t"
                 <<  geocodeinput.folam
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nFirst_line (w.r.t. original_master):  \t"
                 <<  heightinradarsystem.win.linelo
                 << "\nLast_line (w.r.t. original_master):   \t"
                 <<  heightinradarsystem.win.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  heightinradarsystem.win.pixlo
                 << "\nLast_pixel (w.r.t. original_master):  \t"
                 <<  heightinradarsystem.win.pixhi
                 << "\nMultilookfactor_azimuth_direction:    \t"
                 <<  heightinradarsystem.multilookL
                 << "\nMultilookfactor_range_direction:      \t"
                 <<  heightinradarsystem.multilookP
                 << "\n*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_geocoding] << "_NORMAL"
                 << "\n*******************************************************************\n";
  scratchresfile.close();
  PROGRESS.print("finished geocode.");
  } // END geocode

