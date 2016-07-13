/*
 * Copyright (c) 1999-2005 Bert Kampes
 * Copyright (c) 1999-2005 Delft University of Technology, The Netherlands
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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/bk_baseline.hh,v $ *
 * $Revision: 1.7 $ *
 * $Date: 2005/10/06 11:09:20 $ *
 * $Author: kampes $ *
 *
 * definition of BASELINE class *
 * 
 * During computations these should be used with care, they
 * may be a bit approximate.  Better is to do what Doris does,
 * i.e., evaluate positions and relate that to timing.
 * But i guess it is OK.  Use this for ultra fast geocoding...
 *
 #%// Bert Kampes, 06-Apr-2005
 ****************************************************************/

#ifndef BK_BASELINE_H
#define BK_BASELINE_H
using namespace std;

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <iostream>    // cout
#include <cstdlib>     // exit()

#include "constants.hh"   // types etc.
#include "bk_messages.hh" // info etc.
#include "slcimage.hh"    // my slc image class
#include "productinfo.hh" // my products class
#include "orbitbk.hh"     // my orbit class
#include "utilities.hh"   // eye(); normalize(); rad2deg()


// BASELINE is a new type (class) that is initialized using orbits
// and then either models the baseline parameters such as Bperp
// or can give them exact.
// Usage in the programs is something like in main do:
// BASELINE baseline;
// baseline.init(orbit1,orbit2,product?master?);
// and pass baseline to subprograms, there use baseline.get_bperp(x,y) etc.
// Bert Kampes, 31-Mar-2005
// For stability of normalmatrix, internally data are normalized
// using line/1024  pixel/1024  height/1024
// probably it is better to only do this in model_param
// but I also did it in eval_param because no time to check for correctness.
// Bert Kampes, 02-Aug-2005


// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
class BASELINE
  {
  // ______ Private data of BASELINE class ______
  private:
    bool initialized;//          
    real8 master_wavelength;//   tmp for now used for h_amb
    real8 nearrange;// range=nearrange+drange_dp*pixel
    real8 drange_dp;// range=nearrange+drange_dp*pixel
    real8 orbit_convergence;//   tmp for now constant
    real8 orbit_heading;//       tmp for now NOT USED
    real8 L_min;//               for normalization
    real8 L_max;//               for normalization
    real8 P_min;//               for normalization
    real8 P_max;//               for normalization
    real8 H_min;//    height at which parameters are modeled
    real8 H_max;//    to model phase=f(line,pix,hei) and hei=g(line,pix,phase)
    // --- B(l,p,h) = a000 + 
    //                a100*l   + a010*p   + a001*h   +
    //                a110*l*p + a101*l*h + a011*p*h +
    //                a200*l^2 + a020*p^2 + a002*h^2
    int32 N_coeffs;//  ==10 degree of 3D-poly to model B(l,p,h)
    // --- Coefficients ---
    matrix<real8> BPERP_cf;// perpendicular baseline
    matrix<real8> BPAR_cf;// parallel baseline
    matrix<real8> THETA_cf;// viewing angle
    matrix<real8> THETA_INC_cf;// incidence angle to satellite
    //real8 avg_height_ambiguity;//center height ambiguity

    // copy constructor; copy matrices etc as public.
    BASELINE(const BASELINE& A)
      {};// prevent copy constructor usage by privatizing



  /****************************************************************
    // --- B(l,p,h) = a000 + 
    //                a100*l   + a010*p   + a001*h   +
    //                a110*l*p + a101*l*h + a011*p*h +
    //                a200*l^2 + a020*p^2 + a002*h^2
   ****************************************************************/
  real8 polyval(
    const matrix<real8> &C,
    const real8 line, 
    const real8 pixel,
    const real8 height) const
    {
    #ifdef __DEBUG
    TRACE_FUNCTION("BASELINE::polyval (BK 05-Mar-2005)")
    #endif
    if (C.size()!=10) throw("error");
    return
      C(0,0) +
      C(1,0)*line       + C(2,0)*pixel       + C(3,0)*height +
      C(4,0)*line*pixel + C(5,0)*line*height + C(6,0)*pixel*height +
      C(7,0)*sqr(line)  + C(8,0)*sqr(pixel)  + C(9,0)*sqr(height);
    }// END polyval



  /****************************************************************
   * Return baselineparameters                                    *
   ****************************************************************/
  void BBparBperpTheta(real8 &B, real8 &Bpar, real8 &Bperp, real8 &theta,
             const cn Master, const cn Point, const cn Slave) const
    {
    #ifdef __DEBUG
    TRACE_FUNCTION("BBparBperpTheta (BK 05-Mar-2005)")
    #endif
    B                  = Master.dist(Slave);// baseline. abs. value (in plane M,P,S)
    const real8 range1 = Master.dist(Point);
    const real8 range2 = Slave.dist(Point);
    Bpar               = range1-range2;                 // parallel baseline, sign ok
    const cn r1        = Master.min(Point);// points from P to M
    const cn r2        = Slave.min(Point);
    theta              = Master.angle(r1);// viewing angle
    Bperp              = sqr(B)-sqr(Bpar);
    if (Bperp < 0.0) Bperp=0.0;
    else Bperp = (theta > Master.angle(r2)) ?     // perpendicular baseline, sign ok
       sqrt(Bperp) : -sqrt(Bperp);
    }// END BBparBperpTheta



  /****************************************************************
   * returns incidence angle in radians based on coordinate of    *
   * point P on ellips and point M in orbit                       *
   #%// Bert Kampes, 06-Apr-2005
   ****************************************************************/
  real8 IncidenceAngle(const cn Master, const cn Point) const
    {
    #ifdef __DEBUG
    TRACE_FUNCTION("IncidenceAngle (BK 05-Mar-2005)")
    #endif
    const cn r1     = Master.min(Point);// points from P to M
    const real8 inc = Point.angle(r1);// incidence angle (assume P on ellips)
    return inc;
    }; // END BBparBperpTheta


  // ______ Public function of bk_message class ______
  public:
    // ___Constructor/Destructor___
    BASELINE()
      {
      TRACE_FUNCTION("BASELINE() (BK 06-Mar-2005)");
      initialized       = false;//          
      N_coeffs          = 10;
      L_min             =    0.0;//      for normalization
      L_max             = 25000.0;//     for normalization
      P_min             =    0.0;//      for normalization
      P_max             = 5000.0;//      for normalization
      H_min             =    0.0;//      height at which baseline is computed.
      H_max             = 5000.0;//      height at which baseline is computed.
      master_wavelength = 0.0;
      orbit_convergence = 0.0;//   tmp for now constant
      orbit_heading     = 0.0;//       tmp for now NOT USED
      }

    // ___Constructor/Destructor___
    ~BASELINE()                 //{;}// nothing to destruct
      {
      TRACE_FUNCTION("~BASELINE() (BK 06-Mar-2005)");
      ;// nothing to destruct ?
      }// dealloc


    // ___Helper functions___
    void model_parameters(
        const slcimage &master,
        const slcimage &slave,
        orbit &masterorbit,
        orbit &slaveorbit, 
        const input_ell &ellips)
      {
      TRACE_FUNCTION("model_parameters (BK 06-Mar-2005)");
      if (masterorbit.is_initialized() == false)
        {
        DEBUG.print("Baseline cannot be computed, master orbit not initialized.");
        return;
        }
      if (slaveorbit.is_initialized()  == false)
        {
        DEBUG.print("Baseline cannot be computed, slave orbit not initialized.");
        return;
        }

      // --- Get on with it ---------------------------------------
      if (initialized==true) 
        {
        WARNING.print("baseline already initialized??? (returning)");
        return;
        }
      initialized       = true;//          
      master_wavelength = master.wavelength;//
      // ______ Model r=nearrange+drange_dp*p, p starts at 1 ______
      nearrange         = master.pix2range(1.0);
      drange_dp         = master.pix2range(2.0)-master.pix2range(1.0);
      nearrange        -= drange_dp;// (p starts at 1) 

      // ______ Set min/max for normalization ______
      L_min = master.currentwindow.linelo;// also used during polyval
      L_max = master.currentwindow.linehi;// also used during polyval
      P_min = master.currentwindow.pixlo;// also used during polyval
      P_max = master.currentwindow.pixhi;// also used during polyval
      H_min = 0.0;// also used during polyval
      H_max = 5000.0;// also used during polyval

      // ______ Loop counters ______
      register int32 cnt      = 0;// matrix index
      const int N_pointsL     = 10;// every 10km in azimuth
      const int N_pointsP     = 10;// every 10km ground range
      const int N_heights     = 4;// one more level than required for poly
      const real8 deltapixels = master.currentwindow.pixels() / N_pointsP;
      const real8 deltalines  = master.currentwindow.lines()  / N_pointsL;
      const real8 deltaheight = (H_max-H_min)                 / N_heights;
      
      // ______ Matrices for modeling Bperp (BK 21-mar-01) ______
      // --- For stability of normalmatrix, fill AMATRIX with normalized line, etc.
      matrix<real8> BPERP(N_pointsL*N_pointsP*N_heights,1);// perpendicular baseline
      matrix<real8> BPAR(N_pointsL*N_pointsP*N_heights,1);// parallel baseline
      matrix<real8> THETA(N_pointsL*N_pointsP*N_heights,1);// viewing angle
      matrix<real8> THETA_INC(N_pointsL*N_pointsP*N_heights,1);// inc. angle
      matrix<real8> AMATRIX(N_pointsL*N_pointsP*N_heights,N_coeffs);// design matrix

      // ______ Loop over heights, lines, pixels to compute baseline param. ______
      for (register int32 k=0; k<N_heights; ++k) // height levels
        {
        const real8 HEIGHT = H_min + k*deltaheight;
        input_ell ELLIPS(ellips.a+HEIGHT, ellips.b+HEIGHT);
        for (register int32 i=0; i<N_pointsL; ++i) // azimuthlines
          {
          const real8 line = master.currentwindow.linelo + i*deltalines;
          cn P;          // point, returned by lp2xyz
          real8 s_tazi;  // returned by xyz2t
          real8 s_trange;// returned by xyz2t
          const int32 MAXITER   = 10;
          const real8 CRITERPOS = 1e-6;
          const real8 CRITERTIM = 1e-10;
          // ______ Azimuth time for this line ______
          const real8 m_tazi = master.line2ta(line);
          // ______ xyz for master satellite from time ______
          const cn M         = masterorbit.getxyz(m_tazi);
          // ______ Loop over a pixels to compute baseline param. ______
          for (register int32 j=0; j<N_pointsP; ++j) // rangepixels
            {
            const real8 pixel = master.currentwindow.pixlo + j*deltapixels;
            // ______ Range time for this pixel ______
            //const real8 m_trange = master.pix2tr(pixel);
            lp2xyz(line,pixel,ELLIPS,master,masterorbit,P,MAXITER,CRITERPOS);
            // ______ Compute xyz for slave satellite from P ______
            xyz2t(s_tazi,s_trange,slave,slaveorbit,P,MAXITER,CRITERTIM);
            // ______ Slave position ______
            const cn S              = slaveorbit.getxyz(s_tazi);
            // ______ Compute angle between near parallel orbits ______
            const cn Mdot           = masterorbit.getxyzdot(m_tazi);
            const cn Sdot           = slaveorbit.getxyzdot(s_tazi);
            const real8 angleorbits = Mdot.angle(Sdot);
            TRACE << "Angle between orbits master-slave (at l,p="
                  << line << "," << pixel << ") = "
                  << rad2deg(angleorbits) << " [deg]";
            TRACE.print();
            orbit_convergence = angleorbits;// assume constant; store in member
            //const real8 heading = angle(Mdot,[1 0 0])?
            //orbit_heading = 0.0;// not yet used

            // ====== The baseline parameters, derived from the positions (x,y,z) ======
            // ______ alpha is angle counterclockwize(B, vlak met normal=rho1=rho2)
            // ______ theta is angle counterclockwize(rho1=M, r1=M-P, r2=S-P)
            real8 B,Bpar,Bperp,theta;
            BBparBperpTheta(B,Bpar,Bperp,theta,M,P,S);// return B etc.
            const real8 theta_inc = IncidenceAngle(M,P);// [rad]

            // ______ Modelling of Bperp(l,p) = a00 + a10*l + a01*p ______
            BPERP(cnt,0)     = Bperp;
            BPAR(cnt,0)      = Bpar;
            THETA(cnt,0)     = theta;
            THETA_INC(cnt,0) = theta_inc;
            // --- B(l,p,h) = a000 + 
            //                a100*l   + a010*p   + a001*h   +
            //                a110*l*p + a101*l*h + a011*p*h +
            //                a200*l^2 + a020*p^2 + a002*h^2
            AMATRIX(cnt,0)   =  1.0;
            AMATRIX(cnt,1)   = normalize(line,L_min,L_max);
            AMATRIX(cnt,2)   = normalize(pixel,P_min,P_max);
            AMATRIX(cnt,3)   = normalize(HEIGHT,H_min,H_max);
            AMATRIX(cnt,4)   = normalize(line,L_min,L_max) *normalize(pixel,P_min,P_max);
            AMATRIX(cnt,5)   = normalize(line,L_min,L_max) *normalize(HEIGHT,H_min,H_max);
            AMATRIX(cnt,6)   = normalize(pixel,P_min,P_max)*normalize(HEIGHT,H_min,H_max);
            AMATRIX(cnt,7)   = sqr(normalize(line,L_min,L_max));
            AMATRIX(cnt,8)   = sqr(normalize(pixel,P_min,P_max));
            AMATRIX(cnt,9)   = sqr(normalize(HEIGHT,H_min,H_max));
            cnt++;
            // ______ B/alpha representation of baseline ______
            const real8 alpha  = (Bpar==0 && Bperp==0) ?
              NaN :
              theta - atan2(Bpar,Bperp);            // sign ok atan2
            // ______ hor/vert representation of baseline ______
            const real8 Bh    = B * cos(alpha);                        // sign ok
            const real8 Bv    = B * sin(alpha);                        // sign ok
            // ______ Height ambiguity: [h] = -lambda/4pi * (r1sin(theta)/Bperp) * phi==2pi ______
            const real8 hambiguity = (Bperp==0) ? 
              Inf :
              -master.wavelength*(M.min(P)).norm()*sin(theta)/(2.0*Bperp);
      
            // ______ Some extra info if in debug mode ______
            // ====== Write output to screen as TRACE ======
            TRACE << "The baseline parameters for (l,p,h) = " 
                  << line << ", " << pixel << ", " << HEIGHT;
            TRACE.print();
            TRACE << "\talpha (deg), BASELINE: \t" << rad2deg(alpha) << " \t" << B;
            TRACE.print();
            TRACE << "\tBpar, Bperp:      \t" << Bpar << " \t" << Bperp;
            TRACE.print();
            TRACE << "\tBh, Bv:           \t" << Bh << " \t" << Bv;
            TRACE.print();
            TRACE << "\tHeight ambiguity: \t" << hambiguity;
            TRACE.print();
            TRACE << "\ttheta (deg):      \t" << rad2deg(theta);
            TRACE.print();
            TRACE << "\ttheta_inc (deg):  \t" << rad2deg(theta_inc);
            TRACE.print();
            TRACE.precision(10);
            TRACE.width(11);
            TRACE.rewind();
            TRACE << "\tM (x,y,z) = " << M.x << ", " << M.y << ", " << M.z;
            TRACE.print();
            TRACE << "\tS (x,y,z) = " << S.x << ", " << S.y << ", " << S.z;
            TRACE.print();
            TRACE << "\tP (x,y,z) = " << P.x << ", " << P.y << ", " << P.z;
            TRACE.print();
            TRACE.reset();// reset width/precision/pointer
            } // loop pixels
          } // loop lines
        } // loop heights

      // ====== Model the Bperp as 2d polynomial of degree 1 ======
      matrix<real8> N        = matTxmat(AMATRIX,AMATRIX);
      matrix<real8> rhsBperp = matTxmat(AMATRIX,BPERP);
      matrix<real8> rhsBpar  = matTxmat(AMATRIX,BPAR);
      matrix<real8> rhsT     = matTxmat(AMATRIX,THETA);
      matrix<real8> rhsT_INC = matTxmat(AMATRIX,THETA_INC);
      matrix<real8> Qx_hat   = N;
      choles(Qx_hat);               // Cholesky factorisation normalmatrix
      solvechol(Qx_hat,rhsBperp);   // Solution Bperp coefficients in rhsB
      solvechol(Qx_hat,rhsBpar);    // Solution Theta coefficients in rhsT
      solvechol(Qx_hat,rhsT);       // Solution Theta coefficients in rhsT
      solvechol(Qx_hat,rhsT_INC);   // Solution Theta_inc coefficients in rhsT_INC
      invertchol(Qx_hat);           // Covariance matrix of normalized unknowns

      // ______Some other stuff, normalization is ok______
      //matrix<real8> Qy_hat = AMATRIX * (matxmatT(Qx_hat,AMATRIX));
      matrix<real8> y_hatBperp  = AMATRIX * rhsBperp;
      matrix<real8> e_hatBperp  = BPERP - y_hatBperp;
      //matrix<real8> Qe_hat  = Qy - Qy_hat;
      //matrix<real8> y_hatT  = AMATRIX * rhsT;
      //matrix<real8> e_hatT  = THETA - y_hatT;

      // === Copy estimated coefficients to private members ===
      BPERP_cf     = rhsBperp;//
      BPAR_cf      = rhsBpar;//
      THETA_cf     = rhsT;//
      THETA_INC_cf = rhsT_INC;//
    
      // ______Test inverse______
      for (uint i=0; i<Qx_hat.lines(); i++)
        for (uint j=0; j<i; j++)
          Qx_hat(j,i) = Qx_hat(i,j);// repair matrix
      const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
      DEBUG << "BASELINE: max(abs(N*inv(N)-I)) = " << maxdev;
      DEBUG.print();
      if (maxdev > .01)
        {
        WARNING << "BASELINE: max. deviation N*inv(N) from unity = "
           << maxdev << ". This is larger than .01: do not use this!";
        WARNING.print();// no error, no problem
        }
      else if (maxdev > .001)
        {
        WARNING << "BASELINE: max. deviation N*inv(N) from unity = "
           << maxdev << ". This is between 0.01 and 0.001 (maybe not use it)";
        WARNING.print();
        }


      // ______ Output solution and give max error ______
      // --- B(l,p,h) = a000 + 
      //                a100*l   + a010*p   + a001*h   +
      //                a110*l*p + a101*l*h + a011*p*h +
      //                a200*l^2 + a020*p^2 + a002*h^2
      DEBUG.print("--------------------");
      DEBUG.print("Result of modeling: Bperp(l,p) = a000 + a100*l + a010*p + a001*h + ");
      DEBUG.print(" a110*l*p + a101*l*h + a011*p*h + a200*l^2 + a020*p^2 + a002*h^2");
      DEBUG.print("l,p,h in normalized coordinates [-2:2].");
      DEBUG << "Bperp_a000 = " << rhsBperp(0,0); DEBUG.print();
      DEBUG << "Bperp_a100 = " << rhsBperp(1,0); DEBUG.print();
      DEBUG << "Bperp_a010 = " << rhsBperp(2,0); DEBUG.print();
      DEBUG << "Bperp_a001 = " << rhsBperp(3,0); DEBUG.print();
      DEBUG << "Bperp_a110 = " << rhsBperp(4,0); DEBUG.print();
      DEBUG << "Bperp_a101 = " << rhsBperp(5,0); DEBUG.print();
      DEBUG << "Bperp_a011 = " << rhsBperp(6,0); DEBUG.print();
      DEBUG << "Bperp_a200 = " << rhsBperp(7,0); DEBUG.print();
      DEBUG << "Bperp_a020 = " << rhsBperp(8,0); DEBUG.print();
      DEBUG << "Bperp_a002 = " << rhsBperp(9,0); DEBUG.print();
      real8 maxerr = max(abs(e_hatBperp));
      if (maxerr > 2.00)// 
        {
        WARNING << "Max. error bperp modeling at 3D datapoints: " << maxerr << "m";
        WARNING.print();
        }
      else
        {
        INFO << "Max. error bperp modeling at 3D datapoints: " << maxerr << "m";
        INFO.print();
        }
      DEBUG.print("");
      DEBUG.print("--------------------");
      DEBUG.print("Range: r(p) = r0 + dr*p");
      DEBUG.print("l and p in unnormalized, absolute, coordinates (1:N).");
      const real8 range1    = master.pix2range(1.0);
      const real8 range5000 = master.pix2range(5000.0);
      const real8 drange    = (range5000-range1)/5000.0;
      DEBUG << "range = " << range1-drange << " + " << drange << "*p";
      DEBUG.print();
      // ====== Tidy up ======
      }; // END model_parameters()



    // ___ Return RANGE to user ___
    inline real8 get_range(const real8 pixel) const
      {
      return nearrange + drange_dp*pixel;
      };// END get_bperp()


    // === Polyval modeled quantities ===
    // --- B(l,p,h) = a000 + 
    //                a100*l   + a010*p   + a001*h   +
    //                a110*l*p + a101*l*h + a011*p*h +
    //                a200*l^2 + a020*p^2 + a002*h^2
    //
    // l,p,h coefficients take normalized input
    // Bert Kampes, 25-Aug-2005


    // ___ Return BPERP to user ___
    inline real8 get_bperp(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      return polyval(BPERP_cf, 
        normalize(line,L_min,L_max), 
        normalize(pixel,P_min,P_max), 
        normalize(height,H_min,H_max));
      }
    // ___ Return BPAR to user ___
    inline real8 get_bpar(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      return polyval(BPAR_cf, 
        normalize(line,L_min,L_max), 
        normalize(pixel,P_min,P_max), 
        normalize(height,H_min,H_max));
      }
    // ___ Return THETA to user ___
    inline real8 get_theta(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      return polyval(THETA_cf, 
        normalize(line,L_min,L_max), 
        normalize(pixel,P_min,P_max), 
        normalize(height,H_min,H_max));
      }
    // ___ Return THETA_INC to user ___
    inline real8 get_theta_inc(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      return polyval(THETA_INC_cf, 
        normalize(line,L_min,L_max), 
        normalize(pixel,P_min,P_max), 
        normalize(height,H_min,H_max));
      }



    // === Derived quantities: do not normalize these! ===
    // ___ Return B to user ___
    inline real8 get_b(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      return sqrt(sqr(get_bpar(line,pixel,height))+sqr(get_bperp(line,pixel,height)));
      };// END get_b()

    // ___ Return alpha baseline orientation to user ___
    inline real8 get_alpha(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      const real8 Bperp     = get_bperp(line,pixel,height);
      const real8 Bpar      = get_bpar(line,pixel,height);
      const real8 theta     = get_theta(line,pixel,height);
      const real8 alpha     = (Bpar==0 && Bperp==0) ?
        NaN :
        theta-atan2(Bpar,Bperp);            // sign ok atan2
      return alpha;// sign ok
      };// END get_bhor()

    // ___ Return Bh to user ___
    inline real8 get_bhor(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      const real8 B         = get_b(line,pixel,height);
      const real8 alpha     = get_alpha(line,pixel,height);
      return B*cos(alpha);// sign ok
      };// END get_bhor()

    // ___ Return Bv to user ___
    inline real8 get_bvert(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      const real8 B         = get_b(line,pixel,height);
      const real8 alpha     = get_alpha(line,pixel,height);
      return B*sin(alpha);// sign ok
      };// END get_bvert()

    // ___ Return Height ambiguity to user ___
    inline real8 get_hamb(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      //const real8 theta     =  get_theta(line,pixel,height);
      const real8 theta_inc =  get_theta_inc(line,pixel,height);
      const real8 Bperp     =  get_bperp(line,pixel,height);
      const real8 range_MP  =  get_range(pixel);// >
      const real8 h_amb     = (Bperp==0) ?
         Inf : // inf
        -master_wavelength*range_MP*sin(theta_inc)/(2.0*Bperp);// this is wrt local
        //-master_wavelength*range_MP*sin(theta)/(2.0*Bperp);// this is wrt
      return h_amb;
      };// END get_hamb()

    // ___ Return orbit convergence to user ___
    inline real8 get_orb_conv(const real8 line, const real8 pixel, const real8 height=0.0) const
      {
      // do not use l,p..
      return orbit_convergence;
      };// END get_orb_conv()


    // --- Dump overview of all ---
    void dump(const real8 line, const real8 pixel, const real8 height=0.0)
      {
      if (initialized==false)
        {
        DEBUG.print("Exiting dumpbaseline, not initialized.");
        return;
        }
        // ______ Modeled quantities ______
        const real8 Bperp     = get_bperp(line,pixel,height);
        const real8 Bpar      = get_bpar(line,pixel,height);
        const real8 theta     = get_theta(line,pixel,height);
        const real8 theta_inc = get_theta_inc(line,pixel,height);
        // ______ Derived quantities ______
        const real8 B         = get_b(line,pixel,height);
        const real8 alpha     = get_alpha(line,pixel,height);
        const real8 Bh        = get_bhor(line,pixel,height);
        const real8 Bv        = get_bvert(line,pixel,height);
        const real8 h_amb     = get_hamb(line,pixel,height);
        // ______ Height ambiguity: [h] = -lambda/4pi * (r1sin(theta)/Bperp) * phi==2pi ______
        // ====== Write output to screen as INFO ======
        INFO << "The baseline parameters for (l,p,h) = " 
             << line << ", " << pixel << ", " << height;
        INFO.print();
        INFO << "\tBpar, Bperp:      \t" << Bpar << " \t" << Bperp;
        INFO.print();
        DEBUG << "\tB, alpha (deg):  \t" << B << " \t" << rad2deg(alpha);
        DEBUG.print();
        DEBUG << "\tBh, Bv:          \t" << Bh << " \t" << Bv;
        DEBUG.print();
        INFO << "\tHeight ambiguity: \t" << h_amb;
        INFO.print();
        INFO << "\tLook angle (deg): \t" << rad2deg(theta);
        INFO.print();
        DEBUG << "\tIncidence angle (deg): \t" << rad2deg(theta_inc);
        DEBUG.print();
      }// END dump()

  };// eof BASELINE class

#endif // BK_BASELINE_H

