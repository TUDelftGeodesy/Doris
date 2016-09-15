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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/slcimage.hh,v $   *
 * $Revision: 3.15 $                                            *
 * $Date: 2005/08/24 10:03:18 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * The slcimage class contains the definition of the data and   *
 * functions for slc images.                                    *
 * Data mainly public because this used to be a struct and I    *
 * did not want to change the other code.                       *
 * functions for updating the data like prf, reading data into  *
 * a matrix, etc.                                               *
 #%// BK 14-Aug-2000                                            *
 ****************************************************************/


#ifndef SLCIMAGE_H
#define SLCIMAGE_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // typedefs


// ====== Define template functions (no member no friend) ======
// ______ (matrix class is declared way below) ______
template <class Type> class matrix;




// ====== Struct slcimage: information on master/slave ======
class slcimage                          // info per image
  {
  private:
    // copy constructor; copy matrices etc as public.
    slcimage(const slcimage& img)
      {};// prevent copy constructor usage by privatizing

  public:
    char        file[EIGHTY];           // current filename
    char        utc1[25];               // string for getorb
    int16       sensor;                 // flag for ers/asar/rsat/jers [BK 21.08.2002]
    int16       sar_processor;          // VMP or ATLANTIS or TUDELFT [BK 04.2005]
    int16       formatflag;             // flag for current file format
    cn          approxcentreoriginal;   //  wgs84? sphere x,y,z
    window      originalwindow;         // linelo etc. for normalization (master only).
    window      currentwindow;          // linelo etc. of current on file.
    real8       prf;                    // pulse repetition frequency (Hz)
    real8       abw;                    // azimuth band width (Hz)
    real8       rsr2x;                  // 2 times range sampling rate (Hz)
    real8       rbw;                    // range band width (Hz)
    real8       t_azi1;                 // sec. of day of first pixel azimuth time
    real8       t_range1;               // one way time (s) to first range pixel
    real8       wavelength;             // meters
    // real8 hamming_azi;               // weighting function designator/alpha
    // real8 hamming_range;             // weighting function
    // ______ Xtrack f_DC coefficients from leader, not used for now, ______
    // ______ better use matrix for polyval?
    real8       f_DC_a0;                // constant term Hz
    real8       f_DC_a1;                // linear term Hz/s
    real8       f_DC_a2;                // quadratic term Hz/s/s
    
    // ________________  TOPS  ONLY__________________
    real8       f_DC_t_ref_az;             // DC_reference_azimuth_time
    real8       f_DC_t_ref_rn;              // range time reference
    
// ________________  TOPS  ONLY__________________
    // FM polynomial
    real8       FM_t_ref_az;             // azimuth time reference for frequency modulation rate
    real8       FM_t_ref_rn;             // azimuth time reference for frequency modulation rate
    
    real8       FM_a0;                // constant term (Hz) for polynomial of FM rate
    real8       FM_a1;                // linear term (Hz/s) for polynomial of FM rate
    real8       FM_a2;                // quadratic term (Hz/s/s) for polynomial of FM rate
    
    real8       Ks;                   // azimuth steering rate
    real8       dt_az;                 // Azimuth time interval , for most systems this is the same as 1/PRF
    
    // ______ offset = X(l,p) - X(L,P) ______
    // ______ Where l,p are in the local slave coordinate system and ______
    // ______ where L,P are in the local master coordinate system ______
    // ______ These variables are stored in the slaveinfo variable only ______
    int32       coarseoffsetL;          // offset in line (azimuth) direction
    int32       coarseoffsetP;          // offset in pixel (range) direction
    int32       coarseorbitoffsetL;     // orbit offset in line (azimuth) direction [FvL]
    int32       coarseorbitoffsetP;     // orbit offset in pixel (range) direction [FvL]
    
    real4         slopeP;            // initial slope pixels
    real4         slopeL;            // initial slope lines
    real4         realoffsetL;            // initial offset pixels
    real4         realoffsetP;            // initial offset lines
    
    int32       ovs_az;                 // oversampling of SLC
    //real8       ovs_az;                 // oversampling of SLC, multilook test [TODO]
    int32       ovs_rg;                 // oversampling of SLC

    int32       az_timing_error;        // [FvL]
                                        // azimuth timing error in
                                        // lines with respect to
                                        // master (therefore used for slave only)
                                        // [MA] additionally used for master during master timing error 
    int32       r_timing_error;         //[FvL]
                                        // range timing error in
                                        // pixel with respect to
                                        // master (therefore used for slave only)
    bool        timingerror_flag;       //[MA] 0 if master time is not updated 
                                        //     1 if updated
    slavewindow slavemasteroffsets;     // window in master coordinates [FvL]

  
    // ______ Public functions of class ______
    // ___Constructor/Destructor___
    slcimage();// set defaults
    ~slcimage()
      {
      TRACE_FUNCTION("~slcimage() (BK 06-Mar-2005)");
      ;// nothing to destruct ?
      }// dealloc

    // ---- Helper ----
    inline void showdata() const
      {
      DEBUG << "current file:  " << file; DEBUG.print();
      DEBUG << "formatflag:    " << formatflag; DEBUG.print();
      DEBUG << "current win:   "
            << currentwindow.linelo << " " << currentwindow.linehi << " "
            << currentwindow.pixlo  << " " << currentwindow.pixhi; 
      DEBUG.print();
      DEBUG << "original win:  "
            << originalwindow.linelo << " " << originalwindow.linehi << " "
            << originalwindow.pixlo  << " " << originalwindow.pixhi; 
      DEBUG.print();
      DEBUG << "approxcentreoriginal: " << approxcentreoriginal.x
            << " " << approxcentreoriginal.y << " " << approxcentreoriginal.z;
      DEBUG.print();
      DEBUG << "utc1:            " << utc1; DEBUG.print();
      DEBUG << "sensor:          " << sensor; DEBUG.print();
      DEBUG << "sar_processor:   " << sar_processor; DEBUG.print();
      DEBUG << "prf:             " << prf; DEBUG.print();
      DEBUG << "abw:             " << abw;  DEBUG.print();
      DEBUG << "rsr:             " << rsr2x/2.0; DEBUG.print();
      DEBUG << "rbw:             " << rbw;                DEBUG.print();
      DEBUG << "t_azi1:          " << t_azi1;        DEBUG.print();
      DEBUG << "t_range1:        " << t_range1; DEBUG.print();
      DEBUG << "wavelength:      " << wavelength; DEBUG.print();
      DEBUG << "f_DC_a0:         " << f_DC_a0; DEBUG.print();
      DEBUG << "f_DC_a1:         " << f_DC_a1; DEBUG.print();
      DEBUG << "f_DC_a2:         " << f_DC_a2; DEBUG.print();
      DEBUG << "coarseoffsetL:   " << coarseoffsetL; DEBUG.print();
      DEBUG << "coarseoffsetP:   " << coarseoffsetP; DEBUG.print();
      DEBUG << "coarseorbitoffsetL:   " << coarseoffsetL; DEBUG.print(); //[FvL]
      DEBUG << "coarseorbitoffsetP:   " << coarseoffsetP; DEBUG.print(); //[FvL]
      DEBUG << "ovs_rg:          " << ovs_rg; DEBUG.print();
      DEBUG << "ovs_az:          " << ovs_az; DEBUG.print();
      DEBUG << "az_timing_error: " << az_timing_error; DEBUG.print(); //[FvL]
      DEBUG << "r_timing_error:  " << r_timing_error; DEBUG.print(); //[FvL]
      DEBUG << "slave_master_offsets: " 
            << slavemasteroffsets.l00 << " " << slavemasteroffsets.p00 << " "
            << slavemasteroffsets.l0N << " " << slavemasteroffsets.p0N << " "
            << slavemasteroffsets.lN0 << " " << slavemasteroffsets.pN0 << " "
            << slavemasteroffsets.lNN << " " << slavemasteroffsets.pNN; //[FvL]
      }

    // ______ Add a bias 1 way time range ______
    inline void add_rg_t_error(const real8 dt) {t_range1+=dt;}
    inline void add_az_t_error(const real8 dt) {t_azi1+=dt;}
    inline void add_offsetL2az_t_error(const int32 dL) {t_azi1+=real8(dL/prf);}   // [MA]
    inline void add_offsetP2rg_t_error(const int32 dP) {t_range1+=real8(dP/rsr2x);}

    // ______ Add slave master offsets [FvL] _____
    inline void add_offsetl00(const real8 dl) {slavemasteroffsets.l00+=dl;}
    inline void add_offsetp00(const real8 dp) {slavemasteroffsets.p00+=dp;}
    inline void add_offsetl0N(const real8 dl) {slavemasteroffsets.l0N+=dl;}
    inline void add_offsetp0N(const real8 dp) {slavemasteroffsets.p0N+=dp;}
    inline void add_offsetlN0(const real8 dl) {slavemasteroffsets.lN0+=dl;}
    inline void add_offsetpN0(const real8 dp) {slavemasteroffsets.pN0+=dp;}
    inline void add_offsetlNN(const real8 dl) {slavemasteroffsets.lNN+=dl;}
    inline void add_offsetpNN(const real8 dp) {slavemasteroffsets.pNN+=dp;}
 
    // ______ Read readfiles section from resultfile ______
    void fillslcimage(const char *file);
    void updateslcimage(const char *file, const char *iden);    
//    //____RaffaeleNutricato START MODIFICATION SECTION 1
//    void updateslcimageML(const char *file, const char *iden);        
//    //____RaffaeleNutricato END MODIFICATION SECTION 1

    // ______ Read matrix from file (all FORMATS) ______
    //void readdata(matrix<complr4> &Result,window win); 
    //void readdata(matrix<complr4> &Result,window win) const; 
    matrix<complr4> readdata(window win) const; 

    // ______ Convert line number to azimuth time (1 is first line) ______
    inline real8 line2ta(const real8 line) const
      {return t_azi1+(line-1.0)/prf;}

    // ______ Convert pixel number to range time (1 is first pixel) ______
    inline real8 pix2tr(const real8 pixel) const
      {return t_range1+(pixel-1.0)/rsr2x;}

    // ______ Convert pixel number to range (1 is first pixel) ______
    inline real8 pix2range(const real8 pixel) const
      {return SOL*pix2tr(pixel);}

    // ______ Convert azimuth time to line number (1 is first line) ______
    inline real8 ta2line(const real8 azitime) const
      {return 1.0+prf*(azitime-t_azi1);}

    // ______ Convert range time to pixel number (1 is first pixel) ______
    inline real8 tr2pix(const real8 rangetime) const
      {return 1.0+rsr2x*(rangetime-t_range1);}

    // ______ Convert range pixel to fDC (1 is first pixel, can be ovs) ______
    // Bert Kampes, 03-Mar-2005
    inline real8 pix2fdc(const real8 pixel) const
      {const real8 tau=(pixel-1.0)/(rsr2x/2.0);// two-way time
       return f_DC_a0+f_DC_a1*tau+f_DC_a2*sqr(tau);}


  }; // END class slcimage


#endif // SLCIMAGE_H


