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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/slcimage.cc,v $   *
 * $Revision: 3.31 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of slc image class.                           *
 * - data filling/updating.                                     * 
 * - reading in matrices.                                       *
 * - etc.                                                       *
 #%// BK 14-Aug-2000                                            *
 ****************************************************************/


#include "constants.hh"                 // typedefs
#include "slcimage.hh"                  // declarations, constants.h
#include "ioroutines.hh"                // ?
#include "utilities.hh"                 // deg2rad()
#include "ioroutines.hh"                // filesize
#include "matrixbk.hh"                  // matrix
#include "exceptions.hh"                 // my exceptions class

#include <cstdio>                       // some compilers, remove function
#include <strstream>                    // for memory stream
#include <cstring>                      // strcat?
#include <ctime>                        // time

#ifdef WIN32
// Jia defined this here, I put it in constants.hh
// Bert Kampes, 24-Aug-2005
//#define max  _MAX
//#define min  _MIN
#endif

char *strptime(const char *s, const char  *format,  struct tm *tm);


/****************************************************************
 *    slcimage::constructor                                     *
 #%// Bert Kampes, 08-Apr-2005
 ****************************************************************/
slcimage::slcimage()
  {
  TRACE_FUNCTION("slcimage() (BK 14-Aug-2000)")
  DEBUG.print("initializing slcimage with ERS defaults");
  // ______Initialize______
  strcpy(file,"unknown");
  strcpy(utc1,"unknown");
  sensor                  = SLC_ERS;            // default (vs. SLC_ASAR, JERS, RSAT)
  sar_processor           = SARPR_VMP;          // (VMP (esa paf) or ATLANTIS or TUDELFT) // TODO PGS update?
  approxcentreoriginal.x  = 0.0;
  approxcentreoriginal.y  = 0.0;
  approxcentreoriginal.z  = 0.0;
  wavelength              = 0.0565646;          // [m] default ERS2
  t_azi1                  = 0.0;                // [s] sec of day
  t_range1                = 5.5458330/2.0e3;    // [s] one way, default ERS2
  prf                     = 1679.902;           // [Hz] default ERS2
  abw                     = 1378.0;             // [Hz] default ERS2
  f_DC_a0                 = 0.0;                // [Hz] default ERS2
  f_DC_a1                 = 0.0;
  f_DC_a2                 = 0.0;
  
  f_DC_t_ref_az            = 0.0;             // DC_reference_azimuth_time
  f_DC_t_ref_rn           = 0.0;              // DC_reference range time
// ________________  TOPS  ONLY__________________
    // FM polynomial
  FM_t_ref_az            = 0.0 ;        // azimuth time reference for frequency modulation rate 
  FM_t_ref_rn            = 0.0 ;        // range time reference for frequency modulation rate 
  FM_a0                  = 0.0 ;        // constant term Hz
  FM_a1                  = 0.0 ;        // linear term Hz/s
  FM_a2                  = 0.0 ;        // quadratic term Hz/s/s
  Ks                     = 0.0 ;        // azimuth steering rate [deg/s]
  dt_az                  = 0.0 ;        //Azimuth_time_interval [s]]
  
  rsr2x                   = 18.9624680*2.0e6;   // [Hz] default ERS2
  rbw                     = 15.55e6;            // [Hz] default ERS2
  currentwindow.linelo    = 1;
  currentwindow.linehi    = 25000;              // set for normalization?
  currentwindow.pixlo     = 1;
  currentwindow.pixhi     = 5000;               // set for normalization?
  formatflag              = 0;                  // format of file on disk
  originalwindow.linelo   = 1;                  // by default
  originalwindow.linehi   = 25000;              // set for normalization?
  originalwindow.pixlo    = 1;                  // by default
  originalwindow.pixhi    = 5000;               // set for normalization?
  coarseoffsetL           = 0;                  // by default
  coarseoffsetP           = 0;                  // by default
  coarseorbitoffsetL      = 0;                  // by default, [FvL]
  coarseorbitoffsetP      = 0;                  // by default, [FvL]
  ovs_rg                  = 1;                  // by default
  ovs_az                  = 1;                  // by default
  az_timing_error         = 0;                  // by default, [FvL] unit lines
  r_timing_error          = 0;                  // by default, [FvL] unit pixels see slcimage.hh
  timingerror_flag        = false;              // by default, [MA]
  slavemasteroffsets.l00  = 0;                  // window in master coordinates [FvL]
  slavemasteroffsets.p00  = 0;
  slavemasteroffsets.l0N  = 0;
  slavemasteroffsets.p0N  = 0;
  slavemasteroffsets.lN0  = 0;
  slavemasteroffsets.pN0  = 0;
  slavemasteroffsets.lNN  = 0;
  slavemasteroffsets.pNN  = 0;
 
    } // END slcimage::slcimage()




/****************************************************************
 *    slcimage::fillslcimage                                    *
 * Fills data of slscimage object, read from resultfile after   *
 * step readfiles.                                              *
 * Expects certain strings in outpfile.                         *
 *  Break loop and warn if not found but do not exit.           *
 *  Look in first 100 lines for strings, or return sooner.      *
 *  Use SI units (Hz,s,m) for internal representation.          *
 #%// BK 14-Aug-2000                                            *
 * added INFO for fDC.                                          *
 #%// BK 06-Nov-2000                                            *
 * added sarprocessor and check for NORMAL termination          *
 #%// Bert Kampes, 29-Jul-2005                                  *
 ****************************************************************/
void slcimage::fillslcimage(const char* file)
  {
  TRACE_FUNCTION("fillslcimage (BK 14-Aug-2000)")
  INFO << "Reading parameters from: " << file;
  INFO.print();

  // ______Open file______
  ifstream resfile(file, ios::in);
  bk_assert(resfile,file,__FILE__,__LINE__);

  // ______ Lookfor identifiers (strings) ______
  char  dummyline[6*ONE27];
  real8 latitude;
  real8 longitude;
  char  word[6*ONE27]     = " ";
  int32 linecounter     = 0;
  // ______ check all fields are found _______________
  bool found_lat        = false;
  bool found_lon        = false;
  bool found_wavelength = false;
  bool found_t_range1   = false;
  bool found_t_azi1     = false;
  bool found_prf        = false;
  bool found_rsr        = false;
  bool found_nlines     = false;
  bool found_npixels    = false;
  bool found_abw        = false;
  bool found_rbw        = false;
  bool found_fdc_a0     = false;
  bool found_fdc_a1     = false;
  bool found_fdc_a2     = false;
  bool found_f_DC_t_ref  = false;     
  bool found_f_DC_t_ref_rn = false;  
  bool found_FM_t_ref    = false;        // azimuth time reference for frequency modulation rate 
  bool found_FM_t_ref_rn = false;        // range time reference for frequency modulation rate 
  bool found_FM_a0      = false;        // constant term Hz
  bool found_FM_a1      = false;        // linear term Hz/s
  bool found_FM_a2      = false;        // quadratic term Hz/s/s
  bool found_Ks         = false;        // azimuth steering rate [deg/s]
  bool found_dt_az      = false;        //Azimuth_time_interval [s]]
  
  bool found_sarp       = false;
  bool found_product    = false;

  // ______ BK test if this help for DEBUG etc. ______
  DEBUG << setiosflags(ios::fixed)
       << setprecision(8)
       << setiosflags(ios::showpoint);
  while(resfile && linecounter<1000)
    {
    ++linecounter;
    resfile.getline(dummyline,4*ONE27,'\n');              // next line
    TRACE << "dummyline: " << dummyline;
    TRACE.print();
    // ___ Check if we are at end of section already ____
    char *pchq;
    pchq = strstr(dummyline,":_NORMAL");// section terminator
    if (pchq != NULL) 
      {
      char *pchqq;
      pchqq = strstr(dummyline,"leader_datapoints");// but ignore this one...
      if (pchqq != NULL) 
        {
        DEBUG.print("Ignoring leader data points terminator.");
        }
      else
        {
        DEBUG.print("Section terminator found (string \":_NORMAL\").");
        break;// section ends
        }
      }
    // ______ Check if all parameters are read ______
 //   if (found_lat        == true &&
 //       found_lon        == true &&
 //       found_wavelength == true &&
 //       found_t_range1   == true &&
 //       found_t_azi1     == true &&
 //       found_prf        == true &&
 //       found_rsr        == true &&
 //       found_nlines     == true &&
  //      found_npixels    == true &&
 //       found_abw        == true &&
 //       found_rbw        == true &&
 //       found_fdc_a0     == true &&
 //       found_fdc_a1     == true &&
 //       found_fdc_a2     == true &&
 //       found_sarp       == true &&
 //       found_product    == true)     break;

    // ___ read parameters ___
    resfile  >> word;                                   // read word
    if      (!strcmp(word,"Scene_centre_latitude:"))
      {
      found_lat = true;
      resfile >> latitude;
      DEBUG << "[" << linecounter
           << "]: string: \"Scene_centre_latitude:\",      (derived) value: "
           << latitude;
      DEBUG.print();
      if (abs(latitude) > 90)
        WARNING.print("latitude larger than 90 degees?");
      }

    else if (!strcmp(word,"Scene_centre_longitude:"))
      {
      found_lon = true;
      resfile >> longitude;
      DEBUG << "[" << linecounter
           << "]: string: \"Scene_centre_longitude:\",     (derived) value: "
           << longitude;
      DEBUG.print();
      if (abs(longitude) > 360)
        WARNING.print("longitude larger than 360 degees?");
      }

    else if (!strcmp(word,"Radar_wavelength"))          // (m):
      {
      found_wavelength = true;
      resfile >> word >> wavelength;
      DEBUG << "[" << linecounter
           << "]: string: \"Radar_wavelength\",            (derived) value: "
           << wavelength;
      DEBUG.print();
      if (wavelength > 0.50)
        WARNING.print("wavelength is larger than 50 centimeter?");
      if (wavelength < 0.01)
        WARNING.print("wavelength is smaller than 1 centimeter?");
      }

    else if (!strcmp(word,"Range_time_to_first_pixel"))
      {
      found_t_range1 = true;
      resfile >> word >> word >> t_range1;        // (2way) (ms):
      t_range1 /= 2000.;                          // 1 way, s
      DEBUG << "[" << linecounter
           << "]: string: \"Range_time_to_first_pixel\",   (derived) value: "
           << t_range1;
      DEBUG.print();
      }

    else if (!strcmp(word,"First_pixel_azimuth_time"))
      {
      found_t_azi1 = true;
      struct tm tijdstart;
      char c12tijd0[20];
      char c12tijd0_tmp[20];// allow for .123456 ms ASAR in reading
      resfile >> word >> utc1 >> c12tijd0_tmp;  // (UTC): 26-JUL-1995 09:49:23.394
      // ______ utc1 should be including time ______
      // ______ make sure not to put in ms since what are the consequenses
      // ______ for getorb, strptime, etc. ______
      strncpy(c12tijd0,c12tijd0_tmp,12);
      c12tijd0[12]='\0';
      DEBUG << "Extracted first azimuth time string: " << c12tijd0;
      DEBUG.print();
      strcat(utc1," ");
      strcat(utc1,c12tijd0);
      // ______ Compute sec. of day for this date ______
      strptime(c12tijd0,"%T",&tijdstart);
      char c12frac0[20]="0.";
      register int32 i;
      for (i=0; i<20; i++)
        if (c12tijd0_tmp[i]=='.') break;
      int32 j = 2;
      while (c12tijd0_tmp[i] != '\0')                       // use old value of i
        {
        i++;
        c12frac0[j]=c12tijd0_tmp[i];
        j++;
        }
      t_azi1 = tijdstart.tm_sec +
                     atof(c12frac0) +
                     60 * tijdstart.tm_min +
                     3600 * tijdstart.tm_hour;

      DEBUG << "[" << linecounter
           << "]: string: \"First_pixel_azimuth_time\",    (derived) value: "
           << t_azi1;
      DEBUG.print();
      
      }

    else if (!strcmp(word,"Pulse_Repetition_Frequency"))
      {
      found_prf = true;
      resfile >> word >> word >> prf;             // (computed, Hz):
      DEBUG << "[" << linecounter
           << "]: string: \"Pulse_Repetition_Frequency\",  (derived) value: "
           << prf;
      DEBUG.print();
      }

    else if (!strcmp(word,"Range_sampling_rate"))
      {
      found_rsr = true;
      resfile >> word >> word >> rsr2x;           // (computed, MHz):
      rsr2x *= 2000000.;                          // 2 times in Hz
      DEBUG << "[" << linecounter
           << "]: string: \"Range_sampling_rate\",         (derived) value: "
           << rsr2x;
      DEBUG.print();
      }

    else if (!strcmp(word,"Number_of_lines_original:"))
      {
      found_nlines = true;
      resfile >> originalwindow.linehi;
      DEBUG << "[" << linecounter
           << "]: string: \"Number_of_lines_original:\",   (derived) value: "
           << originalwindow.linehi;
      DEBUG.print();
      if (originalwindow.linehi > 200000)
        WARNING.print("more than 200000 lines, very long-strip?"); // [MA] 100k --> 200k 
      }

    else if (!strcmp(word,"Number_of_pixels_original:"))
      {
      found_npixels = true;
      resfile >> originalwindow.pixhi;
      DEBUG << "[" << linecounter
           << "]: string: \"Number_of_pixels_original:\",  (derived) value: "
           << originalwindow.pixhi;
      DEBUG.print();
      if (originalwindow.pixhi > 100000)
        WARNING.print("more than 100000 range pixels?");
      }

    else if (!strcmp(word,"Total_azimuth_band_width"))  // (Hz):
      {
      found_abw = true;
      resfile >> word >> abw;
      DEBUG << "[" << linecounter
           << "]: string: \"Total_azimuth_band_width\",    (derived) value: "
           << abw;
      DEBUG.print();
      }
    else if (!strcmp(word,"Total_range_band_width"))    // (MHz):
      {
      found_rbw = true;
      resfile >> word >> rbw;
      rbw *= 1000000.;                            // (Hz)
      DEBUG << "[" << linecounter
           << "]: string: \"Total_range_band_width\",      (derived) value: "
           << rbw;
      DEBUG.print();
      }
    else if (!strcmp(word,"Xtrack_f_DC_constant"))      // (Hz, early edge):
      {
      found_fdc_a0 = true;
      resfile >> word >> word >> word >> f_DC_a0;// also works for ATL: 1.0E+02
      
      DEBUG << "[" << linecounter
           << "]: string: \"Xtrack_f_DC_constant\",        (derived) value: "
           << f_DC_a0;
      DEBUG.print();
      }
    
    else if (!strcmp(word,"Azimuth_steering_rate"))     // (TOPS only):
      {
      found_Ks= true;
      resfile >> word >>word;
      Ks = atof(word);
      
      DEBUG << "[" << linecounter
           << "]: string: \"Azimuth_steering_rate\",       (derived) value: "
           << word << ", "
           << Ks;
      DEBUG.print();
      }  
    
    else if (!strcmp(word,"Xtrack_f_DC_linear"))        // (Hz/s, early edge):
      {
      found_fdc_a1 = true;
      //resfile >> word >> word >> word >> f_DC_a1;
      resfile >> word >> word >> word >>word;
          
      f_DC_a1 = atof(word);
      if (abs(f_DC_a1)<1.0)
        {
        WARNING << "Strange value for f_DC_a1: " << f_DC_a1 
                << "; setting f_DC_a1=0 (expected ~1e7)";
        WARNING.print();
        WARNING.print("other definition of ATLANTIS? or Focused SLC?");
        if (!found_Ks)//do not change anything for TOPS
          f_DC_a1 = 0.0;
        }
      DEBUG << "[" << linecounter
           << "]: string: \"Xtrack_f_DC_linear\",          (derived) value: "
           << f_DC_a1;
      DEBUG.print();
      }
    else if (!strcmp(word,"Xtrack_f_DC_quadratic"))     // (Hz/s/s, early edge):
      {
      found_fdc_a2 = true;
      resfile >> word >> word >> word >> f_DC_a2;
   
      
      if (abs(f_DC_a2)<1.0)
        {
        WARNING << "strange value for f_DC_a2: " << f_DC_a2 
                << "; setting f_DC_a2=0 (expected ~1e12)";
        WARNING.print();
        WARNING.print("other definition of ATLANTIS? or Focused SLC?");
        if (!found_Ks)//do not change anything for TOPS
          f_DC_a2 = 0.0;
        }
      DEBUG << "[" << linecounter
           << "]: string: \"Xtrack_f_DC_quadratic\",       (derived) value: "
           << f_DC_a2;
      DEBUG.print();
      }
    
    else if (!strcmp(word,"DC_reference_azimuth_time:"))     // (TOPS only):
      {
      found_f_DC_t_ref = true;
      struct tm tijdstart;
      char c12tijd0[20];
      char c12tijd0_tmp[20];// allow for .123456 ms ASAR in reading
      resfile >> utc1 >> c12tijd0_tmp;  // (UTC): 26-JUL-1995 09:49:23.394
      // ______ utc1 should be including time ______
      // ______ make sure not to put in ms since what are the consequences
      // ______ for getorb, strptime, etc. ______
      strncpy(c12tijd0,c12tijd0_tmp,12);
      c12tijd0[12]='\0';
      
      
      strcat(utc1," ");
      strcat(utc1,c12tijd0);
      // ______ Compute sec. of day for this date ______
      strptime(c12tijd0,"%T",&tijdstart);
      char c12frac0[20]="0.";
      register int32 i;
      for (i=0; i<20; i++)
        if (c12tijd0_tmp[i]=='.') break;
      int32 j = 2;
      while (c12tijd0_tmp[i] != '\0')                       // use old value of i
        {
        i++;
        c12frac0[j]=c12tijd0_tmp[i];
        j++;
        }
       f_DC_t_ref_az = tijdstart.tm_sec +
                     atof(c12frac0) +
                     60 * tijdstart.tm_min +
                     3600 * tijdstart.tm_hour;
      
      DEBUG << "[" << linecounter
           << "]: string: \"DC_reference_azimuth_time:\",       (derived) value: "
           << f_DC_t_ref_az;
      DEBUG.print();
      }
     
    else if (!strcmp(word,"DC_reference_range_time:"))     // (TOPS only):
      {
      found_f_DC_t_ref_rn = true;
      resfile >>   f_DC_t_ref_rn;
      
      DEBUG << "[" << linecounter
           << "]: string: \"DC_reference_range_time:\",       (derived) value: "
           << f_DC_t_ref_rn;
      DEBUG.print();
      }
    
    else if (!strcmp(word,"FM_reference_azimuth_time:"))     // (TOPS only):
      {
      found_FM_t_ref = true;
      resfile >>  FM_t_ref_az;
      
      DEBUG << "[" << linecounter
           << "]: string: \"FM_reference_azimuth_time\",       (derived) value: "
           << FM_t_ref_az;
      DEBUG.print();
      }  
    
    else if (!strcmp(word,"FM_reference_range_time:"))     // (TOPS only):
      {
      found_FM_t_ref_rn = true;
      resfile >>  FM_t_ref_rn;
      
      DEBUG << "[" << linecounter
           << "]: string: \"FM_reference_range_time\",       (derived) value: "
           << FM_t_ref_rn;
      DEBUG.print();
      }  
    
    else if (!strcmp(word,"FM_polynomial_constant_coeff"))     // (TOPS only):
      {
      found_FM_a0= true;
      resfile >> word >> word >> word >> FM_a0;
      
      DEBUG << "[" << linecounter
           << "]: string: \"FM_polynomial_constant_coeff\",       (derived) value: "
           << FM_a0;
      DEBUG.print();
      }  
    
    
    else if (!strcmp(word,"FM_polynomial_linear_coeff"))     // (TOPS only):
      {
      found_FM_a1= true;
      resfile >> word >> word >> word >> FM_a1;
      
      DEBUG << "[" << linecounter
           << "]: string: \"FM_polynomial_linear_coeff\",       (derived) value: "
           << FM_a1;
      DEBUG.print();
      }  
    else if (!strcmp(word,"FM_polynomial_quadratic_coeff"))     // (TOPS only):
      {
      found_FM_a2= true;
      resfile >> word >> word >> word >> FM_a2;
      
      DEBUG << "[" << linecounter
           << "]: string: \"FM_polynomial_quadratic_coeff\",       (derived) value: "
           << FM_a2;
      DEBUG.print();
      }  
     
      else if (!strcmp(word,"Azimuth_time_interval"))     // (TOPS only):
      {
      found_dt_az= true;
      resfile >> word >> dt_az;
      
      DEBUG << "[" << linecounter
           << "]: string: \"Azimuth_time_interval\",       (derived) value: "
           << dt_az;
      DEBUG.print();
      }  
    
    
//TODO UPDATE and ADD the FM rate 
    // ___ SAR_PROCESSOR key added, written in readleader() ___
    // ___ if you want to force behavior, change this in master.res ___
    else if (!strcmp(word,"SAR_PROCESSOR:"))     // VMP ATLANTIS TUD
      {
      found_sarp = true;
      resfile >> word;
      DEBUG << "[" << linecounter
           << "]: string: \"SAR_PROCESSOR:\",              (derived) value: "
           << word;
      DEBUG.print();
      // ______ default processor ______
      if (!strcmp(word,"VMP"))
        {
        DEBUG.print("SAR_PROCESSOR: ESA VMP identified.");
        sar_processor = SARPR_VMP;
        }
      // ______ envisat, written by envisat_dump_header2doris.csh ______
      else if (!strcmp(word,"ASAR"))
        {
        DEBUG.print("SAR_PROCESSOR: ESA ASAR identified.");
        sar_processor = SARPR_VMP;
        }
      // ______ Data focussed by user using Atlantis ______
      else if (!strcmp(word,"ATLANTIS"))
        {
        sar_processor = SARPR_ATL;
        DEBUG.print("SAR_PROCESSOR: ATLANTIS identified.");
        }
      else if (!strcmp(word,"TUD"))
        {
        sar_processor = SARPR_TUD;
        DEBUG.print("SAR_PROCESSOR: TUD identified.");
        }
      // Modified by LG for reading ALOS Fine
      else if (!strcmp(word,"JAX"))
        {
          sar_processor = SARPR_JAX;
          DEBUG.print("SAR_PROCESSOR: JAX identified.");
        }
    else if (!strcmp(word,"TSX"))
      {
        sar_processor = SARPR_TSX;
        DEBUG.print("SAR_PROCESSOR: TSX identified.");
      }
    else if ( strstr(word,"TX") != NULL ) // maybe not supported
      {
        sar_processor = SARPR_TSX;
        DEBUG.print("SAR_PROCESSOR: TSX identified.");
      }
    else if ( strstr(word,"TDX") != NULL ) // TandemX
      {
        sar_processor = SARPR_TSX;
        DEBUG.print("SAR_PROCESSOR: TSX identified.");
      }
    else if (!strcmp(word,"RN"))
      {
        sar_processor = SARPR_RS2;
        DEBUG.print("SAR_PROCESSOR: RADARSAT-2 identified.");
      }
    else if (!strcmp(word,"CSK"))   // CSK verify
      {
        sar_processor = SARPR_CSK;
        DEBUG.print("SAR_PROCESSOR: Cosmo-skymed identified.");
      }
    else if (!strcmp(word, "GAMMA"))
      {
        sar_processor = SARPR_GAM;   // 20100916.BO.Gamma support. 
        DEBUG.print("SAR_PROCESSOR: GAMMA identified.");
      }
      else
        {
        sar_processor = SARPR_VMP;
        WARNING.print("SAR_PROCESSOR: not identified, using VMP (ESA).");
        }
      }

    // ______ 11.08.2002 BK, sensor specifier ers or asar? _____
    // ______ Product type specifier: "PRODUCT:ERS-2.SAR.SLC"
    // ______ Product type specifier: "ASAR"  (forced by me)
    // ______ Product type specifier: "JERS"  (???)
    // ______ Product type specifier: "RSAT"  (???)
    else if (!strcmp(word,"Product"))
      {
      found_product = true;
      resfile >> word >> word >> word;
      DEBUG << "[" << linecounter
           << "]: string: \"Product\",                     (derived) value: "
           << word;
      DEBUG.print();
      // Best find substring ERS, but don't know if strstr is standard.
      // We force ASAR as sensor string for this goal.
      // but now with JERS, ERS, RSAT, ASAR, we search a bit...  assume caps
      char *pch;
      pch = strstr (word,"PRODUCT");// Atlantis seems to put "PRODUCT:" or "PRODUCT"
      if (pch != NULL) 
        {
        sensor = SLC_ERS;// Atlantis processor; also ERS... assume  this
        DEBUG.print("Substring \"PRODUCT\" (radarsat?) found in Product type specifier.");
        char next_word[4*ONE27] = " ";

        // [PM]: for ALOS something goes wrong here
        resfile >> next_word;
        //      strcpy(next_word,word);

        char *next_pch;
        next_pch = strstr (next_word,"ERS");
        if (next_pch != NULL) 
          {
          DEBUG.print("Substring \"PRODUCT ERS\" found in Product type specifier.");
          sensor = SLC_ERS;// Atlantis processor; also ERS...
          // assume reader for cropping is not as for RSAT...
          }

        // for ALOS Fine [PM]
        next_pch = strstr (next_word,"ALOS"); // without this
                                              // line,pixs not dumped
        if (next_pch != NULL) 
          {
            DEBUG.print("Substring \"PRODUCT ALOS\" found in Product type specifier.");
            sensor = SLC_ALOS;
          }

        next_pch = strstr (next_word,"RSAT");
        if (next_pch != NULL) 
          {
          DEBUG.print("Substring \"PRODUCT RSAT\" found in Product type specifier.");
          sensor        = SLC_RSAT;// Atlantis processor; also does RSAT
          WARNING.print("Assuming Atlantis processor was used/format of SLC.");
          WARNING.print("If you don;t want this, change line in res file to:");
          WARNING.print("Product type specifier:                         ERS");
          }
        }
//      pch = strstr (word,"RSAT"); // MA moved below TODO remove
//      if (pch != NULL) 
//        {
//        DEBUG.print("Substring \"RSAT\" found in Product type specifier.");
//        sensor = SLC_RSAT;
//        WARNING.print("Assuming Atlantis processor was used/format of SLC.");
//        WARNING.print("If you don;t want this, change line in res file to:");
//        WARNING.print("Product type specifier:                         ERS");
//        }
//      pch = strstr (word,"RADARSAT");
//      if (pch != NULL) 
//        {
//        DEBUG.print("Substring \"RADARSAT\" found in Product type specifier.");
//        sensor        = SLC_RSAT;
//        WARNING.print("Assuming Atlantis processor was used/format of SLC.");
//        WARNING.print("If you don;t want this, change line in res file to:");
//        WARNING.print("Product type specifier:                         ERS");
//        }
      pch = strstr (word,"ES");// envisat official abbrev?
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"ES\" found in Product type specifier.");
        sensor = SLC_ASAR;
        }
      pch = strstr (word,"ASAR");
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"ASAR\" found in Product type specifier.");
        sensor = SLC_ASAR;
        }
      pch = strstr (word,"ENVISAT");
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"ENVISAT\" found in Product type specifier.");
        sensor = SLC_ASAR;
        }
      pch = strstr (word,"ERS");// first ERS since ERS is substring of JERS
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"ERS\" found in Product type specifier.");
        sensor = SLC_ERS;
        }
      pch = strstr (word,"JERS");
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"JERS\" found in Product type specifier.");
        sensor = SLC_JERS;
        }
      pch = strstr (word,"ERS_SAR");// ERS in N1 envisat format
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"ERS_SAR\" found in Product type specifier.");
        sensor = SLC_ERS_N1;
        }
      pch = strstr (word,"TSX");
      if (pch != NULL) 
        {
          DEBUG.print("Substring \"TSX\" found in Product type specifier.");
          sensor = SLC_TSX;
        }
      pch = strstr (word,"TDX");  // TandemX
      if (pch != NULL) 
        {
          DEBUG.print("Substring \"TDX\" found in Product type specifier.");
          sensor = SLC_TSX;
        }
      pch = strstr (word,"RADARSAT-2");
      if (pch != NULL) 
        {
          DEBUG.print("Substring \"Radarsat-2\" found in Product type specifier.");
          sensor = SLC_RS2;
          continue;
        }
      pch = strstr (word,"S1A");
      if (pch != NULL) 
        {
          DEBUG.print("Substring \"Sentinel-1A\" found in Product type specifier.");
          sensor = SLC_S1A;
          continue;
        }
      pch = strstr (word,"RSAT");
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"RSAT\" found in Product type specifier.");
        sensor = SLC_RSAT;
        WARNING.print("Assuming Atlantis processor was used/format of SLC.");
        WARNING.print("If you don;t want this, change line in res file to:");
        WARNING.print("Product type specifier:                         ERS");
        }
      pch = strstr (word,"RADARSAT");
      if (pch != NULL) 
        {
        DEBUG.print("Substring \"RADARSAT\" found in Product type specifier.");
        sensor        = SLC_RSAT;
        WARNING.print("Assuming Atlantis processor was used/format of SLC.");
        WARNING.print("If you don;t want this, change line in res file to:");
        WARNING.print("Product type specifier:                         ERS");
        }
      pch = strstr (word,"CSK"); // CSK verify
      if (pch != NULL) 
        {
          DEBUG.print("Substring \"CSK\" found in Product type specifier.");
          sensor = SLC_CSK;
        }
      DEBUG << "Product/sensor: " << sensor;
      DEBUG.print();
      }
    } // while file and < 1000 lines
  resfile.close();

  // --- Check if all parameters are found ---
  INFO << "Checking file: " << file;
  INFO.print();
  if (found_lat        == false) WARNING.print("found_lat==false");
  if (found_lon        == false) WARNING.print("found_lon==false");
  if (found_wavelength == false) WARNING.print("found_wavelength==false");
  if (found_t_range1   == false) WARNING.print("found_t_range1==false");
  if (found_t_azi1     == false) WARNING.print("found_t_azi1==false");
  if (found_prf        == false) WARNING.print("found_prf==false");
  if (found_rsr        == false) WARNING.print("found_rsr==false");
  if (found_nlines     == false) WARNING.print("found_nlines==false");
  if (found_npixels    == false) WARNING.print("found_npixels==false");
  if (found_abw        == false) WARNING.print("found_abw==false");
  if (found_rbw        == false) WARNING.print("found_rbw==false");
  if (found_fdc_a0     == false) WARNING.print("found_fdc_a0==false");
  if (found_fdc_a1     == false) WARNING.print("found_fdc_a1==false");
  if (found_fdc_a2     == false) WARNING.print("found_fdc_a2==false");
  if (found_sarp       == false) WARNING.print("found_sarp==false");
  if (found_product    == false) WARNING.print("found_product==false");

  // ____ ____
  DEBUG << "sensor id: " << sensor;
  DEBUG.print();
  if (sensor==SLC_ERS || sensor==SLC_ERS_N1 )
    {
    DEBUG.print("ERS 1 or 2 sensor detected.");
    if (abs(wavelength-0.0565646) > 0.01)
      WARNING.print("wavelength seems to deviates more than 1 cm from ERS2 nominal.");
    if (abs(prf - 1680) > 100.0)
      WARNING.print("prf deviates more than 100 Hz from ERS2 nominal.");
    if (abs(rsr2x - 18.9624680*2.0e6) > 100000.0)
      WARNING.print("rsr deviates more than 0.1 MHz from ERS2 nominal.");
    if (abs(abw-1378.0) > 100)
      WARNING.print("ABW deviates more than 100 Hz from ERS2 nominal?");
    if (abs(rbw-15.55e6) > 1000)
      WARNING.print("RBW deviates more than 1 kHz from ERS2 nominal?");
    }
  if (sensor==SLC_ASAR) 
    {
    INFO.print("Yeah, EnviSAT!");
    if (abs(wavelength-0.0562356) > 0.01)
      WARNING.print("wavelength seems to deviates more than 1 cm from ASAR nominal.");
    if (abs(prf-1652.0) > 100.0)
      WARNING.print("prf deviates more than 100 Hz from ASAR nominal.");
    if (abs(rsr2x - 19.207680*2.0e6) > 100000.0)
      WARNING.print("rsr deviates more than 0.1 MHz from ASAR nominal.");
    if (abs(abw-1316.0) > 100)
      WARNING.print("ABW deviates more than 100 Hz from ASAR nominal?");
    if (abs(rbw-16.0e6) > 1000)
      WARNING.print("RBW deviates more than 1 kHz from ASAR nominal?");
    }
  if (sensor==SLC_JERS)
    {
    WARNING.print("JERS detected: still in testing phase.");
    if (abs(wavelength-0.23076920) > 0.01)
      WARNING.print("wavelength seems to deviates more than 1 cm from JERS nominal.");
    if (abs(prf-1555.0) > 100.0)
      WARNING.print("prf deviates more than 100 Hz from JERS nominal.");
    if (abs(rsr2x - 17.075998*2.0e6) > 100000.0)
      WARNING.print("rsr deviates more than 0.1 MHz from JERS nominal.");
    if (abs(abw-1555.0) > 100)
      WARNING.print("ABW deviates more than 100 Hz from JERS nominal?");
    if (abs(rbw-17.075998*2.0e6) > 1000)
      WARNING.print("RBW deviates more than 1 kHz from JERS nominal?");
    }
  if (sensor==SLC_RSAT)
    {
    WARNING.print("RSAT detected: this is still in development.");
    WARNING.print(" +main difficulty with orbit in earth fixed system; doppler.");
    WARNING.print(" +but it seems to work with polyfit(max) for orbit interp.");
    }
  if (sensor==SLC_ALOS)
    {
      INFO.print("Yeah, ALOS!");
      WARNING.print("ALOS detected: this is still in development.");
      WARNING.print(" +main difficulty coregistering actual doppler focused data.");
    }
  if (sensor==SLC_TSX)
    {
      INFO.print("Yeah, TerraSAR-X!");
      if (abs(wavelength-0.031) > 0.01)
        WARNING.print("wavelength seems to deviates more than 1 cm from TSX nominal.");
    }
  if (sensor==SLC_RS2)
    {
      INFO.print("Yeah, Radarsat-2!");
    }
  if (sensor==SLC_CSK)
    {
      INFO.print("Really!, Cosmo-SkyMed!");
      if (abs(wavelength-0.031) > 0.01)
        WARNING.print("wavelength seems to deviates more than 1 cm from CSK nominal.");
    }
  // ______ Some info for range ______
  INFO << setiosflags(ios::fixed)
       << setiosflags(ios::showpoint)
       << setprecision(3);
  INFO << "Range to pixel 1 =    " << pix2range(1.0) << " m";
  INFO.print();
  INFO << "Range to pixel " << originalwindow.pixhi << " = " 
       << pix2range(originalwindow.pixhi) << " m";
  INFO.print();


  // ______ Give some info on Doppler centroid frequency fDC ______
  // ______ x=pixelnumber/RSR [s]; fDC=a0+a1x+a2x^2 [Hz] ______
  const real8 earlyedgefDC = f_DC_a0;// tau=0.0
  // ______ extreme of: y=ax^2+bx+c, dy/dx=0 at -b/2a ______
  const real8 extremepix   = (abs(f_DC_a2)<1e-9) ? 0.0: f_DC_a1/(-2.0*f_DC_a2);
  const real8 extremefDC   = f_DC_a0 +
        extremepix  * f_DC_a1 +
    sqr(extremepix) * f_DC_a2;
  const real8 faredgefDC   = f_DC_a0 + 
        real8(originalwindow.pixhi-1)/(rsr2x/2.0)  * f_DC_a1 +
    sqr(real8(originalwindow.pixhi-1)/(rsr2x/2.0)) * f_DC_a2;
  DEBUG << setiosflags(ios::scientific)
           << setprecision(8)
           << setiosflags(ios::showpos);
  DEBUG << "fDC(col) = " 
       << f_DC_a0
       << f_DC_a1 << "*(col-1)/"
       << rsr2x/2
       << f_DC_a2 << "*((col-1)/"
       << rsr2x/2 << ")^2" << ends;
  DEBUG.print();
  DEBUG << "fDC(col) = "
       << f_DC_a0
       << f_DC_a1/(rsr2x/2)    << "*(col-1)"
       << f_DC_a2/sqr(rsr2x/2) << "*(col-1)^2" << ends;
  DEBUG.print();

  // --- ---
  INFO << "Doppler centroid frequency: at pixel 1 (early edge): \t\t"
       << earlyedgefDC << " Hz";
  INFO.print();
  INFO << "Doppler centroid frequency: at pixel "
       << originalwindow.pixhi << " (far edge): \t"
       << faredgefDC << " Hz";
  INFO.print();
  INFO << "Doppler centroid frequency: at pixel " << extremepix*rsr2x/2.;
  // better if b<0 then max or somthinfg like that
  (extremefDC >= earlyedgefDC) ? 
    INFO << " (maximum):  " : INFO << " (minimum):  ";
  INFO << "\t"
       << extremefDC << " Hz";
  INFO.print();
  //const real8 minfDC = min(min(earlyedgefDC,faredgefDC),);
  //mean    fDC = ... (integral, decide max or min, etc., reshape to rectangle)

  INFO << "sensor: " << sensor;
  INFO.print();

  // ====== Tidy up ======
  //pol2xyz(approxcentreoriginal,                   // convert lat/lon to xyz
  //        deg2rad(latitude),deg2rad(longitude),0.);     //  on sphere
  input_ell wgs84ellips;
  approxcentreoriginal = wgs84ellips.ell2xyz(
    deg2rad(latitude),deg2rad(longitude),0.0);

  currentwindow.linelo = originalwindow.linelo;
  currentwindow.linehi = originalwindow.linehi;
  currentwindow.pixlo  = originalwindow.pixlo;
  currentwindow.pixhi  = originalwindow.pixhi;

#ifdef __DEBUG
  DEBUG.print("Finished fillslcimage");
  DEBUG.print("content of struct:");
  showdata();
#endif
  DEBUG.print("");
  DEBUG.reset();
  INFO.reset();
  } // END slcimage::fillslcimage



/****************************************************************
 *    updateslcimage                                            *
 *                                                              *
 *  Fills struct slcimage after step identified by iden         *
 *  and read information upto string ":_NORMAL" (i.e., only     *
 *  read appropriate section                                    *
 *                                                              *
 * input:                                                       *
 *  - resultfilename                                            *
 *  - identifier                                                *
 *  - struct slcimage                                           *
 * output:                                                      *
 *  - (updated) struct slcimage                                 *
 *      current size of image                                   *
 *                                                              *
 *    Bert Kampes, 22-Dec-1998                                  *
 * put in class.                                                *
#%// BK 14-Aug-2000                                             *
 * added check for "Multilookfactor_azimuth_direction"          *
 * added check for "Multilookfactor_range_direction"            *
 * make sure this does not happen recursively                   *
 * added reading upto ":_NORMAL" termination substring          *
#%// Bert Kampes, 29-Jul-2005                                   *
#%// Mahmut Arikan, 19-May-2009                                 *
 * minor update to file size check                              *
 ****************************************************************/
void slcimage::updateslcimage(
        const char *resultfile,
        const char *iden)
  {
  TRACE_FUNCTION("updateslcimage (BK 14-Aug-2000)")
  char                  word[4*ONE27]=" ";
  char                  dummyline[4*ONE27];

  // ______Open file______
  ifstream resfile(resultfile, ios::in);
  bk_assert(resfile,resultfile,__FILE__,__LINE__);

  bool foundiden = false;
  while(resfile)
    {
    if (strcmp(word,iden))                              // Lookfor identifier
      {
      resfile.getline(dummyline,4*ONE27,'\n');            // next line
      resfile >> word;                                  // read word
      }
    else
      {
      foundiden = true;
      DEBUG << "updateslcimage: found identifier: \"" << iden
           << "\" in file: \"" << resultfile << "\".";
      DEBUG.print();
      break;
      }
    }

  // ______Check if section has been found______
  if (!foundiden)
    {
    ERROR << "(updateslcimage). identifier: "
         << iden
         << " not found in file: " << resultfile;
    PRINT_ERROR(ERROR.get_str())
    throw(file_error);
    }

  // ======Extract info (filename and window size)======
  bool found_file   = false;
  bool found_format = false;
  bool found_line1  = false;
  bool found_lineN  = false;
  bool found_pixel1 = false;
  bool found_pixelN = false;
  bool found_mlL    = false;// optional for oversampled data: adapt prf etc. if present.
  bool found_mlP    = false;// optional for oversampled data: adapt prf etc. if present.
  float tmpmultilookL = 1.0;// default
  float tmpmultilookP = 1.0;// default
  while(resfile)                                        // i.e. rest of file
    {
    resfile.getline(dummyline,4*ONE27,'\n');              // next line
    TRACE << "dummyline: " << dummyline;
    TRACE.print();
    // ___ Check if we are at end of section already ____
    char *pch;
    pch = strstr(dummyline,":_NORMAL");// section terminator
    if (pch != NULL) 
      {
      DEBUG.print("Section terminator found (string \":_NORMAL\").");
      break;// section ends
      }
    // ___ Check if all parameters are found ______
    if (found_file   == true &&
        found_format == true &&
        found_line1  == true &&
        found_lineN  == true &&
        found_pixel1 == true &&
        found_pixelN == true &&
        found_mlL    == true &&
        found_mlP    == true)    break;

    // ___ Check for strings with parameters ___
    resfile  >> word;                                   // read word
    if (!strcmp(word,"Data_output_file:"))
      {
      found_file = true;
      resfile >> file;
      DEBUG << "string: \"Data_output_file:\",        \t\tvalue: "
           << file;
      DEBUG.print();
      }

    else if (!strcmp(word,"Data_output_format:"))
      {
      found_format = true;
      resfile >> word;
      if (!strcmp(word,"complex_short"))
        formatflag = FORMATCI2;
      else if (!strcmp(word,"complex_real4"))
        formatflag = FORMATCR4;
      else
        {
        PRINT_ERROR("wrong format specifier (impossible?)")
        throw(file_error);
        }
      DEBUG << "string: \"Data_output_format:\",      \t\tvalue: "
           << formatflag;
      DEBUG.print();
      }

    else if (!strcmp(word,"First_line"))                // (w.r.t. original):
      {
      found_line1 = true;
      resfile >> word >> word >> currentwindow.linelo;
      DEBUG << "string: \"First_line\",               \t\tvalue: "
           << currentwindow.linelo;
      DEBUG.print();
      }

    else if (!strcmp(word,"Last_line"))                 // (w.r.t. original):
      {
      found_lineN = true;
      resfile >> word >> word >> currentwindow.linehi ;
      DEBUG << "string: \"Last_line\",                \t\tvalue: "
           << currentwindow.linehi;
      DEBUG.print();
      }

    else if (!strcmp(word,"First_pixel"))               // (w.r.t. original):
      {
      found_pixel1 = true;
      resfile >> word >> word >> currentwindow.pixlo ;
      DEBUG << "string: \"First_pixel\",              \t\tvalue: "
           << currentwindow.pixlo;
      DEBUG.print();
      }

    else if (!strcmp(word,"Last_pixel"))                // (w.r.t. original):
      {
      found_pixelN = true;
      resfile >> word >> word >> currentwindow.pixhi ;
      DEBUG << "string: \"Last_pixel\",               \t\tvalue: "
           << currentwindow.pixhi;
      DEBUG.print();
      }

    // --- Only present in OVS section ---
    else if (!strcmp(word,"Multilookfactor_azimuth_direction:"))
      {
      resfile >> tmpmultilookL; 
      found_mlL = true;
      ovs_az    = int32(1.0/tmpmultilookL+0.5);// round to integer
      // ovs_az    = (tmpmultilookL > 1.0) ? (1.0/tmpmultilookL) : int32(1.0/tmpmultilookL+0.5);// round to integer
      // comment below error msg and edit slcimage.hh to define real8 ovs_az  [TODO] multilook test
      DEBUG << "String: \"Multilookfactor_azimuth_direction:\", value: "
            << tmpmultilookL;
      DEBUG.print();
      if (tmpmultilookL > 1.0)
        {
        PRINT_ERROR("Decimation of slc along azimuth direction not implemented yet.");
        throw(file_error);
        }
      }

    // --- Only present in OVS section ---
    else if (!strcmp(word,"Multilookfactor_range_direction:"))
      {
      resfile >> tmpmultilookP;
      found_mlP = true;
      ovs_rg    = int32(1.0/tmpmultilookP+0.5);// round to integer
      DEBUG << "String: \"Multilookfactor_range_direction:\", \tvalue: "
            << tmpmultilookP;
      DEBUG.print();
      if (tmpmultilookP > 1.0)
        {
        PRINT_ERROR("Decimation of slc along range direction not implemented yet.");
        throw(file_error);
        }
      }
    } // while file

  // ___ Check if all required strings are found ______
  INFO << "Checking file: " << file;
  INFO.print();
  showdata();
  if (found_file   == false) {ERROR.print("found_file==false");   throw(some_error);}
  if (found_format == false) {ERROR.print("found_format==false"); throw(some_error);}
  if (found_line1  == false) {ERROR.print("found_line1==false");  throw(some_error);}
  if (found_lineN  == false) {ERROR.print("found_lineN==false");  throw(some_error);}
  if (found_pixel1 == false) {ERROR.print("found_pixel1==false"); throw(some_error);}
  if (found_pixelN == false) {ERROR.print("found_pixelN==false"); throw(some_error);}

  // ___ Check oversampling factor of SLC image ___
  if (found_mlL    == false)  
    DEBUG.print("found_mlL==false (OK if not oversampled)");
  if (found_mlP    == false)
    DEBUG.print("found_mlP==false (OK if not oversampled)");
  // ___ Make changes to image parameters ___
  // ___ Be careful not to do this recursive, e.g., if new sections ___
  // ___ also report this parameter ___
  //
  // [BK]: seems incorrect to change window parameters;
  // [BK]: Doris always stays with original master coordinates.
  //       How is now after azimuth filter the settings, etc. ??
  //       How is function like get_fdc(pix) ?
  //       I suspect that for now the interferogram is bigger
  //       because it says ifg.win=master.win after this ?
  //       this gives incorrect info in ifg.res file, i.e.
  //       last line is 50000 not 25000, i.e., problems?
  //       but correct number of files must be known to 
  //       read data from file, etc.
  //       For now I will leave it. (the code by RN)
  //
  // -- only present in one section; then update for following steps --
  if (found_mlL==true && found_mlP==true)//prevent recursion
    {
    prf                   *= ovs_az;// adapt az sampling freq.
    rsr2x                 *= ovs_rg;// adapt rg sampling freq.
    // ___ originalwindow only used for normalization ___
    originalwindow.linelo  = ovs_az*(originalwindow.linelo-1) + 1;
    originalwindow.linehi *= ovs_az;// for normalization?  
    originalwindow.pixlo   = ovs_rg*(originalwindow.pixlo-1)  + 1;
    originalwindow.pixhi  *= ovs_rg;// for normalization?
    //// ___ currentwindow passed and changed ___
    ////currentwindow.linelo   = ovs_az*(currentwindow.linelo-1) + 1;
    ////currentwindow.linehi  *= ovs_az;
    ////currentwindow.pixlo    = ovs_rg*(currentwindow.pixlo-1)  + 1;
    ////currentwindow.pixhi   *= ovs_rg;
    // ___ give info ___
    DEBUG.print("Do not update abw and rbw: OK.");
    DEBUG.print("Do not update currentwindow, annotated is larger size");
    //DEBUG.print("Do not update originalwindow; normalization unstable for small crops");
    INFO << "Updated the following parameters: ";
    INFO.print();
    INFO << "Pulse_Repetition_Frequency (actual, Hz):  " << prf;
    INFO.print();
    INFO << "Range_sampling_rate (leaderfile, MHz):    " << rsr2x/2.0;
    INFO.print();
    INFO << "Number_of_lines_original:                 " << originalwindow.linehi;
    INFO.print();
    INFO << "Number_of_pixels_original:                " << originalwindow.pixhi;
    INFO.print();
    //INFO << "First_line (w.r.t. original_image):       " << currentwindow.linelo;
    //INFO.print();
    //INFO << "Last_line (w.r.t. original_image):        " << currentwindow.linehi;
    //INFO.print();
    //INFO << "First_pixel (w.r.t. original_image):      " << currentwindow.pixlo;
    //INFO.print();
    //INFO << "Last_pixel (w.r.t. original_image):       " << currentwindow.pixhi;
    //INFO.print();
    }


  // ______ Check filesize with format/dimensions (BK 08-Mar-2001) ______
  ifstream tmpfile(file, ios::in);   // ex: files is name.cint ... etc
  if (tmpfile) 
    {
    tmpfile.seekg(0,ios::end); // internal filesize, normal one exits if not exists
    // uint filesizetrue  = tmpfile.tellg();
    const streamoff &filesizetrue  = tmpfile.tellg();   // [MA] file > 4GB upport, this fix eliminates wrong warning
    // ___ not ok, error in bk_assert: uint filesizetrue  = filesize(file); ___
    int32 bytesperelem = 4;
    if (formatflag==FORMATCI2)          bytesperelem=4;
    if (formatflag==FORMATCR4)          bytesperelem=8;
    if (formatflag==FORMATR4)           bytesperelem=4;
    if (formatflag==FORMATI2)           bytesperelem=2;
    if (formatflag==FORMATI2_BIGENDIAN) bytesperelem=2;
    if (formatflag==FORMATR8)           bytesperelem=8;
    if (formatflag==FORMATHGT)          bytesperelem=8;
    //int32 filesizecomp = currentwindow.lines()*currentwindow.pixels()*bytesperelem;
    uint64 filesizecomp = (uint64)currentwindow.lines()  *
                                  currentwindow.pixels() *
                                  bytesperelem;
    DEBUG << "Checking format/dimensions file=" << file;
    DEBUG.print();
    if (filesizecomp != filesizetrue)
      {
      WARNING << "File: \'" << file << "\' has wrong fileformat or dimensions"   // [MA] reorganized + add file
           << ": bytesperpix="    << bytesperelem
           << ", #l="             << currentwindow.lines()
           << ", #p="             << currentwindow.pixels()
           << "; size_computed="  << filesizecomp << "B" << " v."
           << " size_ondisk="     << filesizetrue << "B";
      WARNING.print();
      }
    else
      {
      DEBUG.print("SLC Fileformat and dimensions are checked and ok.");
      }
    tmpfile.close();
    } // if file stream is OK
  else
    {
    WARNING << "File: " << file
         << " does not seem to exist (may not be a problem).";
    WARNING.print();
    }


  // ______Tidy up______
  DEBUG.print("");
  resfile.close();
#ifdef __DEBUG
  DEBUG.print("Finished updateslcimage");
  DEBUG.print("content of struct:");
  showdata();
#endif
  } // END updateslcimage





/****************************************************************
 * slcimage::readdata                                          *
 *  read data from file in a complr4 matrix                     *
 *  file may be complex short or complex real4 format           *
 * stored in row major order, pixel interleaved                 *
 %// BK 14-Aug-2000                                             *
 ****************************************************************/
matrix<complr4> slcimage::readdata(
        window win) const               // window to be read
  {
  // --- Log debug info ---
  TRACE_FUNCTION("readdata (BK 14-Aug-2000)")
  DEBUG << "Reading file:     " << file;
  DEBUG.print();
  DEBUG << "Formatflag:       " << formatflag;
  DEBUG.print();
  DEBUG << "Currentwindow:    ";  currentwindow.disp();// appends to DEBUG and prints
  DEBUG << "Window from file: ";  win.disp();// appends to DEBUG and prints
  //
  // --- Read data in complex real4 matrix ---
  matrix<complr4> Result(win.lines(),win.pixels());
  switch (formatflag)
    {
    case FORMATCI2:
      {
      DEBUG.print("Format complex short integer (pixel interleaved).");
      fileci2tomatcr4(Result,file,currentwindow.lines(),win,currentwindow);
      break;
      }
    case FORMATCR4:
      {
      DEBUG.print("Format complex float (pixel interleaved).");
      readfile(Result,file,currentwindow.lines(),win,currentwindow);
      break;
      }
    default:
      PRINT_ERROR("readdata::not correct format on file.")
      throw(file_error);
    }
  return Result;
  } // END readdata


