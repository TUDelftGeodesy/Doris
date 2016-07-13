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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readinput.hh,v $
 * $Revision: 3.24 $
 * $Date: 2006/05/18 10:03:18 $
 * $Author: kampes $
 *
 * Some constants for switching processing steps
 * Structs where input parameters are stored
 * Definition of readinput function.
 ****************************************************************/


#ifndef READINPUT_H
#define READINPUT_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                   // typedefs
// ______ only declare here ______
#include "matrixbk.hh"                   // my matrix class



//____RaffaeleNutricato START MODIFICATION SECTION 1
// ______index for: input_gen: process[NUMPROCESSES] (arbitrary index)______
// const int16 // [MA] use enum instead
// With Doris v4, we want to keep SLC processing controls the same for both master and slave result files. [MA]
enum {
                pr_m_readfiles    = 0,  // 0   flag for reading leader etc master
                pr_m_crop         ,     // 1   write slc data to raster format
                pr_m_oversample   ,     // 2   oversample master image //____RaffaeleNutricato added this line
                pr_m_porbits      ,     // 3   calling getorb master
                pr_m_morbits      ,     // 4   modify orbits by correction parameters [HB]
                pr_m_simamp       ,     // 5   simulate amplitude for master timing  [MA]
                pr_m_mtiming      ,     // 6   correlate simamp and master for master timing [MA]
                pr_m_filtazi      ,     // 7   azimuth filtering master
                pr_m_filtrange    ,     // 8   range filtering master
                pr_m_resample     ,     // 9   resample fake entry [MA] 2009, see also ioroutines.cc
                pr_m_EXTRA        ,     // 10  for future use

                pr_s_readfiles    ,     // 11   flag for reading leader etc slave
                pr_s_crop         ,     // 12  writing to raw slave
                pr_s_oversample   ,     // 13  oversample slave image   //____RaffaeleNutricato added this line
                pr_s_porbits      ,     // 14  calling getorb slave
                pr_s_morbits      ,     // 15  modify orbits by correction parameters [HB]
                pr_s_simamp       ,     // 16  simulate amplitude, fake entry no processing is defined  [MA]
                pr_s_mtiming      ,     // 17  master timing, fake entry  [MA]
                pr_s_filtazi      ,     // 18  azimuth filtering slave
                pr_s_filtrange    ,     // 19  range filtering slave
                pr_s_resample     ,     // 20  resample slave, actually at slave.
                pr_s_EXTRA        ,     // 21  for future use

                pr_i_coarse       ,     // 22  coarse coregistration: orbits
                pr_i_coarse2      ,     // 23  coarse coregistration: correlation
                pr_i_fine         ,     // 24  fine coregistration
                pr_i_timing       ,     // 25  relative timing error [FvL]
                pr_i_demassist    ,     // 26  dem assist coregistration [FvL]
                pr_i_coregpm      ,     // 27  compute parameters coregistration
                pr_i_interfero    ,     // 28  computation of interferogram
                pr_i_coherence    ,     // 29  computation of coherence image
                pr_i_comprefpha   ,     // 30  ref.phase (flat earth)
                pr_i_subtrrefpha  ,     // 31  interferogram - refpha
                pr_i_comprefdem   ,     // 32  ref.phase (flat earth)
                pr_i_subtrrefdem  ,     // 33  interferogram - refpha
                pr_i_filtphase    ,     // 34  filtering interferogram
                pr_i_unwrap       ,     // 35  unwrapping interferogram
                pr_i_estorbits    ,     // 36  estimate orbit error [HB]
                pr_i_slant2h      ,     // 37  slant to height conversion
                pr_i_geocoding    ,     // 38  geocoding
                pr_i_dinsar       ,     // 39  3 pass differential
                pr_i_EXTRA2       ,     // 40  for future use
                pr_last_one             // 41  not used; to indicate last one (==NUMPROCESSES-1)
};
const int16     NUMPROCESSES      = pr_last_one+1;        // see above...
//____RaffaeleNutricato END MODIFICATION SECTION 1

// ====== global processcontrol array ======
// BK should be one word (no spaces)
// usage: resfile << "\n*_Start_" << processcontrol[pr_i_interfero]
// usage: resfile << "\n* End_" << processcontrol[pr_i_coarse] << "_NORMAL"

//____RaffaeleNutricato START MODIFICATION SECTION 2
const char processcontrol[NUMPROCESSES][ONE27] = {
   "readfiles:",            // 0  pr_m_readfiles
   "crop:",                 // 1  pr_m_crop
   "oversample:",           // 2  pr_m_oversample //____RaffaeleNutricato added this line
   "precise_orbits:",       // 3  pr_m_porbits
   "modify_orbits:",        // 4  pr_m_morbits [HB]
   "sim_amplitude:",        // 5  pr_m_simamp
   "master_timing:",        // 6  pr_m_mtiming + master timing error // was simamp_corr
   "filt_azi:",             // 7  pr_m_filtazi
   "filt_range:",           // 8  pr_m_filtrange
   "resample:",             // 9  pr_m_resample   [MA] 200903, fake entry
   "NOT_USED:",             // 10 pr_m_EXTRA

   "readfiles:",            // 11 pr_s_readfiles
   "crop:",                 // 12 pr_s_crop
   "oversample:",           // 13 pr_s_oversample //____RaffaeleNutricato added this line
   "precise_orbits:",       // 14 pr_s_porbits
   "modify_orbits:",        // 15 pr_s_morbits [HB]
   "sim_amplitude:",        // 16 pr_s_simamp   [MA] 2009, fake entry
   "master_timing:",        // 17 pr_s_mtiming  [MA] 2009, fake entry
   "filt_azi:",             // 18 pr_s_filtazi: must be same as m_
   "filt_range:",           // 19 pr_s_filtrange
   "resample:",             // 20 pr_s_resample !slave
   "NOT_USED:",             // 21 pr_s_EXTRA

   "coarse_orbits:",        // 22 pr_i_coarse
   "coarse_correl:",        // 23 pr_i_coarse2
   "fine_coreg:",           // 24 pr_i_fine
   "timing_error:",         // 25 pr_i_timing [FvL]
   "dem_assist:",           // 26 pr_i_demassist [FvL]
   "comp_coregpm:",         // 27 pr_i_coregpm
   "interfero:",            // 28 pr_i_interfero
   "coherence:",            // 29 pr_i_coherence
   "comp_refphase:",        // 30 pr_i_comprefpha
   "subtr_refphase:",       // 31 pr_i_subtrrefpha
   "comp_refdem:",          // 32 pr_i_comprefdem
   "subtr_refdem:",         // 33 pr_i_subtrrefdem
   "filtphase:",            // 34 pr_i_filtphase
   "unwrap:",               // 35 pr_i_unwrap
   "est_orbits:",           // 36 pr_i_estorbits [HB
   "slant2h:",              // 37 pr_i_slant2h
   "geocoding:",            // 38 pr_i_geocoding
   "dinsar:",               // 39 pr_i_EXTRA
   "NOT_USED2:",            // 40 pr_i_EXTRA2
   "ERROR update this!"};   // 41 pr_last_one

//  strcpy(processcontrol[pr_m_readfiles],"readfiles:");
//  strcpy(processcontrol[pr_m_crop]     , "crop:");
//____RaffaeleNutricato END MODIFICATION SECTION 2




// ====== Constants, variables ======
// ======Method selectors======
//const int16   readfiles_ers   = 151;          // method for readfiles
//const int16   readfiles_asar  = 152;          // method for readfiles
const int16     cc_magfft       = 21;           // method for coarse corr. coreg
const int16     cc_magspace     = 22;           // method for coarse corr. coreg

const int16     fc_cmplxfft     = 31;           // method for fine coreg
const int16     fc_cmplxspace   = 32;           // method for fine coreg
const int16     fc_magfft       = 33;           // method for fine coreg
const int16     fc_magspace     = 34;           // method for fine coreg
const int16     fc_oversample   = 35;           // method oversample signal,not corr.
const int16     fc_coherence    = 36;           // method oversample signal,CCC.
const int16     fc_intensity    = 37;            //same as fc_oversample but with intensity instead of amplitude

const int16     fe_porbits      = 41;           // method for flat earth correction
const int16     fe_method2      = 42;           // method for flat earth correction

const int16     rs_rect         = 102;          // nearest neighbor, 1 or 2 points
const int16     rs_tri          = 202;          // piecewize linear, 2 point
const int16     rs_cc4p         = 304;          // cubic convolution, 4 point
const int16     rs_cc6p         = 306;          // cubic convolution, 6 point
const int16     rs_ts6p         = 406;          // truncated sinc, 6 point
const int16     rs_ts8p         = 408;          // truncated sinc, 8 point
const int16     rs_ts16p        = 416;          // truncated sinc, 16 point
const int16     rs_knab4p       = 504;          // knab window, 4 point
const int16     rs_knab6p       = 506;          // knab window, 6 point
const int16     rs_knab8p       = 508;          // knab window, 8 point
const int16     rs_knab10p      = 510;          // knab window, 10 point
const int16     rs_knab16p      = 516;          // knab window, 12 point
const int16     rs_rc6p         = 606;          // raised cosine, 6 point
const int16     rs_rc12p        = 612;          // raised cosine, 12 point

const int16     int_oldmethod   = 91;           // method no overs. interf. gen.
const int16     int_oversample  = 92;           // method oversample for int. computation
const int16     coh_oldmethod   = 101;          // coherence computation up to refphase
const int16     coh_newmethod   = 102;          // coherence computation including refdem

const int16     fp_goldstein    = 81;           // method goldstein
const int16     fp_spatialconv  = 82;           // method spatial conv. with kernel
const int16     fp_spectral     = 83;           // method spectral kernel
const int16     fp_modgoldstein = 84;           // method modified goldstein

const int16     uw_method1      = 51;           // TREEF  (delft only)
const int16     uw_method2      = 52;           // SNAPHU Curtis Cheng
const int16     uw_method3      = 53;           // method for unwrapping

const int16     s2h_schwabisch  = 61;           // schwabisch95 thesis
const int16     s2h_rodriguez   = 62;           // see paper rodriguez92
const int16     s2h_ambiguity   = 63;           // exact method

const int16     geo_nointerp    = 70;           // no interpolation in geocding
const int16     geo_regular     = 71;           // ... ? not implemented

// const int16  crd_nearest     = 81;           // nearest neighbor, not used anymore [MA]
// const int16  crd_trilinear   = 82;           // trilinear

const int16     srp_polynomial  = 101;          // eval poly from comp_refpha
const int16     srp_exact       = 102;          // compute ref_pha on he fly

const int16     rf_adaptive     = 301;          // range filter adaptive method
const int16     rf_porbits      = 302;          // precise orbits

const int16     eo_noweighting  = 701;          // estorbit no weithting [HB]
const int16     eo_coh          = 702;          // estorbit weighting by coherence [HB]
const int16     eo_coh2         = 703;          // estorbit weighting by coherence^2 [HB]
const int16     eo_pdf          = 704;          // estorbit weighting by pdf(coherence) [HB]
const int16     eo_lsq          = 751;          // estorbit least squares method [HB]
const int16     eo_gridsearch   = 752;          // estorbit gridsearch method [HB]

// ======Define structs for storing/passing input options======
struct input_gen                        // general input
  {
  char          logfile[4*ONE27];
  char          m_resfile[4*ONE27];
  char          s_resfile[4*ONE27];
  char          i_resfile[4*ONE27];
  //uint          memory;                 // available mem. in Bytes. For 32-bit platform above 4200MB (4.2e9) is a problem even if uint [MA]
  real8          memory;                 // available mem. in Bytes
  bool          process[NUMPROCESSES];  // if .[i] != 0 => process step_(i+1)
  bool          interactive;            // if true, pause
  bool          overwrit;               // 0: don't overwrite existing data output files
  int16         orb_interp;             // method for orbit interpolation
  int16         orb_prm;                // orbit parameters: POSITION or POSITION and VELOCITY
  int32         dumpbaselineL;          // #lines to dump baseline param.
  int32         dumpbaselineP;          // #lines to dump baseline param.
  int32         preview;                // generate sunraster preview file
                                        // 0: no; 1: sh files; 2: sh sh_files.
  real4         terrain_height;         // mean terrain height, or of a point.
  cn            tiepoint;               // lat/lon/hei for e.g., integration const.
  };


// ______These structs are input, filled by readinput______
// ______ ______
struct input_readfiles                  // arguments of readinfo --> m_readfilesinput, --> s_readfilesinput
  {
  int16         sensor_id;              // method selector ers/asar/rsat/jers/alos
  int16         sar_processor;          // atlantis, vmp, set in routines.
  int16         fileid;
  char          volfile[4*ONE27];
  char          leaderfile[4*ONE27];
  char          nullfile[4*ONE27];
  char          datfile[4*ONE27];
  real8         rg_timing_error;        // [BK 27-May-2004]
  real8         az_timing_error;        // [BK 27-May-2004]
  };


// ______ ______
struct input_pr_orbits                    // precise orbits
  {
  char          m_orbdir[4*ONE27];
  char          s_orbdir[4*ONE27];
  int32         timeinterval;               // time in sec.
  int32         timebefore;                   // sec before first line.
  real8         dumpmasterorbit;          // dtime in sec.
  real8         dumpslaveorbit;           // dtime in sec.
  };

// ______ ______ [HB]
struct input_morbits
  {
    matrix<real8> coeff;                  // coefficients of correction polynomials
    char          reforbitfile[4*ONE27];  // result file of supermaster
  };


// ______ ______
struct input_crop                       // arguments of m/s_crop
  {
  int16         fileid;
  int16         sar_processor;                // atlantis, vmp, set in routines.
  char          idcrop[EIGHTY];
  char          filein1[4*ONE27];
  char          fileout1[4*ONE27];
  window        dbow;                                 // cut out of original
  window        dbow_geo;                           // lat_0*1e6, lon_0*1e6, height, width
  };

// ______ ______ [BO, FvL, MA]

struct input_simamp                     // arguments for simulation of amplitude from DEM [MA]
  {
  char          firefdem[4*ONE27];      // input filename reference dem
  int16         iformatflag;            // input format [signed short]
  uint          demrows;                // number of
  uint          demcols;                // number of
  real8         demdeltalat;            // radians
  real8         demdeltalon;            // radians
  real8         demlatleftupper;        // radians
  real8         demlonleftupper;        // radians
  real8         demnodata;              // identifier/flag
  char          fodem[4*ONE27];         // flag+name output of cropped dem
  char          fodemlp[4*ONE27];       // output filename DEM in radarcoord. // [MA] the same size as master image
//  char                fodemi[4*ONE27];                // flag+name output of interpolated dem
  char          fosimamp[4*ONE27];      // flag+name output of simulated amp [MA] TODO define common struct for comprefdem demassist sim-amp
  char          fothetalp[4*ONE27];      // output filename for theta incidence angle in radar coordinates. // MA
  };
// ______ ______


struct input_mtiming                    // arguments of correlation for master timing error
  {
  char          ifpositions[4*ONE27];   // input file name for positions
  int16         method;                 // method selector, [MA] rm if not nec.
  uint          Nwin;                   // #windows
  uint          MasksizeL;              // size of correlation window
  uint          MasksizeP;              // size of correlation window
  uint          AccL;                   // #lines to be searched in 1 direction
  uint          AccP;                   // #pixels to be searched in 1 direction
  int32         initoffsetL;            // initial offset lines
  int32         initoffsetP;            // initial offset pixels
  };

// ______ ______ [MA]

//____RaffaeleNutricato START MODIFICATION SECTION 3
struct input_oversample          // arguments of m/s_oversample
  {
  char      fileoutovs[4*ONE27];  // Name of the oversampled image
  int32     OsrRange;            // Oversampling ratio in range.
  int32     OsrAzimuth;          // Oversampling ratio in azimuth.
  int32     FilterSize;          // Length of the interpolation kernel.
  int32     oformatflag;         // Output format [cr4] ci16, I suggest [cr4].
  };
//____RaffaeleNutricato END MODIFICATION SECTION 3


// ______ azimuth filtering ______
struct input_filtazi                    // arguments for azimuth filter
  {
  int16         method;                 // method selector
  int32         fftlength;              // length per buffer
  int32         overlap;                // 0.5overlap each buffer
  real8         hammingalpha;           // alpha for hamming, 1 is no
  char          foname[4*ONE27];                // output filename passed to routine
  char          fomaster[4*ONE27];      // output filename
  char          foslave[4*ONE27];       // output filename
  int16         oformatflag;            // output format [cr4] ci16
  };



// ______ ______
struct input_coarsecorr                 // arguments for correlation
  {
//  char                idcoarsecorr[EIGHTY];
  char          ifpositions[4*ONE27];   // input file name for positions
  int16         method;                 // method selector
  uint          Nwin;                   // #windows
  uint          MasksizeL;              // size of correlation window
  uint          MasksizeP;              // size of correlation window
  uint          AccL;                   // #lines to be searched in 1 direction
  uint          AccP;                   // #pixels to be searched in 1 direction
  int32         initoffsetL;            // initial offset lines
  int32         initoffsetP;            // initial offset pixels
  };

  
//  struct input_deramp                       // arguments of m/s_crop
//  {
//   
 // char          filein1[4*ONE27];  //name input file
 // char          fileout1[4*ONE27];// name output file
//  
//  };


struct input_reramp                       // arguments of m/s_crop
  {
   
 // char          filein1[4*ONE27];  //name input file
 // char          fileout1[4*ONE27];// name output file
  
  };
  
  
// ______ ______
struct input_fine                       // arguments for fine coreg.
  {
  char          ifpositions[4*ONE27];   // input file name for positions
  int16         method;                 // method selector
  uint          Nwin;                   // #windows
  uint          MasksizeL;              // size of correlation window
  uint          MasksizeP;              // size of correlation window
  uint          AccL;                   // #lines to be searched in l direction
  uint          AccP;                   // #pixels to be searched in p direction
  int32         initoffsetL;            // initial offset lines
  int32         initoffsetP;            // initial offset pixels
  uint          osfactor;               // oversampling factor
  bool          plotoffsets;            // call script
  bool          plotmagbg;              // call script
  real4         plotthreshold;          // call script  
  int32          shiftazi;               // [true] shift spectrum to 0
  
  // * added by MCC
  //For cohernece method only 
  char          firefdem[4*ONE27];      // input filename reference dem
    char          forefdem[4*ONE27];   // output filename DEM in radarcoord.
  int16         iformatflag;            // input format [signed short]
  uint          demrows;                // number of
  uint          demcols;                // number of
  real8         demdeltalat;            // radians
  real8         demdeltalon;            // radians
  real8         demlatleftupper;        // radians
  real8         demlonleftupper;        // radians
  real8         demnodata;              // identifier/flag

  // * added by MCC  
  };


// ______ ______
struct input_reltiming                  // arguments for timing [FvL]
  {
  real4         threshold;              // threshold for correlation
  int32         maxiter;                // max. #pnts to remove (wtests)
  real4         k_alpha;                // critical value for automated outlier removal
  };


// ______ ______
struct input_demassist                  // arguments for DEM assisted coregistration [FvL]
  {
  char          firefdem[4*ONE27];      // input filename reference dem
  int16         iformatflag;            // input format [signed short]
  uint          demrows;                // number of
  uint          demcols;                // number of
  real8         demdeltalat;            // radians
  real8         demdeltalon;            // radians
  real8         demlatleftupper;        // radians
  real8         demlonleftupper;        // radians
  real8         demnodata;              // identifier/flag
  char          forefdemhei[4*ONE27];   // output filename DEM in radarcoord.
  char          fodem[4*ONE27];         // flag+name output of cropped dem
  char          fodemi[4*ONE27];                // flag+name output of interpolated dem       
  };

// ______ ______
struct input_coregpm                    // arguments for coregpm.
  {
  char          idcoregpm[EIGHTY];
  real4         threshold;              // threshold for correlation
  int32         degree;                 // degree of polynomial
  int32         weightflag;             // 0: all same weight
                                        // 1: choice1: weigh with correlation ??
  int32         maxiter;                // max. #pnts to remove (wtests)
  real4         k_alpha;                // critical value for automated outlier removal
  bool          dumpmodel;              // create float files with model
  bool          plot;                   // plot e_hat etc.
  bool          plotmagbg;              // plot magnitude in background
  };


// ______ range filtering ______
struct input_filtrange                  // arguments for range filter
  {
  int16         method;                 // method selector
  int32         oversample;             // factor
  bool          doweightcorrel;         // weighting of correlation values
  int32         nlmean;                 // number of lines to take mean of
  int32         fftlength;              // length for adaptive
  int32         overlap;                // half overlap between blocks of fftlength
  real8         hammingalpha;           // alpha for hamming
  real8         SNRthreshold;           // spectral peak estimation
  real8         terrainslope;           // [rad] porbits method only
  char          fomaster[4*ONE27];      // output filename
  char          foslave[4*ONE27];       // output filename
  int16         oformatflag;            // output format [cr4] ci16
  };


// ______ ______
struct input_comprefpha                 // arguments for flatearth correction.
  {
  //char                idflatearth[EIGHTY];
  char          ifpositions[4*ONE27];   // input file name for positions
  int16         method;                 // method selector
  int32         degree;                 // degree of polynomial
  int32         Npoints;                // number of observations
  };



// ______ ______
struct input_resample                   // arguments for resampling slave
  {
  int16         method;                 // method selector (interpolator) (%100 == Npoints)
  char          fileout[4*ONE27];
  int16         oformatflag;            // output format [cr4] ci16
  window        dbow_geo;               // cut out of original master.geo
  window        dbow;                   // cut out of original master.radar
  int32          shiftazi;               // [1] center the spectrum around zero, if 0 do nothing if 2 deramp using Sentinel-1 params
  };


// ______ ______
struct input_interfero                  // arguments for computation interferogram
  {
  int16         method;                 // method selector
  char          focint[4*ONE27];                // optional output filename complex interferogram.
  char          foint[4*ONE27];         //  ~ of interferogram (phase).
//  char        foflatearth[EIGHTY];    //  ~ of correction (flatearth) model (phase)
                                        //  these are flags as well as arguments.
                                        //  one is man (else no output)
  uint          multilookL;             // multilookfactor in line dir.
  uint          multilookP;             // multilookfactor in pixel dir.
  };


// ______ ______
struct input_coherence                  // arguments for computation coherence
  {
  int16         method;                 // method selector
  char          focoh[4*ONE27];         // opt output filename of real coherence image.
  char          foccoh[4*ONE27];                //  ~ of complex coherence image.
                                        //  these are flags as well as arguments.
  uint          multilookL;             // multilookfactor in line dir.
  uint          multilookP;             // multilookfactor in pixel dir.
  uint          cohsizeL;               // size of estimation window coherence
  uint          cohsizeP;               // size of estimation window coherence
  };


// ______ ______
struct input_filtphase                  // arguments for phase filter
  {
  int16         method;                 // method selector
  char          fofiltphase[4*ONE27];   // output filename
  char          fifiltphase[4*ONE27];   // input filename
  uint          finumlines;             // number of lines input
  // ______ method goldstein ______
  real8         alpha;                  // weighting
  int32         blocksize;              // blocksize filtered blocks
  int32         overlap;                // half overlap
  // ______ method goldstein, spatial conv. and spectral ______
  matrix<real4> kernel;                 // e.g. [1 1 1]
  // ______ method spatial conv. and spectral ______
  char          fikernel2d[4*ONE27];    // input filename
  };


// ______ ______
struct input_dinsar                     // 3 pass differntial
  {
  //int16       method;                 // method selector
  char          fodinsar[4*ONE27];      // output filename complex interferogram
  char          foscaleduint[4*ONE27];  // output filename scaled uint
  char          topomasterresfile[4*ONE27];// input filename
  char          toposlaveresfile[4*ONE27];// input filename
  char          topointresfile[4*ONE27];        // input filename
  };


// ______ ______
struct input_subtrrefpha                // arguments for subtract 'flat earth'
  {
  int16         method;                 // method selector
  uint          multilookL;             // multilookfactor in line dir.
  uint          multilookP;             // multilookfactor in pixel dir.
  char          focint[4*ONE27];                // output filename complex interferogram
  char          forefpha[4*ONE27];      // output filename complex refpha
  char          foh2ph[4*ONE27];                // output filename h2ph, added by FvL
  bool          dumponlyrefpha;         // do nothing except dump refpha
  };

// ______ ______
struct input_comprefdem                 // arguments for reference phase from DEM
  {
//  int16       method;                 // method selector
  char          firefdem[4*ONE27];      // input filename reference dem
  int16         iformatflag;            // input format [signed short]
  uint          demrows;                // number of
  uint          demcols;                // number of
  real8         demdeltalat;            // radians
  real8         demdeltalon;            // radians
  real8         demlatleftupper;        // radians
  real8         demlonleftupper;        // radians
  real8         demnodata;              // identifier/flag
//  real8               extradense;             // extra interpolation factor (4)
  char          forefdem[4*ONE27];      // output filename reference phase
  char          foh2ph[4*ONE27];                // output perp. baseline, added by FvL
  char          forefdemhei[4*ONE27];   // output filename DEM in radarcoord.
  bool          includerefpha;          // flag to include_flatearth correction
  char          fodem[4*ONE27];         // flag+name output of cropped dem
  char          fodemi[4*ONE27];                // flag+name output of interpolated dem
  bool          isCCC;                 // if demassist is used for Coherent complex coregistration (CCC), then isCCC must be set to true. Default is false
  
  
  
  };

// ______ ______
struct input_subtrrefdem                // arguments for subtract reference DEM
  {
  int16         method;                 // method selector
  int32         offsetL;                // offset applied before subtraction
  int32         offsetP;                // offset applied before subtraction
  char          focint[4*ONE27];                // output filename complex interferogram
  };

// ______ ______
struct input_unwrap                     // arguments for unwrapping
  {
  char          fouint[4*ONE27];                // output filename
  int16         oformatflag;            // output format [hgt] real4
  char          foregions[4*ONE27];     // output filename
  char          seedfile[4*ONE27];      // input file with seeds
  int32         deltaLseed;             // seed delta line direction
  int32         deltaPseed;             // seed delta pixel direction
  int16         method;                 // method selector
  char          snaphu_mode[12];        // snaphu TOPO DEFO SMOOTH NOSTATCOSTS
  char          snaphu_log[4*ONE27];    // log filename for snaphu
  char          snaphu_coh[4*ONE27];    // coherence filename for snaphu opt
  char          snaphu_verbose[6];      // snaphu TRUE or FALSE
  bool          snaphu_dumponlyconf;    // dump snaphu.conf but do not run [BO]
  char          snaphu_init[4];         // snaphu MST or MCF
  int16         ntilerow;
  int16         ntilecol;
  int16         rowovrlp;
  int16         colovrlp;
  int16         nproc;
  int16         tilecostthresh;
  };

// ______ ______ [HB]
struct input_estorbits                    // arguments for orbit error estimation
{
  int16 	method; 		  // method selector
  char          fiheightmap[4*ONE27];     // file with heightmap in radar coordinates, obtained with CRD_OUT_DEM_LP
  char          foresiduals[4*ONE27];     // outfile with observation pixels and (in case of lsq method) residuals
  int16         weighting;                // weighting scheme
  uint          nobs;                     // approximate number of observations to be used
  char          ifpositions[4*ONE27];     // file with (l,p) for observations to be used
  real4         threshold;                // min. coherence for a pixel  to be used
  uint          maxiter;                  // max. number of iterations for data snooping
  real4         k_alpha;                  // critical value for data snooping
  uint          maxfringesaz;             // maximum expected number of orbital fringes in azimuth
  uint          maxfringesrg;             // maximum expected number of orbital fringes in range
  char          reforbitfile[4*ONE27];    // resfile of reference acquisition
  char          foobsdata[4*ONE27];       // file with dumped observation data, estimation design, ...
  uint          poldegree;                // degree of orbit error polynomial
  bool          constrained;              // flag indicating if solution is constained of not
  matrix<int16> constraints;              // k x 2 matrix defining k constraints.
                                          //   1st column: component. =0 for Bpar, =1 for Bperp
                                          //   2nd column: derivative. =0 for constant, =1 for linear, ...
};

// ______ ______
struct input_slant2h                    // arguments for slant2height conversion
  {
  int16         method;                 // method selector
  char          fohei[4*ONE27];         // output filename height
  char          folam[4*ONE27];         // output filename lambda
  char          fophi[4*ONE27];         // output filename phi
  int32         Npoints;                //
  int32         degree1d;               // only {1,2} possible now due to solve33
  int32         degree2d;               //
  int32         Nheights;               //
  };


// ______ ______
struct input_geocode                    // arguments for geocode
  {
//  int16       method;                 // method selector not used right now
  char          fophi[4*ONE27];         // output filename phi
  char          folam[4*ONE27];         // output filename lambda
  };


// ______ ______
struct input_m_EXTRA                    // extra step
  {
  int16         method;                 // method selector
  };

// ______ ______
struct input_s_EXTRA                    // extra step
  {
  int16         method;                 // method selector
  };

// ______ ______
struct input_i_EXTRA2                   // extra step
  {
  int16         method;                 // method selector
  };



// ====== Inline functions ======
// ______ to get/set derfaults ______
inline void setunspecified(char *s)
  {strcpy(s," ");}                      // better use string "unspecified"
inline bool specified(const char *s)
  {return bool(strcmp(s," "));}         // return true if specified

// ====== Prototypes ======
//____RaffaeleNutricato START MODIFICATION SECTION 4
// ______ Reads and interprets "inputoptionsfile" ______
// ______ Fills input structs. ______
void readinput(
      input_gen         &generalinput,
      input_ell         &inputellips,
      input_pr_orbits   &porbitsinput,
      input_readfiles   &readfilesmaster,
      input_morbits     &morbitsinputmaster, // [HB]
      input_crop        &cropinputmaster,
      input_oversample  &oversamplemaster, //____RaffaeleNutricato added this line
      input_simamp      &simampmaster,     // master amplitude simulation, related to master timing [FvL],[MA] TODO: do it flexible if we want for slave also
      input_mtiming     &mtiming,          // [MA]
      input_readfiles   &readfilesslave,
      input_morbits     &morbitsinputslave,  // [HB]
      input_crop        &cropinputslave,
      input_oversample  &oversampleslave,  //____RaffaeleNutricato added this line
      input_filtazi     &filtaziinput,
      input_coarsecorr  &coarsecorrinput,
      input_fine        &fineinput,
      input_reltiming   &timinginput,   //[FvL]
      input_demassist   &demassistinput,   //[FvL]
      input_coregpm     &coregpminput,
      input_resample    &resampleinput,
      input_filtrange   &filtrangeinput,
      input_interfero   &interferoinput,
      input_coherence   &coherenceinput,
      input_comprefpha  &comprefphainput,
      input_subtrrefpha &subtrrefphainput,
      input_comprefdem  &comprefdeminput,
      input_subtrrefdem &subtrrefdeminput,
      input_filtphase   &filtphaseinput,
      input_dinsar      &dinsarinput,
      input_unwrap      &unwrapinput,
      input_estorbits   &estorbitsinput,   // [HB]
      input_slant2h     &slant2hinput,
      input_geocode     &geocodeinput);
//____RaffaeleNutricato END MODIFICATION SECTION 4


#endif // READINPUT_H



