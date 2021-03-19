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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readinput.cc,v $
 * $Revision: 3.36 $
 * $Date: 2009/01/09 11:09:20 $
 * $Author: TUDelft $
 *
 * implementation of readinput.
 ****************************************************************/

#include "matrixbk.hh"
#include "constants.hh"         // global constants
#include "ioroutines.hh"        // ?
#include "utilities.hh"         // ispower2
#include "exceptions.hh"        // my exceptions class

#include <fstream>              // for file streams
#include <strstream>            // for file streams
#include <iomanip>              // for setprecision etc.
#include <cstring>              // for strcmp etc.
#include <cstdlib>              // exit, atoi
#include <cctype>               // isspace



// ______ displevel used in ioroutines.h, changed here ______
char WARNS[6][ONE27];           // remember 6 last warnings in WARNS
int32   beeplevel=1;            // global variable for beeping.
                                // 0:  nobeep; BEEP OFF
                                // -1: beep on error exit ; BEEP ERROR
                                // 1:  beep on error, warnings; BEEP WARNING
                                // 2:  beep on error, warnings, progress;
                                //     BEEP PROGRESS, BEEP [ON]
int32   displevel=30000;        // controls level of screen output
                                // -100 only errors
                                // 0:     warnings and errors
                                // 10000: progress, warn and err
                                // 20000: info, pro, warn and err
                                // 30000: debug, info, pro, warn, err


// ______ Checks input prototypes (see below in this file) ______
void checkgeneral     (input_gen               &generalinput, const int16 onlyprocess);
void checkreadfiles   (const input_readfiles   &readfilesinput, const int16 id);
void checkcrop        (const input_crop         &cropinput, const int16 id);
//____RaffaeleNutricato START MODIFICATION SECTION 1
void checkoversample  (const input_oversample  &oversampleinput, const int16 id);
//____RaffaeleNutricato END MODIFICATION SECTION 1
void checkporbits     (const input_pr_orbits   &porbitsinput, const int16 id);
void checksimamp      (const input_simamp      &simampinput);        //[MA] ?nec. id check?
void checkmtiming     (const input_mtiming     &mtiminginput);    //[MA] ?nec. id check?
void checkslant2h     (const input_slant2h     &slant2hinput);
void checkunwrap      (const input_unwrap      &unwrapinput);
void checkgeocode     (const input_geocode     &geocodeinput);
void checkcoarsecorr  (const input_coarsecorr  &coarsecorrinput);
void checkfine        (const input_fine        &fineinput);
void checkreltiming   (const input_reltiming   &reltiminginput); //[FvL]
void checkdemassist   (const input_demassist   &demassistinput); //[FvL]
void checkcoregpm     (const input_coregpm     &coregpminput);
void checkcomprefpha  (const input_comprefpha  &comprefphainput);
void checksubtrrefpha (const input_subtrrefpha &subtrrefphainput);
void checkresample    (const input_resample    &resampleinput);
void checkinterfero   (const input_interfero   &interferoinput);
void checkcoherence   (const input_coherence   &coherenceinput);
void checkcomprefdem  (const input_comprefdem  &comprefdeminput);
void checksubtrrefdem (const input_subtrrefdem &subtrrefdeminput);
void checkfiltrange   (const input_filtrange   &filtrangeinput);
void checkdinsar      (const input_dinsar      &dinsarinput);
void checkfiltphase   (const input_filtphase   &filtphaseinput);
void checkfiltazi     (const input_filtazi     &filtaziinput, const int16 id);

/****************************************************************
 *    writearg                                                  *
 * echo arg to screen if debug defined                          *
 #%// BK 13-Jul-2000                                            *
 ****************************************************************/
template <class Type>
void writearg(const Type argument)
  {
  TRACE_FUNCTION("writearg");
  DEBUG << "Read argument: " << argument;
  DEBUG.print();
  } // END writearg



/****************************************************************
 *    readinput                                                 *
 *                                                              *
 * Read and interpret                                           *
 * "inputoptionsfile" (file in variable: logfile).              *
 * Write to logfile (via "scratchreadinput")                    *
 * mandatory input is checked on presence by checksums (cs).    *
 * If there are more methods, process card switches default     *
 *  methodselector (if present) is used later to correct        *
 *  see for example unwrap methodeselector                      *
 * linecounter only used in case of errors                      *
 *                                                              *
 * input:                                                       *
 *  - struct: general input options                             *
 *  - struct: master readfiles input options                    *
 *  - struct: slave readfiles input options                     *
 *  - struct: master crop input options                         *
 *  - struct: slave crop input options                          *
 *  - ...                                                       *
 * output:                                                      *
 *  - (updated input structs)                                   *
 *  - (file copy of input)                                      *
 *                                                              *
 *    Bert Kampes,      11-Dec-1998                             *
 *    Mahmut Arikan,    09-Jan-2009 (code to fix linecnt bug    *
 *    G.J. van Zwieten, 09-Jan-2009  and line extend         )  *
 *    Mahmut Arikan,    29-Aug-2010 Fix: Read only those prms   *
 *                                  for processing              *
 ****************************************************************/
void readinput(
        input_gen         &generalinput, 
        input_ell         &ellipsinput,
        input_pr_orbits   &porbitsinput,
        input_readfiles   &m_readfilesinput,
        input_morbits     &morbitsinputmaster, // [HB]
        input_crop         &m_cropinput,
//____RaffaeleNutricato START MODIFICATION SECTION 2
        input_oversample  &m_oversample, 
        input_simamp      &simampinput,                   // master simamp [FvL],[MA]
        input_mtiming     &mtiminginput,                  // mtiming correl [MA]
//____RaffaeleNutricato END MODIFICATION SECTION 2
        input_readfiles   &s_readfilesinput,
        input_morbits     &morbitsinputslave,  // [HB]
        input_crop         &s_cropinput,
//____RaffaeleNutricato START MODIFICATION SECTION 3
        input_oversample  &s_oversample, 
//____RaffaeleNutricato END MODIFICATION SECTION 3
        input_filtazi     &filtaziinput,
        input_coarsecorr  &coarsecorrinput,
 //       input_deramp      &derampinput, //MCC
        input_fine        &fineinput,
        input_reltiming   &reltiminginput, //[FvL]
        input_demassist   &demassistinput, //[FvL]
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
        input_geocode     &geocodeinput) 
  {
  //TRACE_FUNCTION("readinput (BK 11-Dec-1998)");
  TRACE_FUNCTION("readinput rev.5 (TUDelft 29-Aug-2010)");

  // ______ Set ids ______
  m_readfilesinput.fileid = MASTERID;
  m_cropinput.fileid      = MASTERID;
  s_readfilesinput.fileid = SLAVEID;
  s_cropinput.fileid      = SLAVEID;

  // ______ Misuse generalinput.logfile to store name: open file here! ______
  //ifstream optionsfile(generalinput.logfile, ios::in | ios::nocreate);
  ifstream optionsfile(generalinput.logfile, ios::in);
  bk_assert(optionsfile,generalinput.logfile,__FILE__,__LINE__);
  string inputoptionsfile = generalinput.logfile;   // [MA] keep input filename
                                                    // later variable is updated
                                                    // as logfile name.

  int16                 onlyprocess     = -1;        // flag for ONLYPROCESS card
  int16                 linecnt         =  0;        // counter
  const int16           BASE10          = 10;        // [MA] base 10 defined for strtol
  //char                  keyword[EIGHTY];
  //char                  filename[4*ONE27];            // string for filenames // [MA] changed EIGHTY --> 2*ONE27, due to comments line gets longer
  char                  eachline[4*ONE27];           // assuming maximum char lenght of the line is 4*ONE27. It should be sufficient.

  // ______ Check (multiple) occurrence of cards ______
  bool                  priorscreen     = false;     // no screen card present
  bool                  priormemory     = false;     // check if present for info
  bool                  priorbatch      = false;     // check if present for info
  bool                  prioroverwrite  = false;     // check if present for info
  bool                  priorlistinput  = false;     // check if present for info
  bool                  priorrs_fileout = false;     // 


// ====== Initialization, defaults ======
  input_ell             WGS84;
  bool                  listinput  = true;      // default copy input to log
  bool                  ellipsoid  = false;     // set default if no card present
  ellipsinput                   = WGS84;        // default

  register int32 i;
  for (i=0; i<NUMPROCESSES; i++)
    generalinput.process[i] = 0;                         // default no processing step
  generalinput.interactive      = true;                 // default interactive mode
  generalinput.overwrit         = false;                // default no overwriting
  generalinput.memory           = real8(MEMORY_DEF);           // default memory (500MB), see constants.hh
  strcpy(generalinput.logfile,    "log.out");           // default logfile
  strcpy(generalinput.m_resfile,  "master_result.out"); // default resultfile
  strcpy(generalinput.s_resfile,  "slave_result.out");  // default resultfile
  strcpy(generalinput.i_resfile,  "interferogram.out"); // default interf_out
  generalinput.orb_interp       = ORB_DEFAULT;          // default polyfit
  generalinput.dumpbaselineL    = 0;                    // default no dump
  generalinput.dumpbaselineP    = 0;                    // default no dump
  generalinput.preview          = 0;                    // default no preview
  generalinput.terrain_height   = 0.0;                  // above ellipsoid
  generalinput.tiepoint.x       = 0.0;                  // default none
  generalinput.tiepoint.y       = 0.0;                  // default none
  generalinput.tiepoint.z       = 0.0;                  // default none

  //m_readfilesinput.method     = readfiles_ers;        // default
  m_readfilesinput.sensor_id    = SLC_ERS;              // default
  setunspecified(m_readfilesinput.volfile);             // check later, then set default
  setunspecified(m_readfilesinput.leaderfile);          // check later, then set default
  setunspecified(m_readfilesinput.nullfile);            // check later, then set default
  setunspecified(m_readfilesinput.datfile);             // check later, then set default
  m_readfilesinput.rg_timing_error = 0.0;               // default
  m_readfilesinput.az_timing_error = 0.0;               // default
  //s_readfilesinput.sensor_id          = readfiles_ers;        // default
  s_readfilesinput.sensor_id    = SLC_ERS;              // default
  setunspecified(s_readfilesinput.volfile);             // check later, then set default
  setunspecified(s_readfilesinput.leaderfile);          // check later, then set default
  setunspecified(s_readfilesinput.nullfile);            // check later, then set default
  setunspecified(s_readfilesinput.datfile);             // check later, then set default
  s_readfilesinput.rg_timing_error = 0.0;               // default
  s_readfilesinput.az_timing_error = 0.0;               // default

  setunspecified(porbitsinput.m_orbdir);                // check later, then set default
  setunspecified(porbitsinput.s_orbdir);                // check later, then set default
  porbitsinput.timeinterval     = 1;                    // default time interval
  porbitsinput.timebefore       = 4*porbitsinput.timeinterval;   // default 4 extra datapoints
  porbitsinput.dumpmasterorbit  = -1.;                  // default no dump
  porbitsinput.dumpslaveorbit   = -1.;                  // default no dump

  strcpy(m_cropinput.idcrop,"master step01");           // default identifier
  strcpy(s_cropinput.idcrop,"slave step01");            // default identifier
  strcpy(m_cropinput.fileout1,"master.raw");            // default output filename
  strcpy(s_cropinput.fileout1,"slave.raw");             // default output filename
  m_cropinput.dbow.linelo       = 0;                    // min. line coord. initialization
  m_cropinput.dbow.linehi       = 0;                    // max. line coord. initialization
  m_cropinput.dbow.pixlo        = 0;                    // min. pixel coord. initialization
  m_cropinput.dbow.pixhi        = 0;                    // max. pixel coord. initialization
  s_cropinput.dbow.linelo       = 0;                    // min. line coord. initialization
  s_cropinput.dbow.linehi       = 0;                    // max. line coord. initialization
  s_cropinput.dbow.pixlo        = 0;                    // min. pixel coord. initialization
  s_cropinput.dbow.pixhi        = 0;                    // max. pixel coord. initialization
  m_cropinput.dbow_geo.linelo   = 0;                    // min. line coord. initialization
  m_cropinput.dbow_geo.linehi   = 0;                    // max. line coord. initialization
  m_cropinput.dbow_geo.pixlo    = 0;                    // min. pixel coord. initialization
  m_cropinput.dbow_geo.pixhi    = 0;                    // max. pixel coord. initialization
  s_cropinput.dbow_geo.linelo   = 0;                    // min. line coord. initialization
  s_cropinput.dbow_geo.linehi   = 0;                    // max. line coord. initialization
  s_cropinput.dbow_geo.pixlo    = 0;                    // min. pixel coord. initialization
  s_cropinput.dbow_geo.pixhi    = 0;                    // max. pixel coord. initialization

//____RaffaeleNutricato START MODIFICATION SECTION 4
  m_oversample.OsrRange          = 1;                    // Default for oversampling ratio in range.
  m_oversample.OsrAzimuth        = 1;                    // Default for oversampling ratio in azimuth.
  m_oversample.FilterSize        = 16;                   // Default for length of the interpolation kernel in range. 
  m_oversample.oformatflag       = FORMATCI2;            // Default for output format.
  strcpy(m_oversample.fileoutovs,"master_ovs.raw");      // Default output filename.
  s_oversample.OsrRange          = 1;                    // Default for oversampling ratio in range.
  s_oversample.OsrAzimuth        = 1;                    // Default for oversampling ratio in azimuth.
  s_oversample.FilterSize        = 16;                   // Default for length of the interpolation kernel in range. 
  s_oversample.oformatflag       = FORMATCI2;            // Default for output format.
  strcpy(s_oversample.fileoutovs,"slave_ovs.raw");       // Default output filename.
//____RaffaeleNutricato END MODIFICATION SECTION 4

  //____ added by MA ____

  setunspecified(simampinput.firefdem);                                // check later, mandatory
  //setunspecified(simampinput.fodemi);                                // check later, then set default
  setunspecified(simampinput.fodemlp);                           // check later, then set default  // MA SPT
  setunspecified(simampinput.fothetalp);                                // check later, then set default  // MA SPT
  strcpy(simampinput.fodem,"demcrop_sam.raw");                         // default name
  setunspecified(simampinput.fosimamp);                                // check later, then set default
  simampinput.iformatflag     =  FORMATI2;                             // default gtopo30
  simampinput.demrows         =  6000;                                 // default gtopo30
  simampinput.demcols         =  4800;                                 // default gtopo30
  simampinput.demnodata       = -9999;                                 // default gtopo30
  simampinput.demdeltalat     = deg2rad(0.00833333333333333333);       // default gtopo30
  simampinput.demdeltalon     = deg2rad(0.00833333333333333333);       // default gtopo30
  simampinput.demlatleftupper = deg2rad(89.995833333333333333);        // w020n90.DEM
  simampinput.demlonleftupper = deg2rad(-19.995833333333333333);       // w020n90.DEM
  strcpy(simampinput.fosimamp,"master.sam");                           // default name


  const int32 def_mte_nwin   = 16;                                      // default #windows           
  setunspecified(mtiminginput.ifpositions);                             // check later, then set default
  mtiminginput.MasksizeL     = 256;                                     // default correlation size
  //mtiminginput.MasksizeP     = mtiminginput.MasksizeL;                // default correlation size
  mtiminginput.MasksizeP     = 128;                                     // default correlation size
  mtiminginput.AccL          = 32;                                      // default searching limit
  mtiminginput.AccP          = mtiminginput.AccL;                       // default searching limit
  mtiminginput.initoffsetL   = 0;                                       // default initial offset
  mtiminginput.initoffsetP   = 0;                                       // default initial offset

  // ____ end added by MA ____

  // ____ start added by HB ____
  morbitsinputmaster.coeff        = matrix<real8>(10,2);  // 
  morbitsinputmaster.coeff        = 0;                    // default: no modification
  setunspecified(morbitsinputmaster.reforbitfile);        //
  morbitsinputslave.coeff         = matrix<real8>(10,2);  // 
  morbitsinputslave.coeff         = 0;                    // default: no modification
  setunspecified(morbitsinputslave.reforbitfile);         //
  // ____ end added by HB ____

  filtaziinput.oformatflag      = FORMATCR4;            // default
  filtaziinput.fftlength        = 1024;                 // default
  filtaziinput.overlap          = -1;                   // default to fftlength/8
  filtaziinput.hammingalpha     = 0.75;                 // default (slc)
  strcpy(filtaziinput.fomaster,"master.afilter");       // default
  strcpy(filtaziinput.foslave,"slave.afilter");         // default

  const int32 def_cc_nwin       = 11;                   // default #windows
  setunspecified(coarsecorrinput.ifpositions);          // check later, then set default
  coarsecorrinput.MasksizeL     = 64;                   // default correlation size
  coarsecorrinput.MasksizeP     = coarsecorrinput.MasksizeL;// default correlation size
  coarsecorrinput.AccL          = 8;                    // default searching limit
  coarsecorrinput.AccP          = coarsecorrinput.AccL; // default searching limit
  coarsecorrinput.initoffsetL   = 0;                    // default initial offset
  coarsecorrinput.initoffsetP   = 0;                    // default initial offset

  //deramp.filein1
  setunspecified(fineinput.ifpositions);                // check later, then set default
  const int32 def_fc_nwin       = 601;                  // default #windows
  fineinput.MasksizeL           = 64;                   // default correlation size
  fineinput.MasksizeP           = fineinput.MasksizeL;  // default correlation size
  fineinput.AccL                = 8;                    // default searching limit
  fineinput.AccP                = fineinput.AccL;       // default searching limit
  fineinput.initoffsetL         = 0;                    // default initial offset
  fineinput.initoffsetP         = 0;                    // default initial offset
  fineinput.osfactor            = 32;                   // default oversampling factor
  fineinput.plotoffsets         = false;                // default no plotting
  fineinput.plotmagbg           = false;                // default no plotting
  fineinput.plotthreshold       = 0.3;                  // default no plotting
  fineinput.shiftazi            = 1;                 // [1] shift spectrum to 0
  
  //ADD by MCC for fine CCCoregistration
  
  setunspecified(fineinput.forefdem);
  setunspecified(fineinput.firefdem);
  fineinput.iformatflag   = FORMATI2;              // default gtopo30
  fineinput.demrows       = 6000;                  // default gtopo30
  fineinput.demcols       = 4800;                  // default gtopo30
  fineinput.demnodata     = -9999;                 // default gtopo30
  fineinput.demdeltalat   = deg2rad(0.00833333333333333333);// default gtopo30
  fineinput.demdeltalon   = deg2rad(0.00833333333333333333);// default gtopo30
  fineinput.demlatleftupper = deg2rad(89.995833333333333333);      // w020n90.DEM
  fineinput.demlonleftupper = deg2rad(-19.995833333333333333);// w020n90.DEM
 
  //added by MCC for fine CCCoregistration
  
 
  
  
  //____ added by FvL ____
  reltiminginput.threshold         = 0.4;               // default threshold data
  reltiminginput.maxiter           = 10000;             // default max. 10000 outliers removed
  reltiminginput.k_alpha           = 1.97;              // critical value for outliers

  setunspecified(demassistinput.firefdem);              // check later, mandatory
  setunspecified(demassistinput.fodemi);                // check later, then set default
  setunspecified(demassistinput.forefdemhei);           // check later, then set default
  strcpy(demassistinput.fodem,"demcrop.raw");           // default name
  demassistinput.iformatflag   = FORMATI2;              // default gtopo30
  demassistinput.demrows       = 6000;                  // default gtopo30
  demassistinput.demcols       = 4800;                  // default gtopo30
  demassistinput.demnodata     = -9999;                 // default gtopo30
  demassistinput.demdeltalat   = deg2rad(0.00833333333333333333);// default gtopo30
  demassistinput.demdeltalon   = deg2rad(0.00833333333333333333);// default gtopo30
  demassistinput.demlatleftupper = deg2rad(89.995833333333333333);      // w020n90.DEM
  demassistinput.demlonleftupper = deg2rad(-19.995833333333333333);// w020n90.DEM
  // ____ end added by FvL ____

  coregpminput.threshold        = 0.2;                  // default threshold data
  coregpminput.degree           = 1;                    // default degree polynomial
  coregpminput.weightflag       = 0;                    // default no weighting of data
  coregpminput.maxiter          = 10000;                // default max. 10000 outliers removed, changed by [FvL]
  coregpminput.k_alpha          = 1.97;                 // critical value for outliers
  coregpminput.dumpmodel        = false;                // default no files
  coregpminput.plot             = false;                // default no plots
  coregpminput.plotmagbg        = false;                // default no plots

  filtrangeinput.method         = rf_adaptive;          // default
  filtrangeinput.terrainslope   = 0.0;                  // default porbits
  filtrangeinput.fftlength      = -999;                 // set default later
  filtrangeinput.overlap        = 0;                    // default no overlap
  filtrangeinput.nlmean         = 15;                   // default
  filtrangeinput.hammingalpha   = 0.75;                 // default
  filtrangeinput.SNRthreshold   = 5.0;                  // default
  filtrangeinput.oversample     = 2;                    // default
  filtrangeinput.doweightcorrel = false;                // default
  strcpy(filtrangeinput.fomaster,"master.rfilter");     // default
  strcpy(filtrangeinput.foslave,"slave.rfilter");       // default
  filtrangeinput.oformatflag    = FORMATCR4;            // default

  setunspecified(comprefphainput.ifpositions);          // check later, then set default
  const int32 def_fe_Npoints    = 501;                  // default
  const int32 def_fe_degree     = 5;                    // default

  strcpy(resampleinput.fileout,"s_resampled.raw");      // default
  resampleinput.oformatflag     = FORMATCR4;            // default
  resampleinput.dbow.linelo     = 0;                    // min. line coord. initialization
  resampleinput.dbow.linehi     = 0;                    // max. line coord. initialization
  resampleinput.dbow.pixlo      = 0;                    // min. pixel coord. initialization
  resampleinput.dbow.pixhi      = 0;                    // max. pixel coord. initialization
  resampleinput.dbow_geo.linelo = 0;                    // min. line coord. initialization
  resampleinput.dbow_geo.linehi = 0;                    // max. line coord. initialization
  resampleinput.dbow_geo.pixlo  = 0;                    // min. pixel coord. initialization
  resampleinput.dbow_geo.pixhi  = 0;                    // max. pixel coord. initialization
  resampleinput.shiftazi        = 1;                    // default [1] apply shift

  interferoinput.method         = int_oldmethod;        // default method
  setunspecified(interferoinput.focint);                // check later, then set default
  setunspecified(interferoinput.foint);                 // use later if specified
  //setunspecified(interferoinput.foflatearth);         // check later, then set default
  interferoinput.multilookL     = 5;                    // default multilookfactor
  interferoinput.multilookP     = 1;                    // default multilookfactor

  coherenceinput.method         = coh_oldmethod;        // default method
  setunspecified(coherenceinput.focoh);                 // check later, then set default
  setunspecified(coherenceinput.foccoh);                // check later, then set default
  coherenceinput.multilookL     = 10;                   // default multilookfactor
  coherenceinput.multilookP     = 2;                    // default multilookfactor
  coherenceinput.cohsizeL       = coherenceinput.multilookL;    // default windowsize. 
  coherenceinput.cohsizeP       = coherenceinput.multilookP;    // default windowsize. 

  subtrrefphainput.method       = srp_polynomial;       // default method
  strcpy(subtrrefphainput.forefpha, "refphase.raw");    // default name
  strcpy(subtrrefphainput.focint,"cint.minrefpha.raw");// default
  subtrrefphainput.multilookL   = 1;                    // default multilookfactor
  subtrrefphainput.multilookP   = 1;                    // default multilookfactor
  subtrrefphainput.dumponlyrefpha  = false;             // default not
  //_____ added by FvL
  setunspecified(subtrrefphainput.foh2ph);              // check later, then set default
  // ____ end added by FvL

  filtphaseinput.method         = fp_goldstein;         // default method
  filtphaseinput.alpha          = 0.2;                  // default
  // ______ 32 blocks default for goldstein, kernel as large as possible, ______
  // ______ 2Dkernel then specify ______
  filtphaseinput.blocksize      = 32;                   // default
  filtphaseinput.overlap        = 3;                    // default
  // set default later with alpha in name
  setunspecified(filtphaseinput.fofiltphase);           // check later, then set default
  setunspecified(filtphaseinput.fifiltphase);           // if specified, use it
  setunspecified(filtphaseinput.fikernel2d);            // if specified, use it
  filtphaseinput.finumlines     = 0;                    // numlines.

  strcpy(dinsarinput.fodinsar,"differentialinterf.raw");// default
  setunspecified(dinsarinput.foscaleduint);             // default no output
  // ______ set mastertopo file to same as master resfile if not specified ______
  setunspecified(dinsarinput.topomasterresfile);        // if specified, use 4 pass
  setunspecified(dinsarinput.toposlaveresfile);         // check later, mandatory
  setunspecified(dinsarinput.topointresfile);           // check later, mandatory

  setunspecified(comprefdeminput.firefdem);             // check later, mandatory
  setunspecified(comprefdeminput.fodemi);               // check later, then set default
  //_____ added by FvL
  setunspecified(comprefdeminput.foh2ph);               // check later, then set default
  // ____ end added by FvL
  setunspecified(comprefdeminput.forefdemhei);          // check later, then set default
  strcpy(comprefdeminput.forefdem,"refdem.raw");        // default name
  strcpy(comprefdeminput.fodem,"demcrop.raw");          // default name [FvL]
  comprefdeminput.iformatflag   = FORMATI2;             // default gtopo30
//  comprefdeminput.method        = crd_trilinear;      // default method
//  comprefdeminput.extradense    = 0.5;                        // default (now interpolated in l,p)
  comprefdeminput.demrows       = 6000;                 // default gtopo30
  comprefdeminput.demcols       = 4800;                 // default gtopo30
  comprefdeminput.demnodata     = -9999;                // default gtopo30
  comprefdeminput.demdeltalat   = deg2rad(0.00833333333333333333);// default gtopo30
  comprefdeminput.demdeltalon   = deg2rad(0.00833333333333333333);// default gtopo30
  comprefdeminput.demlatleftupper = deg2rad(89.995833333333333333);     // w020n90.DEM
  comprefdeminput.demlonleftupper = deg2rad(-19.995833333333333333);// w020n90.DEM
  
  comprefdeminput.isCCC =false; // MCC
  
  strcpy(subtrrefdeminput.focint,"cint.minrefdem.raw"); // default name
  subtrrefdeminput.offsetL      = 0;                    // default no offset
  subtrrefdeminput.offsetP      = 0;                    // default no offset

  strcpy(unwrapinput.fouint,"unwrapped_interferogram.raw");     // default name
  strcpy(unwrapinput.foregions,"regions_unwrapped.raw");        // default name
  setunspecified(unwrapinput.seedfile);                 // check later, then set default
  unwrapinput.deltaLseed        = 100;                  // default 100 pixels;
  unwrapinput.deltaPseed        = unwrapinput.deltaLseed;       // default 100 pixels;
  setunspecified(unwrapinput.snaphu_log);// " "
  setunspecified(unwrapinput.snaphu_coh);// " "
  strcpy(unwrapinput.snaphu_mode,"DEFO");               // default to DEFO from TOPO
  strcpy(unwrapinput.snaphu_init,"MST");                // default method
  strcpy(unwrapinput.snaphu_verbose,"TRUE");            // default verbose
  unwrapinput.oformatflag       = FORMATR4;             // default to
                                                        // REAL4 from FORMATHGT

  // block for Tile.CONTROL: only for SNAPHU; defaults for single CPU
  unwrapinput.ntilerow          = 1; // number of tiles in range
  unwrapinput.ntilecol          = 1; // number of tiles in azimuth
  unwrapinput.rowovrlp          = 0; // overlap between tiles in rng
  unwrapinput.colovrlp          = 0; // overlap between tiles in az
  unwrapinput.nproc             = 1; // no.cpus or nodes on load
                                     // balancing cluster
  unwrapinput.tilecostthresh   = 500; // cost threshold boundaries of reliable regions
 

  slant2hinput.Npoints          = 200;                  // default 100 pixels;
  slant2hinput.degree1d         = 2;                    // default 2 
  slant2hinput.degree2d         = 5;                    // default 5
  slant2hinput.Nheights         = slant2hinput.degree1d+1;// minimum
  strcpy(slant2hinput.fohei,"hei.raw");                 // default
  strcpy(slant2hinput.fophi,"phi.raw");                 // default
  strcpy(slant2hinput.folam,"lam.raw");                 // default

  strcpy(geocodeinput.fophi,"geo_phi.raw");             // default name
  strcpy(geocodeinput.folam,"geo_lambda.raw");          // default name

  // ____ start added by HB ____
  estorbitsinput.method          = eo_lsq;               // default: least squares method
  setunspecified(estorbitsinput.fiheightmap);            // default: use no heights 
  setunspecified(estorbitsinput.foresiduals);            // default: no output
  estorbitsinput.weighting       = eo_noweighting;       // default: unweighted
  const uint def_eo_nobs         = 1000;                 // default value, apply later
  setunspecified(estorbitsinput.ifpositions);            // default: no file input
  estorbitsinput.threshold       = 0;                    // default: ignore coherence estimate
  estorbitsinput.maxiter         = 0;                    // default: no iteration
  estorbitsinput.k_alpha         = 3.29;                 // default: norminv(.9995) \approx tinv(.9995,redund) 
  estorbitsinput.maxfringesaz    = 15;                   //
  estorbitsinput.maxfringesrg    = 15;                   //
  setunspecified(estorbitsinput.reforbitfile);           // default: use master as reference
  setunspecified(estorbitsinput.foobsdata);              // default: do not dump observation data
  estorbitsinput.poldegree       = 1;                    //
  estorbitsinput.constrained     = true;                 // default: constrain BPAR0 BPERP1 (define later)
  // end added by HB ____

// ====== (default) Methods ======
  const int16 def_mte_method     = cc_magspace;
  const int16 def_cc_method     = cc_magfft;
  const int16 def_fc_method     = fc_magfft;
  const int16 def_fe_method     = fe_porbits;
  const int16 def_rs_method     = rs_cc4p;
  const int16 def_uw_method     = uw_method2;// snaphu, system call, see www
  const int16 def_s2h_method    = s2h_ambiguity;

// ______ To check later if default should be used ______
  mtiminginput.method           = def_mte_method - 999; // default method (repair later)
  mtiminginput.Nwin             = def_mte_nwin   + 999; // default #windows
  coarsecorrinput.method        = def_cc_method  - 999; // default method (repair later)
  coarsecorrinput.Nwin          = def_cc_nwin    + 999; // default #windows
  fineinput.method              = def_fc_method  - 999; // default method (repair later)
  fineinput.Nwin                = def_fc_nwin    + 999; // default #windows
  unwrapinput.method            = def_uw_method  - 999; // default method (repair later)
  comprefphainput.method        = def_fe_method  - 999; // default method (repair later)
  comprefphainput.degree        = def_fe_degree  - 999; // default degree polynomial
  comprefphainput.Npoints       = def_fe_Npoints - 999; // default degree polynomial
  resampleinput.method          = def_rs_method  - 999; // default method (repair later)
  slant2hinput.method           = def_s2h_method - 999; // default method (repair later)
  estorbitsinput.nobs            = 0;                    // set default later





// ====== Process "inputoptionsfile" ======
  bool continuereading = true;
  //while (! optionsfile.eof())                               // read file until end of file.
  //while (continuereading || ! optionsfile.eof())              // read file until STOP card or if fails exit at end of the file
  while (continuereading)              // read file until STOP card
    {
    linecnt++;
    optionsfile.getline(eachline,4*ONE27,'\n');                 // get line. [MA] if the line is longer than 4*ONE27 then loops forever.
    //cerr << "line: " << linecnt << " : " << eachline << endl;

    if ( optionsfile.eof() ) // && strcmp(keyword,"STOP") )      // [MA] STOP is missing
      {
      ERROR << "STOP:   \t is missing"
            << " in file '" << inputoptionsfile << "'";
      PRINT_ERROR(ERROR.get_str())
      WARNING.print("Please make sure the inputfile line: STOP has ended with [eol] character");
      throw(keyword_error);
      continuereading = false;                        // break while loop
      }

    const int maxwords = 8;   // [GJ, MA]
    char *word[ maxwords ];

    char *c = eachline;
    for ( int i = 0; i < maxwords; i++ )
    {
      while ( *c == ' ' || *c == '\t' ) // get rid of leadin white space
      {
        c++;                            // next address
      }
      word[ i ] = c;                    // pass first char address of a word

      while ( *c != ' ' && *c != '\t' && *c != '\0' )
      {
       c++;
      }

      if ( *c != '\0' ) // at last char
      {
        *c = '\0';
        c++;
      }
    }

    char *keyword = word[0];    // start off with the card.
                                // word[1] --> first argument
                                // word[2] --> second argument and so on.
    toupper(keyword);

    DEBUG << linecnt << ": Read keyword: " << keyword;
    DEBUG.print();

// *******************************************************************
// *** GENERAL
// *******************************************************************
    if (!strcmp(keyword,"COMMENT") ||
        !strcmp(keyword,"C")) 
      {
      ;                                                                 // comment: no action 
      }
    else if (!strncmp(keyword,"//",2) ||            // assume user comments out
             !strncmp(keyword,"#",1))                     // with // or '#' as delimiter
      {                                                                   // but no blank after keyword
      ;                                                                       // comment: no action
      }
    //else if (!strncmp(keyword,'\0',1))          // empty line? crashes
    else if (!strlen(keyword))                        // empty line
      {                                                
      ;                                                                       // comment: no action
      }

// **********************************************************************
    else if (!strcmp(keyword,"BEEP"))     // level of beep output (THE CARD)
      {                                   // /progress/[warning]/error/ON/OFF
      //keyword =  word[1] ;    // pass keyword                 // argument
      //writearg((char*)keyword);
      keyword = word[1];                  // pass next word       (the argument)
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"ERROR"))
        {
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(0);
        WARNING.bellrings(0);
        ERROR.bellrings(1);
        beeplevel = -1;
        INFO.print("BEEP: \tbeeping enabled at level: \tERROR");
        }
      else if (!strcmp(keyword,"PROGRESS"))
        {
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(1);
        WARNING.bellrings(2);
        ERROR.bellrings(3);
        beeplevel =  2;
        INFO.print("BEEP: \tbeeping enabled at level: \tPROGRESS");
        }
      else if (!strcmp(keyword,"WARNING"))
        {
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(0);
        WARNING.bellrings(1);
        ERROR.bellrings(2);
        beeplevel =  1;
        INFO.print("BEEP: \tbeeping enabled at level: \tWARNING");
        }
      else if (!strcmp(keyword,"OFF"))
        {
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(0);
        WARNING.bellrings(0);
        ERROR.bellrings(0);
        beeplevel =  0;
        INFO.print("BEEP: \tbeeping disabled");
        }
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||      // comment
               !strncmp(keyword,"#",1) ||       // comment
         //!strcmp(keyword,\'\\0\'))
               !(keyword[0] == '\0'))           // no keyword
//             !strcmp(keyword,""))             // no keyword
        {
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(1);
        WARNING.bellrings(2);
        ERROR.bellrings(3);
        beeplevel =  2;
        INFO.print("BEEP: \tbeeping enabled for all levels: \tON");
        }
      else
        {
        beeplevel = 1;
        TRACE.bellrings(0);
        DEBUG.bellrings(0);
        INFO.bellrings(0);
        PROGRESS.bellrings(0);
        WARNING.bellrings(1);
        ERROR.bellrings(2);
        WARNING << "BEEP:   line " << linecnt
             << ": Argument " << keyword 
             << " not recognized."
             << " [error/warning/progress/on/off] I used WARNING.";
        WARNING.print();
        }
      } // BEEP key

// **********************************************************************
    else if (!strcmp(keyword,"SCREEN")) // level of screen output 
      {                                 // debug/info/progress/warning or error
      switch (priorscreen)
        {
        case true:
          WARNING << "SCREEN: line " << linecnt << ": stdout: "
               << " ignored due to prior occurrence.";
          WARNING.print();
          break;

        default:
          priorscreen = true;
          //keyword =  word[1] ;        // pass keyword                         // argument
          keyword = word[1];                    // argument
          writearg(keyword);
          toupper(keyword);
          if (!strcmp(keyword,"INFO"))
            {
            TRACE.doprint(0);
            DEBUG.doprint(0);
            INFO.doprint(1);
            PROGRESS.doprint(1);
            WARNING.doprint(1);
            ERROR.doprint(1);
            displevel = 20000 + displevel%10000;        // for cnt #warnings
            INFO.print("SCREEN: \tverboseness: \t\t\tINFO");
            }
          else if (!strcmp(keyword,"PROGRESS"))
            {
            TRACE.doprint(0);
            DEBUG.doprint(0);
            INFO.doprint(1);
            PROGRESS.doprint(1);
            WARNING.doprint(1);
            ERROR.doprint(1);
            displevel = 10000 + displevel%10000;        // for cnt #warnings
            INFO.print("SCREEN: \tverboseness: \t\t\tPROGRESS");
            }
          else if (!strcmp(keyword,"DEBUG"))
            {
            TRACE.doprint(0);
            DEBUG.doprint(1);
            INFO.doprint(1);
            PROGRESS.doprint(1);
            WARNING.doprint(1);
            ERROR.doprint(1);
            displevel = 30000 + displevel%10000;        // for cnt #warnings
            INFO.print("SCREEN: \tverboseness: \t\t\tDEBUG");
            }
          else if (!strcmp(keyword,"TRACE"))
            {
            TRACE.doprint(1);
            DEBUG.doprint(1);
            INFO.doprint(1);
            PROGRESS.doprint(1);
            WARNING.doprint(1);
            ERROR.doprint(1);
            }
          else if (!strcmp(keyword,"WARNING"))
            {
            TRACE.doprint(0);
            DEBUG.doprint(0);
            INFO.doprint(0);
            PROGRESS.doprint(0);
            WARNING.doprint(1);
            ERROR.doprint(1);
            displevel = 0 + displevel%10000;            // for cnt #warnings;
            INFO.print("SCREEN: \tverboseness: \t\t\tWARNING");
            }
          else if (!strcmp(keyword,"ERROR"))
            {
            TRACE.doprint(0);
            DEBUG.doprint(0);
            INFO.doprint(0);
            PROGRESS.doprint(0);
            WARNING.doprint(0);
            ERROR.doprint(1);
            displevel = -100 + displevel%10000;         // for cnt #warnings
            INFO.print("SCREEN: \tverboseness: \t\t\tERROR");
            }
          else
            {
            TRACE.doprint(0);
            DEBUG.doprint(1);
            INFO.doprint(1);
            PROGRESS.doprint(1);
            WARNING.doprint(1);
            ERROR.doprint(1);
            WARNING << "SCREEN: line " << linecnt
                 << ": Argument " << keyword 
                 << " not recognized."
                 << " [error/warning/progress/info/debug] I used DEBUG.";
            WARNING.print();
            displevel = 30000 + displevel%10000;       // for cnt #warnings
            }
        } // switch
      } // SCREEN key

// **********************************************************************
    else if (!strcmp(keyword,"MEMORY"))         // available mem in MB for processing in buffers
      {
      switch (priormemory)
        {
        case true:
          WARNING << "MEMORY: line " << linecnt
               << ": ignored due to prior occurrence.";
          WARNING.print();
          break;

        default:
          priormemory = true;
          //generalinput.memory =  word[1] ;    // pass keyword
          keyword = word[1];
          char *pLast = NULL;
          //generalinput.memory = uint(strtod( keyword, &pLast )); // [MA] try strtoul( keyword, &pLast, 10 ) for uint 
          generalinput.memory = strtod( keyword, &pLast ); 
          if ( pLast == keyword ) // fail
          {
            ERROR << "memory argument: "  << keyword << " is not valid.";
            PRINT_ERROR(ERROR.get_str())
            throw(keyword_error);
          }
          writearg(generalinput.memory);
          if (generalinput.memory > real8(MEMORY_MAX))
            WARNING << "MEMORY: > " << MEMORY_MAX << " MB seems unlikely.";

          //generalinput.memory *= 1000000;                       // in B
          generalinput.memory *= 1e6;                       // convert Mb --> b
        } // switch
      } // MEMORY card

// **********************************************************************
    else if (!strcmp(keyword,"BATCH"))                  // overrides interactive mode
      {
      switch (priorbatch)
        {
        case true:
          WARNING << "BATCH: line: " << linecnt << ": "
               << "ignored due to prior occurrence.";
          WARNING.print();
          break;
        default:
          priorbatch = true;                            // flag for occurrence
         // keyword =  word[1] ;        // pass keyword                 // argument
          keyword = word[1];                  // pass next word       (the argument)
          writearg(keyword);
    toupper(keyword);
          if (!strcmp(keyword,"OFF"))
            generalinput.interactive = true;
          else if (!strcmp(keyword,"ON")    ||
                   !strncmp(keyword,"//",2) ||          // comment
                   !strncmp(keyword,"#",1)  ||          // comment
                   !(keyword[0] == '\0'))               // no keyword
            generalinput.interactive = false;
          else
            {
            generalinput.interactive = true;
            WARNING << "BATCH: line: " << linecnt << ": "
                 << "argument: " << keyword 
                 << " not recognized, interactive processing.";
            WARNING.print();
            }
        } // switch
      } // BATCH key

// **********************************************************************
    else if (!strcmp(keyword,"OVERWRITE")) 
      {
      switch (prioroverwrite)
        {
        case true:
          WARNING << "OVERWRITE: line: " << linecnt << ": "
               << "ignored due to prior occurrence.";
          WARNING.print();
          break;
        default:
          prioroverwrite = true;                                // flag for occurrence
          keyword = word[1];                  // pass next word       (the argument)
          writearg(keyword);
          toupper(keyword);
          if (!strcmp(keyword,"OFF")) 
            generalinput.overwrit = false;
          else if (!strcmp(keyword,"ON")    ||
                   !strncmp(keyword,"//",2) ||          // comment
                   !strncmp(keyword,"#",1)  ||          // comment
                   !(keyword[0] == '\0'))               // no keyword
            generalinput.overwrit = true;
          else
            {
            generalinput.overwrit=false;                    // don't overwrite files
            WARNING << "OVERWRITE: line " << linecnt 
                 << ": argument: " << keyword 
                 << " not recognized, existing files are not overwritten.";
            WARNING.print();
            }
        } // switch
      } // OVERWRITE key

// **********************************************************************
    else if (!strcmp(keyword,"LISTINPUT")) 
      {
      switch (priorlistinput)
        {
        case true:
          WARNING << "LISTINPUT: line: " << linecnt << ": "
               << "ignored due to prior occurrence.";
          WARNING.print();
          break;
        default:
          priorlistinput = true;                                // flag for occurrence
          keyword = word[1];                  // pass next word       (the argument)
          writearg(keyword);
          toupper(keyword);
          if (!strcmp(keyword,"OFF")) 
            listinput = false;
          else if (!strcmp(keyword,"ON")    ||
                   !strncmp(keyword,"//",2) ||          // comment
                   !strncmp(keyword,"#",1)  ||          // comment
                   !(keyword[0] == '\0'))               // no keyword
            listinput = true;
          else
            {
            listinput = true;                                   // default list input
            WARNING << "LISTINPUT: line " << linecnt 
                 << ": argument: " << keyword 
                 << " not recognized, input will be appended to logfile.";
            WARNING.print();
            }
        } // switch
      } // LISTINPUT key

// **********************************************************************
    else if (!strcmp(keyword,"ONLYPROCESS"))     // process only one step
      {
      // Read argument and set onlyprocess to value for filling
      //  the flag input array 'process[NUMPROCESSES]' after reading reset input
      //  to avoid interference with PROCESS cards (ONLYPROCESS overrides)
      //
      if (onlyprocess == -1)                             // check multiple occurrences
        {
          keyword = word[1];                  // pass next word       (the argument)
        writearg(keyword);
    toupper(keyword);
        if      (!strcmp(keyword,"M_READFILES"))
          onlyprocess=pr_m_readfiles;
        else if (!strcmp(keyword,"M_CROP"))
          onlyprocess=pr_m_crop;
//____RaffaeleNutricato START MODIFICATION SECTION 5
        else if (!strcmp(keyword,"M_OVS"))
          onlyprocess=pr_m_oversample;
//____RaffaeleNutricato END MODIFICATION SECTION 5
        else if (!strcmp(keyword,"M_PORBITS"))
          onlyprocess=pr_m_porbits;
        else if (!strcmp(keyword,"M_MORBITS")) // [HB]
          onlyprocess=pr_m_morbits;
        else if (!strcmp(keyword,"M_SIMAMP")) // [MA]
          onlyprocess=pr_m_simamp;
        else if (!strcmp(keyword,"M_TIMING"))
          onlyprocess=pr_m_mtiming;
        else if (!strcmp(keyword,"M_FILTAZI"))
          onlyprocess=pr_m_filtazi;
        else if (!strcmp(keyword,"FILTRANGE"))
          {
          onlyprocess=pr_m_filtrange;
          onlyprocess=pr_s_filtrange;
          }
        else if (!strcmp(keyword,"M_EXTRA"))
          onlyprocess=pr_m_EXTRA;

        else if (!strcmp(keyword,"S_READFILES"))
          onlyprocess=pr_s_readfiles;
        else if (!strcmp(keyword,"S_MORBITS")) // [HB]
          onlyprocess=pr_s_morbits;
        else if (!strcmp(keyword,"S_CROP"))
          onlyprocess=pr_s_crop;
//____RaffaeleNutricato START MODIFICATION SECTION 6
        else if (!strcmp(keyword,"S_OVS"))
          onlyprocess=pr_s_oversample;
//____RaffaeleNutricato END MODIFICATION SECTION 6
        else if (!strcmp(keyword,"S_PORBITS"))
          onlyprocess=pr_s_porbits;
        else if (!strcmp(keyword,"S_FILTAZI"))
          onlyprocess=pr_s_filtazi;
        else if (!strcmp(keyword,"RESAMPLE"))
          onlyprocess=pr_s_resample;
        else if (!strcmp(keyword,"S_EXTRA"))
          onlyprocess=pr_s_EXTRA;

        else if (!strcmp(keyword,"COARSEORB"))
          onlyprocess=pr_i_coarse;
        else if (!strcmp(keyword,"COARSECORR"))
          onlyprocess=pr_i_coarse2;
        else if (!strcmp(keyword,"FINE"))
          onlyprocess=pr_i_fine;                         // see methodselector
        else if (!strcmp(keyword,"RELTIMING")) // [FvL]
          onlyprocess=pr_i_timing;
        else if (!strcmp(keyword,"DEMASSIST")) // [FvL]
          onlyprocess=pr_i_demassist;
        else if (!strcmp(keyword,"COREGPM"))
          onlyprocess=pr_i_coregpm;
        else if (!strcmp(keyword,"INTERFERO"))
          onlyprocess=pr_i_interfero;
        else if (!strcmp(keyword,"COHERENCE"))
          onlyprocess=pr_i_coherence;
        else if (!strcmp(keyword,"FILTPHASE"))
          onlyprocess=pr_i_filtphase;
        else if (!strcmp(keyword,"COMPREFPHA"))
          onlyprocess=pr_i_comprefpha;
        else if (!strcmp(keyword,"SUBTRREFPHA"))
          onlyprocess=pr_i_subtrrefpha;
        else if (!strcmp(keyword,"COMPREFDEM"))
          onlyprocess=pr_i_comprefdem;
        else if (!strcmp(keyword,"SUBTRREFDEM"))
          onlyprocess=pr_i_subtrrefdem;
        else if (!strcmp(keyword,"UNWRAP"))
          onlyprocess=pr_i_unwrap;                       // see methodselector
        else if (!strcmp(keyword,"ESTORBITS")) // [HB]
          onlyprocess=pr_i_estorbits;         
        else if (!strcmp(keyword,"SLANT2H"))
          onlyprocess=pr_i_slant2h;
        else if (!strcmp(keyword,"GEOCODE"))
          onlyprocess=pr_i_geocoding;
        else if (!strcmp(keyword,"DINSAR"))
          onlyprocess=pr_i_dinsar;
        else if (!strcmp(keyword,"I_EXTRA2"))
          onlyprocess=pr_i_EXTRA2;
        else
          {
          ERROR << "ONLYPROCESS: line " << linecnt << ": Argument " 
                   << keyword << " not recognized.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
          }
        INFO << "ONLYPROCESS: \tonly processing step: \t\t" << keyword;
        INFO.print();
        }
      else
        {
        WARNING << "ONLYPROCESS: more than one occurrence of card, ignored line: "
                   << linecnt << ".";
        WARNING.print();

        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"PROCESS"))                 // which routine to run
      {
      if (onlyprocess+1)                                 // initialized to -1;
        {
        WARNING << "PROCESS card on line " << linecnt 
             << " ignored due to presence of ONLYPROCESS card.";
        WARNING.print();
        }
      else
        {
          keyword = word[1];                  // pass next word       (the argument)
        writearg(keyword);
    toupper(keyword);
        if (!strcmp(keyword,"M_READFILES"))
          generalinput.process[pr_m_readfiles] = 1;
        else if (!strcmp(keyword,"M_CROP"))
          generalinput.process[pr_m_crop]   = 1;
//____RaffaeleNutricato START MODIFICATION SECTION 7
        else if (!strcmp(keyword,"M_OVS"))
          generalinput.process[pr_m_oversample]   = 1; 
//____RaffaeleNutricato END MODIFICATION SECTION 7
        else if (!strcmp(keyword,"M_PORBITS"))
          generalinput.process[pr_m_porbits]   = 1;
        else if (!strcmp(keyword,"M_MORBITS")) // [HB]
          generalinput.process[pr_m_morbits]   = 1;
        else if (!strcmp(keyword,"M_SIMAMP"))         // [MA]
          generalinput.process[pr_m_simamp]    = 1;
        else if (!strcmp(keyword,"M_TIMING"))
          generalinput.process[pr_m_mtiming]   = 1;
        else if (!strcmp(keyword,"M_FILTAZI"))
          generalinput.process[pr_m_filtazi]   = 1;
        else if (!strcmp(keyword,"FILTRANGE"))
          {
          generalinput.process[pr_m_filtrange] = 1;     // use for s_ as well
          generalinput.process[pr_s_filtrange] = 1;
          }
        else if (!strcmp(keyword,"M_EXTRA"))
          generalinput.process[pr_m_EXTRA] = 1;

        else if (!strcmp(keyword,"S_READFILES"))
          generalinput.process[pr_s_readfiles] = 1;
        else if (!strcmp(keyword,"S_MORBITS")) // [HB]
          generalinput.process[pr_s_morbits]   = 1;
        else if (!strcmp(keyword,"S_CROP"))
          generalinput.process[pr_s_crop]   = 1;
//____RaffaeleNutricato START MODIFICATION SECTION 8
        else if (!strcmp(keyword,"S_OVS"))
          generalinput.process[pr_s_oversample]   = 1;
//____RaffaeleNutricato END MODIFICATION SECTION 8
        else if (!strcmp(keyword,"S_PORBITS"))
          generalinput.process[pr_s_porbits]   = 1;
        else if (!strcmp(keyword,"S_FILTAZI"))
          generalinput.process[pr_s_filtazi]   = 1;
        else if (!strcmp(keyword,"RESAMPLE"))
          generalinput.process[pr_s_resample]  = 1;
        else if (!strcmp(keyword,"S_EXTRA"))
          generalinput.process[pr_s_EXTRA] = 1;

        else if (!strcmp(keyword,"COARSEORB"))
          generalinput.process[pr_i_coarse]    = 1;
        else if (!strcmp(keyword,"COARSECORR"))
          generalinput.process[pr_i_coarse2]   = 1;
        else if (!strcmp(keyword,"FINE"))
          generalinput.process[pr_i_fine]      = 1;
        else if (!strcmp(keyword,"RELTIMING"))         // [FvL]
          generalinput.process[pr_i_timing]    = 1;
        else if (!strcmp(keyword,"DEMASSIST"))         // [FvL]
          generalinput.process[pr_i_demassist]    = 1;
        else if (!strcmp(keyword,"COREGPM"))
          generalinput.process[pr_i_coregpm]   = 1;
        else if (!strcmp(keyword,"COMPREFPHA"))
          generalinput.process[pr_i_comprefpha] = 1;
        else if (!strcmp(keyword,"SUBTRREFPHA"))
          generalinput.process[pr_i_subtrrefpha] = 1;
        else if (!strcmp(keyword,"COMPREFDEM"))
          generalinput.process[pr_i_comprefdem] = 1;
        else if (!strcmp(keyword,"SUBTRREFDEM"))
          generalinput.process[pr_i_subtrrefdem] = 1;
        else if (!strcmp(keyword,"INTERFERO"))
          generalinput.process[pr_i_interfero] = 1;
        else if (!strcmp(keyword,"COHERENCE"))
          generalinput.process[pr_i_coherence] = 1;
        else if (!strcmp(keyword,"FILTPHASE"))
          generalinput.process[pr_i_filtphase] = 1;
        else if (!strcmp(keyword,"UNWRAP"))
          generalinput.process[pr_i_unwrap]    = 1;
        else if (!strcmp(keyword,"ESTORBITS"))    // [HB]
          generalinput.process[pr_i_estorbits]    = 1;
        else if (!strcmp(keyword,"SLANT2H"))
          generalinput.process[pr_i_slant2h]   = 1;
        else if (!strcmp(keyword,"GEOCODE"))
          generalinput.process[pr_i_geocoding] = 1;
        else if (!strcmp(keyword,"DINSAR"))
          generalinput.process[pr_i_dinsar] = 1;

        else if (!strcmp(keyword,"I_EXTRA2"))
          generalinput.process[pr_i_EXTRA2] = 1;
        else
          {
          ERROR << "PROCESS: line " << linecnt 
               << ": Argument " << keyword 
               << " not recognized.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
          }
        INFO << "PROCESS: \tI will process step: \t\t" << keyword;
        INFO.print();
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"ELLIPSOID"))            // ref. system
      {                                               //  inputoptionsfile
      ellipsoid = true;                               // use below 
      keyword = word[1];
      char *keyword2 = word[2]; 
      writearg(keyword);
      toupper(keyword);
      writearg(keyword2);
      if (!strcmp(keyword,"WGS84"))
        {
        ellipsinput   = WGS84;// default
        }
      else if (!strcmp(keyword,"GRS80"))
        {
        WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
              input_ell GRS80(6378137.0,6356752.3);
              GRS80.set_name("GRS80");
              ellipsinput = GRS80;// copy
        }
      else if (!strcmp(keyword,"BESSEL"))
        {
        WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
              input_ell BESSEL(6377397.155,6356078.963);
              BESSEL.set_name("BESSEL");
              ellipsinput = BESSEL;// copy
        }
      else if (isdigit(keyword2[0]))                      // likely to be a,b
        {
        WARNING.print("ELLIPS: not ok, sat. ephemerides should be in this system.");
        input_ell ELL_USER_DEFINED(atof(keyword),atof(keyword2));// a,b
              ELL_USER_DEFINED.set_name("user_defined");
        ellipsinput = ELL_USER_DEFINED;// copy
        if (ellipsinput.a<ellipsinput.b || ellipsinput.b<EPS)
          {
          ERROR << "ELLIPSOID keyword (real8A real8B): B==0 or A<B: "
               << ellipsinput.a << "<" << ellipsinput.b 
               << " at line " << linecnt << ".";
                PRINT_ERROR(ERROR.get_str())
                throw(keyword_error);
          }
        }
      else
        {
        PRINT_ERROR("unknown argument for ellipsoid card.")
        throw(keyword_error);
        }
      INFO << "ELLIPSOID: \tsemimajor=" << ellipsinput.a
           << ", semiminor=" << ellipsinput.b << "; line " << linecnt << ".";
      INFO.print();
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_RESFILE"))            // resultfile filename
      {                                               
      //generalinput.m_resfile =  word[1] ;     // pass keyword
      keyword = word[1];
      // strcpy(generalinput.m_resfile, keyword);
      strcpy(generalinput.m_resfile, keyword);
      writearg(generalinput.m_resfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_RESFILE"))            // resultfile filename
      {                                               
      keyword = word[1];
      strcpy(generalinput.s_resfile, keyword);
      writearg(generalinput.s_resfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"LOGFILE"))              // logfile filename
      {                                               
      keyword = word[1];
      strcpy(generalinput.logfile, keyword);
      writearg(generalinput.logfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"I_RESFILE"))            // interferogram.out
      {                                               
      keyword = word[1];
      strcpy(generalinput.i_resfile, keyword);
      writearg(generalinput.i_resfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"ORB_INTERP"))           // orbit
      {                                               
      //keyword =  word[1] ;    // pass keyword                       // argument
      keyword = word[1]; 
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"POLYFIT"))
        {
        INFO.print("ORB_INTERP:  polynomial fit for interpolation");
        generalinput.orb_interp = ORB_DEFAULT;// depends on number of points
        // ___ Check second optional argument with degree ___
        // keyword =  word[1] ;  // pass keyword
        char *keyword2= word[2];
        int32 degree = atoi(keyword2); 
        if (degree > 0)// atoi returns 0 if not convertible
          {
          generalinput.orb_interp = degree;
          INFO << "ORB_INTERP:  second argument read: degree = " << degree;
          INFO.print();
          }
        }
      else if (!strcmp(keyword,"SPLINE"))
        {
        INFO.print("ORB_INTERP:  natural cubic splines used fit interpolation");
        generalinput.orb_interp = ORB_SPLINE;// natural cubic splines
        }
      else
        {
        WARNING.print("argument ORB_INTERP not recognized, using polyfit");
        generalinput.orb_interp = ORB_DEFAULT;// depends on number of points
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"ORB_PRM"))           // orbit parameters
      {                                               
      keyword = word[1];    // pass keyword                       // argument 
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"POSITION") || !strcmp(keyword,"POS") )
        {
        INFO.print("ORB_PRM:  identified, using position state vectors for orbit interpolation");
        generalinput.orb_prm = ORB_PRM_POS;// user positional state vectors only
        }
      else if (!strcmp(keyword,"POSVEL"))
        {
        INFO.print("ORB_PRM:  identified, using positions and velocities for orbit interpolation");
        generalinput.orb_prm = ORB_PRM_VEL;// include velocity vectors too
        }
      else
        {
        WARNING.print("argument ORB_PRM not identified, using position state vectors only (this ignores any velocity column.)");
        generalinput.orb_prm = ORB_PRM_POS;// default positional.
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"DUMPBASELINE"))         // eval on grid
      {                                               
//      generalinput.dumpbaselineL >> generalinput.dumpbaselineP =  word[1] ;   // pass keyword
          keyword = word[1];
          char *keyword2= word[2]; // in local scope
          generalinput.dumpbaselineL = atoi( keyword ) ;
          generalinput.dumpbaselineP = atoi( keyword2 );
          writearg(generalinput.dumpbaselineL);
          writearg(generalinput.dumpbaselineP);
      if (generalinput.dumpbaselineL==0 || generalinput.dumpbaselineP==0)
        {
        ERROR << "DUMPBASELINE: " << generalinput.dumpbaselineL
             << " " << generalinput.dumpbaselineP << " line: "
             << linecnt << ": ==0.\n";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"PREVIEW"))        // system call to cpxfiddle to get SUNraster
      {
      //keyword =  word[1] ;    // pass keyword                 // argument
      keyword = word[1];                        // argument
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF"))
        {
        generalinput.preview = 0;
        INFO.print("PREVIEW: \tOFF: generation of SUNraster files disabled.");
        }
      else if (!strcmp(keyword,"XV"))
        {
        generalinput.preview = 2;
        INFO.print("PREVIEW: \tON: generation of SUNraster files enabled + XV system call.");
        }
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        {
        generalinput.preview = 1;
        INFO.print("PREVIEW: \tON: generation of SUNraster files enabled.");
        }
      else
        {
        generalinput.preview = 0;
        WARNING << "PREVIEW: line: " << linecnt << ": "
             << "argument: " << keyword 
             << " not recognized, no preview generated.";
        WARNING.print();
        }
      } // PREVIEW key

// **********************************************************************
    else if (!strcmp(keyword,"HEIGHT"))               // mean height or for CROP
      {
      //generalinput.terrain_height =  word[1] ;        // pass keyword
          keyword = word[1];
          char *pLast = NULL;
          generalinput.terrain_height = strtod( keyword, &pLast );
          if ( pLast == keyword ) // fail
          {
            ERROR << "Height argument: "  << keyword << " is not valid.";
            PRINT_ERROR(ERROR.get_str())
            throw(keyword_error);
          }
      writearg(generalinput.terrain_height);
      }

// **********************************************************************
    else if (!strcmp(keyword,"TIEPOINT"))          // lat/lon/hei dec.degrees
      {

        keyword = word[1];
        char *pLast, *pLast2, *pLast3 = NULL;
        generalinput.tiepoint.x = strtod( keyword, &pLast  );
        generalinput.tiepoint.y = strtod( word[2], &pLast2 ); // 2nd arg
                    generalinput.tiepoint.z = strtod( word[3], &pLast3 ); // 3rd arg
        if ( pLast == keyword || pLast2 == word[2] || pLast3 == word[3] ) // fail
        {
          ERROR << "Tiepoints: "  << keyword 
          << " " << word[2] 
          << " " << word[3] << " are not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
       writearg(generalinput.tiepoint.x);
       writearg(generalinput.tiepoint.y);
       writearg(generalinput.tiepoint.z);
      }

// **********************************************************************
    else if (!strcmp(keyword,"STOP"))                 // STOP interpreting input
      {
      INFO.print("STOP:   \tEncountered.");
      DEBUG << "STOP card encountered at line "
           << linecnt << endl;
      DEBUG.print();
      continuereading = false;                        // break while loop
      }

// *******************************************************************
// *** ?_READFILES
// *******************************************************************
    else if (!strcmp(keyword,"M_IN_METHOD"))         // ERS or ASAR ENVISAT
      {                                               
      //filename =  word[1] ;   // pass keyword
      keyword = word[1];
      writearg(keyword);
      toupper(keyword);

      if (!strcmp(keyword,"ERS"))
        m_readfilesinput.sensor_id=SLC_ERS;          // default ers
      else if (!strcmp(keyword,"ERS-1"))             // ers
        m_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS1"))              // ers
        m_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS-2"))             // ers
        m_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS2"))              // ers
        m_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS_N1"))           // ers in envisat file format [MA]; treat as ASAR
        //m_readfilesinput.sensor_id=SLC_ASAR;
        m_readfilesinput.sensor_id=SLC_ERS_N1;
      else if (!strcmp(keyword,"N1"))                // envisat
        m_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ASAR"))              // envisat
        m_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ENVISAT"))           // envisat
        m_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ASAR_AP_HH"))        // envisat Alternating Pol. HH
        m_readfilesinput.sensor_id=SLC_ASAR_AP_HH;
      else if (!strcmp(keyword,"ASAR_AP_VV"))        // envisat Alternating Pol. VV
        m_readfilesinput.sensor_id=SLC_ASAR_AP_VV;
      else if (!strcmp(keyword,"ATLANTIS"))          // radarsat ceos reader
        m_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RSAT"))              // radarsat
        m_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RADARSAT"))          // radarsat
        m_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RADARSAT-1"))        // radarsat
        m_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"JERS"))              // jers
        m_readfilesinput.sensor_id=SLC_JERS;
      else if (!strcmp(keyword,"J-ERS"))             // jers
        m_readfilesinput.sensor_id=SLC_JERS;
      else if (!strcmp(keyword,"ALOS"))              // alos FBS  [DON] 
        m_readfilesinput.sensor_id=SLC_ALOS;
      else if (!strcmp(keyword,"TSX"))               // TSX    [PM]
        m_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"TERRASARX"))         // TSX
        m_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"TERRASAR-X"))        // TSX
        m_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"RS2"))               // RS2  [MA]
        m_readfilesinput.sensor_id=SLC_RS2;
      else if (!strcmp(keyword,"RADARSAT-2"))        // RS2
        m_readfilesinput.sensor_id=SLC_RS2;
      else if (!strcmp(keyword,"CSK"))               // CSK  [PD]
        m_readfilesinput.sensor_id=SLC_CSK;
      else if (!strcmp(keyword,"CSX"))               // CSK  [PD]
        m_readfilesinput.sensor_id=SLC_CSK;
      else if (!strcmp(keyword,"GAMMA"))	     // GAMMA [BO]
        m_readfilesinput.sensor_id=SLC_GAMMA;	     
      else
        {
        ERROR << "M_IN_METHOD: method " <<  keyword
             << " not known for reading input files on line "
             << linecnt << ".";
            PRINT_ERROR(ERROR.get_str())
            throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_IN_VOL"))            // volumefile filename
      {                                               
      //m_readfilesinput.volfile =  word[1] ;   // pass keyword
      strcpy(m_readfilesinput.volfile, word[1]);          // pass keyword
      writearg(m_readfilesinput.volfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_IN_LEA"))            // leaderfile filename
      {                                               
      strcpy(m_readfilesinput.leaderfile, word[1]);          // pass keyword
      writearg(m_readfilesinput.leaderfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_IN_NULL"))           // nullfile filename
      {                                               
      strcpy(m_readfilesinput.nullfile, word[1]);          // pass keyword
      writearg(m_readfilesinput.nullfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_IN_DAT"))            // datafile filename
      {                                               
      strcpy(m_readfilesinput.datfile, word[1]);          // pass keyword
      writearg(m_readfilesinput.datfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_RG_T_ERROR"))            // master timing error
      {                                               
       keyword = word[1];
       char *pLast = NULL;
       m_readfilesinput.rg_timing_error= strtod( keyword, &pLast );
       if ( pLast == keyword ) // fail
        {
          ERROR << "M_RG_T_ERROR argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(m_readfilesinput.rg_timing_error);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_AZ_T_ERROR"))            // master timing error
      {                                               
       keyword = word[1];
       char *pLast = NULL;
       m_readfilesinput.az_timing_error = strtod( keyword, &pLast );
       if ( pLast == keyword ) // fail
        {
          ERROR << "M_AZ_T_ERROR argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(m_readfilesinput.az_timing_error);
      }

// *******************************************************************
    else if (!strcmp(keyword,"S_IN_METHOD"))         // ERS or ASAR ENVISAT
      {                                               
      keyword = word[1];
      writearg(keyword);
      toupper(keyword);

      if (!strcmp(keyword,"ERS"))
        s_readfilesinput.sensor_id=SLC_ERS;         // default
      else if (!strcmp(keyword,"ERS-1"))             // ers
        s_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS1"))              // ers
        s_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS-2"))             // ers
        s_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS2"))              // ers
        s_readfilesinput.sensor_id=SLC_ERS;
      else if (!strcmp(keyword,"ERS_N1"))           // ers in envisat file format [MA]; treat as ASAR
        s_readfilesinput.sensor_id=SLC_ERS_N1;
      else if (!strcmp(keyword,"N1"))                // envisat
        s_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ASAR"))              // envisat
        s_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ENVISAT"))           // envisat
        s_readfilesinput.sensor_id=SLC_ASAR;
      else if (!strcmp(keyword,"ATLANTIS"))          // radarsat ceos reader
        s_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RSAT"))              // radarsat
        s_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RADARSAT"))          // radarsat
        s_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"RADARSAT-1"))        // radarsat
        s_readfilesinput.sensor_id=SLC_RSAT;
      else if (!strcmp(keyword,"JERS"))              // jers
        s_readfilesinput.sensor_id=SLC_JERS;
      else if (!strcmp(keyword,"J-ERS"))             // jers
        s_readfilesinput.sensor_id=SLC_JERS;
      else if (!strcmp(keyword,"ALOS"))              // [DON]
        s_readfilesinput.sensor_id=SLC_ALOS;
      else if (!strcmp(keyword,"TSX"))               // TSX [PM]
        s_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"TERRASARX"))         // TSX
        s_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"TERRASAR-X"))        // TSX
        s_readfilesinput.sensor_id=SLC_TSX;
      else if (!strcmp(keyword,"RS2"))               // RS2  [MA]
        s_readfilesinput.sensor_id=SLC_RS2;
      else if (!strcmp(keyword,"RADARSAT-2"))        // RS2
        s_readfilesinput.sensor_id=SLC_RS2;
      else if (!strcmp(keyword,"CSK"))               // CSK  [PD]
        s_readfilesinput.sensor_id=SLC_CSK;
      else if (!strcmp(keyword,"CSX"))               // CSK  [PD]
        s_readfilesinput.sensor_id=SLC_CSK;
      else if (!strcmp(keyword,"GAMMA"))	     // GAMMA [BO]
        s_readfilesinput.sensor_id=SLC_GAMMA;	     
      else
        {
        ERROR << "S_IN_METHOD: method " <<  keyword
             << " not known for reading input files on line "
             << linecnt << ".";
            PRINT_ERROR(ERROR.get_str())
            throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_IN_VOL"))            // volumefile filename
      {                                               
      strcpy(s_readfilesinput.volfile, word[1]);          // pass keyword
      writearg(s_readfilesinput.volfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_IN_LEA"))            // leaderfile filename
      {                                               
      strcpy(s_readfilesinput.leaderfile, word[1]);          // pass keyword
      writearg(s_readfilesinput.leaderfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_IN_NULL"))           // nullfile filename
      {                                               
      strcpy(s_readfilesinput.nullfile, word[1]);          // pass keyword
      writearg(s_readfilesinput.nullfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_IN_DAT"))            // datfile filename
      {                                               
      strcpy(s_readfilesinput.datfile, word[1]);          // pass keyword
      writearg(s_readfilesinput.datfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_RG_T_ERROR"))            // slave timing error
      {                                               
       keyword = word[1];
       char *pLast = NULL;
       s_readfilesinput.rg_timing_error = strtod( keyword, &pLast );
       if ( pLast == keyword ) // fail
        {
          ERROR << "S_RG_T_ERROR argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(s_readfilesinput.rg_timing_error);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_AZ_T_ERROR"))            // slave timing error
      {                                               
       keyword = word[1];
       char *pLast = NULL;
       s_readfilesinput.az_timing_error = strtod( keyword, &pLast );
       if ( pLast == keyword ) // fail
        {
          ERROR << "S_AZ_T_ERROR argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(s_readfilesinput.az_timing_error);
      }


// **********************************************************************
// *** ?_PORBITS
// **********************************************************************
    else if (!strcmp(keyword,"M_ORBDIR"))             // orbitfile filename
      {                                               
      if (specified(porbitsinput.m_orbdir))
               WARNING.print("Prior occurrence of M_ORBDIR ignored.");
      strcpy(porbitsinput.m_orbdir,  word[1] );         // pass keyword
      writearg(porbitsinput.m_orbdir);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_ORB_INTERVAL") ||      // time interval in sec
                   !strcmp(keyword,"S_ORB_INTERVAL")) 
      {                                               
       keyword = word[1];         // pass keyword
       char *pLast = NULL;
       porbitsinput.timeinterval = strtol( keyword, &pLast , BASE10); // int32
       if ( pLast == keyword ) // fail
        {
          ERROR << "[M|S]_ORB_INTERVAL argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(porbitsinput.timeinterval);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_ORB_EXTRATIME") ||    // time before 1st line in sec
             !strcmp(keyword,"S_ORB_EXTRATIME")) 
      {                                               
       keyword = word[1];         // pass keyword
       char *pLast = NULL;
       porbitsinput.timebefore = strtol( keyword, &pLast , BASE10); // int32 
       if ( pLast == keyword ) // fail
        {
          ERROR << "[M|S]_ORB_EXTRATIME argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(porbitsinput.timebefore);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_ORB_DUMP"))           // dump orbit with dt
      {                                               
       keyword = word[1];         // pass keyword
       char *pLast = NULL;
       porbitsinput.dumpmasterorbit = strtod( keyword, &pLast ); 
       if ( pLast == keyword ) // fail
        {
          ERROR << "M_ORB_DUMP argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(porbitsinput.dumpmasterorbit);
      INFO.print("dumping master orbit to ascii file: masterorbit.dat");
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_ORBDIR"))             // orbitfile filename
      {                                               
      if (specified(porbitsinput.s_orbdir))
              WARNING.print("Prior occurrence of S_ORBDIR ignored.");
      strcpy(porbitsinput.s_orbdir,  word[1] );         // pass keyword
      writearg(porbitsinput.s_orbdir);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_ORB_DUMP"))           // dump orbit with dt
      {                                               
       keyword = word[1];         // pass keyword
       char *pLast = NULL;
       porbitsinput.dumpslaveorbit = strtod( keyword, &pLast ); 
       if ( pLast == keyword ) // fail
        {
          ERROR << "S_ORB_DUMP argument: "  << keyword << " is not valid.";
          PRINT_ERROR(ERROR.get_str())
          throw(keyword_error);
        }
      writearg(porbitsinput.dumpslaveorbit);
      INFO.print("dumping slave orbit to ascii file: slaveorbit.dat");
      }


// ____ start added by HB ____
// *******************************************************************
// *** ?_MORBITS
// *******************************************************************
    else if (!strncmp(keyword,"M_MO_D",6))
      {
	char *keywordpart = &keyword[6];
        int16 component;
 	if (!strncmp(keywordpart,"BH",2)) component = 0;
	else if (!strncmp(keywordpart,"BV",2)) component = 1;
	else
	  {
	    ERROR << "Unknown keyword: \"" << keyword << "\" at line: " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str()) throw(keyword_error);
	  }
	int16 degree = int(keywordpart[2])-48; // char2int
	morbitsinputmaster.coeff(degree,component) = atof(word[1]);
	writearg(morbitsinputmaster.coeff(degree,component));
      }

// *******************************************************************
    else if (!strcmp(keyword,"M_MO_REFORBIT"))
      {
        strcpy(morbitsinputmaster.reforbitfile, word[1]);
	writearg(morbitsinputmaster.reforbitfile);
      }

// *******************************************************************
    else if (!strncmp(keyword,"S_MO_D",6))
      {
	char *keywordpart = &keyword[6];
        int16 component;
 	if (!strncmp(keywordpart,"BH",2)) component = 0;
	else if (!strncmp(keywordpart,"BV",2)) component = 1;
	else
	  {
	    ERROR << "Unknown keyword: \"" << keyword << "\" at line: " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str()) throw(keyword_error);
	  }
	int16 degree = int(keywordpart[2])-48; // char2int
	morbitsinputslave.coeff(degree,component) = atof(word[1]);
	writearg(morbitsinputslave.coeff(degree,component));
      }

// *******************************************************************
    else if (!strcmp(keyword,"S_MO_REFORBIT"))
      {
        strcpy(morbitsinputslave.reforbitfile, word[1]);
	writearg(morbitsinputslave.reforbitfile);
      }
// ____ end added by HB ____
    

// *******************************************************************
// *** ?_CROP
// *******************************************************************
    else if (!strcmp(keyword,"M_CROP_ID"))          // identifier of run
      {                                               
      strcpy(m_cropinput.idcrop,  word[1] );    // pass keyword
      writearg(m_cropinput.idcrop);
      }

// *******************************************************************   CROP
    else if (!strcmp(keyword,"S_CROP_ID"))          // identifier of run
      {                                               
      strcpy(s_cropinput.idcrop,  word[1] );    // pass keyword
      writearg(s_cropinput.idcrop);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_CROP_IN"))              // SLC input image filename
      {                                               
      strcpy(m_cropinput.filein1,  word[1] );   // pass keyword
      writearg(m_cropinput.filein1);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_CROP_IN"))              // SLC input image filename
      {                                               
      strcpy(s_cropinput.filein1,  word[1] );   // pass keyword
      writearg(s_cropinput.filein1);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_CROP_OUT"))             // SLC output image filename
      {                                               
      strcpy(m_cropinput.fileout1,  word[1] );  // pass keyword
      writearg(m_cropinput.fileout1);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_CROP_OUT"))             // SLC output image filename
      {                                               
      strcpy(s_cropinput.fileout1,  word[1] );  // pass keyword
      writearg(s_cropinput.fileout1);
      }

// **********************************************************************
    else if (!strcmp(keyword,"M_DBOW"))               // min. line coord.
      {                                               
       m_cropinput.dbow.linelo = atoi(word[1]);  // pass keywords
       m_cropinput.dbow.linehi = atoi(word[2]);
       m_cropinput.dbow.pixlo  = atoi(word[3]);
       m_cropinput.dbow.pixhi  = atoi(word[4]);

      writearg(m_cropinput.dbow.linelo);
      writearg(m_cropinput.dbow.linehi);
      writearg(m_cropinput.dbow.pixlo);
      writearg(m_cropinput.dbow.pixhi);

// ______initial check, later checked to #lines of image______
      if (m_cropinput.dbow.linelo <= 0 || 
          m_cropinput.dbow.pixlo <= 0 ||
          m_cropinput.dbow.linelo  > m_cropinput.dbow.linehi ||
          m_cropinput.dbow.pixlo > m_cropinput.dbow.pixhi)
        {
        ERROR << "code 300: Arguments of M_DBOW card on line " << linecnt
             << " missing. [DBOW  min_line  max_line  min_pixel  max_pixel].";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_DBOW"))               // min. line coord.
      {                                               
       s_cropinput.dbow.linelo = atoi(word[1]);  // pass keywords
       s_cropinput.dbow.linehi = atoi(word[2]);
       s_cropinput.dbow.pixlo  = atoi(word[3]);
       s_cropinput.dbow.pixhi  = atoi(word[4]);

      writearg(s_cropinput.dbow.linelo);
      writearg(s_cropinput.dbow.linehi);
      writearg(s_cropinput.dbow.pixlo);
      writearg(s_cropinput.dbow.pixhi);
      // ______initial check, later checked to #lines of image______
      if (s_cropinput.dbow.linelo <= 0 || 
          s_cropinput.dbow.pixlo <= 0 ||
          s_cropinput.dbow.linelo  > s_cropinput.dbow.linehi ||
          s_cropinput.dbow.pixlo > s_cropinput.dbow.pixhi)
        {
        ERROR << "code 300: Arguments of S_DBOW card on line " << linecnt
             << " missing. [DBOW  min_line  max_line  min_pixel  max_pixel].";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    // BK 15-Dec-2003
    else if (!strcmp(keyword,"M_DBOW_GEO")) // 1e6*lat_0 [deg], lon_0, height, width [pix]
      {                                               
      // use dbow to store tmp, compute later in crop if not zero
      real8 tmp_lat_0, tmp_lon_0, tmp_height, tmp_width;
      tmp_lat_0  =  atof(word[1]) ;     // pass keyword
      tmp_lon_0  =  atof(word[2]) ;     // pass keyword
      tmp_height =  atof(word[3]) ;     // pass keyword
       tmp_width =  atof(word[4]) ;     // pass keyword
      writearg(tmp_lat_0);
      writearg(tmp_lon_0);
      writearg(tmp_height);
      writearg(tmp_width);
      m_cropinput.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
      m_cropinput.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
      m_cropinput.dbow_geo.pixlo  = uint(tmp_height);
      m_cropinput.dbow_geo.pixhi  = uint(tmp_width);
      }

// **********************************************************************
    // BK 15-Dec-2003
    else if (!strcmp(keyword,"S_DBOW_GEO")) // 1e6*lat_0 [deg], lon_0, height, width [pix]
      {                                               
      // use dbow to store tmp, compute later in crop if not zero
      real8 tmp_lat_0, tmp_lon_0, tmp_height, tmp_width;
      tmp_lat_0  =  atof(word[1]) ;     // pass keyword
      tmp_lon_0  =  atof(word[2]) ;     // pass keyword
      tmp_height =  atof(word[3]) ;     // pass keyword
       tmp_width =  atof(word[4]) ;     // pass keyword
      writearg(tmp_lat_0);
      writearg(tmp_lon_0);
      writearg(tmp_height);
      writearg(tmp_width);
      s_cropinput.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
      s_cropinput.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
      s_cropinput.dbow_geo.pixlo  = uint(tmp_height);
      s_cropinput.dbow_geo.pixhi  = uint(tmp_width);
      }

//____RaffaeleNutricato START MODIFICATION SECTION 9
// *******************************************************************
// *** ?_OVS
// **********************************************************************
    else if (!strcmp(keyword,"M_OVS_OUT"))             // oversampled SLC output image filename
      {                                               
      strcpy(m_oversample.fileoutovs,  word[1] );       // pass keyword
      writearg(m_oversample.fileoutovs);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S_OVS_OUT"))             // oversampled SLC output image filename
      {                                               
      strcpy(s_oversample.fileoutovs,  word[1] );       // pass keyword
      writearg(s_oversample.fileoutovs);
      }

// *******************************************************************
    else if (!strcmp(keyword,"M_OVS_FACT_RNG"))               // Oversampling ratio in the master range direction.
      {                                                  
      m_oversample.OsrRange =  strtol(word[1],NULL, BASE10);      // pass keyword               
      writearg(m_oversample.OsrRange);                    
      if (m_oversample.OsrRange < 1)
              {
        PRINT_ERROR(
        "(M_OVS_FACT_RNG < 1) Oversampling ratio must be at least 1.");
              throw(keyword_error);
              }
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"S_OVS_FACT_RNG"))               // Oversampling ratio in the slave range direction.
      {                                                       
      s_oversample.OsrRange =  strtol(word[1],NULL, BASE10) ;  // pass keyword               
      writearg(s_oversample.OsrRange);                    
      if (s_oversample.OsrRange < 1)
        {
        PRINT_ERROR(
        "(S_OVS_FACT_RNG < 1) Oversampling ratio must be at least 1.");
        throw(keyword_error);
        }
      }                                                  

// *******************************************************************
    else if (!strcmp(keyword,"M_OVS_FACT_AZI"))               // Oversampling ratio in the master azimuth direction.
      {                                                  
      m_oversample.OsrAzimuth = strtol(word[1],NULL, BASE10) ;         // pass keyword               
      writearg(m_oversample.OsrAzimuth);                    
      if (m_oversample.OsrAzimuth < 1)
        {
        PRINT_ERROR(
        "(M_OVS_FACT_AZI < 1) Oversampling ratio must be at least 1.");
        throw(keyword_error);
        }
      if (m_oversample.OsrAzimuth > 2)
        {
        PRINT_ERROR(
        "(M_OVS_FACT_AZI > 2) Not implemented!");
        throw(keyword_error);
        }
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"S_OVS_FACT_AZI"))               // Oversampling ratio in the slave azimuth direction.
      {                                                       
      s_oversample.OsrAzimuth =  strtol(word[1],NULL, BASE10) ;        // pass keyword               
      writearg(s_oversample.OsrAzimuth);                    
      if (s_oversample.OsrAzimuth < 1)
        {
        PRINT_ERROR("(S_OVS_FACT_AZI < 1) Oversampling ratio must be at least 1.");
        throw(keyword_error);
        }
      if (s_oversample.OsrAzimuth > 2)
        {
        PRINT_ERROR(
        "(S_OVS_FACT_AZI > 2) Not implemented!");
        throw(keyword_error);
        }
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"M_OVS_KERNELSIZE"))     // Length of the interpolation kernel. 
      {                                                  
      m_oversample.FilterSize =  strtol(word[1],NULL, BASE10) ;        // pass keyword        
      writearg(m_oversample.FilterSize);             
      if (m_oversample.FilterSize < 2)
        {
         PRINT_ERROR("(M_OVS_KERNELSIZE < 2) Interpolation kernel length must be > 1.");
        throw(keyword_error);
        }
      if (m_oversample.FilterSize % 2)
        {
        PRINT_ERROR("(M_OVS_KERNELSIZE not even) Range Interpolation kernel length must be an even number.");
        throw(keyword_error);
        }
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"S_OVS_KERNELSIZE"))     // Length of the interpolation kernel. 
      {                                                  
      s_oversample.FilterSize =  strtol(word[1],NULL, BASE10) ;        // pass keyword        
      writearg(s_oversample.FilterSize);             
      if (s_oversample.FilterSize < 2)
        {
        PRINT_ERROR("code ???: (S_OVS_KERNELSIZE < 2) Interpolation kernel length must be > 1.");
        throw(keyword_error);
        }
      if (s_oversample.FilterSize % 2)
        {
        PRINT_ERROR(
        "code ???: (S_OVS_KERNELSIZE not even) Interpolation kernel length must be an even number.");
        throw(keyword_error);
        }
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"M_OVS_OUT_FORMAT"))         // Output format [cr4] ci16, I suggest [cr4].
      {                                                  
      keyword =  word[1] ;      // pass keyword                           
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CR4"))
        m_oversample.oformatflag = FORMATCR4;         // default
      else if (!strcmp(keyword,"CI2"))
        m_oversample.oformatflag = FORMATCI2;
      else
        {  
        ERROR << "M_OVS_OUT_FORMAT: output format "
             <<  keyword
             << " not known for master range oversampling. line "
             << linecnt << ".";
          PRINT_ERROR(ERROR.get_str())      
        throw(keyword_error);
        }             
      }                                                  

// ********************************************************************** 
    else if (!strcmp(keyword,"S_OVS_OUT_FORMAT"))         // Output format [cr4] ci16, I suggest [cr4].
      {                                                  
      keyword =  word[1] ;      // pass keyword                           
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CR4"))
        s_oversample.oformatflag = FORMATCR4;         // default
      else if (!strcmp(keyword,"CI2"))
        s_oversample.oformatflag = FORMATCI2;
      else
        {  
        ERROR << "S_OVS_OUT_FORMAT: output format "
             <<  keyword
             << " not known for slave range oversampling. line "
             << linecnt << ".";
          PRINT_ERROR(ERROR.get_str()) 
        throw(keyword_error);
        }             
      }               
//____RaffaeleNutricato END MODIFICATION SECTION 9

// ____ start added by MA ____


// **********************************************************************
// *** SIMULATE AMPLITUDE FOR MASTER
// **********************************************************************

      else if (!strcmp(keyword,"SAM_IN_DEM") && 
                 (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )          // input file
                                                                                                 // [MA] read parm only if there is (ONLY)PROCESS flag
                                                                                                 // default onlyprocess is -1
      {
      strcpy(simampinput.firefdem,  word[1] );  // pass keyword
      writearg(simampinput.firefdem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_IN_FORMAT") && 
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )       //  format input file
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
        simampinput.iformatflag = FORMATR4;
      else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
        simampinput.iformatflag = FORMATI2;     // default
      else if (!strcmp(keyword,"I2_BIGENDIAN") || 
               !strcmp(keyword,"SHORT_BIGENDIAN"))
        simampinput.iformatflag = FORMATI2_BIGENDIAN;   // default
      else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
        simampinput.iformatflag = FORMATR8;
      else
        {
        ERROR << "SAM_IN_FORMAT: input format "
             <<  keyword
             << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line "
             << linecnt << ".";
              PRINT_ERROR(ERROR.get_str())
              throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_IN_SIZE") &&
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )         // nrow ncols (lat lon)
      {
      char *pLast1, *pLast2 = NULL;
      simampinput.demrows = strtoul(word[1], &pLast1, BASE10); 
      simampinput.demcols = strtoul(word[2], &pLast2, BASE10);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "SAM_IN_SIZE: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(simampinput.demrows);
      writearg(simampinput.demcols);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_IN_DELTA") && 
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )        // degrees delta lat lon
      {
      simampinput.demdeltalat  = atof(word[1]);           // pass keyword
      keyword                  =  word[2];              // pass keyword
      writearg(simampinput.demdeltalat);
      writearg(keyword);
      if (isdigit(keyword[0]) || keyword[0]=='.')    // likely to be 2 numbers
        simampinput.demdeltalon = atof(keyword);
      else // default same gridsize
        simampinput.demdeltalon = simampinput.demdeltalat;

      // ______ Store as radians ______
      simampinput.demdeltalat = deg2rad(simampinput.demdeltalat);
      simampinput.demdeltalon = deg2rad(simampinput.demdeltalon);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_IN_UL")  && 
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )           // upperleft coordinates
      {
      char *pLast1, *pLast2 = NULL;
      simampinput.demlatleftupper = strtod(word[1], &pLast1); 
      simampinput.demlonleftupper = strtod(word[2], &pLast2);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "SAM_IN_UL: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(simampinput.demlatleftupper);
      writearg(simampinput.demlonleftupper);
      simampinput.demlatleftupper = deg2rad(simampinput.demlatleftupper);
      simampinput.demlonleftupper = deg2rad(simampinput.demlonleftupper);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_IN_NODATA")  &&
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )       // flag for no data
      {
      char *pLast = NULL;
      simampinput.demnodata = strtod(word[1], &pLast); 
      if ( pLast == word[1]  ) // fails to convert to double.
       {
        ERROR << "SAM_IN_NODATA: "  << word[1] << " is not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(simampinput.demnodata);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_OUT_DEM")  &&
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )          // name of output file
      {
      strcpy(simampinput.fodem,  word[1] );     // pass keyword
      writearg(simampinput.fodem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SAM_OUT_FILE") &&  
               (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )          // name of output file
      {
      strcpy(simampinput.fosimamp,  word[1] );  // pass keyword
      writearg(simampinput.fosimamp);
      }

// // **********************************************************************
//     else if (!strcmp(keyword,"SAM_OUT_DEMI"))          // name of output file
//       {
//       strcpy(demassistinput.fodemi,  word[1] );      // pass keyword
//       writearg(demassistinput.fodemi);
//       }
// 
// // **********************************************************************
     else if (!strcmp(keyword,"SAM_OUT_DEM_LP") && 
                (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )        // name of output file  // MA enabled for SPT tracking.
       {
       strcpy(simampinput.fodemlp,  word[1] );         // pass keyword
       writearg(simampinput.fodemlp);
       }
     else if (!strcmp(keyword,"SAM_OUT_THETA_LP") && 
                (  generalinput.process[pr_m_simamp] || onlyprocess == pr_m_simamp ) )        // name of output file  // MA enabled for SPT tracking.
       {
       strcpy(simampinput.fothetalp,  word[1] );         // pass keyword
       writearg(simampinput.fothetalp);
       }


// **********************************************************************
// *** MASTER TIMING ERROR ESTIMATION using COREGISTRATION                             
// **********************************************************************
    else if (!strcmp(keyword,"MTE_METHOD"))            // method selector simamp coreg
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"MAGFFT"))
        {
         mtiminginput.method=cc_magfft;                // default MTE_magfft
        }
      else if (!strcmp(keyword,"MAGSPACE"))
        {
         mtiminginput.method=cc_magspace;              // MTE_magspace
        }
      else
        {
         ERROR << "MTE_METHOD: method " <<  keyword
               << " not known for simamp correlation coregistration on line "
               << linecnt << ".";
         PRINT_ERROR(ERROR.get_str())
               throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"MTE_NWIN"))              // #windows for simamp correlation
      {                                               
      char *pLast = NULL;
      mtiminginput.Nwin = strtoul(word[1], &pLast, BASE10); 
      if ( pLast == word[1]  ) // fails to convert to double.
       {
        ERROR << "MTE_NWIN: "  << word[1] << " is not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(mtiminginput.Nwin);
      if (mtiminginput.Nwin > 10000)
             {
        PRINT_ERROR("Too many windows requested (MTE_NWIN > 10000).")
        throw(keyword_error);
             }
      }

// **********************************************************************
    else if (!strcmp(keyword,"MTE_IN_POS"))            // file with #windows positions
      {                                               
      strcpy(mtiminginput.ifpositions,  word[1] );      // pass keyword
      writearg(mtiminginput.ifpositions);
      }

// **********************************************************************
    else if (!strcmp(keyword,"MTE_WINSIZE"))           // windowsize for simamp correlation
      {                                               
      char *pLast1, *pLast2 = NULL;
      mtiminginput.MasksizeL = strtoul(word[1], &pLast1, BASE10); 
      mtiminginput.MasksizeP = strtoul(word[2], &pLast2, BASE10);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "MTE_WINSIZE: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(mtiminginput.MasksizeL);
      writearg(mtiminginput.MasksizeP);
      if (mtiminginput.MasksizeL > 4096 || mtiminginput.MasksizeP > 4096)
              {
        PRINT_ERROR("Too large correlation window (MTE_WINSIZE > 4096).");
        throw(keyword_error);
              }
      }

// **********************************************************************
    else if (!strcmp(keyword,"MTE_ACC"))               // Searchdomain for correlation
      {                                               //  simamp correlation
      char *pLast1, *pLast2 = NULL;
      mtiminginput.AccL  = strtoul(word[1], &pLast1, BASE10); 
      mtiminginput.AccP  = strtoul(word[2], &pLast2, BASE10);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "MTE_ACC: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(mtiminginput.AccL);
      writearg(mtiminginput.AccP);
      if (mtiminginput.AccL > 1000 || mtiminginput.AccP > 1000)
              {
        PRINT_ERROR("Too large searchwindow (MTE_ACC > 1000).");
        throw(keyword_error);
              }
      if (mtiminginput.AccL == 0 || mtiminginput.AccP == 0)
              {
        PRINT_ERROR("Acc = 0 ?(MTE_ACC).");
        throw(keyword_error);
              }
      }

// **********************************************************************
    else if (!strcmp(keyword,"MTE_INITOFF"))            // Initial offset
      {
      keyword   =  word[1] ;    // pass keyword
      char *keyword2  =  word[2] ;      // pass keyword
      writearg(keyword);
      writearg(keyword2);
      if (isdigit(keyword2[0]) || keyword2[0]=='-') // likely to be 2 numbers
                                                    // BK19/1/00 thanks to m.goos for '-'
        {
        mtiminginput.initoffsetL = atoi(keyword);
        mtiminginput.initoffsetP = atoi(keyword2);
        }
      else
        {
        ERROR << "MTE_INITOFF: unknown input: "
             <<  keyword << ", " << keyword2
             << " on line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }


// ____ end added by MA ____


// **********************************************************************
// *** AZIMUTH FILTERING
// **********************************************************************
    else if (!strcmp(keyword,"AF_BLOCKSIZE"))
      {
      char *pLast = NULL;
      filtaziinput.fftlength = strtol(word[1], &pLast, BASE10); // int32 
      if ( pLast == word[1]  ) // fails to convert to double.
       {
        ERROR << "AF_BLOCKSIZE: "  << word[1] << " is not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(filtaziinput.fftlength);
      }

// **********************************************************************
    else if (!strcmp(keyword,"AF_HAMMING"))
      {
      char *pLast = NULL;
      filtaziinput.hammingalpha = strtod(word[1], &pLast); 
      if ( pLast == word[1]  ) // fails to convert to double.
       {
        ERROR << "AF_HAMMING: "  << word[1] << " is not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(filtaziinput.hammingalpha);
      }

// **********************************************************************
    else if (!strcmp(keyword,"AF_OVERLAP"))
      {
      char *pLast = NULL;
      filtaziinput.overlap = strtol(word[1], &pLast, BASE10); // int32 
      if ( pLast == word[1]  ) // fails to convert to double.
       {
        ERROR << "AF_OVERLAP: "  << word[1] << " is not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(filtaziinput.overlap);
      }

// **********************************************************************
    else if (!strcmp(keyword,"AF_OUT_MASTER"))
      {
      strcpy(filtaziinput.fomaster,  word[1] );         // pass keyword
      writearg(filtaziinput.fomaster);
      }

// **********************************************************************
    else if (!strcmp(keyword,"AF_OUT_SLAVE"))
      {
      strcpy(filtaziinput.foslave,  word[1] );  // pass keyword
      writearg(filtaziinput.foslave);
      }

// **********************************************************************
    else if (!strcmp(keyword,"AF_OUT_FORMAT"))          // output format
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CR4"))
        filtaziinput.oformatflag = FORMATCR4;         // default
      else if (!strcmp(keyword,"CI2"))
        filtaziinput.oformatflag = FORMATCI2;
      else
        {
        ERROR << "AF_OUT_FORMAT: output format "
             <<  keyword
             << " not known for azimuth filtering. line "
             << linecnt << ".";
              PRINT_ERROR(ERROR.get_str())
              throw(keyword_error);
        }
      }



// **********************************************************************
// *** COARSE CORR COREGISTRATION
// **********************************************************************
    else if (!strcmp(keyword,"CC_METHOD"))            // method selector coarse coreg
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"MAGFFT"))
        coarsecorrinput.method=cc_magfft;             // default
      else if (!strcmp(keyword,"MAGSPACE"))
        coarsecorrinput.method=cc_magspace;
      else
        {
        ERROR << "CC_METHOD: method " <<  keyword
             << " not known for coarse correlation coregistration on line "
             << linecnt << ".";
              PRINT_ERROR(ERROR.get_str())
              throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CC_NWIN"))              // #windows for coarse correlation
      {                                               
      coarsecorrinput.Nwin =  atoi(word[1]) ;   // pass keyword
      writearg(coarsecorrinput.Nwin);
      if (coarsecorrinput.Nwin > 10000)
        {
        PRINT_ERROR("Too many windows requested (CC_NWIN > 10000).")
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CC_IN_POS"))            // file with #windows positions
      {                                               
      strcpy(coarsecorrinput.ifpositions,  word[1] );   // pass keyword
      writearg(coarsecorrinput.ifpositions);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CC_WINSIZE"))           // windowsize for coarse correlation
      {                                               
      coarsecorrinput.MasksizeL  =  atoi(word[1]) ;     // pass keyword
      coarsecorrinput.MasksizeP  =  atoi(word[2]) ;     // pass keyword
      writearg(coarsecorrinput.MasksizeL);
      writearg(coarsecorrinput.MasksizeP);
      if (coarsecorrinput.MasksizeL > 10240 || coarsecorrinput.MasksizeP > 4096)   // [TODO] move MAX constant values to constant.hh
              {
        PRINT_ERROR("Too large correlation window (CC_WINSIZE(l,p) > [10240,4096] ).");
        throw(keyword_error);
              }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CC_ACC"))               // Searchdomain for correlation
      {                                               //  coarse correlation
      coarsecorrinput.AccL  =  atoi(word[1]) ;  // pass keyword
      coarsecorrinput.AccP  =  atoi(word[2]) ;  // pass keyword
      writearg(coarsecorrinput.AccL);
      writearg(coarsecorrinput.AccP);
      if (coarsecorrinput.AccL > 1000 || coarsecorrinput.AccP > 1000)
              {
        PRINT_ERROR("Too large searchwindow (CC_ACC > 1000).");
        throw(keyword_error);
              }
      if (coarsecorrinput.AccL == 0 || coarsecorrinput.AccP == 0)
              {
        PRINT_ERROR("Acc = 0 ?(CC_ACC).");
        throw(keyword_error);
              }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CC_INITOFF"))             // Initial offset
      {
      keyword  =  word[1] ;     // pass keyword
      char *keyword2 =  word[2] ;       // pass keyword
      writearg(keyword);
      writearg(keyword2);
      if (!strcmp(keyword,"ORBIT") ||                   // use value of precise orbits
          !strcmp(keyword,"orbit"))
        {
        coarsecorrinput.initoffsetL = NaN;              // flag used in main, see processor.cc
        coarsecorrinput.initoffsetP = NaN;              // flag used in main
        INFO << "CC_INITOFF: \tInitial offsets from COARSEORB: " << generalinput.i_resfile;
        INFO.print();
        }
      else if (isdigit(keyword2[0]) || keyword2[0]=='-') // likely to be 2 numbers
                                                                           // BK19/1/00 thanks to m.goos for '-'
        {
        //coarsecorrinput.initoffsetL = atof(keyword);
        //coarsecorrinput.initoffsetP = atof(keyword2);
        coarsecorrinput.initoffsetL = atoi(keyword);
        coarsecorrinput.initoffsetP = atoi(keyword2);
        }
      else
        {
        ERROR << "CC_INITOFF: unknown input: "
             <<  keyword << ", " << keyword2
             << " on line "
             << linecnt << ".";
              PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
// *** FINE COREGISTRATION
// **********************************************************************
    else if (!strcmp(keyword,"FC_NWIN"))              // #windows for fine correlation
      {                                               
      fineinput.Nwin =  atoi(word[1]) ;         // pass keyword
      writearg(fineinput.Nwin);
      if (fineinput.Nwin > 100000)              // [MA] change this to a higher value for SPT
              {
        PRINT_ERROR("Too many windows requested (FC_NWIN).")
        throw(keyword_error);
              }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_IN_POS"))            // file with #windows positions
      {                                               
      strcpy(fineinput.ifpositions,  word[1] );         // pass keyword
      writearg(fineinput.ifpositions);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_WINSIZE"))           // windowsize for fine correlation
      {                                               
      fineinput.MasksizeL  =  atoi(word[1]) ;   // pass keyword
      fineinput.MasksizeP  =  atoi(word[2]) ;   // pass keyword
      writearg(fineinput.MasksizeL);
      writearg(fineinput.MasksizeP);
      if (fineinput.MasksizeL > 1024 || fineinput.MasksizeP > 1024)
        {
        PRINT_ERROR("Too large correlation window (FC_WINSIZE).")
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_ACC"))               // Searchdomain for correlation
      {                                               //  fine correlation
      fineinput.AccL =  atoi(word[1]);  // pass keyword
      fineinput.AccP =  atoi(word[2]);  // pass keyword
      writearg(fineinput.AccL);
      writearg(fineinput.AccP);
      if (fineinput.AccL > 1000 || fineinput.AccP > 1000)
        {
        PRINT_ERROR("Too large searchwindow (FC_ACC).")
        throw(keyword_error);
        }
      if (fineinput.AccL == 0 || fineinput.AccP == 0)
        {
        PRINT_ERROR("Acc = 0 ?(FC_ACC).")
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_INITOFF"))             // Initial offset
      {
      keyword        =  word[1] ;       // pass keyword
      char *keyword2 =  word[2] ;       // pass keyword
      writearg(keyword);
      writearg(keyword2);
      if (!strcmp(keyword,"COARSECORR") 
       || !strcmp(keyword,"coarsecorr"))                // use value of coarse correlation
        {
        fineinput.initoffsetL = NaN;                    // flag used in main
        fineinput.initoffsetP = NaN;                    // flag used in main
        INFO << "FC_INITOFF: \tInitial offset from COARSECORR: "
        << generalinput.i_resfile;
        INFO.print();
        }
      else if (isdigit(keyword2[0]) || keyword2[0]=='-') // likely to be 2 numbers
                                                // BK19/1/00 thanks to M.Goos for '-'
        {
        fineinput.initoffsetL = atoi(keyword);
        fineinput.initoffsetP = atoi(keyword2);
        }
      else
        {
        ERROR << "FC_INITOFF: unknown input: "
             <<  keyword << ", " << keyword2
             << " on line "
             << linecnt << ".";
             PRINT_ERROR(ERROR.get_str())
             throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_METHOD"))            // method selector fine coreg
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CMPLXFFT")) 
       {
        //fineinput.method=fc_cmplxfft;                 // default
        PRINT_ERROR("CMPLXFFT not implemented in v1.0 of Doris.")
        throw(keyword_error);
       }
      else if (!strcmp(keyword,"CMPLXSPACE"))
       {
        //fineinput.method=fc_cmplxspace;
       PRINT_ERROR("CMPLXSPACE not implemented in v1.0 of Doris.")
       throw(keyword_error);
       }
      else if (!strcmp(keyword,"MAGFFT"))
        fineinput.method=fc_magfft;
      else if (!strcmp(keyword,"MAGSPACE"))
        fineinput.method=fc_magspace;
      else if (!strcmp(keyword,"OVERSAMPLE"))
        fineinput.method=fc_oversample;
      else if (!strcmp(keyword,"COHERENCE"))
        fineinput.method=fc_coherence;
      else if (!strcmp(keyword,"INTENSITY"))
        fineinput.method=fc_intensity;
      else
        {
        ERROR << "FC_METHOD: method " <<  keyword
             << " not known for fine coregistration on line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_OSFACTOR"))          // oversampling factor
      {                                               
      fineinput.osfactor =  atoi(word[1]) ;     // pass keyword
      writearg(fineinput.osfactor);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_PLOT"))              // plotting results
      {                                               
      fineinput.plotoffsets = true;
      fineinput.plotthreshold =  atof(word[1]) ;        // pass keyword
      keyword =  word[2] ;      // pass keyword
      writearg(fineinput.plotthreshold);
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"BG"))
        fineinput.plotmagbg = true;
      else if (!strcmp(keyword,"NOBG"))
        fineinput.plotmagbg = false;
      else if (!strcmp(keyword,"ON"))           // actually not allowed...
        fineinput.plotoffsets = true;
      else if (!strcmp(keyword,"OFF"))          // actually not allowed...
        fineinput.plotoffsets = false;
      else
        WARNING.print("FC_PLOT: missing argument(s). (default: 0.4 NOBG)");
      }
    
      else if (!strcmp(keyword,"FC_DEM"))          // input file
      {
      strcpy(fineinput.firefdem,  word[1] );       // pass keyword
      writearg(fineinput.firefdem);
      }

    // **********************************************************************
    else if (!strcmp(keyword,"FC_DEM_FORMAT"))       //  format input file
      {
      keyword =  word[1];       // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
        fineinput.iformatflag = FORMATR4;
      else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
        fineinput.iformatflag = FORMATI2;  // default
      else if (!strcmp(keyword,"I2_BIGENDIAN") || 
               !strcmp(keyword,"SHORT_BIGENDIAN"))
        fineinput.iformatflag = FORMATI2_BIGENDIAN;        // default
      else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
        fineinput.iformatflag = FORMATR8;
      else
        {
        ERROR << "FC_IN_FORMAT: input format "
             <<  keyword
             << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_DEM_SIZE"))         // nrow ncols (lat lon)
      {
      char *pLast1, *pLast2 = NULL;
      fineinput.demrows = strtoul(word[1], &pLast1, BASE10); 
      fineinput.demcols = strtoul(word[2], &pLast2, BASE10);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "FC_DEM_SIZE: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(fineinput.demrows);
      writearg(fineinput.demcols);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_DEM_DELTA"))        // degrees delta lat lon
      {
      fineinput.demdeltalat =  atof(word[1]) ;     // pass keyword
      keyword                    =  word[2] ;         // update keyword
      writearg(fineinput.demdeltalat);
      writearg(keyword);  // lon
      if (isdigit(keyword[0]) || keyword[0]=='.')    // likely to be 2 numbers
        fineinput.demdeltalon = atof(keyword);
      else // default same gridsize
        fineinput.demdeltalon = fineinput.demdeltalat;

      // ______ Store as radians ______
      fineinput.demdeltalat = deg2rad(fineinput.demdeltalat);
      fineinput.demdeltalon = deg2rad(fineinput.demdeltalon);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_DEM_UL"))           // upperleft coordinates
      {
      char *pLast1, *pLast2 = NULL;
      fineinput.demlatleftupper = strtod(word[1], &pLast1);
                  fineinput.demlonleftupper = strtod(word[2], &pLast2);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert
       {
        ERROR << "FC_DEM_UL: "  << word[1] << " : "
                 << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }

      writearg(fineinput.demlatleftupper);
      writearg(fineinput.demlonleftupper);
      fineinput.demlatleftupper = deg2rad(fineinput.demlatleftupper);
      fineinput.demlonleftupper = deg2rad(fineinput.demlonleftupper);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FC_DEM_NODATA"))       // flag for no data
      {
      fineinput.demnodata =  atof(word[1]) ;       // pass keyword
      writearg(fineinput.demnodata);
      }
   

    //added by MCC for fine CCCoregistration
  
// **********************************************************************
    else if (!strcmp(keyword,"FC_SHIFTAZI"))            // true: shift before rs.
      {
      keyword =  word[1] ;      // pass keyword      
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"ON") ||                   // consistent with previous versions
          (keyword[0] == '\0')  || // no keyword
          !strcmp(keyword,"DC"))                        // Doppler centroid polynomial
        fineinput.shiftazi = 1;
      else if (!strcmp(keyword,"DERAMP")) //Only for TOPS
        fineinput.shiftazi = 2;
      else if (!strcmp(keyword,"OFF")   ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)   ||  //comment
                !strcmp(keyword,"NONE"))                 // just in case                
        fineinput.shiftazi = 0;
      else 
        {
        fineinput.shiftazi = 1;
        WARNING << "FC_SHIFTAZI: line: " << linecnt << ": unknown argument: "
             << keyword << "; Set to ON (do shift azimuth spectrum).";
        WARNING.print();
        }        
      } 
    
//  ***************************************
    
    
// ____ start added by FvL ____

// **********************************************************************
// *** RELATIVE TIMING ERROR
// **********************************************************************

    else if (!strcmp(keyword,"RTE_THRESHOLD"))         // treshhold value
      {                                               
      reltiminginput.threshold =  atof(word[1]) ;       // pass keyword
      writearg(reltiminginput.threshold);
      if (reltiminginput.threshold > 1)
        {
        PRINT_ERROR("RTE_THRESHOLD: threshold > 1.")
              throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RTE_MAXITER"))     // max number of offsets to reject
      {                                               
      reltiminginput.maxiter =  atoi(word[1]) ;         // pass keyword
      writearg(reltiminginput.maxiter);
      if (reltiminginput.maxiter < 0)
        {
        WARNING.print("RTE_MAXITER: max. number of points to remove < 0? (using 0)");
        reltiminginput.maxiter = 0;
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RTE_K_ALPHA"))     // critical value for outlier removal
      {                                               
      reltiminginput.k_alpha =  atof(word[1]) ;         // pass keyword
      writearg(reltiminginput.k_alpha);
      if (reltiminginput.k_alpha < 0)
        {
        WARNING.print("RTE_K_ALPHA: critical value < 0.0?");
        reltiminginput.k_alpha = 1.97;
        }
      }


// **********************************************************************
// *** DEM ASSISTED COREGISTRATION
// **********************************************************************

    else if (!strcmp(keyword,"DAC_IN_DEM"))          // input file
      {
      strcpy(demassistinput.firefdem,  word[1] );       // pass keyword
      writearg(demassistinput.firefdem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_IN_FORMAT"))       //  format input file
      {
      keyword =  word[1];       // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
        demassistinput.iformatflag = FORMATR4;
      else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
        demassistinput.iformatflag = FORMATI2;  // default
      else if (!strcmp(keyword,"I2_BIGENDIAN") || 
               !strcmp(keyword,"SHORT_BIGENDIAN"))
        demassistinput.iformatflag = FORMATI2_BIGENDIAN;        // default
      else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
        demassistinput.iformatflag = FORMATR8;
      else
        {
        ERROR << "DAC_IN_FORMAT: input format "
             <<  keyword
             << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_IN_SIZE"))         // nrow ncols (lat lon)
      {
      char *pLast1, *pLast2 = NULL;
      demassistinput.demrows = strtoul(word[1], &pLast1, BASE10); 
      demassistinput.demcols = strtoul(word[2], &pLast2, BASE10);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "DAC_IN_SIZE: "  << word[1] << " : " 
              << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(demassistinput.demrows);
      writearg(demassistinput.demcols);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_IN_DELTA"))        // degrees delta lat lon
      {
      demassistinput.demdeltalat =  atof(word[1]) ;     // pass keyword
      keyword                    =  word[2] ;         // update keyword
      writearg(demassistinput.demdeltalat);
      writearg(keyword);  // lon
      if (isdigit(keyword[0]) || keyword[0]=='.')    // likely to be 2 numbers
        demassistinput.demdeltalon = atof(keyword);
      else // default same gridsize
        demassistinput.demdeltalon = demassistinput.demdeltalat;

      // ______ Store as radians ______
      demassistinput.demdeltalat = deg2rad(demassistinput.demdeltalat);
      demassistinput.demdeltalon = deg2rad(demassistinput.demdeltalon);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_IN_UL"))           // upperleft coordinates
      {
      char *pLast1, *pLast2 = NULL;
      demassistinput.demlatleftupper = strtod(word[1], &pLast1);
                  demassistinput.demlonleftupper = strtod(word[2], &pLast2);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert
       {
        ERROR << "DAC_IN_UL: "  << word[1] << " : "
                 << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }

      writearg(demassistinput.demlatleftupper);
      writearg(demassistinput.demlonleftupper);
      demassistinput.demlatleftupper = deg2rad(demassistinput.demlatleftupper);
      demassistinput.demlonleftupper = deg2rad(demassistinput.demlonleftupper);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_IN_NODATA"))       // flag for no data
      {
      demassistinput.demnodata =  atof(word[1]) ;       // pass keyword
      writearg(demassistinput.demnodata);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_OUT_DEM"))          // name of output file
      {
      strcpy(demassistinput.fodem,  word[1] );  // pass keyword
      writearg(demassistinput.fodem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_OUT_DEMI"))          // name of output file
      {
      strcpy(demassistinput.fodemi,  word[1] );         // pass keyword
      writearg(demassistinput.fodemi);
      }

// **********************************************************************
    else if (!strcmp(keyword,"DAC_OUT_DEM_LP"))        // name of output file
      {
      strcpy(demassistinput.forefdemhei,  word[1] );    // pass keyword
      writearg(demassistinput.forefdemhei);
      }

// ____ end added by FvL ____


// **********************************************************************
// *** COMPUTATION OF COREGISTRATION PARAMETERS
// **********************************************************************
    else if (!strcmp(keyword,"CPM_THRESHOLD"))         // treshhold value
      {                                               
      coregpminput.threshold =  atof(word[1]) ;         // pass keyword
      writearg(coregpminput.threshold);
      if (coregpminput.threshold > 1)
        {
        PRINT_ERROR("CPM_THRESHOLD: threshold > 1.")
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_DEGREE"))           // degree of polynomial
      {                                               
      coregpminput.degree =  atoi(word[1]) ;    // pass keyword
      writearg(coregpminput.degree);
      if (coregpminput.degree > 4)
        WARNING.print("CPM_DEGREE: degree > 4 dangerous at edges?");
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_WEIGHT"))           // weightmatrix
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"NONE"))
        coregpminput.weightflag = 0;          // default
      else if (!strcmp(keyword,"LINEAR"))
        coregpminput.weightflag = 1;
      else if (!strcmp(keyword,"QUADRATIC"))
        coregpminput.weightflag = 2;
      else if (!strcmp(keyword,"BAMLER"))
        coregpminput.weightflag = 3;
      else
        {
        ERROR << "CPM_WEIGHT: data weighting option: "
             <<  keyword
             << " not known for computation of coregistration parameters on line "
             << linecnt << ".";
            PRINT_ERROR(ERROR.get_str())
            throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_MAXITER"))     // max number of offsets to reject
      {                                               
      coregpminput.maxiter =  atoi(word[1]) ;   // pass keyword
      writearg(coregpminput.maxiter);
      if (coregpminput.maxiter < 0)
        {
        WARNING.print("CPM_MAXITER: max. number of points to remove < 0? (using 0)");
        coregpminput.maxiter = 0;
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_K_ALPHA"))     // critical value for outlier removal
      {                                               
      coregpminput.k_alpha =  atof(word[1]) ;   // pass keyword
      writearg(coregpminput.k_alpha);
      if (coregpminput.k_alpha < 0)
        {
        WARNING.print("CPM_K_ALPHA: critical value < 0.0?");
        coregpminput.k_alpha = 1.97;
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_DUMP"))               // boolean
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF")) 
        coregpminput.dumpmodel = false;
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        coregpminput.dumpmodel = true;
      else
        {
        WARNING << "CPM_DUMP: line " << linecnt 
             << ": argument: " << keyword 
             << " not recognized, no dumping to files.";
        WARNING.print();
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CPM_PLOT"))           // plotting results
      {                                               
      coregpminput.plot = true;
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"BG"))
        coregpminput.plotmagbg = true;
      else if (!strcmp(keyword,"NOBG"))
        coregpminput.plotmagbg = false;
      else if (!strcmp(keyword,"ON"))           // actually not allowed...
        coregpminput.plot = true;
      else if (!strcmp(keyword,"OFF"))          // actually not allowed...
        coregpminput.plot = false;
      else
        WARNING.print("CPM_PLOT: missing argument. (default NOBG magnitude background)");
      }


// **********************************************************************
// *** RANGE FILTERING
// **********************************************************************
    else if (!strcmp(keyword,"RF_METHOD"))            // method range filtering
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"ADAPTIVE"))
        filtrangeinput.method = rf_adaptive;
      else if (!strcmp(keyword,"PORBITS"))
        filtrangeinput.method = rf_porbits;
      else
        {
        ERROR << "RF_METHOD: method "
             <<  keyword
             << " not known for range filtering. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_SLOPE"))               // terrain method porbits
      {
      filtrangeinput.terrainslope =  atof(word[1]) ;    // pass keyword
      writearg(filtrangeinput.terrainslope);
      deg2rad(filtrangeinput.terrainslope);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_THRESHOLD"))           // threshhold value
      {
      filtrangeinput.SNRthreshold =  atof(word[1]) ;    // pass keyword
      writearg(filtrangeinput.SNRthreshold);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_HAMMING"))             // alpha
      {
      filtrangeinput.hammingalpha =  atof(word[1]) ;    // pass keyword
      writearg(filtrangeinput.hammingalpha);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_WEIGHTCORR"))          // boolean
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF")) 
        filtrangeinput.doweightcorrel = false;
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        filtrangeinput.doweightcorrel = true;
      else
        {
        filtrangeinput.doweightcorrel = false;          // default already...
        WARNING << "RF_WEIGHTCORREL: line " << linecnt 
             << ": argument: " << keyword 
             << " not recognized, weighting correlation set to OFF.";
        WARNING.print();
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_OVERSAMPLE"))         // int
      {
      filtrangeinput.oversample =  atoi(word[1]) ;      // pass keyword
      writearg(filtrangeinput.oversample);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_NLMEAN"))              // take mean over nlmean lines
      {
      filtrangeinput.nlmean =  atoi(word[1]) ;         // pass keyword
      if (!isodd(filtrangeinput.nlmean))               // check if ODD value  [MA]                   
        {
        PRINT_ERROR("RF_NLMEAN value has to be odd.")
        throw(argument_error);
        }
      writearg(filtrangeinput.nlmean);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_FFTLENGTH"))           // adaptive length
      {
      filtrangeinput.fftlength =  atoi(word[1]) ;       // pass keyword
      writearg(filtrangeinput.fftlength);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_OVERLAP"))             // overlap blocks
      {
      filtrangeinput.overlap =  atoi(word[1]) ;         // pass keyword
      writearg(filtrangeinput.overlap);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_OUT_MASTER"))          // filename
      {
      strcpy(filtrangeinput.fomaster,  word[1] );       // pass keyword
      writearg(filtrangeinput.fomaster);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_OUT_SLAVE"))           // filename
      {
      strcpy(filtrangeinput.foslave,  word[1] );        // pass keyword
      writearg(filtrangeinput.foslave);
      }

// **********************************************************************
    else if (!strcmp(keyword,"RF_OUT_FORMAT"))          // output format
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CR4"))
        filtrangeinput.oformatflag = FORMATCR4;         // default
      else if (!strcmp(keyword,"CI2"))
        filtrangeinput.oformatflag = FORMATCI2;
      else
        {
        ERROR << "RF_OUT_FORMAT: output format "
             <<  keyword
             << " not known for range filtering. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }


// **********************************************************************
// *** FLAT EARTH CORRECTION == compute reference phase since Feb-2000
// **********************************************************************
    else if (!strcmp(keyword,"FE_METHOD"))            // method selector flatearth
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"PORBITS"))
        comprefphainput.method=fe_porbits;
      else if (!strcmp(keyword,"METHOD2"))
        comprefphainput.method=fe_method2;
      else
        {
        ERROR << "FE_METHOD: argument: " << keyword << " not recognized.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"FE_DEGREE"))            // degree for flat earth correction
      {                                               
      comprefphainput.degree =  atoi(word[1]) ;         // pass keyword
      writearg(comprefphainput.degree);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FE_NPOINTS"))           // number of points used 
      {                                               //  for estimation of polynomial
      comprefphainput.Npoints =  atoi(word[1]) ;        // pass keyword          //  flat earth correction.
      writearg(comprefphainput.Npoints);
      }

// **********************************************************************
    else if (!strcmp(keyword,"FE_IN_POS"))            // file with #windows positions
      {                                               
      strcpy(comprefphainput.ifpositions,  word[1] );   // pass keyword
      writearg(comprefphainput.ifpositions);
      }


// **********************************************************************
// *** RESAMPLING (SLAVE)
// **********************************************************************
    else if (!strcmp(keyword,"RS_METHOD"))            // method selector resampling
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CC4P"))
        resampleinput.method = rs_cc4p;                       // default
      else if (!strcmp(keyword,"CC6P"))
        resampleinput.method = rs_cc6p;
      else if (!strcmp(keyword,"TS6P"))
        resampleinput.method = rs_ts6p;
      else if (!strcmp(keyword,"TS8P"))
        resampleinput.method = rs_ts8p;
      else if (!strcmp(keyword,"TS16P"))
        resampleinput.method = rs_ts16p;
      else if (!strcmp(keyword,"KNAB4P"))
        resampleinput.method = rs_knab4p;
      else if (!strcmp(keyword,"KNAB6P"))
        resampleinput.method = rs_knab6p;
      else if (!strcmp(keyword,"KNAB8P"))
        resampleinput.method = rs_knab8p;
      else if (!strcmp(keyword,"KNAB10P"))
        resampleinput.method = rs_knab10p;
      else if (!strcmp(keyword,"KNAB16P"))
        resampleinput.method = rs_knab16p;
      else if (!strcmp(keyword,"RC6P"))
        resampleinput.method = rs_rc6p;
      else if (!strcmp(keyword,"RC12P"))
        resampleinput.method = rs_rc12p;
      else if (!strcmp(keyword,"RECT"))
        resampleinput.method = rs_rect;
      else if (!strcmp(keyword,"TRI"))
        resampleinput.method = rs_tri;
      else
        {
        ERROR << "RS_METHOD: method "
             <<  keyword
             << " not known for resampling. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RS_DBOW"))              // database output window
      {
      resampleinput.dbow.linelo  = atoi(word[1]);
      resampleinput.dbow.linehi  = atoi(word[2]);
      resampleinput.dbow.pixlo   = atoi(word[3]);
      resampleinput.dbow.pixhi   = atoi(word[4]);
      writearg(resampleinput.dbow.linelo);
      writearg(resampleinput.dbow.linehi);
      writearg(resampleinput.dbow.pixlo);
      writearg(resampleinput.dbow.pixhi);
      if (resampleinput.dbow.linelo <= 0 || 
          resampleinput.dbow.pixlo <= 0 ||
          resampleinput.dbow.linelo  > resampleinput.dbow.linehi ||
          resampleinput.dbow.pixlo > resampleinput.dbow.pixhi)
        {
        ERROR << "code 300: Arguments of RS_DBOW card on line " << linecnt
              << " missing. [RS_DBOW  min_line  max_line  min_pixel  max_pixel].";
              PRINT_ERROR(ERROR.get_str())
              throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RS_DBOW_GEO"))  // database output
                                              // window on GEO.LOCATION
      {
       // use dbow to store temp, compute later
       real8 tmp_lat_0;
       real8 tmp_lon_0; 
       real8 tmp_height; // total height, not half.window
       real8 tmp_width;  // total width, not half.windo
        
       char *pLast1, *pLast2, *pLast3, *pLast4 = NULL;
       tmp_lat_0  = strtod(word[1], &pLast1); 
       tmp_lon_0  = strtod(word[2], &pLast2);
       tmp_height = strtod(word[3], &pLast2);
       tmp_width  = strtod(word[4], &pLast2);
      if ( pLast1 == word[1] || pLast2 == word[2] ||
           pLast3 == word[3] || pLast4 == word[4]    ) // fails to convert one of them to double.
       {
        ERROR << "RS_DBOW_GEO: "  << word[1] << " : " 
              << word[2] << " : " << word[3] << " : " << word[4] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
        writearg(tmp_lat_0);
        writearg(tmp_lon_0);
        writearg(tmp_height);
        writearg(tmp_width);
        
        resampleinput.dbow_geo.linelo = uint((360.0+tmp_lat_0)*1e6);
        resampleinput.dbow_geo.linehi = uint((360.0+tmp_lon_0)*1e6);
        resampleinput.dbow_geo.pixlo  = uint(tmp_height);
        resampleinput.dbow_geo.pixhi  = uint(tmp_width);

      }

// **********************************************************************
    else if (!strcmp(keyword,"RS_OUT_FILE"))           // name of output file
      {
      switch (priorrs_fileout)
        {
        case true:
          WARNING << "RS_OUT_FILE: line: " << linecnt << ": "
              << "ignored due to prior occurrence.";
          WARNING.print();
          break;
        default:
          priorrs_fileout = true;
          strcpy(resampleinput.fileout,  word[1] );     // pass keyword (filename)
          writearg(resampleinput.fileout);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"RS_OUT_FORMAT"))          // output format
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"CR4"))
        resampleinput.oformatflag = FORMATCR4;          // default
      else if (!strcmp(keyword,"CI2"))
        resampleinput.oformatflag = FORMATCI2;
      else
        {
        ERROR << "RS_OUT_FORMAT: output format "
             <<  keyword
             << " not known for resampling. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
   
     else if (!strcmp(keyword,"RS_SHIFTAZI"))            // true: shift before rs.
      {
      keyword =  word[1] ;      // pass keyword      
      writearg(keyword);
      toupper(keyword);
      
      if (!strcmp(keyword,"OFF")   ||
          !strncmp(keyword,"//",2) ||              // comment
          !strncmp(keyword,"#",1)  ||              //comment
          !strcmp(keyword,"NONE"))                 // just in case                
        resampleinput.shiftazi = 0;
      else if (!strcmp(keyword,"DERAMP")) //Only for TOPS
        resampleinput.shiftazi = 2;
      else if (!strcmp(keyword,"ON")  ||                   // consistent with previous versions          
               !strcmp(keyword,"DC") ||                   // Doppler centroid polynomial
               (keyword[0] == '\0'))                        // no keyword
        resampleinput.shiftazi = 1;            
      else 
        {
        resampleinput.shiftazi = 1;
        WARNING << "RS_SHIFTAZI: line: " << linecnt << ": unknown argument: "
             << keyword << "; Set to ON (do shift azimuth spectrum).";
        WARNING.print();
        } 
      } 

// **********************************************************************
// *** COMPUTATION OF INTERFEROGRAM
// **********************************************************************
    else if (!strcmp(keyword,"INT_METHOD"))            // method selector interfero
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OLD"))
        interferoinput.method = int_oldmethod;
      else if (!strcmp(keyword,"OVERSAMPLE"))
        interferoinput.method = int_oversample;
      else
        {
        ERROR << "INT_METHOD: method "
             <<  keyword
             << " not known for interfero. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"INT_OUT_INT"))           // name of output file
      {
      strcpy(interferoinput.foint,  word[1] );  // pass keyword
      writearg(interferoinput.foint);
      }

// **********************************************************************
    else if (!strcmp(keyword,"INT_OUT_CINT"))          // name of output file
      {
      strcpy(interferoinput.focint,  word[1] );         // pass keyword
      writearg(interferoinput.focint);
      }

// **********************************************************************
//    else if (!strcmp(keyword,"INT_OUT_FE"))         // name of output file
//      interferoinput.foflatearth =  word[1] ;         // pass keyword

// **********************************************************************
    else if (!strcmp(keyword,"INT_MULTILOOK"))        // multilookfactors
      {
      interferoinput.multilookL  = atoi(word[1]);       // pass keyword
      interferoinput.multilookP  = atoi(word[2]);       // pass keyword
      writearg(interferoinput.multilookL);
      writearg(interferoinput.multilookP);
      }

// **********************************************************************
// *** COMPUTATION OF COHERENCE
// **********************************************************************
    else if (!strcmp(keyword,"COH_METHOD"))            // method selector coherence
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"REFPHASE_ONLY"))                    // OLD
        coherenceinput.method = coh_oldmethod;
      else if (!strcmp(keyword,"INCLUDE_REFDEM"))              // NEW
        coherenceinput.method = coh_newmethod;
      else
        {
        ERROR << "COH_METHOD: method "
             <<  keyword
             << " not known for coherence. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"COH_OUT_COH"))          // name of output file
      {
      strcpy(coherenceinput.focoh,  word[1] );  // pass keyword
      writearg(coherenceinput.focoh);
      }

// **********************************************************************
    else if (!strcmp(keyword,"COH_OUT_CCOH"))          // name of output file
      {
      strcpy(coherenceinput.foccoh,  word[1] );         // pass keyword
      writearg(coherenceinput.foccoh);
      }

// **********************************************************************
    else if (!strcmp(keyword,"COH_MULTILOOK"))        // multilookfactors
      {
      coherenceinput.multilookL =  atoi(word[1]);       // pass keyword
      coherenceinput.multilookP =  atoi(word[2]);       // pass keyword
      writearg(coherenceinput.multilookL);
      writearg(coherenceinput.multilookP);
      }

// **********************************************************************
    else if (!strcmp(keyword,"COH_WINSIZE"))          // estimator winsize
      {
      coherenceinput.cohsizeL =  atoi(word[1]);         // pass keyword
      coherenceinput.cohsizeP =  atoi(word[2]);         // pass keyword
      writearg(coherenceinput.cohsizeL);
      writearg(coherenceinput.cohsizeP);
      }

// **********************************************************************
// *** SUBTRACTION OF REFERENCE PHASE
// **********************************************************************
    else if (!strcmp(keyword,"SRP_METHOD"))           // method selector ref. phase
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"POLYNOMIAL"))
        subtrrefphainput.method = srp_polynomial;
      else if (!strcmp(keyword,"EXACT"))
        subtrrefphainput.method = srp_exact;
      else
        {
        ERROR << "SRP_METHOD: method "
             <<  keyword
             << " not known for subtraction of ref. phase. line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"SRP_OUT_CINT"))          // name of output file
      {
      strcpy(subtrrefphainput.focint,  word[1] );       // pass keyword
      writearg(subtrrefphainput.focint);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SRP_OUT_REFPHA"))       // name of output file
      {
      strcpy(subtrrefphainput.forefpha,  word[1] );     // pass keyword
      writearg(subtrrefphainput.forefpha);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SRP_DUMPREFPHA"))      // true: dump ref.pha
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"ON"))
        subtrrefphainput.dumponlyrefpha = true;
      else if (!strcmp(keyword,"OFF")   ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        subtrrefphainput.dumponlyrefpha = false;
      else
        {
        subtrrefphainput.dumponlyrefpha = false;
        WARNING << "SRP_DUMPREFPHA: line: " << linecnt << ": unknown argument: "
             << keyword << "; Set to OFF (no dump).";
        WARNING.print();
        }
      }

// **********************************************************************

// ___________ added by FvL      
    else if (!strcmp(keyword,"SRP_OUT_H2PH"))          // name of output file
      {
      strcpy(subtrrefphainput.foh2ph,  word[1] );       // pass keyword
      writearg(subtrrefphainput.foh2ph);
      }
// ___________ end added by FvL


// **********************************************************************
    else if (!strcmp(keyword,"SRP_MULTILOOK"))        // multilookfactors
      {
      subtrrefphainput.multilookL  =  atoi(word[1]) ;   // pass keyword
      keyword                      =  word[2] ;               // pass keyword
      writearg(subtrrefphainput.multilookL);
      writearg(keyword);
      if (isdigit(keyword[0]))
        subtrrefphainput.multilookP = atoi(keyword);
      else                                              // default same factor
        subtrrefphainput.multilookP = subtrrefphainput.multilookL;
      }

// **********************************************************************
// *** PHASE FILTER
// **********************************************************************
    else if (!strcmp(keyword,"PF_METHOD"))           // method selector phase filtering
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"GOLDSTEIN"))
        filtphaseinput.method = fp_goldstein;
      else if (!strcmp(keyword,"MODGOLDSTEIN"))
        filtphaseinput.method = fp_modgoldstein;
      else if (!strcmp(keyword,"SPATIALCONV"))
        filtphaseinput.method = fp_spatialconv;
      else if (!strcmp(keyword,"SPECTRAL"))
        filtphaseinput.method = fp_spectral;
      else
        {
        ERROR << "PF_METHOD: method "
             <<  keyword
             << " not known for phase filtering. line "
             << linecnt << ".";
             PRINT_ERROR(ERROR.get_str())
             throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_OUT_FILE"))          // filename
      {
      strcpy(filtphaseinput.fofiltphase,  word[1] );    // pass keyword
      writearg(filtphaseinput.fofiltphase);
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_IN_FILE"))           // filename
      {
      strcpy(filtphaseinput.fifiltphase ,  word[1] );   // pass keyword
      char *pLast = NULL;
      filtphaseinput.finumlines = strtoul(word[1], &pLast, BASE10); // pass numoflines
      if ( pLast == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "PF_IN_FILE (numoflines): "  <<  word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(filtphaseinput.fifiltphase);
      writearg(filtphaseinput.finumlines);
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_BLOCKSIZE"))         // buffersize
      {
      filtphaseinput.blocksize =  atoi(word[1]) ;       // pass keyword
      writearg(filtphaseinput.blocksize);
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_ALPHA"))             // alpha
      {
      filtphaseinput.alpha =  atof(word[1]) ;   // pass keyword
      writearg(filtphaseinput.alpha);
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_OVERLAP"))           // overlap
      {
      filtphaseinput.overlap =  atoi(word[1]) ;         // pass keyword
      writearg(filtphaseinput.overlap);
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_KERNEL"))            // conv. kernel e.g. 3 1 1 1
      {
      keyword =  word[1] ;      // pass keyword
      const int32 sizekernel = atoi(keyword);
      writearg(sizekernel);
      if (!(isodd(sizekernel)))
        {
         PRINT_ERROR("PF_KERNEL: size must be odd! (add 0, center around midpix)")
         throw(keyword_error);
        }
      filtphaseinput.kernel.resize(1,sizekernel);
      real4 sum=0.;
      for (int32 argnum=0; argnum<sizekernel; ++argnum)
        {
         keyword =  word[1] ;   // pass keyword
         if (!(isdigit(keyword[0])))
         WARNING.print("kernel seems to be wrong?");
         const real4 in = atof(keyword);
         writearg(in);
         filtphaseinput.kernel(0,argnum) = in;
         sum += abs(in);                                        // kernel -1 1
        }
        if (sum!=1)
        filtphaseinput.kernel /= sum;                   // normalize
        INFO << "PF_KERNEL: Input kernel normalized by: " <<  sum;
        INFO.print();
      }

// **********************************************************************
    else if (!strcmp(keyword,"PF_IN_KERNEL2D"))       // filename for ascii 2d kernel
      {
      strcpy(filtphaseinput.fikernel2d,  word[1] );     // pass keyword
      writearg(filtphaseinput.fikernel2d);
      }

// *******************************************************************
// *** DINSAR
// *******************************************************************
    else if (!strcmp(keyword,"DI_OUT_FILE"))            // output filename
      {                                               
      strcpy(dinsarinput.fodinsar,  word[1] );  // pass keyword
      writearg(dinsarinput.fodinsar);
      }

    else if (!strcmp(keyword,"DI_OUT_SCALED"))          // output filename
      {                                               
      strcpy(dinsarinput.foscaleduint,  word[1] );      // pass keyword
      writearg(dinsarinput.foscaleduint);
      }

    else if (!strcmp(keyword,"DI_IN_TOPOMASTER"))       // input resultfilename
      {                                               
      strcpy(dinsarinput.topomasterresfile,  word[1] );         // pass keyword
      writearg(dinsarinput.fodinsar);
      }

    else if (!strcmp(keyword,"DI_IN_TOPOSLAVE"))        // input resultfilename
      {                                               
      strcpy(dinsarinput.toposlaveresfile,  word[1] );  // pass keyword
      writearg(dinsarinput.fodinsar);
      }

    else if (!strcmp(keyword,"DI_IN_TOPOINT"))          // input resultfilename
      {                                               
      strcpy(dinsarinput.topointresfile,  word[1] );    // pass keyword
      writearg(dinsarinput.fodinsar);
      }


// **********************************************************************
// *** COMPUTATION OF REFERENCE DEM (phase)
// **********************************************************************
//    else if (!strcmp(keyword,"CRD_METHOD"))          // name of output file
//      {
//      filename =  word[1] ;   // pass keyword
//      writearg(keyword2);
//      toupper(keyword2);
//      if (!strcmp(keyword2,"TRILINEAR") || !strcmp(keyword2,"TRI_LINEAR") ||
//          !strcmp(keyword2,"TRILIN") || !strcmp(keyword2,"TRI_LIN"))
//        comprefdeminput.method = crd_trilinear;
//      else if (!strcmp(keyword2,"NN") || !strcmp(keyword2,"NEAREST") ||
//               !strcmp(keyword2,"NEAREST_NEIGHBOR" ))
//        comprefdeminput.method = crd_nearest;
//      else
//        {
//        ERROR << "CRD_METHOD: "
//             <<  filename
//             << " not known (use TRILINEAR or NEAREST); line "
//             << linecnt << ".";
//      PRINT_ERROR(ERROR.get_str())
//      throw(keyword_error);
//        }
//      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_DEM"))          // input file
      {
      strcpy(comprefdeminput.firefdem,  word[1] );      // pass keyword
      writearg(comprefdeminput.firefdem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_FORMAT"))       //  format input file
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
        comprefdeminput.iformatflag = FORMATR4;
      else if (!strcmp(keyword,"I2") || !strcmp(keyword,"SHORT"))
        comprefdeminput.iformatflag = FORMATI2; // default
      else if (!strcmp(keyword,"I2_BIGENDIAN") || 
               !strcmp(keyword,"SHORT_BIGENDIAN"))
        comprefdeminput.iformatflag = FORMATI2_BIGENDIAN;       // default
      else if (!strcmp(keyword,"R8") || !strcmp(keyword,"REAL8"))
        comprefdeminput.iformatflag = FORMATR8;
      else
        {
        ERROR << "CRD_IN_FORMAT: input format "
             <<  keyword
             << " not known (R4 R8 I2 (native) SHORT_BIGENDIAN); line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_SIZE"))         // nrow ncols (lat lon)
      {
      comprefdeminput.demrows =  strtoul(word[1], NULL, BASE10);         // pass keyword [MA] atoi instead strtoul 
      comprefdeminput.demcols =  strtoul(word[2], NULL, BASE10);         // pass keyword
      writearg(comprefdeminput.demrows);
      writearg(comprefdeminput.demcols);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_DELTA"))        // degrees delta lat lon
      {
      comprefdeminput.demdeltalat =  atof(word[1]) ;    // pass keyword
      keyword                     =  word[2] ;        // update keyword
      writearg(comprefdeminput.demdeltalat);
      writearg(keyword);
      if (isdigit(keyword[0]) || keyword[0]=='.')    // likely to be 2 numbers
        comprefdeminput.demdeltalon = atof(keyword);
      else // default same gridsize
        comprefdeminput.demdeltalon = comprefdeminput.demdeltalat;

      // ______ Store as radians ______
      comprefdeminput.demdeltalat = deg2rad(comprefdeminput.demdeltalat);
      comprefdeminput.demdeltalon = deg2rad(comprefdeminput.demdeltalon);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_UL"))           // upperleft coordinates
      {
      char *pLast1, *pLast2 = NULL;
      comprefdeminput.demlatleftupper = strtod(word[1], &pLast1); 
                  comprefdeminput.demlonleftupper = strtod(word[2], &pLast2);
      if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
       {
        ERROR << "CRD_IN_UL: "  << word[1] << " : " << word[2] << " are not valid.";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
       }
      writearg(comprefdeminput.demlatleftupper);
      writearg(comprefdeminput.demlonleftupper);
      comprefdeminput.demlatleftupper = deg2rad(comprefdeminput.demlatleftupper);
      comprefdeminput.demlonleftupper = deg2rad(comprefdeminput.demlonleftupper);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_IN_NODATA"))       // flag for no data
      {
      comprefdeminput.demnodata =  atof(word[1]) ;      // pass keyword
      writearg(comprefdeminput.demnodata);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_INCLUDE_FE"))      // true: ref.pha incl. flat earth
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF"))
        comprefdeminput.includerefpha = false;
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        comprefdeminput.includerefpha = true;
      else
        {
        comprefdeminput.includerefpha = false;
        WARNING << "CRD_INCLUDE_FE: line: " << linecnt << ": unknown argument: "
             << keyword 
             << "; Set to OFF (computing pure topo phase w.r.t. flat earth).";
       WARNING.print();
       }
      }

// **********************************************************************
//    else if (!strcmp(keyword,"CRD_DENSE"))            // factor
//      {
//      comprefdeminput.extradense =  word[1] ;         // pass keyword
//      writearg(comprefdeminput.extradense);
//      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_OUT_DEM"))          // name of output file
      {
      strcpy(comprefdeminput.fodem,  word[1] );         // pass keyword
      writearg(comprefdeminput.fodem);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_OUT_DEMI"))          // name of output file
      {
      strcpy(comprefdeminput.fodemi,  word[1] );        // pass keyword
      writearg(comprefdeminput.fodemi);
      }

// **********************************************************************
    else if (!strcmp(keyword,"CRD_OUT_FILE"))          // name of output file
      {
      strcpy(comprefdeminput.forefdem,  word[1] );      // pass keyword
      writearg(comprefdeminput.forefdem);
      }

// ___________ added by FvL      
// **********************************************************************
    else if (!strcmp(keyword,"CRD_OUT_H2PH"))          // name of output file
      {
      strcpy(comprefdeminput.foh2ph,  word[1] );        // pass keyword
      writearg(comprefdeminput.foh2ph);
      }
// ___________ end added by FvL

// **********************************************************************
    else if (!strcmp(keyword,"CRD_OUT_DEM_LP"))        // name of output file
      {
      strcpy(comprefdeminput.forefdemhei,  word[1] );   // pass keyword
      writearg(comprefdeminput.forefdemhei);
      }

 
// **********************************************************************
// *** SUBTRACTION OF REFERENCE DEM (phase)
// **********************************************************************
    else if (!strcmp(keyword,"SRD_OUT_CINT"))          // name of output file
      {
      strcpy(subtrrefdeminput.focint,  word[1] );       // pass keyword
      writearg(subtrrefdeminput.focint);
      }

// **********************************************************************
    else if (!strcmp(keyword,"SRD_OFFSET"))          // Line Pixel
      {
      subtrrefdeminput.offsetL =  atoi(word[1]);        // pass keyword
      subtrrefdeminput.offsetP =  atoi(word[2]);        // pass keyword
      writearg(subtrrefdeminput.offsetL);
      writearg(subtrrefdeminput.offsetP);
      }

//  // **********************************************************************
//  // *** ADDITION OF REFERENCE DEM (phase)
//  // **********************************************************************
//      else if (!strcmp(keyword,"ARD_OUT_CINT"))          // name of output file
//        {
//        addrefdeminput.focint =  word[1] ;    // pass keyword
//        }

// **********************************************************************
// *** UNWRAPPING
// **********************************************************************
    else if (!strcmp(keyword,"UW_METHOD"))  // method selector unwrapping
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"RAMON"))
        unwrapinput.method = uw_method1;
      else if (!strcmp(keyword,"SNAPHU"))
        unwrapinput.method = uw_method2;//   default
      else if (!strcmp(keyword,"MCF_DLR"))
        unwrapinput.method = uw_method3;
      else
        {
        ERROR << "UW_METHOD: method "
             <<  keyword
             << " not known for unwrapping on line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SEEDS"))             // position of seeds
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      if (isdigit(keyword[0]))
        {
         unwrapinput.deltaLseed = atoi(keyword);
         keyword =  word[2] ;   // update keyword
         writearg(keyword);
         if (isdigit(keyword[0]))
            unwrapinput.deltaPseed = atoi(keyword);
         else
            unwrapinput.deltaPseed = unwrapinput.deltaLseed;
        }
      else      // assume no numbers but filename with seeds
        {
         strcpy(unwrapinput.seedfile,keyword);
        }
       }

// **********************************************************************
    else if (!strcmp(keyword,"UW_OUT_FILE"))          // filename output unwrapped int.
      {
      strcpy(unwrapinput.fouint,  word[1] );    // pass keyword
      writearg(unwrapinput.fouint);
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_OUT_FORMAT"))       // output format
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"R4") || !strcmp(keyword,"REAL4"))
        unwrapinput.oformatflag = FORMATR4;
      else if (!strcmp(keyword,"CR4") || !strcmp(keyword,"COMPLEXR4"))
        {
        WARNING.print("UW_OUT_FORMAT = CR4 --> Using hgt format");
        unwrapinput.oformatflag = FORMATHGT;// default
        }
      else if (!strcmp(keyword,"HGT"))
        unwrapinput.oformatflag = FORMATHGT;// default
      else
        {
        unwrapinput.oformatflag = FORMATHGT;// default
        WARNING << "UW_OUT_FORMAT: output format "
             <<  keyword
             << " not known (R4 or HGT, using HGT); line "
             << linecnt << ".";
        WARNING.print();
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_OUT_REGIONS"))       // filename output regions
      {
      strcpy(unwrapinput.foregions,  word[1] );         // pass keyword
      writearg(unwrapinput.foregions);
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SNAPHU_MODE"))       // filename output regions
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"TOPO"))
        strcpy(unwrapinput.snaphu_mode,"TOPO");        // default TOPO
      else if (!strcmp(keyword,"DEFO"))
        strcpy(unwrapinput.snaphu_mode,"DEFO");
      else if (!strcmp(keyword,"SMOOTH"))
        strcpy(unwrapinput.snaphu_mode,"SMOOTH");
      else if (!strcmp(keyword,"NOSTATCOSTS"))
        strcpy(unwrapinput.snaphu_mode,"NOSTATCOSTS");
      else
        {
        ERROR << "UW_SNAPHU_MODE: "
             <<  keyword
             << " not known for unwrapping on line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SNAPHU_INIT"))
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"MST"))
        strcpy(unwrapinput.snaphu_init,"MST");     // default mst
      else if (!strcmp(keyword,"MCF"))
        strcpy(unwrapinput.snaphu_init,"MCF");
      else
        {
        WARNING << "UW_SNAPHU_INIT: "
             <<  keyword
             << " not known for unwrapping on line "
             << linecnt << " (using MST).";
        WARNING.print();
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SNAPHU_LOG"))
      {
      strcpy(unwrapinput.snaphu_log,  word[1] );        // pass keyword
      writearg(unwrapinput.snaphu_log);
      //strcpy(unwrapinput.snaphu_log,"-l ");
      //strcat(unwrapinput.snaphu_log,keyword);
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SNAPHU_COH"))
      {
      strcpy(unwrapinput.snaphu_coh,  word[1] );        // pass keyword
      writearg(unwrapinput.snaphu_coh);
      }

// **********************************************************************
    else if (!strcmp(keyword,"UW_SNAPHU_VERBOSE"))
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF"))
        strcpy(unwrapinput.snaphu_verbose,"FALSE");
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        strcpy(unwrapinput.snaphu_verbose,"TRUE");
      else
        {
        strcpy(unwrapinput.snaphu_verbose,"TRUE");
        WARNING << "UW_SNAPHU_VERBOSE: line: " 
             << linecnt << ": unknown argument: " << keyword 
             << "; Set to ON.";
       WARNING.print();
        }
      }

    else if (!strcmp(keyword,"UW_SNAPHU_NTILEROW"))
      {
             unwrapinput.ntilerow =  atoi(word[1]) ;    // pass keyword
             writearg(unwrapinput.ntilerow);
       if (unwrapinput.ntilerow > 50)
         {
          WARNING.print("UW_SNAPHU_NTILEROW > 100 tiles, is this okay? ");
         } 
      INFO << "UW_SNAPHU_NTILEROW: \t " 
           << unwrapinput.ntilerow;
      INFO.print();
      }

    else if (!strcmp(keyword,"UW_SNAPHU_NTILECOL"))
      {
       unwrapinput.ntilecol =  atoi(word[1]) ;  // pass keyword
       writearg(unwrapinput.ntilecol);
       if (unwrapinput.ntilecol > 50)
       {
         WARNING.print("UW_SNAPHU_NTILECOL > 100 tiles, is this okay? ");
       }
     
       INFO << "UW_SNAPHU_NTILECOL: \t "
            << unwrapinput.ntilecol;
       INFO.print();
      }   

    else if (!strcmp(keyword,"UW_SNAPHU_ROWOVRLP"))
      {
        unwrapinput.rowovrlp =  atoi(word[1]) ;         // pass keyword
        writearg(unwrapinput.rowovrlp);
        INFO << "UW_SNAPHU_ROWOVRLP: \t "
             << unwrapinput.rowovrlp;
        INFO.print();
      }   

    else if (!strcmp(keyword,"UW_SNAPHU_COLOVRLP"))
      {
        unwrapinput.colovrlp =  atoi(word[1]) ;         // pass keyword
        writearg(unwrapinput.colovrlp);
        INFO << "UW_SNAPHU_COLOVRLP: \t "
             << unwrapinput.colovrlp;
        INFO.print();
      }   

    else if (!strcmp(keyword,"UW_SNAPHU_NPROC"))
      {
        unwrapinput.nproc =  atoi(word[1]) ;    // pass keyword
        writearg(unwrapinput.nproc);
        if (unwrapinput.ntilecol > 2)
        {
          WARNING.print("UW_SNAPHU_NPROC > 2CPUs, do you have a cluster?");
        }
      
        INFO << "UW_SNAPHU_NPROC: \t "
             << unwrapinput.ntilecol;
        INFO.print();
      }   

    else if (!strcmp(keyword,"UW_SNAPHU_TILECOSTTHRESH"))
      {
        unwrapinput.tilecostthresh =  atoi(word[1]) ;   // pass keyword
        writearg(unwrapinput.tilecostthresh);
        if (unwrapinput.ntilecol > 500)
        {
          WARNING.print("UW_SNAPHU_TILECOSTTHRESH > 500, do you have a cluster?");
        }
      
        INFO << "UW_SNAPHU_TILECOSTTHRESH: \t "
             << unwrapinput.tilecostthresh;
        INFO.print();
      }   
    else if (!strcmp(keyword,"UW_SNAPHU_DUMPONLYCONF")) //
      {
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"OFF"))
        unwrapinput.snaphu_dumponlyconf=false;
      else if (!strcmp(keyword,"ON")    ||
               !strncmp(keyword,"//",2) ||              // comment
               !strncmp(keyword,"#",1)  ||              // comment
               !(keyword[0] == '\0'))                   // no keyword
        unwrapinput.snaphu_dumponlyconf=true;
      else 
        {
        unwrapinput.snaphu_dumponlyconf=true;
        WARNING << "UW_SNAPHU_DUMPONLYCONF: line: " 
             << linecnt << ": unknown argument: " << keyword 
             << "; Set to ON.";
       WARNING.print();
        }
      }

// **********************************************************************
// *** SLANT to HEIGHT CONVERSION
// **********************************************************************
    else if (!strcmp(keyword,"S2H_METHOD"))           // method selector slant2height
      {                                               
      keyword =  word[1] ;      // pass keyword
      writearg(keyword);
      toupper(keyword);
      if (!strcmp(keyword,"SCHWABISCH"))
        slant2hinput.method = s2h_schwabisch;
      else if (!strcmp(keyword,"AMBIGUITY"))
        slant2hinput.method = s2h_ambiguity;
      else if (!strcmp(keyword,"RODRIGUEZ"))
        slant2hinput.method = s2h_rodriguez;
      else
        {
        ERROR << "S2H_METHOD: method "
             <<  keyword
             << " not known for slant to height conversion, line "
             << linecnt << ".";
        PRINT_ERROR(ERROR.get_str())
        throw(keyword_error);
        }
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_NPOINTS"))          // number of points to use
      {
      slant2hinput.Npoints =  atoi(word[1]) ;   // pass keyword
      writearg(slant2hinput.Npoints);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_DEGREE1D"))         // degree of 1d polynomial
      {
      slant2hinput.degree1d =  atoi(word[1]) ;  // pass keyword
      writearg(slant2hinput.degree1d);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_DEGREE2D"))          // degree of 2d polynomial
      {
      slant2hinput.degree2d =  atoi(word[1]) ;  // pass keyword
      writearg(slant2hinput.degree2d);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_NHEIGHTS"))         // #heights to evaluate ref.pha
      {
      slant2hinput.Nheights =  atoi(word[1]) ;  // pass keyword
      writearg(slant2hinput.Nheights);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_OUT_HEI"))          // filename output height
      {
      strcpy(slant2hinput.fohei,  word[1] );    // pass keyword
      writearg(slant2hinput.fohei);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_OUT_PHI"))          // filename output latitude
      {
      strcpy(slant2hinput.fophi,  word[1] );    // pass keyword
      writearg(slant2hinput.fophi);
      }

// **********************************************************************
    else if (!strcmp(keyword,"S2H_OUT_LAM"))          // filename output longitude
      {
      strcpy(slant2hinput.folam,  word[1] );    // pass keyword
      writearg(slant2hinput.folam);
      }

// **********************************************************************
// *** GEOCODING
// **********************************************************************
//      else if (!strcmp(keyword,"GEO_METHOD"))
//        {
//        geocodeinput.method =  word[1] ;      // pass keyword
//        }

// **********************************************************************
    else if (!strcmp(keyword,"GEO_OUT_PHI"))          // output file latitude
      {
      strcpy(geocodeinput.fophi,  word[1] );    // pass keyword
      writearg(geocodeinput.fophi);
      }

// **********************************************************************
    else if (!strcmp(keyword,"GEO_OUT_LAM"))          // output file longitude
      {
      strcpy(geocodeinput.folam,  word[1] );    // pass keyword
      writearg(geocodeinput.folam);
      }


// ____ start added by HB ____
// **********************************************************************
// *** ESTORBITS 
// **********************************************************************
    else if (!strcmp(keyword,"EO_METHOD"))
      {
	keyword = word[1];
	writearg(keyword);
	toupper(keyword);
	if (!strcmp(keyword,"LSQ"))
	  estorbitsinput.method = eo_lsq;
	else if (!strcmp(keyword,"GRIDSEARCH"))
	  estorbitsinput.method = eo_gridsearch;
	else
	  {
	    ERROR << "EO_METHOD: unknown method " << keyword 
		  << ", line " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str())
	      throw(keyword_error);
	  }
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_IN_DEM_LP"))
      {
        strcpy(estorbitsinput.fiheightmap, word[1]);
	writearg(estorbitsinput.fiheightmap);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_OUT_RES"))
      {
        strcpy(estorbitsinput.foresiduals, word[1]);
	writearg(estorbitsinput.foresiduals);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_WEIGHTING"))
      {
	keyword = word[1];
	writearg(keyword);
	toupper(keyword);
	if (!strcmp(keyword,"NONE"))
	  estorbitsinput.weighting = eo_noweighting;
	else if (!strcmp(keyword,"COH"))
	  estorbitsinput.weighting = eo_coh;
	else if (!strcmp(keyword,"COH2"))
	  estorbitsinput.weighting = eo_coh2;
	else if (!strcmp(keyword,"PDF"))
	  estorbitsinput.weighting = eo_pdf;
	else
	  {
	    ERROR << "EO_WEIGHTING: unknown weighting scheme " << keyword
		  << ", line " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str())
	      throw(keyword_error);
	  }
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_NPOINTS"))
      {
	estorbitsinput.nobs = atoi(word[1]);
	writearg(estorbitsinput.nobs);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_IN_POS"))
      {
	strcpy(estorbitsinput.ifpositions,word[1] ); 
	writearg(estorbitsinput.ifpositions);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_THRESHOLD"))
      {
	estorbitsinput.threshold = atof(word[1]);
	writearg(estorbitsinput.threshold);
	if (estorbitsinput.threshold > 1)
	  {
	    ERROR << "EO_THRESHOLD: threshold > 1, line " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str())
	      throw(keyword_error);
	  }
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_MAXITER"))
      {
	estorbitsinput.maxiter = atoi(word[1]);
	writearg(estorbitsinput.maxiter);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_K_ALPHA"))
      {
	estorbitsinput.k_alpha = atof(word[1]);
	writearg(estorbitsinput.k_alpha);
	if (estorbitsinput.k_alpha < 0)
	  {
	    ERROR << "EO_K_ALPHA: critical value < 0, line " << linecnt << ".";
	    PRINT_ERROR(ERROR.get_str())
	      throw(keyword_error);
	  }
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_SEARCHSPACE"))
      {                                               
	char *pLast1, *pLast2 = NULL;
	estorbitsinput.maxfringesaz = strtoul(word[1], &pLast1, BASE10); 
	estorbitsinput.maxfringesrg = strtoul(word[2], &pLast2, BASE10);
	if ( pLast1 == word[1] || pLast2 == word[2] ) // fails to convert one of them to double.
	  {
	    ERROR << "EO_SEARCHSPACE: "  << word[1] << " : " 
		  << word[2] << " are not valid.";
	    PRINT_ERROR(ERROR.get_str())
	      throw(keyword_error);
	  }
	writearg(estorbitsinput.maxfringesaz);
	writearg(estorbitsinput.maxfringesrg);
	if (estorbitsinput.maxfringesaz > 500 || estorbitsinput.maxfringesrg > 500)
	  {
	    WARNING << "EO_SEARCHSPACE: Parameters (" << estorbitsinput.maxfringesaz
		    << "," << estorbitsinput.maxfringesrg << ") appear to be very high.";
	    WARNING.print();
	  }
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_REFORBIT"))
      {
        strcpy(estorbitsinput.reforbitfile, word[1]);
	writearg(estorbitsinput.reforbitfile);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_OUT_DATA"))
      {
        strcpy(estorbitsinput.foobsdata, word[1]);
	writearg(estorbitsinput.foobsdata);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_DEGREE"))
      {
	estorbitsinput.poldegree = atoi(word[1]);
	writearg(estorbitsinput.poldegree);
      }

// **********************************************************************
    else if (!strcmp(keyword,"EO_CONSTRAIN"))
      {
	int16 npar = atoi(word[1]);
	writearg(npar);
	if (npar==0)
	  estorbitsinput.constrained = false;
	else
	  {
	    estorbitsinput.constraints = matrix<int16>(npar,2);
	    for (uint i=0; i<npar; i++)
	      {
		keyword = word[2+i];
		toupper(keyword);
		int16 component;
		int16 deriv = -1;
		if (!strncmp(keyword, "BPAR", 4))
		  {
		    component = 0;
		    deriv = int(keyword[4])-48;  // char2int
		  }
		else if (!strncmp(keyword, "BPERP", 5))
		  {
		    component = 1;
		    deriv = int(keyword[5])-48;  // char2int
		  }
		if (deriv<0)
		  {
		    ERROR << "EO_CONSTRAIN: parameter "
			  << keyword
			  << " not interpretable for orbit estimation, line "
			  << linecnt << ".";
		    PRINT_ERROR(ERROR.get_str())
		      throw(keyword_error);
		  }
		else if (deriv>15)
		  {
		    ERROR << "EO_CONSTRAIN: degree of parameter "
			  << keyword
			  << " is too high, line "
			  << linecnt << ".";
		    PRINT_ERROR(ERROR.get_str())
		      throw(keyword_error);
		  }
		else
		  {
		    estorbitsinput.constraints(i,0) = component;
		    estorbitsinput.constraints(i,1) = deriv;
		  }
	      }
	  }
      }
// ____ end added by HB ____    


// **********************************************************************
// Assume wrong keyword, but if it starts with //, a c, or comment
// then continue with a warning, no blank between key and arg
// **********************************************************************
    else if (!strncmp(keyword,"COMMENT",7) ||   // assume user wanted to comment out
             !strncmp(keyword,"C",1))           // but forgot a space
      {
      WARNING << "Obsolete or unknown keyword: \"" << keyword 
           << "\" at line: " << linecnt
           << ". (Interpreted as comment, ignored)";
      WARNING.print();
      }
    // ______ Really cannot make anything from your input ______
    else
      {
      //ERROR << "Unknown keyword: \"" << keyword 
      //     << "\" at line: " << linecnt << ".";
      //PRINT_ERROR(ERROR.get_str())
      //throw(keyword_error);
      WARNING << "Unknown keyword: \"" << keyword 
           << "\" at line: " << linecnt << ".";
      WARNING.print();
      }
    //optionsfile.getline(eachline,4*ONE27,'\n');               // goto next line. [MA] needs debug if the line is long loops forever.
    } // end while (file)






// ====== Check input/ info to screen ======
  checkgeneral(generalinput, onlyprocess);              // also repair onlyprocess
  INFO << "LISTINPUT: \tAppend input to logfile: \t" << listinput;
  INFO.print();


// ______ info ELLIPSOID card, compute and fill e2, e2b ______
  ellipsinput.showdata();


// ====== Check input of step reading info from files ======
  if (generalinput.process[pr_m_readfiles])
    checkreadfiles(m_readfilesinput, MASTERID);         // check mandatory cards

  if (generalinput.process[pr_s_readfiles])
    checkreadfiles(s_readfilesinput, SLAVEID);          // check mandatory cards


// ====== Check input of step conversion slc to raw ======
  if (generalinput.process[pr_m_crop])
    {
    checkcrop(m_cropinput, MASTERID);
    if (!strcmp(s_cropinput.fileout1,m_cropinput.fileout1))
      {
      PRINT_ERROR(
      "code 301: same name outputfile CROP_OUT for master and slave not allowed.")
      throw(keyword_error);
      }
    if (generalinput.overwrit)
      if (existed(m_cropinput.fileout1))
        {
        INFO << "OVERWRIT: file " << m_cropinput.fileout1 
             << " will be overwritten.";
        INFO.print();
        }
    }

  if (generalinput.process[pr_s_crop])
    {
    checkcrop(s_cropinput, SLAVEID);
    if (!strcmp(s_cropinput.fileout1,m_cropinput.fileout1))
      {
      PRINT_ERROR(
      "code 301: same name outputfile CROP_OUT for master and slave not allowed.")
      throw(keyword_error);
      }
    if (generalinput.overwrit)
      if (existed(s_cropinput.fileout1))
        {
        INFO << "OVERWRIT: file " << s_cropinput.fileout1 
             << " will be overwritten.";
        INFO.print();
        }
    }

//____RaffaeleNutricato START MODIFICATION SECTION 10
// ====== Check input of step oversample master======
  if (generalinput.process[pr_m_oversample])
    {
    checkoversample(m_oversample, MASTERID);
    }
// ====== Check input of step oversampling slave ======
  if (generalinput.process[pr_s_oversample])
    {
    checkoversample(s_oversample, SLAVEID);
    }
//____RaffaeleNutricato END MODIFICATION SECTION 10


// ====== Check input of step porbits ======
  if (generalinput.process[pr_m_porbits])
    checkporbits(porbitsinput, MASTERID);
  if (generalinput.process[pr_s_porbits])
    checkporbits(porbitsinput, SLAVEID);


  // ____ start added by HB ____
  // ====== Check input of steps m_modorb and s_modorb ======
  if (generalinput.process[pr_m_morbits])
    {
      // ______ determine polynomial degree => number of parameters ______
      uint npar=2;
      for (int16 i=9; i>=0; i--)
	{
	if (morbitsinputmaster.coeff(i,0)!=0 || morbitsinputmaster.coeff(i,1)!=0)
	  {
	    npar = 2*(i+1);
	    break;
	  }
	}

      // ______ resize coefficient matrices ______
      morbitsinputmaster.coeff.setdata(npar/2,0,morbitsinputmaster.coeff.getdata(window(0,npar/2-1,1,1)));
      matrix<real8> temp = morbitsinputmaster.coeff.getdata(window(0,npar-1,0,0));
      morbitsinputmaster.coeff.resize(npar,1);
      morbitsinputmaster.coeff = temp; 
      for (int p=0; p<npar; p++)
	{
	  INFO << "MODORB master coefficient (" << p << "): " << morbitsinputmaster.coeff(p,0);
	  INFO.print();
	}
    }
  if (generalinput.process[pr_s_morbits])
    {
      // ______ determine polynomial degree => number of parameters ______
      uint npar=2;
      for (int16 i=9; i>=0; i--)
	if (morbitsinputslave.coeff(i,0)!=0 || morbitsinputslave.coeff(i,1)!=0)
	  {
	    npar = 2*(i+1);
	    break;
	  }

      // ______ resize coefficient matrices ______
      morbitsinputslave.coeff.setdata(npar/2,0,morbitsinputslave.coeff.getdata(window(0,npar/2-1,1,1)));
      matrix<real8> temp = morbitsinputslave.coeff.getdata(window(0,npar-1,0,0));
      morbitsinputslave.coeff.resize(npar,1);
      morbitsinputslave.coeff = temp; 
      for (int p=0; p<npar; p++)
	{
	  INFO << "MODORB slave coefficient (" << p << "): " << morbitsinputslave.coeff(p,0);
	  INFO.print();
	}
    }
  // ____ end added by HB ____


// ====== Check input of SIMAMP step master====== [MA]
if (generalinput.process[pr_m_simamp])
  {
  checksimamp(simampinput);
  }

// ====== Check input of MTIMING step master====== [MA] rev.2
if (generalinput.process[pr_m_mtiming])
  {
    if (mtiminginput.method == def_mte_method-999) // see at defaults, no method card used
      {
      mtiminginput.method = def_mte_method;           // default method
      INFO.print("Default method will be used for MASTER TIMING ERROR estimation.");
      }
    if (specified(mtiminginput.ifpositions))         // filename specified
      {
      if (mtiminginput.Nwin != def_mte_nwin+999)      // something is specified
        {
        WARNING << "MTE_NWIN: \t" << mtiminginput.Nwin
             << " ignored due to existence of input file (MTE_IN_POS) "
             << mtiminginput.ifpositions;
        WARNING.print();
        }
      mtiminginput.Nwin = filelines(mtiminginput.ifpositions);
      }
    else if (mtiminginput.Nwin == def_mte_nwin+999)   // no inputfile, default
      {
      mtiminginput.Nwin = def_mte_nwin;
      INFO.print("Default number of windows will be used for MASTER TIMING ERROR estimation.");
      }
  checkmtiming(mtiminginput);
  }

// ====== Check input of step azimuth filtering ======
  if (generalinput.process[pr_m_filtazi] ||
      generalinput.process[pr_s_filtazi])
    {
    // ______ Set defaults ______
    if (filtaziinput.overlap == -1)
      filtaziinput.overlap = filtaziinput.fftlength/8;  // default
    // ______ Check input ______
    if (generalinput.process[pr_m_filtazi] &&
        generalinput.process[pr_s_filtazi])
      checkfiltazi(filtaziinput,MASTERID+SLAVEID);
    else if (generalinput.process[pr_m_filtazi])
      checkfiltazi(filtaziinput,MASTERID);
    else if (generalinput.process[pr_s_filtazi])
      checkfiltazi(filtaziinput,SLAVEID);
    else 
      {
      PRINT_ERROR("PANIC, this cannot be")
      throw(keyword_error);
      }
    }

// ====== Check input of step range filtering ======
  if (generalinput.process[pr_m_filtrange])
    {
    // ______ Check defaults ______
    if (filtrangeinput.fftlength==-999) // not specified
      {
      if (filtrangeinput.method==rf_adaptive)
        filtrangeinput.fftlength=64;            // default
      else if (filtrangeinput.method==rf_porbits)
        filtrangeinput.fftlength=1024;          // default
      }
    checkfiltrange(filtrangeinput);
    }

// ====== Check input of reserved EXTRA step master/slave ======
  if (generalinput.process[pr_m_EXTRA])
    {
    PRINT_ERROR("extra step master not implemented.")
    throw(keyword_error);
    }
  if (generalinput.process[pr_s_EXTRA])
    {
    PRINT_ERROR("extra step slave not implemented.")
    throw(keyword_error);
    }


// ====== Check input of step coarse (orbits) ======
// ______ no checks required. (no cards as well)
  if (generalinput.process[pr_i_coarse])
    {
    ; // do nothing
    }


// ====== Check coarse coregistration based on correlation ======
// ______ Check + repair method selector coarse correlation ______
  if (generalinput.process[pr_i_coarse2])               // correlation
    {
    if (coarsecorrinput.method == def_cc_method-999) // see at defaults, no method card used
      {
      coarsecorrinput.method = def_cc_method;           // default method
      INFO.print("Default method will be used for coarse coregistration.");
      }
    if (specified(coarsecorrinput.ifpositions))         // filename specified
      {
      if (coarsecorrinput.Nwin != def_cc_nwin+999)      // something is specified
        {
        WARNING << "CC_NWIN: \t" << coarsecorrinput.Nwin
             << " ignored due to existence of input file (CC_IN_POS) "
             << coarsecorrinput.ifpositions;
        WARNING.print();
        }
      coarsecorrinput.Nwin = filelines(coarsecorrinput.ifpositions);
      }
    else if (coarsecorrinput.Nwin == def_cc_nwin+999)   // no inputfile, default
      {
      coarsecorrinput.Nwin = def_cc_nwin;
      INFO.print("Default number of windows will be used for coarse coregistration.");
      }
    checkcoarsecorr(coarsecorrinput);
    }


// ====== Check input of step fine ======
  if (generalinput.process[pr_i_fine])
    {
    if (fineinput.method == def_fc_method-999)  // see at defaults, no method card used
      {
      fineinput.method = def_fc_method;         // default method
      INFO.print("Default method will be used for fine coregistration.");
      }
    if (specified(fineinput.ifpositions))         // filename specified
      {
      if (fineinput.Nwin != def_fc_nwin+999)    // something is specified
        {
        WARNING << "FC_NWIN: \t" << fineinput.Nwin
             << " ignored due to existence of input file (FC_IN_POS) "
             << fineinput.ifpositions;
        WARNING.print();
        }
      fineinput.Nwin = filelines(fineinput.ifpositions);
      }
    else if (fineinput.Nwin == def_fc_nwin+999) // no inputfile, default
      {
      fineinput.Nwin = def_fc_nwin;
      INFO.print("Default number of windows will be used for fine coregistration.");
      }
    checkfine(fineinput);
    }

  // ====== Check input of TIMING step interferogram ====== [FvL]
  if (generalinput.process[pr_i_timing])
    {
    checkreltiming(reltiminginput);
    }

  // ====== Check input of DEMASSIST step interferogram ====== [FvL]
  if (generalinput.process[pr_i_demassist])
    {
    checkdemassist(demassistinput);
    }

  // ====== Check input step coregpm ======
  if (generalinput.process[pr_i_coregpm])
    {
    checkcoregpm(coregpminput);
    }

  // ====== Check + repair method selector resampling ======
  if (generalinput.process[pr_s_resample])              // request for process resample
    {
    if (resampleinput.method == def_rs_method-999)      // see at defaults, no method card used
      {
      resampleinput.method = def_rs_method;             // default method
      INFO.print("RS_METHOD: Using default.");
      }
    if (!priorrs_fileout)
      INFO.print("RS_OUT_FILE: Using default.");
    checkresample(resampleinput);
    }

  // ====== Check + repair method selector flatearth correction ======
  if (generalinput.process[pr_i_comprefpha])
    {
    if (comprefphainput.method == def_fe_method-999)     // see at defaults, no method card used
      {
      comprefphainput.method = def_fe_method;            // default method
      INFO.print("FE_METHOD: Using default.");
      }
    if (comprefphainput.Npoints == def_fe_Npoints-999)   // default applies
      {
      comprefphainput.Npoints = def_fe_Npoints;          // default value
      // INFO.print("FE_NPOINTS: Using default.");
      }
    // ______ if read from file, count number of points ______
    if (specified(comprefphainput.ifpositions))         // filename specified
      {
      INFO.print("FE_IN_POS: Using file to read positions.");
      comprefphainput.Npoints = filelines(comprefphainput.ifpositions);
      }
    if (comprefphainput.degree == def_fe_degree-999)     // default
      {
      comprefphainput.degree = def_fe_degree;            // default
      INFO.print("FE_DEGREE: Using default.");
      }
    checkcomprefpha(comprefphainput);
    }

  // ====== Check input of step interfero ======
  if (generalinput.process[pr_i_interfero])
    {
    checkinterfero(interferoinput);
    }

  // ====== Check input of SUBTRREFPHA step interferogram ======
  if (generalinput.process[pr_i_subtrrefpha])
    checksubtrrefpha(subtrrefphainput);

  // ====== Check input of COMPREFDEM step interferogram ======
  if (generalinput.process[pr_i_comprefdem])
    {
    checkcomprefdem(comprefdeminput);
    }

  // ====== Check input of SUBTRDEM step interferogram ======
  if (generalinput.process[pr_i_subtrrefdem])
    {
    checksubtrrefdem(subtrrefdeminput);
    }

  // ====== Check input of step coherence ======
  if (generalinput.process[pr_i_coherence])
    {
    checkcoherence(coherenceinput);
    }

  // ====== Check input of step phase filtering ======
  if (generalinput.process[pr_i_filtphase])
    {
    if (!specified(filtphaseinput.fofiltphase))  // not specified, use default
      {
      if (filtphaseinput.method==fp_goldstein || filtphaseinput.method==fp_modgoldstein)
        {
        char dummy127[ONE27];
        ostrstream omemfo(dummy127,ONE27);
        omemfo << "cint." << filtphaseinput.alpha << "filtered" << ends;
        strcpy(filtphaseinput.fofiltphase,dummy127);
        if (filtphaseinput.kernel.size()==0)
          {
          INFO.print("PF_KERNEL: Using default kernel [1 2 3 2 1].");
          filtphaseinput.kernel.resize(1,5);
          filtphaseinput.kernel(0,0) = 1./9.;
          filtphaseinput.kernel(0,1) = 2./9.;
          filtphaseinput.kernel(0,2) = 3./9.;
          filtphaseinput.kernel(0,3) = 2./9.;
          filtphaseinput.kernel(0,4) = 1./9.;
          }
        }
      else if (filtphaseinput.method==fp_spatialconv)
        {
        strcpy(filtphaseinput.fofiltphase,"cint.filtered");
        // ______ Set default kernel ______
        if (!specified(filtphaseinput.fikernel2d))      // no input file
          {
          if (filtphaseinput.kernel.size()==0)
            {
            INFO.print("PF_KERNEL: Using default kernel [1 1 1].");
            filtphaseinput.kernel.resize(1,3);
            filtphaseinput.kernel(0,0) = 1./3.;
            filtphaseinput.kernel(0,1) = 1./3.;
            filtphaseinput.kernel(0,2) = 1./3.;
            }
          }
        }
      else if (filtphaseinput.method==fp_spectral)
        {
        strcpy(filtphaseinput.fofiltphase,"cint.filtered");
        }
      else
        {
        PRINT_ERROR("Method filtphase not known")
        throw(keyword_error);
        }
      }
    checkfiltphase(filtphaseinput);
    }


  // ====== Check input of step unwrapping ======
  if (generalinput.process[pr_i_unwrap])                 // request for process unwrapping
    {
    if (unwrapinput.method == def_uw_method-999)        // see at defaults, no method card used
      {
      unwrapinput.method = def_uw_method;               // default method
      INFO.print("Default method will be used for unwrapping.");
      }
    checkunwrap(unwrapinput);
    }


  // ____ start added by HB ____
  // ====== Check input of step estorbits ======
  if (generalinput.process[pr_i_estorbits])
    { 
      // ______ check constraints______
      bool setDefaultConstraints = false;
      if(estorbitsinput.method==eo_gridsearch)
	{
	  if (estorbitsinput.poldegree != 1)
	    {
	      estorbitsinput.poldegree = 1;
	      WARNING.print("Degree of orbit error polynomial changed to 1 for method GRIDSEARCH");
	    }
	  setDefaultConstraints = true;
 	  INFO.print("Using constraints for method GRIDSEARCH: BPAR0 BPERP1");
	}
      else if (estorbitsinput.constrained==true && estorbitsinput.constraints.lines()==0) 
	{
	  setDefaultConstraints = true;
 	  INFO.print("Using default orbit constraints: BPAR0 BPERP1");
	}
      if (setDefaultConstraints)
	{
	  estorbitsinput.constraints = matrix<int16>(2,2);
	  estorbitsinput.constraints(0,0) = 0;
	  estorbitsinput.constraints(0,1) = 0;
	  estorbitsinput.constraints(1,0) = 1;
	  estorbitsinput.constraints(1,1) = 1;
	}

      // ______ check selection method of observation points ______
      if (specified(estorbitsinput.ifpositions))         // filename specified
      {
	if (estorbitsinput.nobs != 0)      // something is specified
	  {
	    WARNING << "EO_NOBS: \t" << estorbitsinput.nobs
		    << " ignored due to existence of input file (EO_IN_POS) "
		    << estorbitsinput.ifpositions;
	    WARNING.print();
	  }
	estorbitsinput.nobs = filelines(estorbitsinput.ifpositions);
      }
      else if (estorbitsinput.nobs == 0)   // no inputfile, default
	{
	  estorbitsinput.nobs = def_eo_nobs;
	  INFO << "Default number of " << def_eo_nobs 
	       << " observations will be used for orbit error estimation.";
	  INFO.print();
	}
    }


  // ====== Check input of step slant2height ======
  if (generalinput.process[pr_i_slant2h])
    {
    if (slant2hinput.method == def_s2h_method-999)
      {
      slant2hinput.method = def_s2h_method;               // default method
      INFO.print("Default method will be used for slant2h.");
      }
    if (slant2hinput.Nheights <= slant2hinput.degree1d)
      {
      WARNING << "S2H_NHEIGHTS: \t" << slant2hinput.Nheights 
           << " is too large because S2H_DEGREE1D="
           << slant2hinput.degree1d 
           << "; set S2H_NHEIGHTS = " << slant2hinput.degree1d+1
           << " (minimum)";
      WARNING.print();
      slant2hinput.Nheights = slant2hinput.degree1d + 1;
      }
    checkslant2h(slant2hinput);
    }


  // ====== Check input of step geocoding ======
  if (generalinput.process[pr_i_geocoding])
    {
    checkgeocode(geocodeinput);
    }

  // ====== Check input of dinsar step interferogram ======
  if (generalinput.process[pr_i_dinsar])
    {
    if (!strcmp(dinsarinput.topomasterresfile,generalinput.m_resfile))
      setunspecified(dinsarinput.topomasterresfile);    // used in check...
    checkdinsar(dinsarinput);
    if (!specified(dinsarinput.topomasterresfile))
      strcpy(dinsarinput.topomasterresfile,generalinput.m_resfile);
    }

  // ====== Check input of reserved EXTRA step interferogram ======
  if (generalinput.process[pr_i_EXTRA2])
    {
    PRINT_ERROR("extra step2 interferogram not implemented.")
    throw(keyword_error);
    }



// ====== Copy input to logfile (via scratchfile, update in main) ======
//   to avoid problems with not existing file it is always opened here
  ofstream scratchreadinput("scratchreadinput", ios::out | ios::trunc);
  bk_assert(scratchreadinput,"readinput: scratchreadinput",__FILE__,__LINE__);

  if (listinput)
    {
    // ______Start writing______
    linecnt=0;
    scratchreadinput << "\n----------------------------------------------------\n"
                     << "      BEGIN: copy of input (non-interpretive).\n"
                     << "-line---keyword----argument-------------comment-----\n";
    optionsfile.seekg(0,ios::beg);
    optionsfile.clear(); // clear eofbit
    optionsfile.getline(eachline,4*ONE27,'\n');
    while (optionsfile)
      {
      linecnt++;
      scratchreadinput << setw(3) << linecnt << ": " << eachline << endl;
      optionsfile.getline(eachline,4*ONE27,'\n');
      }
    scratchreadinput << "\n----------------------------------------------------\n"
                     << "      END: copy of input.\n"
                     << "----------------------------------------------------\n\n";
    DEBUG.print("Finished writing to scratchreadinput.");
    }

// ______Tidy up______
  PROGRESS.print("Interpretation inputoptionsfile finished.");
  scratchreadinput.close();
  optionsfile.close();
  } // END READINPUT



/****************************************************************
 *    checkgeneral                                              *
 *                                                              *
 * Checks general cards.                                        *
 * Repair process card                                          *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void checkgeneral(
        input_gen &generalinput,
        const int16 onlyprocess)
  {
  TRACE_FUNCTION("checkgeneral (BK 06-Sep-1999)")
  register int32 i;
  int32 cs_processcard = 0;

// ______ Repair process control if ONLYPROCESS card was present ______
  for (i=0; i<NUMPROCESSES; i++)
    cs_processcard += generalinput.process[i];
  if (onlyprocess != -1)                                     // initialized to -1;
    {
    generalinput.interactive=false;
    INFO.print("ONLYPROCESS card present, batch processing.");
    if (cs_processcard == 1)
      WARNING.print("PROCESS card ignored due to presence ONLYPROCESS card.");
    if (cs_processcard >  1)
      WARNING.print("PROCESS cards ignored due to presence ONLYPROCESS card.");
    for (i=0; i<NUMPROCESSES; i++)                      // do not process anything...
      generalinput.process[i]=0;
    generalinput.process[onlyprocess]=1;                  // exept this one.
    }

  else                                                  // check process card presence
    {
    if (cs_processcard == 0)                            // no cards
      {
      PRINT_ERROR("code 303: No (ONLY)PROCESS card present, exiting.")
      throw(keyword_error);
      }
    }

  INFO.print("\n\t*** General input cards ***");
  INFO << "MEMORY: \tAvailable to Doris [MB]: \t"
       << generalinput.memory/1e6;
  INFO.print();

  INFO << "M_RESFILE: \tResultfile for master: \t\t"
       << generalinput.m_resfile;
  INFO.print();
  INFO << "S_RESFILE: \tResultfile for slave: \t\t"
       << generalinput.s_resfile;
  INFO.print();
  INFO << "I_RESFILE: \tResultfile for products: \t"
       << generalinput.i_resfile;
  INFO.print();
  INFO << "LOGFILE: \tOut file for logging: \t\t"
       << generalinput.logfile;
  INFO.print();
  INFO << "ORB_INTERP: \tmethod selector value:\t\t"
       << generalinput.orb_interp;
  INFO.print();
  INFO << "ORB_PRM:   \torbit parameters selection value:\t\t"
       << generalinput.orb_prm;
  INFO.print();
  INFO << "DUMPBASELINE: evaluation grid for baseline: \t"
       << generalinput.dumpbaselineL << " lines x "
       << generalinput.dumpbaselineP << " pixels: ";
  INFO.print();
  INFO << "HEIGHT: \taverage terrain height:\t\t"
       << generalinput.terrain_height;
  INFO.print();
  INFO << "TIEPOINT: \tlat/lon/hei: "
       << generalinput.tiepoint.x << " " 
       << generalinput.tiepoint.y << " "
       << generalinput.tiepoint.z;
  INFO.print();

  if (!strcmp(generalinput.m_resfile,generalinput.s_resfile))
    {
    PRINT_ERROR("same name master and slave resultfile not allowed.")
    throw(keyword_error);
    }
  if (!strcmp(generalinput.m_resfile,generalinput.i_resfile))
    {
    PRINT_ERROR("same name master and interferogram resultfile not allowed.");
    throw(keyword_error);
    }
  if (!strcmp(generalinput.m_resfile,generalinput.logfile))
    {
    PRINT_ERROR("same name master resultfile and logfile not allowed.");
    throw(keyword_error);
    }
  if (!strcmp(generalinput.i_resfile,generalinput.logfile))
    {
    PRINT_ERROR("same name interferogram resultfile and logfile not allowed.");
    throw(keyword_error);
    }
  } // END checkgeneral


/****************************************************************
 *    checkreadfiles                                            *
 *                                                              *
 * Checks cards for step readfiles.                             *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void checkreadfiles(
        const input_readfiles &readfilesinput,
        const int16 id)
  {
  TRACE_FUNCTION("checkreadfiles (BK 06-Sep-1999)")
  switch (id)
    {
    case MASTERID:
      INFO.print("\n\t*** Input for step M_READFILES (master) ***");
      INFO << "M_IN_METHOD: \tmethod selected for master: \t\t"
           << readfilesinput.sensor_id;
      INFO.print();
      INFO << "M_IN_VOL:    \tVolumefile of master: \t\t"
           << readfilesinput.volfile;
      INFO.print();
      INFO << "M_IN_LEA:    \tLeaderfile of master: \t\t"
           << readfilesinput.leaderfile;
      INFO.print();
      INFO << "M_IN_NULL:   \tNullfile of master:  \t\t"
           << readfilesinput.nullfile;
      INFO.print();
      INFO << "M_IN_DAT:    \tDatfile of master:    \t\t"
           << readfilesinput.datfile;
      INFO.print();
      if (readfilesinput.sensor_id == SLC_ERS)
        {
        if (!specified(readfilesinput.volfile))
          {
          PRINT_ERROR("M_IN_VOL not defined");
          throw(keyword_error);
          }
        if (!specified(readfilesinput.leaderfile))
          {
          PRINT_ERROR("M_IN_LEA not defined");
          throw(keyword_error);
          }
        if (!specified(readfilesinput.nullfile))
          WARNING.print("M_IN_NULL not defined");
        }
      // ___ use datfile for asar input file ___
      if (!specified(readfilesinput.datfile))
        {
        PRINT_ERROR("M_IN_DAT not defined");
        throw(keyword_error);
        }
      if (readfilesinput.sensor_id == SLC_TSX)      // [MA] TSX 
        {
        if (!specified(readfilesinput.leaderfile))  // if changed, update also processor for leaderfile
          {
          PRINT_ERROR("M_IN_LEA not defined");
          throw(keyword_error);
          }
        }
      if (readfilesinput.sensor_id == SLC_RS2 || readfilesinput.sensor_id == SLC_RS2_QUAD )      // [MA] RS2 
        {
        if (!specified(readfilesinput.leaderfile))  // if changed, update also processor for leaderfile
          {
          PRINT_ERROR("M_IN_LEA not defined");
          throw(keyword_error);
          }
        }
      if (readfilesinput.sensor_id == SLC_CSK || readfilesinput.sensor_id == SLC_CSK_POL )      // [MA] CSK 
        {
        //if (!specified(readfilesinput.leaderfile))  // if changed, update also processor for leaderfile
        if (!specified(readfilesinput.datfile))  // if changed, update also processor for leaderfile
          {
          //PRINT_ERROR("M_IN_LEA not defined");
          PRINT_ERROR("M_IN_DAT not defined");
          throw(keyword_error);
          }
        }
      break;

    case SLAVEID:
      INFO.print("\n\t*** Input for step S_READFILES (slave) ***");
      INFO << "S_IN_METHOD: \tmethod selected for slave: \t\t"
           << readfilesinput.sensor_id;
      INFO.print();
      INFO << "S_IN_VOL:    \tvolumefile of slave:  \t\t"
           << readfilesinput.volfile;
      INFO.print();
      INFO << "S_IN_LEA:    \tleaderfile of slave:  \t\t"
           << readfilesinput.leaderfile;
      INFO.print();
      INFO << "S_IN_NULL:   \tnullfile of slave:   \t\t"
           << readfilesinput.nullfile;
      INFO.print();
      INFO << "S_IN_DAT:    \tdatfile of slave:     \t\t"
           << readfilesinput.datfile;
      INFO.print();
      if (readfilesinput.sensor_id == SLC_ERS)
        {
        if (!specified(readfilesinput.volfile))
          {
          PRINT_ERROR("S_IN_VOL not defined");
          throw(keyword_error);
          }
        if (!specified(readfilesinput.leaderfile))
          {
          PRINT_ERROR("S_IN_LEA not defined");
          throw(keyword_error);
          }
        if (!specified(readfilesinput.nullfile))
          WARNING.print("S_IN_NULL not defined");
        }
      // ___ use datfile for asar input file ___
      if (!specified(readfilesinput.datfile))
        {
        PRINT_ERROR("S_IN_DAT not defined");
        throw(keyword_error);
        }
      if (readfilesinput.sensor_id == SLC_TSX)      // [MA] TSX 
        {
        if (!specified(readfilesinput.leaderfile))  // if changed, update also processor.cc for the [leader]file
          {
          PRINT_ERROR("S_IN_LEA not defined");
          throw(keyword_error);
          }
        }
      if (readfilesinput.sensor_id == SLC_RS2 || readfilesinput.sensor_id == SLC_RS2_QUAD )      // [MA] RS2 
        {
        if (!specified(readfilesinput.leaderfile))  // if changed, update also processor for leaderfile
          {
          PRINT_ERROR("S_IN_LEA not defined");
          throw(keyword_error);
          }
        }
      if (readfilesinput.sensor_id == SLC_CSK || readfilesinput.sensor_id == SLC_CSK_POL )      // [MA] CSK 
        {
        //if (!specified(readfilesinput.leaderfile))  // if changed, update also processor for leaderfile
        if (!specified(readfilesinput.datfile))  // if changed, update also processor for leaderfile
          {
          PRINT_ERROR("S_IN_DAT not defined");
          throw(keyword_error);
          }
        }
      break;

    default:
      PRINT_ERROR("panic: impossible");
      throw(keyword_error);
    } // switch
  if (readfilesinput.sensor_id == SLC_ERS)
    {
    if (!strcmp(readfilesinput.volfile,readfilesinput.leaderfile))
      {
      PRINT_ERROR("same file name volume and leader file not allowed.");
      throw(keyword_error);
      }
    if (!strcmp(readfilesinput.volfile,readfilesinput.nullfile))
      {
      PRINT_ERROR("same file name volume and null file not allowed.");
      throw(keyword_error);
      }
    if (!strcmp(readfilesinput.volfile,readfilesinput.datfile))
      {
      PRINT_ERROR("same file name volume and data file not allowed.");
      throw(keyword_error);
      }
    if (!strcmp(readfilesinput.leaderfile,readfilesinput.nullfile))
      {
      PRINT_ERROR("same file name leader and null file not allowed.");
      throw(keyword_error);
      }
    if (!strcmp(readfilesinput.leaderfile,readfilesinput.datfile))
      {
      PRINT_ERROR("same file name leader and data file not allowed.");
      throw(keyword_error);
      }
    if (!strcmp(readfilesinput.nullfile,readfilesinput.datfile))
      {
      PRINT_ERROR("same file name null and data file not allowed.");
      throw(keyword_error);
      }
    }
  } // END checkreadfiles


/****************************************************************
 *    checkcrop                                                 *
 *                                                              *
 * Checks cards for step crop.                          *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void checkcrop(
        const input_crop &cropinput,
        const int16 id)
  {
  TRACE_FUNCTION("checkcrop (BK 06-Sep-1999)")
  switch (id)
    {
    case MASTERID:
      INFO.print("\n\t*** Input for step M_CROP (master) ***");
      INFO << "M_IDCROP: \tidentifier master write slc data to raster format: \t"
           << cropinput.idcrop;
      INFO.print();
      INFO << "M_CROP_IN: \tslc data inputfile for master:   \t"
           << cropinput.filein1;
      INFO.print();
      INFO << "M_CROP_OUT: \traw data outputfile for master: \t"
           << cropinput.fileout1;
      INFO.print();
      INFO << "M_DBOW: \tProcessing master line "
           << cropinput.dbow.linelo << " to "
           << cropinput.dbow.linehi << ". pixel "
           << cropinput.dbow.pixlo << " to "
           << cropinput.dbow.pixhi << ".";
      INFO.print();
      if (cropinput.dbow_geo.pixhi != 0) 
        {
        INFO.print("M_DBOW_GEO: overrides M_DBOW! processing: ");
        INFO << "center latitude " << cropinput.dbow_geo.linelo/1e6-360.0 <<
                "; center longitude " << cropinput.dbow_geo.linehi/1e6-360.0 <<
                "; height, width: " << cropinput.dbow_geo.pixlo << ", " << 
                cropinput.dbow_geo.pixhi;
        INFO.print();
        }
      break;

    case SLAVEID:
      INFO.print("\n\t*** Input for step S_CROP (slave) ***");
      INFO << "S_IDCROP: \tidentifier slave write slc data to raster format: \t"
           << cropinput.idcrop;
      INFO.print();
      INFO << "S_CROP_IN: \tslc data inputfile for slave:   \t"
           << cropinput.filein1;
      INFO.print();
      INFO << "S_CROP_OUT: \traw data outputfile for slave: \t"
           << cropinput.fileout1;
      INFO.print();
      INFO << "S_DBOW: \tProcessing slave line "
           << cropinput.dbow.linelo << " to "
           << cropinput.dbow.linehi << ". pixel "
           << cropinput.dbow.pixlo << " to "
           << cropinput.dbow.pixhi << ".";
      INFO.print();
      if (cropinput.dbow_geo.pixhi != 0) 
        {
        INFO.print("S_DBOW_GEO: overrides S_DBOW! processing: ");
        INFO << "center latitude " << cropinput.dbow_geo.linelo/1e6-360.0 <<
                "; center longitude " << cropinput.dbow_geo.linehi/1e6-360.0 <<
                "; height, width: " << cropinput.dbow_geo.pixlo << ", " << 
                cropinput.dbow_geo.pixhi;
        INFO.print();
        }
      break;

    default:
      PRINT_ERROR("panic: impossible");
      throw(keyword_error);
    }
  } // END checkcrop


//____RaffaeleNutricato START MODIFICATION SECTION 11
/****************************************************************
 *    checkoversample                                           *
 *                                                              *
 * Checks cards for step oversample.                            *
 *                                                              *
 *    Raffaele Nutricato, 12-Jan-2004                           *
 ****************************************************************/
void checkoversample(
        const input_oversample  &oversampleinput,
        const int16 id)
  {
  TRACE_FUNCTION("checkoversample (Raffaele Nutricato 12-Jan-2004)")
  switch (id)
    {
    case MASTERID:
      INFO.print("\n\t*** Input for step M_OVS (master) ***");
      INFO << "M_OVS_OUT: \t\tData output file for ovs master: "
           << oversampleinput.fileoutovs; // Output file for the oversampled master
      INFO.print(); 
      INFO << "M_OVS_FACT_RNG: \tOversampling ratio in the master range direction: "   
           << oversampleinput.OsrRange;  // Oversampling ratio in the range direction.
      INFO.print(); 
      INFO << "M_OVS_FACT_AZI: \tOversampling ratio in the master azimuth direction: "   
           << oversampleinput.OsrAzimuth;  // Oversampling ratio in the azimuth direction.
      INFO.print(); 
      INFO << "M_OVS_KERNELSIZE: \tKernel length for the master oversampling: "   
           << oversampleinput.FilterSize;  // Length of the interpolation kernel in range. 
      INFO.print(); 
      INFO << "M_OVS_OUT_FORMAT: \tOutput data format for the oversampled master: " ;  
      if (oversampleinput.oformatflag==FORMATCR4)
        INFO << "complex_real4";
      if (oversampleinput.oformatflag==FORMATCI2)
        INFO << "complex_short";
      INFO.print();
      break;

    case SLAVEID:
      INFO.print("\n\t*** Input for step S_OVS (slave) ***");
      INFO << "S_OVS_OUT: \t\tData output file for ovs slave: "
           << oversampleinput.fileoutovs; // Output file for the oversampled slave
      INFO.print();
      INFO << "S_OVS_FACT_RNG: \tOversampling ratio in the slave range direction: "   
           << oversampleinput.OsrRange;  // Oversampling ratio in the range direction.
      INFO.print(); 
      INFO << "S_OVS_FACT_AZI: \tOversampling ratio in the slave azimuth direction: "   
           << oversampleinput.OsrAzimuth;  // Oversampling ratio in the azimuth direction.
      INFO.print(); 
      INFO << "S_OVS_KERNELSIZE: \tKernel length for the slave oversampling: "   
           << oversampleinput.FilterSize;  // Length of the interpolation kernel in range. 
      INFO.print(); 
      INFO << "S_OVS_OUT_FORMAT: \tOutput data format for the oversampled slave: " ;  
      if (oversampleinput.oformatflag==FORMATCR4)
        INFO << "complex_real4";
      if (oversampleinput.oformatflag==FORMATCI2)
        INFO << "complex_short";
      INFO.print();
      break;

    default:
      PRINT_ERROR("panic: impossible");
      throw(keyword_error);
    }
  } // END checkoversample

//____RaffaeleNutricato END MODIFICATION SECTION 11


/****************************************************************
 *    checkporbits                                              *
 *                                                              *
 * Checks cards for step porbits.                               *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void checkporbits(
        const input_pr_orbits &porbitsinput,
        const int16 id)
  {
  TRACE_FUNCTION("checkporbits (BK 06-Sep-1999)")
  switch (id)
    {
    case MASTERID:
      INFO.print("\n\t*** Input for step M_PORBITS (master) ***");
      INFO << "M_ORBDIR: \tPrecise orbits master in: "
           << porbitsinput.m_orbdir;
      INFO.print();
      if (!specified(porbitsinput.m_orbdir))
        {
        PRINT_ERROR("M_ORBDIR: no directory specified.");
        throw(keyword_error);
        }
      break;

    case SLAVEID:
      INFO.print("\n\t*** Input for step S_PORBITS (slave) ***");
      INFO << "S_ORBDIR: \tPrecise orbits slave in: "
           << porbitsinput.s_orbdir;
      INFO.print();
      if (!specified(porbitsinput.s_orbdir))
        {
        PRINT_ERROR("S_ORBDIR: no directory specified.");
        throw(keyword_error);
        }
      break;

    default:
      PRINT_ERROR("panic: impossible");
      throw(keyword_error);
    }
  INFO << "ORB_INTERVAL: \ttime between ephemerides: \t"
       << porbitsinput.timeinterval;
  INFO.print();
  INFO << "ORB_EXTRATIME: \ttime before first, after last line: \t"
       << porbitsinput.timebefore;
  INFO.print();
  if (!porbitsinput.dumpmasterorbit<0)
    INFO.print("dumping masterorbit to file.");
  if (!porbitsinput.dumpslaveorbit<0)
    INFO.print("dumping slaveorbit to file.");
  } // END checkporbits


/****************************************************************
 *    checkslant2h                                              *
 *                                                              *
 * Checks cards for step slant2h.                               *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkslant2h(
        const input_slant2h &slant2hinput)
  {
  TRACE_FUNCTION("checkslant2h (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step SLANT2H ***");
  switch (slant2hinput.method)
    {
    case s2h_schwabisch:
      INFO.print("Method schwabisch is used for slant2height conversion.");

      INFO << "S2H_NPOINTS: \tNumber of points used in computation:   "
           << slant2hinput.Npoints;
      INFO.print();
      INFO << "S2H_DEGREE1D: \tDegree of 1d polynomial at each point: "
           << slant2hinput.degree1d;
      INFO.print();
      INFO << "S2H_DEGREE2D: \tDegree of 2d polynomial to compute     "
           << "coefficients of 1D polynomial: "
           << slant2hinput.degree2d;
      INFO.print();
      INFO << "S2H_NHEIGHTS: \t#heights evaluation ref. phase for 1d polynomial: "
           << slant2hinput.Nheights;
      INFO.print();
      break;

    case s2h_ambiguity:
      INFO.print("Method exact is used for slant2height conversion.");
      INFO << "S2H_OUT_LAM: \tData output file for lambda: "
           << slant2hinput.folam;
      INFO.print();
      INFO << "S2H_OUT_PHI: \tData output file for phi:    "
           << slant2hinput.fophi;
      INFO.print();
      if (!strcmp(slant2hinput.folam,slant2hinput.fophi))
        {
        PRINT_ERROR("Same filename S2H_OUT_LAM and S2H_OUT_PHI not allowed.");
        throw(keyword_error);
        }
      if (!strcmp(slant2hinput.fohei,slant2hinput.fophi))
        {
        PRINT_ERROR("Same filename S2H_OUT_HEI and S2H_OUT_PHI not allowed.");
        throw(keyword_error);
        }
      if (!strcmp(slant2hinput.folam,slant2hinput.fohei))
        {
        PRINT_ERROR("Same filename S2H_OUT_LAM and S2H_OUT_HEI not allowed.");
        throw(keyword_error);
        }
      break;

    case s2h_rodriguez:
      DEBUG.print("do some checks here as well.");
      break;

    default:
      PRINT_ERROR("impossible, unknown method s2h");
      throw(keyword_error);
    }

    INFO << "S2H_OUT_HEI: Data output file for height:   "
         << slant2hinput.fohei;
    INFO.print();
  } // END checkslant2h


/****************************************************************
 *    checkunwrap                                               *
 *                                                              *
 * Checks cards for step unwrap.                                *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkunwrap(
        const input_unwrap &unwrapinput)
  {
  TRACE_FUNCTION("checkunwrap (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step UNWRAP ***");
  INFO << "UW_OUT_FILE: \tOutput data file for unwrapped interferogram: "
       << unwrapinput.fouint;
  INFO.print();
  INFO << "UW_OUT_FORMAT: \tOutput data format for unwrapped interferogram: "
       << unwrapinput.oformatflag;
  INFO.print();
  switch (unwrapinput.method)
    {
    case uw_method1:
      INFO.print("Method 1 is used for unwrapping (treef_ramon unix calls).");
      INFO << "UW_SEEDS: ";
      if (isspace(unwrapinput.seedfile[0]))     // no file specified
        INFO << "Delta seed in line/pixel direction: "
             << unwrapinput.deltaLseed << " " << unwrapinput.deltaPseed;
      else
        INFO << " Input file with seeds for unwrapping: "
             << unwrapinput.seedfile;
      INFO.print();
      INFO << "UW_OUT_REGIONS: \tOutput data file with region numbers: "
           << unwrapinput.foregions;
      INFO.print();
      break;

    case uw_method2:
      INFO.print("Method 2: SNAPHU is used for unwrapping.");
      INFO.print("Please make sure snaphu is installed.  check results carefully.");
      INFO << "UW_SNAPHU_LOG: \tOutput log file of snaphu: "
           << unwrapinput.snaphu_log;
      INFO.print();
      INFO << "UW_SNAPHU_INIT: \tinitialization method of snaphu: "
           << unwrapinput.snaphu_init;
      INFO.print();
      INFO << "UW_SNAPHU_COH: \tinput coherence file for snaphu: "
           << unwrapinput.snaphu_coh;
      INFO.print();
      INFO << "UW_SNAPHU_MODE: \tsnaphu mode: "
           << unwrapinput.snaphu_mode;
      INFO.print();
      break;

    case uw_method3:
      INFO.print("Method 3 is used for unwrapping (?).");
      PRINT_ERROR("only for delft, standalone application.");
      throw(keyword_error);
      break;

    default:
      PRINT_ERROR("probably forgot to update this, should not be possible.");
      throw(keyword_error);
    }
  } // END checkunwrap


/****************************************************************
 *    checkgeocode                                              *
 *                                                              *
 * Checks cards for step geocode.                               *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkgeocode(
        const input_geocode &geocodeinput)
  {
  TRACE_FUNCTION("checkgeocode (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step GEOCODE ***");
  INFO << "GEO_OUT_PHI: \tOutput file for latitude: \t"
       << geocodeinput.fophi;
  INFO.print();
  INFO << "GEO_OUT_PHI: \tOutput file for longitude: \t"
       << geocodeinput.folam;
  INFO.print();
  if (!strcmp(geocodeinput.fophi,geocodeinput.folam))
    {
    ERROR << "Same file name GEO_OUT_PHI and GEO_OUT_LAM not allowed.";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  } // END checkgeocode



/****************************************************************
 *    checkcoarsecorr                                           *
 *                                                              *
 * Checks cards for step coarsecorr.                            *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkcoarsecorr(
        const input_coarsecorr &coarsecorrinput)
  {
  TRACE_FUNCTION("checkcoarsecorr (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step COARSE_CORR ***");
  switch (coarsecorrinput.method)
    {
    case cc_magfft:
      INFO.print("CC_METHOD: \tMAGFFT is used for coarse correlation."); break;
    case cc_magspace:
      INFO.print("CC_METHOD: \tMAGSPACE is used for coarse correlation."); break;
    default:
      PRINT_ERROR("panic.");
      throw(keyword_error);
    }
  if (specified(coarsecorrinput.ifpositions))
    {
    INFO << "CC_IN_POS: \tFile: " <<  coarsecorrinput.ifpositions 
         << " containing " << coarsecorrinput.Nwin
         << " positions is used.";
    INFO.print();
    }
  else                                          // no input file
    {
    INFO << "CC_NWIN: \tDistributing "
         <<  coarsecorrinput.Nwin 
         << " windows for coarse coregistration based on correlation.";
    INFO.print();
    }
  } // END checkcoarsecorr


/****************************************************************
 *    checkfine                                                 *
 *                                                              *
 * Checks cards for step fine.                                  *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkfine(
        const input_fine &fineinput)
  {
  TRACE_FUNCTION("checkfine (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step FINE ***");
  INFO << fineinput.method;
  INFO.print();
  switch (fineinput.method)
    {
//    case fc_cmplxfft:
//      INFO.print("FC_METHOD: \tCMPLXFFT is used for fine correlation.");
//      break;
//    case fc_cmplxspace:
//      INFO.print("FC_METHOD: \tCMPLXSPACE is used for fine correlation.");
//      break;
    case fc_magfft:
      INFO.print("FC_METHOD: \tMAGFFT is used for fine correlation.");
      break;
    case fc_magspace:
      INFO.print("FC_METHOD: \tMAGSPACE is used for fine correlation.");
      break;
    case fc_oversample:
      INFO.print("FC_METHOD: \tOVERSAMPLE is used for fine correlation.");
      break;
    case fc_coherence:
      INFO.print("FC_METHOD: \tCOHERENCE is used for fine correlation.");
      break;
       case fc_intensity:
     INFO.print("FC_METHOD: \tINTENSITY is used for fine correlation.");
      break;
    default:
      PRINT_ERROR("panic.");
      throw(keyword_error);
    }

  if (specified(fineinput.ifpositions))
    {
    INFO << "FC_IN_POS: File: "
         <<  fineinput.ifpositions 
         << " containing "
         << fineinput.Nwin << " positions is used.";
    INFO.print();
    }
  else                                          // no input file
    {
    INFO << "FC_NWIN: \tDistributing "
         <<  fineinput.Nwin 
         << " windows for fine coregistration.";
    INFO.print();
    }
  if (fineinput.plotoffsets)
    {
    INFO << "FC_PLOT " << fineinput.plotthreshold 
         << ": Plotting estimated offsets with correlation > "
         << fineinput.plotthreshold;
    INFO.print();
    }
  else
    {
    WARNING.print("It is highly recommended to use FC_PLOT to plot the computed FINE offset vectors.");
    }
  if (fineinput.plotmagbg)
    INFO.print("FC_PLOT: magnitude will be in background of plot.");

  INFO << "FC_WINSIZE for correlation: (l,p) = ("
       <<  fineinput.MasksizeL << ", "
       <<  fineinput.MasksizeP << ").";
  INFO.print();
  INFO.print("Max. estimable offset is half the window size");
  INFO.print("If the offset varies a lot over the image, use larger winsize");
  INFO << "FC_WINSIZE for oversampling (2*Acc): (l,p) = ("
       <<  2*fineinput.AccL << ", "
       <<  2*fineinput.AccP << ").";
  INFO.print();
  INFO << "FC_OSFACTOR: "
       <<  fineinput.osfactor
       <<  ". Oversampling factor for fine coregistration.";
  INFO.print();
  if (!ispower2(fineinput.osfactor))
    {
    ERROR << " no power of 2.";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  } // END checkfine


/****************************************************************
 *    checksimamp                                               *
 *                                                              *
 * Checks cards for step simamp.                                *
 *                                                              *
 *    Mahmut Arikan, 30-Oct-2008                                *
 ****************************************************************/
void checksimamp(
        const input_simamp &simampinput)
  {
  TRACE_FUNCTION("checksimamp (MA 30-oct-2008)")
  INFO.print("\n\t*** Input for step SIMAMP ***");
  INFO << "SAM_IN_DEM:    \t" << simampinput.firefdem;
  INFO.print();
  if (specified(simampinput.fothetalp))                              // MA SPT
    {
    INFO << "SAM_OUT_THETA_LP: \t" << simampinput.fothetalp
       << "; output requested of incidence angle THETA [rad] in radar coordinates.";
    INFO.print();
    }
  if (specified(simampinput.fodemlp))                              // MA SPT
    {
    INFO << "SAM_OUT_DEM_LP: \t" << simampinput.fodemlp
       << "; output requested of DEM [m] in radar coordinates.";
    INFO.print();
    }
  if (specified(simampinput.fodem))
    {
    INFO << "SAM_OUT_DEM:   \t" << simampinput.fodem
       << "; output requested of input DEM.";
    INFO.print();
    if (!strcmp(simampinput.fosimamp,simampinput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
      throw(keyword_error);
      }
    }
  if (specified(simampinput.fosimamp))
    {
    INFO << "SAM_OUT_FILE:   \t" << simampinput.fosimamp
         << "; output requested of simulated amplitude.";
    INFO.print();
    if (!strcmp(simampinput.fosimamp,simampinput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
      throw(keyword_error);
      }
    }
  INFO << "SAM_IN_SIZE:   \t" << simampinput.demrows
       << " " << simampinput.demcols
       << "; number of rows (latitude), columns (lon) in DEM.";
  INFO.print();
  INFO << "SAM_IN_UL:     \t" 
       << rad2deg(simampinput.demlatleftupper) << " "
       << rad2deg(simampinput.demlonleftupper)
       << "; coordinates of upper left corner (first row/col).";
  INFO.print();
  INFO << "SAM_IN_DELTA:  \t"
       << rad2deg(simampinput.demdeltalat) << " "
       << rad2deg(simampinput.demdeltalon);
  INFO.print();
  INFO << "SAM_IN_NODATA:  \t" << simampinput.demnodata
       << "; this number in DEM will be set to 0 reference phase.";
  INFO.print();

  INFO << "SAM_IN_FORMAT: \tinput format DEM: \t";
  switch (simampinput.iformatflag)
    {
    case FORMATR4:
      INFO << "real4.";
      break;
    case FORMATR8:
      INFO << "real8.";
      break;
    case FORMATI2:
      INFO << "signed short.";
      break;
    case FORMATI2_BIGENDIAN:
      INFO << "signed short big endian.";
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  INFO.print();



// ______ Check some things ______
  if (!existed(simampinput.firefdem))
    {
    ERROR << "SAM_IN_DEM:   \t" << simampinput.firefdem
         << " can not be opened.";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(simampinput.demdeltalat)<.0002) // [MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
    {
    WARNING << "SAM_IN_DELTA: \t" << rad2deg(simampinput.demdeltalat)
         << " [deg] seems very small (input in decimal degrees).";
    WARNING.print();
    }
  if (simampinput.demrows<1 || simampinput.demrows>100000)
    {
    WARNING << "SAM_DEM_SIZE: numrows: \t" << simampinput.demrows
         << " seems to be wrong.";
    WARNING.print();
    }
  if (simampinput.demcols<1 || simampinput.demcols>100000)
    {
    WARNING << "SAM_DEM_SIZE: numcols: \t" << simampinput.demcols
         << " seems to be wrong.";
    WARNING.print();
    }
  if (rad2deg(simampinput.demdeltalon)<.0002)  
    {
    WARNING << "SAM_IN_DELTA: \t" << simampinput.demdeltalon
         << " seems very small (it should be in decimal degrees).";
    WARNING.print();
    }
  if (rad2deg(simampinput.demlatleftupper) < -90. ||
      rad2deg(simampinput.demlatleftupper) >  90.   )
    {
    ERROR << "SAM_IN_LU:    \t" << rad2deg(simampinput.demlatleftupper)
         << " out of range (-90:90).";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(simampinput.demlonleftupper) < -180. ||
      rad2deg(simampinput.demlonleftupper) >  180.   )
    {
    WARNING << "SAM_IN_LU:    \t" << rad2deg(simampinput.demlonleftupper)
         << " out of range (-180:180).";
    WARNING.print();
    }
  } // END simamp  [MA]



/****************************************************************
 *    checkmtiming                                              *
 *                                                              *
 * Checks cards for step simamp corr.                           *
 *                                                              *
 *    Mahmut Arikan, 11-Nov-2008                                *
 ****************************************************************/
void checkmtiming(
        const input_mtiming &mtiminginput)
  {
  TRACE_FUNCTION("checkmtiming (MA, BO 11-Nov-2008)")
  INFO.print("\n\t*** Input for step MASTER TIMING ***");
  switch (mtiminginput.method)
    {
    case cc_magfft:
      INFO.print("MTE_METHOD: \tMAGFFT is used for MASTER TIMING ERROR estimation."); break;
    case cc_magspace:
      INFO.print("MTE_METHOD: \tMAGSPACE is used for MASTER TIMING ERROR estimation."); break;
    default:
      PRINT_ERROR("panic.");
      throw(keyword_error);
    }
  if (specified(mtiminginput.ifpositions))
    {
    INFO << "MTE_IN_POS: \tFile: " <<  mtiminginput.ifpositions 
         << " containing " << mtiminginput.Nwin
         << " positions is used.";
    INFO.print();
    }
  else                                          // no input file
    {
    INFO << "MTE_NWIN: \tDistributing "
         <<  mtiminginput.Nwin 
         << " windows for master timing error estimation (correlation).";
    INFO.print();
    }
  } // END checkmtiming


// ____end added by MA ____


// ____start added by FvL ____

/****************************************************************
 *    checkreltiming                                            *
 *                                                              *
 * Checks cards for step timing                                 *
 *                                                              *
 *    Freek van Leijen, 22-Aug-2007                             *
 ****************************************************************/
void checkreltiming(
        const input_reltiming &timinginput)
  {
  TRACE_FUNCTION("checkreltiming (FvL 22-aug-2007)")
  INFO.print("\n\t*** Input for step REL_TIMING ***");
  INFO << "RTE_THRESHOLD: \tThreshold correlation for model: \t"
       <<  timinginput.threshold;
  INFO.print();
  INFO << "RTE_MAXITER: \tNumber of points to remove (max): \t"
       <<  timinginput.maxiter;
  INFO.print();
  INFO << "RTE_K_ALPHA: \tCritical value for outlier test: \t"
       <<  timinginput.k_alpha;
  INFO.print();
  } // END checkreltiming


/****************************************************************
 *    checkdemassist                                            *
 *                                                              *
 * Checks cards for step demassist.                             *
 *                                                              *
 *    Freek van Leijen, 22-Aug-2007                             *
 ****************************************************************/
void checkdemassist(
        const input_demassist &demassistinput)
  {
  TRACE_FUNCTION("checkdemassist (FvL 22-aug-2007)")
  INFO.print("\n\t*** Input for step DEMASSIST ***");
  INFO << "DAC_IN_DEM:    \t" << demassistinput.firefdem;
  INFO.print();
  if (specified(demassistinput.forefdemhei))
    {
    INFO << "DAC_OUT_DEM_LP: \t" << demassistinput.forefdemhei
       << "; output requested of DEM [m] in radar coordinates.";
    INFO.print();
    }
  if (specified(demassistinput.fodem))
    {
    INFO << "DAC_OUT_DEM:   \t" << demassistinput.fodem
       << "; output requested of input DEM.";
    INFO.print();
    if (!strcmp(demassistinput.fodemi,demassistinput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
      throw(keyword_error);
      }
    }
  if (specified(demassistinput.fodemi))
    {
    INFO << "DAC_OUT_DEMI:   \t" << demassistinput.fodemi
         << "; output requested of interpolated DEM.";
    INFO.print();
    if (!strcmp(demassistinput.fodemi,demassistinput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
      throw(keyword_error);
      }
    }
  INFO << "DAC_IN_SIZE:   \t" << demassistinput.demrows
       << " " << demassistinput.demcols
       << "; number of rows (latitude), columns (lon) in DEM.";
  INFO.print();
  INFO << "DAC_IN_UL:     \t" 
       << rad2deg(demassistinput.demlatleftupper) << " "
       << rad2deg(demassistinput.demlonleftupper)
       << "; coordinates of upper left corner (first row/col).";
  INFO.print();
  INFO << "DAC_IN_DELTA:  \t"
       << rad2deg(demassistinput.demdeltalat) << " "
       << rad2deg(demassistinput.demdeltalon);
  INFO.print();
  INFO << "DAC_IN_NODATA:  \t" << demassistinput.demnodata
       << "; this number in DEM will be set to 0 reference phase.";
  INFO.print();

  INFO << "DAC_IN_FORMAT: \tinput format DEM: \t";
  switch (demassistinput.iformatflag)
    {
    case FORMATR4:
      INFO << "real4.";
      break;
    case FORMATR8:
      INFO << "real8.";
      break;
    case FORMATI2:
      INFO << "signed short.";
      break;
    case FORMATI2_BIGENDIAN:
      INFO << "signed short big endian.";
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  INFO.print();



// ______ Check some things ______
  if (!existed(demassistinput.firefdem))
    {
    ERROR << "DAC_IN_DEM:   \t" << demassistinput.firefdem
         << " can not be opened.";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(demassistinput.demdeltalat)<.0002) //[MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
    {
    WARNING << "DAC_IN_DELTA: \t" << rad2deg(demassistinput.demdeltalat)
         << " [deg] seems very small (input in decimal degrees).";
    WARNING.print();
    }
  if (demassistinput.demrows<1 || demassistinput.demrows>100000)
    {
    WARNING << "DAC_DEM_SIZE: numrows: \t" << demassistinput.demrows
         << " seems to be wrong.";
    WARNING.print();
    }
  if (demassistinput.demcols<1 || demassistinput.demcols>100000)
    {
    WARNING << "DAC_DEM_SIZE: numcols: \t" << demassistinput.demcols
         << " seems to be wrong.";
    WARNING.print();
    }
  if (rad2deg(demassistinput.demdeltalon)<.0002)
    {
    WARNING << "DAC_IN_DELTA: \t" << demassistinput.demdeltalon
         << " seems very small (it should be in decimal degrees).";
    WARNING.print();
    }
  if (rad2deg(demassistinput.demlatleftupper) < -90. ||
      rad2deg(demassistinput.demlatleftupper) >  90.   )
    {
    ERROR << "DAC_IN_LU:    \t" << rad2deg(demassistinput.demlatleftupper)
         << " out of range (-90:90).";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(demassistinput.demlonleftupper) < -180. ||
      rad2deg(demassistinput.demlonleftupper) >  180.   )
    {
    WARNING << "DAC_IN_LU:    \t" << rad2deg(demassistinput.demlonleftupper)
         << " out of range (-180:180).";
    WARNING.print();
    }
  } // END demassist

// ____end added by FvL ____


/****************************************************************
 *    checkcoregpm                                              *
 *                                                              *
 * Checks cards for step coregpm.                               *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkcoregpm(
        const input_coregpm &coregpminput)
  {
  TRACE_FUNCTION("checkcoregpm (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step COREGPM ***");
  INFO << "CPM_THRESHOLD: \tThreshold correlation for model: \t"
       <<  coregpminput.threshold;
  INFO.print();
  INFO << "CPM_DEGREE: \tPolynomial for coregistration: \t"
       <<  coregpminput.degree;
  INFO.print();
  INFO << "CPM_WEIGHT: \tData weighting option: \t"
       <<  coregpminput.weightflag
       << " (0 none, 1 linear, 2 quadratic, 3 bamler paper) is used.";
  INFO.print();
  INFO << "CPM_MAXITER: \tNumber of points to remove (max): \t"
       <<  coregpminput.maxiter;
  INFO.print();
  INFO << "CPM_K_ALPHA: \tCritical value for outlier test: \t"
       <<  coregpminput.k_alpha;
  INFO.print();
  if (coregpminput.dumpmodel)
    INFO.print("CPM_DUMP: dumping model to files (see INFO).");
  if (coregpminput.plot)
    INFO.print("CPM_PLOT: error vectors etc. will be plotted.");
  else
    WARNING.print("It is higly recommended to use CPM_PLOT to review the estimated model.");
  if (coregpminput.plotmagbg)
    INFO.print("CPM_PLOT: magnitude will be in background.");
  } // END checkcoregpm


/****************************************************************
 *    checkcomprefpha                                           *
 *                                                              *
 * Checks cards for step comprefpha.                            *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkcomprefpha(
        const input_comprefpha &comprefphainput)
  {
  TRACE_FUNCTION("checkcomprefpha (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step COMPREFPHA ***");
  switch (comprefphainput.method)
    {
    case fe_porbits:
      INFO.print("FE_METHOD: method for flatearth correction: PORBITS");
      // ______ Check points from file or random distributed ______
      if (specified(comprefphainput.ifpositions))
        {
        INFO << "FE_IN_POS: file: "
             <<  comprefphainput.ifpositions 
             << " containing "
             << comprefphainput.Npoints << " positions used for ref. phase estimation.";
        INFO.print();
        }
      else                                              // no input file
        {
        INFO << "FE_NPOINTS: Using " << comprefphainput.Npoints
             << " (random like) distributed points for estimation of refpha polynomial.";
        INFO.print();
        }
      if (comprefphainput.Npoints > 5000)
        WARNING.print("FE_NPOINTS: Too many points requested?");
      INFO << "FE_DEGREE: Using " << comprefphainput.degree
           << " order polynomial for flat Earth correction.";
      INFO.print();
      if (comprefphainput.degree > 10)
        WARNING.print("FE_DEGREE: degree > 10?");
      break;
    case fe_method2:
      INFO.print("FE_METHOD: method for flatearth correction: method2 :NOT IMPLEMENTED");
      break;
    default:
      PRINT_ERROR("impossible, method name is checked while reading cards.");
      throw(keyword_error);
    }
  } // END checkcomprefpha


/****************************************************************
 *    checksubtrrefpha                                          *
 *                                                              *
 * Checks cards for step subtrrefpha.                           *
 *                                                              *
 *    Bert Kampes, 09-Feb-2000                                  *
 ****************************************************************/
void checksubtrrefpha(
        const input_subtrrefpha &subtrrefphainput)
  {
  TRACE_FUNCTION("checksubtrrefpha (BK 09-Feb-2000)")
  INFO.print("\n\t*** Input for step SUBTRREFPHA ***");
  // ______ Print info ______
  switch (subtrrefphainput.method)
    {
    case srp_polynomial:
      INFO.print("SRP_METHOD: \tpolynomial: \tPolynomial from COMPREFPHA used.");
      break;
    case srp_exact:
      INFO.print("SRP_METHOD: \texact:      \treference phase computed foreach pixel.");
      break;
    default:
      PRINT_ERROR("impossible, checked above.");
      throw(keyword_error);
    }
  if (subtrrefphainput.dumponlyrefpha==false)
    {
    INFO << "SRP_OUT_CINT: \toutputfile complex interferogram: \t"
       << subtrrefphainput.focint;
    INFO.print();
    }
  INFO << "SRP_MULTILOOK: \tFactors (line pixel): \t" 
       << subtrrefphainput.multilookL << " "
       << subtrrefphainput.multilookP;
  INFO.print();
  if (subtrrefphainput.dumponlyrefpha==true)
    {
    WARNING.print("SRP_DUMPREFPHA: Only dumping reference phase, no subtraction.");
    INFO << "SRP_DUMPREFPHA: only dump refphase: "
         << subtrrefphainput.dumponlyrefpha;
    INFO.print();
    INFO << "SRP_OUT_REFPHA: Output file reference phase: "
         << subtrrefphainput.forefpha;
    INFO.print();
    }
    // __________ added by FvL
  if (specified(subtrrefphainput.foh2ph))
    {
    INFO << "SRP_OUT_H2PH:   \t" << subtrrefphainput.foh2ph
         << "; output requested of height-to-phase constant per resolution cell.";
    INFO.print();
    }
  // ____________ end added by FvL

  // ______ Check ______ // [MA] increased from 100 to 1000
  if (subtrrefphainput.multilookL > 1000 || subtrrefphainput.multilookP > 1000 ||
      subtrrefphainput.multilookL < 1   || subtrrefphainput.multilookP < 1 )
    WARNING.print("SRP_MULTILOOK: multilookfactor seems very high.");
  } // END checksubtrrefpha


/****************************************************************
 *    checkresample                                             *
 *                                                              *
 * Checks cards for step resample.                              *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkresample(
        const input_resample &resampleinput)
  {
  TRACE_FUNCTION("checkresample (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step RESAMPLE ***");
  INFO << "RS_OUT_FILE: \toutput filename: \t\t" 
       << resampleinput.fileout;
  INFO.print();
  INFO << "RS_OUT_FORMAT: output format: \t\t";
  switch (resampleinput.oformatflag)
    {
    case FORMATCR4:
      INFO << "complex_real4.";
      break;
    case FORMATCI2:
      INFO << "complex_short.";
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  INFO.print();
  if (resampleinput.dbow.linehi != 0 || resampleinput.dbow.pixhi != 0)
    {
    INFO << "RS_DBOW: \tOutput window: \t\t\t" 
         << resampleinput.dbow.linelo << " " 
         << resampleinput.dbow.linehi << " " 
         << resampleinput.dbow.pixlo  << " " 
         << resampleinput.dbow.pixhi ;
    INFO.print();
    }
  INFO << "RS_SHIFTAZI: \tshift azimuth spectrum: \t" 
       << resampleinput.shiftazi ;
  INFO.print();
  } // END checkresample


/****************************************************************
 *    checkinterfero                                            *
 *                                                              *
 * Checks cards for step interfero.                             *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkinterfero(
        const input_interfero &interferoinput)
  {
  TRACE_FUNCTION("checkinterfero (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step INTERFERO ***");
  bool filespecified = false;
  bool samefilename  = false;
  if (specified(interferoinput.foint))
    {
    filespecified = true;
    INFO << "INT_OUT_INT: \tOutputfile interferogram: \t"
         <<  interferoinput.foint;
    INFO.print();
    if (!(strcmp(interferoinput.foint,interferoinput.focint)))
      samefilename = true;
    }
  if (specified(interferoinput.focint))
    {
    filespecified = true;
    INFO << "INT_OUT_CINT: Outfile complex interferogram: \t"
         <<  interferoinput.focint;
    INFO.print();
    if (!(strcmp(interferoinput.focint,interferoinput.foint)))
      samefilename = true;
    }
//  if (strcmp(interferoinput.foflatearth," "))         // something is specified
//    {
//    filespecified = true;
//    INFO << "INT_OUT_FE: data outputfile reference phase: \t"
//         <<  interferoinput.foflatearth;
//    INFO.print();
//    if (!(strcmp(interferoinput.foflatearth,interferoinput.foint))  ||
//        !(strcmp(interferoinput.foflatearth,interferoinput.focint)))
//      samefilename = true;
//    }

  INFO << "INT_MULTILOOK: \tFactor (line pixel): \t"
       <<  interferoinput.multilookL << " " 
       << interferoinput.multilookP;
  INFO.print();
  // [MA] increased from 100 to 1000
  if (interferoinput.multilookL > 1000 || interferoinput.multilookP > 1000)
    {
    PRINT_ERROR("code ???: INT_MULTILOOK: > 1000.");
    throw(keyword_error);
    }
  if (interferoinput.multilookL < 1 || interferoinput.multilookP < 1)
    {
    PRINT_ERROR("code ???: INT_MULTILOOK: < 1.");
    throw(keyword_error);
    }

  if (!filespecified)
    {
    PRINT_ERROR("code ???: INT_OUT_*: at least one output file must be specified.");
    throw(keyword_error);
    }
  if (samefilename)
    {
    PRINT_ERROR("code ???: INT_OUT_*: same name output files.");
    throw(keyword_error);
    }
  } // END checkinterfero


/****************************************************************
 *    checkcoherence                                            *
 *                                                              *
 * Checks cards for step coherence.                             *
 *                                                              *
 *    Bert Kampes, 29-Sep-1999                                  *
 ****************************************************************/
void checkcoherence(
        const input_coherence &coherenceinput)
  {
  TRACE_FUNCTION("checkcoherence (BK 29-Sep-1999)")
  INFO.print("\n\t*** Input for step COHERENCE ***");
  bool filespecified = false;
  bool samefilename  = false;
  if (specified(coherenceinput.foccoh))
    {
    filespecified = true;
    INFO << "COH_OUT_CCOH: \tOutfile complex coherence: \t"
         <<  coherenceinput.foccoh;
    INFO.print();
    if (!(strcmp(coherenceinput.foccoh,coherenceinput.focoh)))
      samefilename = true;
    }

  if (specified(coherenceinput.focoh))
    {
    filespecified = true;
    INFO << "COH_OUT_COH: \tOutputfile coherence image: "
         <<  coherenceinput.focoh;
    INFO.print();
    if (!(strcmp(coherenceinput.focoh,coherenceinput.foccoh)))
      samefilename = true;
    }

  INFO << "COH_MULTILOOK: \tFactor (line pixel): \t"
       <<  coherenceinput.multilookL << " "
       << coherenceinput.multilookP;
  INFO.print();
  INFO << "COH_WINSIZE: \t window size coh. estimation (l/p): \t"
       <<  coherenceinput.cohsizeL << " " << coherenceinput.cohsizeP;
  INFO.print();
   // [MA] increased from 100 to 1000
  if (coherenceinput.multilookL > 1000 || coherenceinput.multilookP > 1000)
    {
    PRINT_ERROR("code ???: COH_MULTILOOK: > 1000.");
    throw(keyword_error);
    }
  if (coherenceinput.multilookL < 1 || coherenceinput.multilookP < 1)
    {
    PRINT_ERROR("code ???: COH_MULTILOOK: < 1.");
    throw(keyword_error);
    }
  if (coherenceinput.cohsizeL > 500 || coherenceinput.cohsizeP > 500)
    {
    PRINT_ERROR("code ???: COH_WINSIZE: > 500.");
    throw(keyword_error);
    }
  if (coherenceinput.cohsizeL < 1 || coherenceinput.cohsizeP < 1)
    {
    PRINT_ERROR("code ???: COH_WINSIZE: < 1.");
    throw(keyword_error);
    }

  if (!filespecified)
    {
    PRINT_ERROR("code ???: COH_OUT_*: at least one output file must be specified.");
    throw(keyword_error);
    }
  if (samefilename)
    {
    PRINT_ERROR("code ???: COH_OUT_*: same name output files.");
    throw(keyword_error);
    }
  } // END checkcoherence



/****************************************************************
 *    checkcomprefdem                                           *
 *                                                              *
 * Checks cards for step comprefdem.                            *
 *                                                              *
 *    Bert Kampes, 14-Feb-2000                                  *
 ****************************************************************/
void checkcomprefdem(
        const input_comprefdem &comprefdeminput)
  {
  TRACE_FUNCTION("checkcomprefdem (BK 14-Feb-2000)")
  INFO.print("\n\t*** Input for step COMPREFDEM ***");
//  switch (comprefdeminput.method)
//    {
//    case crd_nearest:
//      INFO.print("NEAREST_NEIGHBOR, use DENSE=2 or so.");
//      break;
//    case crd_trilinear:
//      INFO.print("TRI_LINEAR; use DENSE=0.2 for speed.");
//      break;
//    default:
//      PRINT_ERROR("totally impossible, checked input.");
//      throw(keyword_error);
//    }
  INFO << "CRD_IN_DEM:    \t" << comprefdeminput.firefdem;
  INFO.print();
  INFO << "CRD_OUT_FILE:  \t" << comprefdeminput.forefdem;
  INFO.print();
  if (specified(comprefdeminput.forefdemhei))
    {
    INFO << "CRD_OUT_DEM_LP: \t" << comprefdeminput.forefdemhei
       << "; output requested of DEM [m] in radar coordinates.";
    INFO.print();
    if (!strcmp(comprefdeminput.forefdem,comprefdeminput.forefdemhei))
      {
      PRINT_ERROR("CRD_OUT_FILE, CRD_OUT_DEM_LP: Same output file name.");
      throw(keyword_error);
      }
    }
  if (specified(comprefdeminput.fodem))
    {
    INFO << "CRD_OUT_DEM:   \t" << comprefdeminput.fodem
       << "; output requested of input DEM.";
    INFO.print();
    if (!strcmp(comprefdeminput.fodemi,comprefdeminput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
      throw(keyword_error);
      }
    }
  if (specified(comprefdeminput.fodemi))
    {
    INFO << "CRD_OUT_DEMI:   \t" << comprefdeminput.fodemi
         << "; output requested of interpolated DEM.";
    INFO.print();
    if (!strcmp(comprefdeminput.fodemi,comprefdeminput.fodem))
      {
      PRINT_ERROR("OUT_DEM, OUT_DEMI: Same output file name.");
      throw(keyword_error);
      }
    }
    // __________ added by FvL
  if (specified(comprefdeminput.foh2ph))
    {
    INFO << "CRD_OUT_H2PH:   \t" << comprefdeminput.foh2ph
         << "; output requested of height-to-phase constant per resolution cell.";
    INFO.print();
    }
  // ____________ end added by FvL
  INFO << "CRD_IN_SIZE:   \t" << comprefdeminput.demrows
       << " " << comprefdeminput.demcols
       << "; number of rows (latitude), columns (lon) in DEM.";
  INFO.print();
  INFO << "CRD_IN_UL:     \t" 
       << rad2deg(comprefdeminput.demlatleftupper) << " "
       << rad2deg(comprefdeminput.demlonleftupper)
       << "; coordinates of upper left corner (first row/col).";
  INFO.print();
  INFO << "CRD_IN_DELTA:  \t"
       << rad2deg(comprefdeminput.demdeltalat) << " "
       << rad2deg(comprefdeminput.demdeltalon);
  INFO.print();
  INFO << "CRD_IN_NODATA:  \t" << comprefdeminput.demnodata
       << "; this number in DEM will be set to 0 reference phase.";
  INFO.print();
//  INFO << "CRD_DENSE:      \t" << comprefdeminput.extradense
//       << "; this is the factor for oversampling DEM more than minimum.";
//  INFO.print();
  if (comprefdeminput.includerefpha)
    INFO << "CRD_INCLUDE_FE: \tref. dem is computed including flat earth.";
  else
    INFO << "CRD_INCLUDE_FE: \tref. dem is computed w.r.t. ellipsoid (topo only).";
  INFO.print();

  INFO << "CRD_IN_FORMAT: \tinput format DEM: \t";
  switch (comprefdeminput.iformatflag)
    {
    case FORMATR4:
      INFO << "real4.";
      break;
    case FORMATR8:
      INFO << "real8.";
      break;
    case FORMATI2:
      INFO << "signed short.";
      break;
    case FORMATI2_BIGENDIAN:
      INFO << "signed short big endian.";
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  INFO.print();



// ______ Check some things ______
  if (!existed(comprefdeminput.firefdem))
    {
    ERROR << "CRD_IN_DEM:   \t" << comprefdeminput.firefdem
         << " can not be opened.";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(comprefdeminput.demdeltalat)<.0002) //[MA] 1 arc secs = 1/3600=0.00027 : allow finer resolutions
    {
    WARNING << "CRD_IN_DELTA: \t" << rad2deg(comprefdeminput.demdeltalat)
         << " [deg] seems very small (input in decimal degrees).";
    WARNING.print();
    }
  if (comprefdeminput.demrows<1 || comprefdeminput.demrows>100000)
    {
    WARNING << "CRD_DEM_SIZE: numrows: \t" << comprefdeminput.demrows
         << " seems to be wrong.";
    WARNING.print();
    }
  if (comprefdeminput.demcols<1 || comprefdeminput.demcols>100000)
    {
    WARNING << "CRD_DEM_SIZE: numcols: \t" << comprefdeminput.demcols
         << " seems to be wrong.";
    WARNING.print();
    }
//  if (comprefdeminput.extradense>5.0)
//    {
//    WARNING << "CRD_DENSE:    \t" << comprefdeminput.extradense
//         << " seems to be quite large.";
//    WARNING.print();
//    }
//  if (comprefdeminput.extradense<0.2)
//    {
//    WARNING << "CRD_DENSE:    \t" << comprefdeminput.extradense
//         << " seems too small.";
//    WARNING.print();
//    }
  if (rad2deg(comprefdeminput.demdeltalon)<.0002)
    {
    WARNING << "CRD_IN_DELTA: \t" << comprefdeminput.demdeltalon
         << " seems very small (it should be in decimal degrees).";
    WARNING.print();
    }
  if (rad2deg(comprefdeminput.demlatleftupper) < -90. ||
      rad2deg(comprefdeminput.demlatleftupper) >  90.   )
    {
    ERROR << "CRD_IN_LU:    \t" << rad2deg(comprefdeminput.demlatleftupper)
         << " out of range (-90:90).";
    PRINT_ERROR(ERROR.get_str())
    throw(keyword_error);
    }
  if (rad2deg(comprefdeminput.demlonleftupper) < -180. ||
      rad2deg(comprefdeminput.demlonleftupper) >  180.   )
    {
    WARNING << "CRD_IN_LU:    \t" << rad2deg(comprefdeminput.demlonleftupper)
         << " out of range (-180:180).";
    WARNING.print();
    }
  if (!strcmp(comprefdeminput.fodem,comprefdeminput.forefdem))
    {
    PRINT_ERROR("OUT_DEM, OUT_FILE: Same output file name.");
    throw(keyword_error);
    }
  if (!strcmp(comprefdeminput.firefdem,comprefdeminput.forefdem))
    {
    PRINT_ERROR("OUT_FILE, IN_DEM: Same file name.");
    throw(keyword_error);
    }
  } // END comprefdem



/****************************************************************
 *    checksubtrrefdem                                          *
 *                                                              *
 * Checks cards for step subtrrefdem.                           *
 *                                                              *
 *    Bert Kampes, 14-Feb-2000                                  *
 ****************************************************************/
void checksubtrrefdem(
        const input_subtrrefdem &subtrrefdeminput)
  {
  TRACE_FUNCTION("checksubtrrefdem (BK 14-Feb-2000)")
  INFO.print("\n\t*** Input for step SUBTRREFDEM ***");
  INFO << "SRD_OUT_CINT:   \t" << subtrrefdeminput.focint;
  INFO.print();
  INFO << "SRD_OFFSET:     \t" << subtrrefdeminput.offsetL 
       << " " << subtrrefdeminput.offsetP;
  INFO.print();
  if (abs(subtrrefdeminput.offsetL)>5) 
    WARNING.print("Apply offset in azimuth larger than 5 pixels?");
  if (abs(subtrrefdeminput.offsetP)>5) 
    WARNING.print("Apply offset in range larger than 5 pixels?");
  } // END checksubtrrefdem


/****************************************************************
 *    checkfiltrange                                            *
 * Checks cards for step filtrange.                             *
 *    Bert Kampes, 14-Feb-2000                                  *
 ****************************************************************/
void checkfiltrange(
        const input_filtrange &filtrangeinput)
  {
  TRACE_FUNCTION("checkfiltrange (BK 14-Feb-2000)")
  INFO.print("\n\t*** Input for step FILTRANGE ***");
  // ______ Give info ______
  switch (filtrangeinput.method)
    {
    case rf_adaptive:
      INFO.print("RF_METHOD:        ADAPTIVE \t(estimate fringe freq.)");
      INFO << "RF_NLMEAN:     " << filtrangeinput.nlmean;
      INFO.print();
      INFO << "RF_THRESHOLD:  " << filtrangeinput.SNRthreshold;
      INFO.print();
      INFO << "RF_OVERSAMPLE: " << filtrangeinput.oversample;
      INFO.print();
      INFO << "RF_WEIGHTCORR: " << filtrangeinput.doweightcorrel;
      INFO.print();
      INFO << "RF_OVERLAP:    " << filtrangeinput.overlap;
      INFO.print();
      if (filtrangeinput.nlmean > 51)
        WARNING.print("RF_NLMEAN:     mean over more than 51 lines (?)");
      if (filtrangeinput.SNRthreshold<1.99)
        WARNING.print("RF_THRESHOLD:  < 2");
      if (filtrangeinput.SNRthreshold>10.01)
        WARNING.print("RF_THRESHOLD:  > 10 ?");
      if (filtrangeinput.oversample<=1)
        WARNING.print("RF_OVERSAMPLE: no oversampling.");
      if (filtrangeinput.oversample>8)
        WARNING.print("RF_OVERSAMPLE: >8 ?");
      if (!ispower2(filtrangeinput.oversample))
        WARNING.print("RF_OVERSAMPLE: not power of two.");
      if (filtrangeinput.doweightcorrel==true)
        WARNING.print("RF_WEIGHTCORR: weighting, not sure it has effect.");
      if (filtrangeinput.fftlength > 1024)
        WARNING.print("RF_FFTLENGTH:  adaptive filterlength > 1024 ?");
      if (filtrangeinput.SNRthreshold<0.)
        {
        PRINT_ERROR(  "RF_THRESHOLD:  < 0.");
        throw(keyword_error);
        }
      break;
    case rf_porbits:
      INFO.print("RF_METHOD:        PORBITS  \t(based on orbits.)");
      INFO << "RF_SLOPE:      " << rad2deg(filtrangeinput.terrainslope)
           << "\t[deg] terrainslope.";
      INFO.print();
      if (filtrangeinput.fftlength < 256)
        WARNING.print("RF_FFTLENGTH:  porbits filterlength < 256 (?)");
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  // ______ Both methods cards ______
  INFO << "RF_FFTLENGTH:  " << filtrangeinput.fftlength;
  INFO.print();
  INFO << "RF_HAMMING:    " << filtrangeinput.hammingalpha;
  INFO.print();
  INFO << "RF_OUT_MASTER: " << filtrangeinput.fomaster;
  INFO.print();
  INFO << "RF_OUT_SLAVE:  " << filtrangeinput.foslave;
  INFO.print();
  INFO << "RF_OUT_FORMAT: output format: ";
  switch (filtrangeinput.oformatflag)
    {
    case FORMATCR4:
      INFO << "complex_real4.";
      break;
    case FORMATCI2:
      INFO << "complex_short.";
      break;
    default:
      PRINT_ERROR("totally impossible, checked input.");
      throw(keyword_error);
    }
  INFO.print();

  // ______ Check input ______
  if (filtrangeinput.hammingalpha>0.999)
    WARNING.print("RF_HAMMING: no hamming filtering.");
  if (existed(filtrangeinput.fomaster))
    WARNING.print("RF_OUT_MASTER: file exists.");
  if (existed(filtrangeinput.foslave))
    WARNING.print("RF_OUT_SLAVE: file exists.");
  if (!ispower2(filtrangeinput.fftlength))
    {
    PRINT_ERROR(  "RF_FFTLENGTH: not power of 2.");
    throw(keyword_error);
    }
  if (filtrangeinput.overlap >= 0.5*filtrangeinput.fftlength)
    {
    PRINT_ERROR(  "RF_OVERLAP >= 0.5*RF_FFTLENGTH");
    throw(keyword_error);
    }
  if (filtrangeinput.hammingalpha>1. || filtrangeinput.hammingalpha<0.)
    {
    PRINT_ERROR(  "RF_HAMMING: not e[0,1].");
    throw(keyword_error);
    }
  } // END checkfiltrange



/****************************************************************
 *    checkdinsar                                               *
 *                                                              *
 * Checks cards for step dinsar.                                *
 *                                                              *
 #%// BK 25-Sep-2000
 ****************************************************************/
void checkdinsar(
        const input_dinsar &dinsarinput)
  {
  TRACE_FUNCTION("checkdinsar (BK 25-Sep-2000)")
  INFO.print("\n\t*** Input for step DINSAR ***");

  if (!specified(dinsarinput.topomasterresfile))
    {
    INFO.print("Using 3 pass differential (for 4 pass, see DI_IN_TOPOMASTER card).");
    }
  else
    {
    INFO << "DI_IN_TOPOMASTER: \t" << dinsarinput.topomasterresfile
         << " (4 pass)";
    INFO.print();
    }
  INFO << "DI_IN_TOPOSLAVE: \t" << dinsarinput.toposlaveresfile;
  INFO.print();
  INFO << "DI_IN_TOPOINT:   \t" << dinsarinput.topointresfile;
  INFO.print();
  INFO << "DI_OUT_FILE:     \t" << dinsarinput.fodinsar;
  INFO.print();
  if (!specified(dinsarinput.foscaleduint))
    INFO << "DI_OUT_SCALED: \tNo (debug) output requested scaled topography interf.";
  else
    INFO << "DI_OUT_SCALED: \t" << dinsarinput.foscaleduint 
         << "; debug output requested.";
  INFO.print();
  if (!specified(dinsarinput.toposlaveresfile))
    {
    PRINT_ERROR("DI_IN_TOPOSLAVE: result file topo slave not specified.");
    throw(keyword_error);
    }
  if (!specified(dinsarinput.topointresfile))
    {
    PRINT_ERROR("DI_IN_TOPOINT: result file topo interferogram not specified.");
    throw(keyword_error);
    }
  if (!strcmp(dinsarinput.toposlaveresfile,dinsarinput.topointresfile))
    {
    PRINT_ERROR("IN_TOPOSLAVE, IN_TOPOINT: Same input file name.");
    throw(keyword_error);
    }
  } // END checkdinsar



/****************************************************************
 *    checkfiltphase                                            *
 *                                                              *
 * Checks cards for step filtphase.                             *
 *                                                              *
 #%// BK 25-Sep-2000
 ****************************************************************/
void checkfiltphase(
        const input_filtphase &filtphaseinput)
  {
  TRACE_FUNCTION("checkfiltphase (BK 25-Sep-2000)")
  INFO.print("\n\t*** Input for step FILTPHASE ***");
  if (specified(filtphaseinput.fifiltphase))
    {
    INFO << "PF_IN_FILE: \t" << filtphaseinput.fifiltphase
         << " " <<  filtphaseinput.finumlines
         << " (this cr4 file will be filtered)";
    INFO.print();
    if (!existed(filtphaseinput.fifiltphase))
      WARNING.print("Impossible? PF input file does not exist?");
    }
  INFO << "PF_OUT_FILE: \t" << filtphaseinput.fofiltphase 
       << " (output filename).";
  INFO.print();

  // ______ Method goldstein filter ______
  if (filtphaseinput.method==fp_goldstein)
    {
    INFO.print("FILTPHASE: Method goldstein.");
    INFO << "PF_ALPHA: \t" << filtphaseinput.alpha
         << " (weigthing parameter for spectrum).";
    INFO.print();
    INFO << "PF_BLOCKSIZE: " << filtphaseinput.blocksize
         << " (size of block to perform filtering on).";
    INFO.print();
    INFO << "PF_OVERLAP: \t" << filtphaseinput.overlap
         << " (half overlap between consequetive blocks).";
    INFO.print();

    // ______ Use 1d kernel to smooth amplitude, e.g. 12321 ______
    // ______ Which is normalized by me ______
    INFO << "PF_KERNEL: \t";
    for (int32 ii=0; ii<filtphaseinput.kernel.pixels(); ++ii)
      INFO << " " << filtphaseinput.kernel(0,ii);
    INFO << " (smooth |spectrum| with this).";
    INFO.print();
    if (filtphaseinput.kernel.pixels()==1)
      INFO.print("No smoothing of amplitude spectrum!");

    // ______ Check errors _____
    if (filtphaseinput.alpha<0. || filtphaseinput.alpha>1.)
      WARNING.print("PF_ALPHA not 0<a<1");
    if (filtphaseinput.blocksize>64 || filtphaseinput.blocksize<16)
      WARNING.print("PF_BLOCKSIZE very small or large?");
    if (filtphaseinput.kernel.pixels()>11)
      WARNING.print("smoothing kernel > 11:  very large?");
    if (filtphaseinput.overlap<0)
      {
      PRINT_ERROR("PF_OVERLAP < 0");
      throw(keyword_error);
      }
    if (2*filtphaseinput.overlap>filtphaseinput.blocksize)
      {
      PRINT_ERROR("2*PF_OVERLAP > PF_BLOCKSIZE");
      throw(keyword_error);
      }
    if (filtphaseinput.kernel.pixels()>filtphaseinput.blocksize)
      {
      PRINT_ERROR("smoothing kernel > PF_BLOCKSIZE");
      throw(keyword_error);
      }
    if (!ispower2(filtphaseinput.blocksize))
      {
      PRINT_ERROR("PF_BLOCKSIZE not a power of 2");
      throw(keyword_error);
      }
    }

  // ______ Method modified goldstein filter ______
  else if (filtphaseinput.method==fp_modgoldstein)
    {
    INFO.print("FILTPHASE: Method modgoldstein.");
    INFO << "PF_ALPHA: \t" << -1
         << " (Calculated based on coherence).";
    INFO.print();
    INFO << "PF_BLOCKSIZE: " << filtphaseinput.blocksize
         << " (size of block to perform filtering on).";
    INFO.print();
    INFO << "PF_OVERLAP: \t" << filtphaseinput.overlap
         << " (half overlap between consequetive blocks).";
    INFO.print();

    // ______ Use 1d kernel to smooth amplitude, e.g. 12321 ______
    // ______ Which is normalized by me ______
    INFO << "PF_KERNEL: \t";
    for (int32 ii=0; ii<filtphaseinput.kernel.pixels(); ++ii)
      INFO << " " << filtphaseinput.kernel(0,ii);
    INFO << " (smooth |spectrum| with this).";
    INFO.print();
    if (filtphaseinput.kernel.pixels()==1)
      INFO.print("No smoothing of amplitude spectrum!");

    // ______ Check errors _____
    if (filtphaseinput.alpha<0. || filtphaseinput.alpha>1.)
      WARNING.print("PF_ALPHA not 0<a<1");
    if (filtphaseinput.blocksize>64 || filtphaseinput.blocksize<16)
      WARNING.print("PF_BLOCKSIZE very small or large?");
    if (filtphaseinput.kernel.pixels()>11)
      WARNING.print("smoothing kernel > 11:  very large?");
    if (filtphaseinput.overlap<0)
      {
      PRINT_ERROR("PF_OVERLAP < 0");
      throw(keyword_error);
      }
    if (2*filtphaseinput.overlap>filtphaseinput.blocksize)
      {
      PRINT_ERROR("2*PF_OVERLAP > PF_BLOCKSIZE");
      throw(keyword_error);
      }
    if (filtphaseinput.kernel.pixels()>filtphaseinput.blocksize)
      {
      PRINT_ERROR("smoothing kernel > PF_BLOCKSIZE");
      throw(keyword_error);
      }
    if (!ispower2(filtphaseinput.blocksize))
      {
      PRINT_ERROR("PF_BLOCKSIZE not a power of 2");
      throw(keyword_error);
      }
    }

  // ______ Method spatial convolution ______
  else if (filtphaseinput.method==fp_spatialconv)
    {
    INFO.print("FILTPHASE: Method spatial convolution.");
    if (!specified(filtphaseinput.fikernel2d))
      {
      INFO.print("Using 1d kernel for spatial convolution (no PF_IN_KERNEL2D).");
      INFO << "PF_KERNEL: used: \t";
      for (int32 ii=0; ii<filtphaseinput.kernel.pixels(); ++ii)
        INFO << " " << filtphaseinput.kernel(0,ii);
      INFO.print();
      }
    else
      {
      INFO.print("Using 2d kernel for spatial convolution.");
      INFO << "PF_IN_KERNEL2D: \t" << filtphaseinput.fikernel2d
           << " (ascii input file with 2d kernel).";
      INFO.print();
      INFO.print("PF_IN_KERNEL2D: \t(input file has 1 line header: numrows numcols scale");
      if (filtphaseinput.kernel.size()!=0)
        WARNING.print("PF_KERNEL card ignored due to card PF_IN_KERNEL2D.");
      if (!existed(filtphaseinput.fikernel2d))
        WARNING.print("PF_IN_KERNEL2D infile cannot be found.");
      }
    }
  else if (filtphaseinput.method==fp_spectral)
    {
    INFO.print("FILTPHASE: Method spectral filter with 2D kernel.");
    INFO << "PF_BLOCKSIZE: " << filtphaseinput.blocksize
         << " (size of block to perform filtering on).";
    INFO.print();
    INFO << "PF_OVERLAP: \t" << filtphaseinput.overlap
         << " (half overlap between consequetive blocks).";
    INFO.print();
    if (filtphaseinput.kernel.size()!=0)
      WARNING.print("PF_KERNEL card ignored for method spectral.");
    if (filtphaseinput.overlap<0)
      {
      PRINT_ERROR("PF_OVERLAP < 0");
      throw(keyword_error);
      }
    if (2*filtphaseinput.overlap>filtphaseinput.blocksize)
      {
      PRINT_ERROR("2*PF_OVERLAP > PF_BLOCKSIZE");
      throw(keyword_error);
      }
    if (!ispower2(filtphaseinput.blocksize))
      {
      PRINT_ERROR("PF_BLOCKSIZE not a power of 2");
      throw(keyword_error);
      }
    if (!specified(filtphaseinput.fikernel2d))
      {
      PRINT_ERROR("method spectral needs card PF_IN_KERNEL2D");
      throw(keyword_error);
      }
    }
  else
    {
    PRINT_ERROR("Method phasefiltering != {goldstein,modgolstein,spatialconv,spectral}.");
    throw(keyword_error);
    }
  } // END checkfiltphase


/****************************************************************
 *    checkfiltazi                                              *
 * Checks cards for step filtazi.                               *
 *    Bert Kampes, 02-Nov-2000                                  *
 ****************************************************************/
void checkfiltazi(
        const input_filtazi &filtaziinput,
        const int16 id)         // either master, slave, or m+s
  {
  TRACE_FUNCTION("checkfiltazi (BK 02-Nov-2000)")
  INFO.print("\n\t*** Input for step FILTAZI ***");
  INFO << "AF_BLOCKSIZE:   \t" << filtaziinput.fftlength;
  INFO.print();
  INFO << "AF_OVERLAP:     \t" << filtaziinput.overlap;
  INFO.print();
  INFO << "AF_HAMMING:     \t" << filtaziinput.hammingalpha;
  INFO.print();
  if (filtaziinput.oformatflag==FORMATCR4)
    INFO.print("AF_OUT_FORMAT:  \tcomplex_real4");
  else if (filtaziinput.oformatflag==FORMATCI2)
    INFO.print("AF_OUT_FORMAT:  \tcomplex_short");
  else
    {
    PRINT_ERROR("formatflag not ok for output.");
    throw(keyword_error);
    }

  if (id!=SLAVEID)
    {
    INFO << "AF_OUT_MASTER:  \t" << filtaziinput.fomaster;
    INFO.print();
    }
  if (id!=MASTERID)
    {
    INFO << "AF_OUT_SLAVE:   \t" << filtaziinput.foslave;
    INFO.print();
    }

  if (filtaziinput.fftlength<256)
    WARNING.print("AF_BLOCKSIZE < 256 (too small?)");
  if (2*filtaziinput.overlap>0.5*filtaziinput.fftlength)
    WARNING.print("2*AF_OVERLAP > .5*AF_BLOCKSIZE (very large?)");
  if (filtaziinput.hammingalpha<0 || filtaziinput.hammingalpha>1)
    {
    PRINT_ERROR("AF_HAMMING not e[0,1]");
    throw(keyword_error);
    }
  if (filtaziinput.overlap<0)
    {
    PRINT_ERROR("AF_BLOCKSIZE < 0");
    throw(keyword_error);
    }
  if (2*filtaziinput.overlap>filtaziinput.fftlength)
    {
    PRINT_ERROR("AF_BLOCKSIZE < 2*AF_BLOCKSIZE");
    throw(keyword_error);
    }
  if (!ispower2(filtaziinput.fftlength))
    {
    PRINT_ERROR("AF_BLOCKSIZE must be power of 2.");
    throw(keyword_error);
    }
  } // END checkfiltazi


