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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/unwrap.cc,v $     *
 * $Revision: 3.18 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of unwrapping routines                        *
 * -unwrapramon (stanford)                                      *
 * -snaphu_unwrap (Curtis Cheng, SU)                            *
 ****************************************************************/

#include "matrixbk.hh"
#include "constants.hh"                 // global constants
#include "productinfo.hh"               // my 'products' class
#include "exceptions.hh"                 // my exceptions class

#include "unwrap.hh"                    // header file
#include "ioroutines.hh"                // error messages
#include "coregistration.hh"            // distributepoints
#include "readinput.hh"                 // "specified" function
#include "orbitbk.hh"                   // getxyz etc.
#include "utilities.hh"                 // getBhBv etc.

#include <cctype>                       // isspace
#include <cstdio>                       // some compilers, remove function
// #include <cmath>                     // rint function



/****************************************************************
 *    unwraptreeframon                                          *
 *                                                              *
 * unwrap with code howard from ramon:                          *
 * this will not work if executable treef_ramon is not in       *
 * path.                                                        *
 * Place seeds over IFG, and start unwrapping there.  Compose   *
 * a map of totally unwrapped regions.  With matlab program     *
 * we can check these regions on consistency and adapt them.    *
 * Input:                                                       *
 *  -                                                           *
 * Output:                                                      *
 *  -                                                           *
 *                                                              *
 *    Bert Kampes, 07-Jun-1999                                  *
 ****************************************************************/
void unwraptreeframon(
        const input_gen     &generalinput,
        const input_unwrap  &unwrapinput,
        const productinfo   &interferogram)
  {
  TRACE_FUNCTION("unwraptreeframon (BK 07-Jun-1999)")

//  char prog[ONE27]="/home/hanssen/projects/cerroprieto/matlab/manuw/treef_ramon ";
  char prog[ONE27]              = "treef_ramon ";
  char fileouttreeframon[ONE27] = "ramon.uw";
  char dummyline[4*ONE27];


  const int32 dL = unwrapinput.deltaLseed;
  const int32 dP = unwrapinput.deltaPseed;
//  const int32 dL = 100;
//  const int32 dP = dL;
//  char unwrapinput.foregions[ONE27]="REGIONS";
//  char fileoutuint[ONE27]="UWPHASE";


// ______ Number of lines/pixels on file (complex interferogram) ______
  const int32 fileliness = (interferogram.win.linehi - 
                           interferogram.win.linelo + 1) / 
                           interferogram.multilookL;
  const int32 filepixels = (interferogram.win.pixhi  - 
                           interferogram.win.pixlo  + 1) /
                           interferogram.multilookP;

// ______ Allocate matrices ______
// BAND INTERLEAVED !! output of treef_ramon.f (hgt format...)
//  matrix<complr4>     AMP_PHA(fileliness,filepixels);         // read from file ramon.uw
  matrix<real4>         AMP(fileliness,filepixels);             // read from file ramon.uw
  matrix<real4>         PHA(fileliness,filepixels);             // read from file ramon.uw
  matrix<int16>         REGIONS(fileliness,filepixels);         // initialize to 0

DEBUG.print("only for now like this, maak matrix allocation with setdata??");
//  matrix<real4>       PHASE(fileliness,filepixels,NaN);       // initialize to NaN
  matrix<real4>         PHASE(fileliness,filepixels);
  PHASE.setdata(NaN);


// ______ Make commandstring ______
  char basecmdstring[3*ONE27];
  ostrstream basecmdstr(basecmdstring,3*ONE27);                 // to convert int to char
  basecmdstr.seekp(0);
  basecmdstr << prog << interferogram.file << " " 
//           << filepixels << " " << fileliness << " 1 0  ";    // FIRSTLINE = 0 ! (see source)
             << filepixels << " " << fileliness << " 0 0  ";    // no ends!!
  const int32 posbasecmdstr = basecmdstr.tellp();


// ______ Windows for band interleaved file loading, cmp readhgt ______
  const uint dummy = 999999;                    // large to force error if not ok
  window winamp(1, fileliness, 1, filepixels);
  window winpha(1, fileliness, filepixels+1, 2*filepixels);
  window offsetbuffer(1,dummy,1,dummy); // dummy not used in readfile, no offset

  register int32 regionnumber = 0;
  register int32 line, pixel;
  register int32 npoints = 0;
  uint seedL,seedP;
  matrix<uint> SEEDS;
  if (isspace(unwrapinput.seedfile[0]))         // no filename specified
    {
// ______ First determine size of matrix the easy way (just counting) ______
    for (seedL=4; seedL<=fileliness-4; seedL=seedL+dL)
      {
      for (seedP=4; seedP<=filepixels-4; seedP=seedP+dP)
        {
        npoints++;
        }
      }

// ______ Then fill the matrix ______
    SEEDS.resize(npoints,2);
    register int32 cnt = 0;
    for (seedL=4; seedL<=fileliness-4; seedL=seedL+dL)
      {
      for (seedP=4; seedP<=filepixels-4; seedP=seedP+dP)
        {
        SEEDS(cnt,0) = seedL;
        SEEDS(cnt,1) = seedP;
        cnt++;
        }
      }
    }

  else                                  // read in file + check
    {
    INFO.print("Using seedfile for unwrapping seeds.");
    if (!existed(unwrapinput.seedfile))
      {
      PRINT_ERROR("UW_SEED seedfile: does not exist.")
      throw(input_error);
      }
    npoints = filelines(unwrapinput.seedfile);          // number of enters
    INFO << "Number of points in seedfile (should be an enter after last one): "
         << npoints;
    INFO.print();
    ifstream ifile(unwrapinput.seedfile, ios::in);      // not binary
    SEEDS.resize(npoints,2);
#ifdef __DEBUG
    bk_assert(ifile,unwrapinput.seedfile,__FILE__,__LINE__);
#endif

// ______ Read in data ______
    for (line=0; line<npoints; line++)          // read in data
      {
      ifile >> seedL >> seedP;
      SEEDS(line,0) = seedL;
      SEEDS(line,1) = seedP;
      ifile.getline(dummyline,4*ONE27,'\n');
      }

// ______ Check input (wrapped interferogram starts at (1,1)) ______
    if (max(SEEDS.getcolumn(0)) > fileliness)
      {
      PRINT_ERROR("max line in seedfile exceeds number of lines on file (complex interferogram).")
      throw(input_error);
      }
    if (max(SEEDS.getcolumn(1)) > filepixels)
      {
      PRINT_ERROR("max pixel in seedfile exceeds number of pixels on file (complex interferogram).")
      throw(input_error);
      }
    } // seedfile specified


  for (npoints=0; npoints<SEEDS.lines(); npoints++)
    {
    seedL = SEEDS(npoints,0);
    seedP = SEEDS(npoints,1);

// ______ Check if region already unwrapped ______
    INFO << "Seed (l,p): " << seedL << ", " << seedP;
    INFO.print();
    if (REGIONS(seedL,seedP) == 0)            // not yet unwrapped
      {
      regionnumber++;
      basecmdstr.seekp(posbasecmdstr);
      basecmdstr << seedP << " " << seedL << ends;
      INFO << "basecmdstring: " << basecmdstring;
      INFO.print();

      remove(fileouttreeframon);                        // try to remove always
      system(basecmdstring);

// ______ Read in output of treef, band interleaved ______
// ______ Use amp for check if unwrapped ok.
      readfile (AMP, fileouttreeframon, AMP.lines(), winamp, 
      offsetbuffer);
      readfile (PHA, fileouttreeframon, PHA.lines(), winpha,
      offsetbuffer);

// ______ Loop over all points, fill REGIONS and PHASE matrix ______
      for (line=0; line<AMP.lines(); line++)
        {
        for (pixel=0; pixel<AMP.pixels(); pixel++)
          {
          if (AMP(line,pixel) != 0)             // unwrapping ok ?
            {
            REGIONS(line,pixel) = regionnumber;
            PHASE(line,pixel)   = PHA(line,pixel);
            }
          }
        }
      }
    } // all seeds

  if (remove(fileouttreeframon))
    WARNING.print("code 101: could not remove ramon.uw file.");


// ______ Remove regions that are not ok unwrapped (too small) ______
  const int32 SEVENTY = 70;                             // criterium
  const int32 maxregions = regionnumber;
  int32 numberofregions = regionnumber;
  register int32 numberofpixels  = 0;

  for (regionnumber=1; regionnumber<=maxregions; regionnumber++)
    {
  
// ______ Count number of pixels ______
    for (line=0; line<REGIONS.lines(); line++)
      {
      if (numberofpixels >= SEVENTY) 
        break;                                          // break inner loop
      for (pixel=0; pixel<REGIONS.pixels(); pixel++)
        {
        if (REGIONS(line,pixel) == regionnumber)
          {
          numberofpixels++;
          }
        }
      }
  
// ______ If not ok, delete info ______
    if (numberofpixels == 0)
      {
      WARNING.print("this is impossible..., but still i got this message ones???");
      }
  
    else if (numberofpixels < SEVENTY)                  // too small
      {
      numberofregions--;
      for (line=0; line<REGIONS.lines(); line++)
        {
        for (pixel=0; pixel<REGIONS.pixels(); pixel++)
          {
          if (REGIONS(line,pixel) == regionnumber)
            {
            PHASE(line,pixel) = NaN;                    // delete unwrapped phase
            REGIONS(line,pixel) = 0;                    // delete regionnumber
            }
          }
        }
      }
    }   // loop all regions



// ====== Write result to output files ======
  ofstream ofileregions;
  openfstream(ofileregions,unwrapinput.foregions,generalinput.overwrit);
  bk_assert(ofileregions,unwrapinput.foregions,__FILE__,__LINE__);
  ofileregions << REGIONS;
  ofileregions.close();

  ofstream ofileuint;
  openfstream(ofileuint,unwrapinput.fouint,generalinput.overwrit);
  bk_assert(ofileuint,unwrapinput.fouint,__FILE__,__LINE__);
  ofileuint << PHASE;
  ofileuint.close();


// ====== Write info to resultfiles ======
  ofstream scratchlogfile("scratchlogunwrap", ios::out | ios::trunc);
  scratchlogfile << "\n\n*******************************************************************"
                 << "\n**** UNWRAP                                                   *****"
                 << "\n*******************************************************************"
                 << "\nData_output_file: \t\t"
                 <<  unwrapinput.fouint
                 << "\nData_output_format: \t\t"
                 << "real4"
                 << "\nData_output_file_regions: \t"
                 << unwrapinput.foregions 
                 << "\nData_output_format: \t\t"
                 << "short int (2B)"
                 << "\nProgram for unwrappng: \t\t"
                 <<  prog
                 << "\nOutput program for unwrapping: \t"
                 <<  fileouttreeframon;
  if (isspace(unwrapinput.seedfile[0]))         // no seedfile specified
    {
    scratchlogfile << "\nDelta lines for seed: \t\t"
                   <<  dL
                   << "\nDelta pixels for seed: \t\t"
                   <<  dP;
    }
  else
    { 
    scratchlogfile << "\nSeeds where read from file: \t"
                   <<  unwrapinput.seedfile;
    }
  scratchlogfile << "\nNumber of patches used: \t"
                 <<  numberofregions
                 << "\n* End_unwrap:"
                 << "\n*******************************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchresunwrap", ios::out | ios::trunc);
  bk_assert(scratchresfile,"unwrap: scratchresunwrap",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 //<< "\n*_Start_unwrap"
                 << "\n*_Start_" << processcontrol[pr_i_unwrap]
                 << "\n*******************************************************************"
                 << "\nData_output_file:                     \t"
                 <<  unwrapinput.fouint
                 << "\nData_output_format:                   \t"
                 << "real4"
                 << "\nData_output_file_regions:             \t"
                 << unwrapinput.foregions 
                 << "\nData_output_format:                   \t"
                 << "short int (2B)"
                 << "\nFirst_line (w.r.t. original_master):  \t"
                 <<  interferogram.win.linelo
                 << "\nLast_line (w.r.t. original_master):   \t"
                 <<  interferogram.win.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  interferogram.win.pixlo
                 << "\nLast_pixel (w.r.t. original_master):  \t"
                 <<  interferogram.win.pixhi
                 << "\nMultilookfactor_azimuth_direction:    \t"
                 <<  interferogram.multilookL
                 << "\nMultilookfactor_range_direction:      \t"
                 <<  interferogram.multilookP

                 << "\nProgram for unwrappng:                \t"
                 <<  prog
                 << "\nOutput program for unwrapping:        \t"
                 <<  fileouttreeframon;
  if (isspace(unwrapinput.seedfile[0]))         // no seedfile specified
    {
    scratchresfile
                 << "\nDelta lines for seed:                 \t"
                 <<  dL
                 << "\nDelta pixels for seed:                \t"
                 <<  dP;
    }
  else
    { 
    scratchresfile << "\nSeeds where read from file:      \t"
                   <<  unwrapinput.seedfile;
    }
  scratchresfile << "\nNumber of patches used:               \t"
                 <<  numberofregions
                 << "\n*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_unwrap] << "_NORMAL"
                 << "\n*******************************************************************\n";
  scratchresfile.close();

  // ______ Tidy up ______
  PROGRESS.print("finished unwraptreeframon.");
  } // END unwraptreeframon



/****************************************************************
 *    snaphu_unwrap                                             *
 *                                                              *
 * unwrap with snaphu public domain code, system call           *
 * this will not work if executable snaphu is not in path       *
 *                                                              *
 * Input:                                                       *
 *  wrapped ifg                                                 *
 * Output:                                                      *
 *  unwrapped ifg                                               *
 *                                                              *
 #%// BK 01-Nov-2002
#%// PM 26-Feb-2007: multi.CPU support
 ****************************************************************/
void snaphu_unwrap(
        const input_gen     &generalinput,
        const input_unwrap  &unwrapinput,
        const productinfo   &interferogram,
        const slcimage      &master,
        const slcimage      &slave,
              orbit         &masterorbit,
              orbit         &slaveorbit,
        const input_ell     &ellips)
  {
  TRACE_FUNCTION("snaphu_unwrap (Bert Kampes 01-Nov-2002)")
  char prog[ONE27]       = "snaphu ";// run this executable
  char configfile[ONE27] = "snaphu.conf";// create this file
  INFO << "Configuration file for snaphu: " << configfile;
  INFO.print();


  // ______ Number of lines/pixels on file (complex interferogram) ______
  const int32 filelines  = (interferogram.win.linehi - 
                           interferogram.win.linelo + 1) / 
                           interferogram.multilookL;
  const int32 filepixels = (interferogram.win.pixhi  - 
                           interferogram.win.pixlo  + 1) /
                           interferogram.multilookP;

  // ______ Make commandstring ______
  char basecmdstring[3*ONE27];
  ostrstream basecmdstr(basecmdstring,3*ONE27);// to convert int to char
  basecmdstr.seekp(0);
  basecmdstr << prog << " -f " << configfile << " "
             << interferogram.file << " " 
             << filepixels << ends;


  // _____ Compute some parameters for snaphu here ______
  real8 NEARRANGE = master.pix2range(interferogram.win.pixlo);
    //pix2range(interferogram.win.pixlo,master.t_range1,master.rsr2x);
  real8 DR        = real8(interferogram.multilookP) *
    (master.pix2range(interferogram.win.pixlo+1.0) - NEARRANGE);
    //(pix2range(interferogram.win.pixlo+1,master.t_range1,master.rsr2x) - NEARRANGE);
  real8 RANGERES  = ((master.rsr2x/2.)/master.rbw) * 
                    (DR/real8(interferogram.multilookP));

  cn P;   // point at ellips for mid image, returned by lp2xyz
  real8 line  = 0.5*real8(interferogram.win.linelo+interferogram.win.linehi);
  real8 pixel = 0.5*real8(interferogram.win.pixlo+interferogram.win.pixhi);
  lp2xyz(line,pixel,ellips,master,masterorbit,P);// returns P
  real8 EARTHRADIUS = P.norm();
  // ______ Compute xyz for master satellite ______
  real8 time_azi;                               // returned
  real8 time_range;                             // returned
  xyz2t(time_azi,time_range,master, masterorbit, P);
  cn M = masterorbit.getxyz(time_azi);          // master satellite position
  xyz2t(time_azi,time_range,slave, slaveorbit, P);
  cn S = slaveorbit.getxyz(time_azi);           // slave satellite position
  real8 ORBITRADIUS = M.norm();
  // ______ Baseline parametrizations returned ______
  real8 Bh, Bv, Bpar, Bperp, theta;
  real8 BASELINE;
  real8 BASELINEANGLE_RAD;
  BalphaBhBvBparBperpTheta(
    BASELINE, BASELINEANGLE_RAD,
    Bh, Bv, Bpar, Bperp, theta, M, P, S);
  real8 BASELINEANGLE_DEG = rad2deg(BASELINEANGLE_RAD);
  ;// dA via orbit...
  lp2xyz(line+1.,pixel,ellips,master,masterorbit,P);// returns P
  xyz2t(time_azi,time_range,master, masterorbit, P);// returns time
  cn M1 = masterorbit.getxyz(time_azi);
  real8 DA        = abs(real8(interferogram.multilookL)*M.dist(M1));
  real8 AZRES     = ((master.prf)/master.abw) * 
                    (DA/real8(interferogram.multilookL));

  // ______ Create config file ______
  ofstream snaphuconfigfile(configfile, ios::out | ios::trunc);
  snaphuconfigfile
    << "# snaphu configuration file\n"
    << "#\n"
    << "# Lines with fewer than two fields and lines whose first non-whitespace\n"
    << "# characters are not alphnumeric are ignored.  For the remaining lines,\n"
    << "# anything after the first two fields (delimited by whitespace) is\n"
    << "# also ignored.  Inputs are converted in the order they appear in the file;\n"
    << "# if multiple assignments are made to the same parameter, the last one\n"
    << "# given is the one used.  Parameters in this file will be superseded by\n"
    << "# parameters given on the command line after the -f flag specifying this\n"
    << "# file.  Multiple configuration files may be given on the command line.\n"
    << "\n"
    << "# CONFIG FOR SNAPHU\n"
    << "# -----------------------------------------------------------\n"
    << "# Created by Doris software#\n"
    << "# Has been run from within Doris with command:\n\n"
    << "#    " << basecmdstring << "\n"
    << "# -----------------------------------------------------------\n"
    << "\n"
    << "# Statistical-cost mode (TOPO, DEFO, SMOOTH, or NOSTATCOSTS)\n"
    << "STATCOSTMODE    " << unwrapinput.snaphu_mode << "\n"
    << "\n"
    << "# Output file name\n"
    << "OUTFILE         " << unwrapinput.fouint << "\n"
    << "\n"
    << "# Correlation file name\n";
  if (specified(unwrapinput.snaphu_coh))
    {
    DEBUG.print("Using coherence of Doris for snaphu:");
    DEBUG.print("Please make sure window sizes used in config file are correct");
    DEBUG.print("Maybe you will have to re-run snaphu by hand? (see notes in snaphu.conf)");
    snaphuconfigfile
      << "CORRFILE      " << unwrapinput.snaphu_coh << "\n"
      << "\n"
      << "# Correlation file format\n"
      << "CORRFILEFORMAT        FLOAT_DATA\n";
    }
  else
    { 
    snaphuconfigfile
      << "# CORRFILE      " << "doris.coh" << "\n"
      << "\n"
      << "# Correlation file format\n"
      << "# CORRFILEFORMAT        FLOAT_DATA\n";
    }
  snaphuconfigfile
    << "# Text file to which runtime parameters will be logged.  The format of\n"
    << "# that file will be suitable so that it can also be used as a\n"
    << "# configuration file.\n";
  if (specified(unwrapinput.snaphu_log))
    {
    snaphuconfigfile
      << "LOGFILE           " << unwrapinput.snaphu_log
      << "\n";
    }
  else
    {
    snaphuconfigfile
      << "# LOGFILE           " << "snaphu.log"
      << "\n";
    }
  snaphuconfigfile
    << "# Algorithm used for initialization of wrapped phase values.  Possible\n"
    << "# values are MST and MCF.\n"
    << "INITMETHOD      " << unwrapinput.snaphu_init << "\n"
    << "\n"
    << "# Verbose-output mode (TRUE or FALSE)\n"
    << "VERBOSE         " << unwrapinput.snaphu_verbose << "\n"
    << "\n"
    << "\n"
    << "################\n"
    << "# File formats #\n"
    << "################\n"
    << "\n"
    << "# Valid data formats:\n"
    << "#\n"
    << "# COMPLEX_DATA:      complex values: real, imag, real, imag\n"
    << "# ALT_LINE_DATA:     real values from different arrays, alternating by line\n"
    << "# ALT_SAMPLE_DATA:   real values from different arrays, alternating by sample\n"
    << "# FLOAT_DATA:        single array of floating-point data\n"
    << "#\n"
    << "\n"
    << "# Input file format\n"
    << "INFILEFORMAT          COMPLEX_DATA\n"
    << "\n"
    << "# Output file format (this is almost hgt, BK)\n";
  if (unwrapinput.oformatflag == FORMATHGT)
    {
    snaphuconfigfile
      << "OUTFILEFORMAT         ALT_LINE_DATA\n"
      << "# OUTFILEFORMAT         FLOAT_DATA\n";
    }
  else if (unwrapinput.oformatflag == FORMATR4)
    {
    snaphuconfigfile
      << "# OUTFILEFORMAT         ALT_LINE_DATA\n"
      << "OUTFILEFORMAT         FLOAT_DATA\n";
    }
  else
    {
    PRINT_ERROR("unknown format specified for output with SNAPHU")
    throw(unhandled_case_error);
    }
  snaphuconfigfile
    << "\n"
    << "\n"
    << "###############################\n"
    << "# SAR and geometry parameters #\n"
    << "###############################\n"
    << "\n"
    << "# Orbital radius (double, meters) or altitude (double, meters).  The\n"
    << "# radius should be the local radius if the orbit is not circular.  The\n"
    << "# altitude is just defined as the orbit radius minus the earth radius.\n"
    << "# Only one of these two parameters should be given.\n"
    << "#ORBITRADIUS             7153000.0 (example)\n"
    << "ORBITRADIUS             "
    << setiosflags(ios::fixed) << setw(10) << setprecision(2)
    << ORBITRADIUS << "\n"
    << "#ALTITUDE               775000.0\n"
    << "\n"
    << "# Local earth radius (double, meters).  A spherical-earth model is\n"
    << "# used.\n"
    << "#EARTHRADIUS             6378000.0 (example)\n"
    << "EARTHRADIUS             " 
    << setiosflags(ios::fixed) << setw(10) << setprecision(2)
    << EARTHRADIUS << "\n"
    << "\n"
    << "# The baseline parameters are not used in deformation mode, but they\n"
    << "# are very important in topography mode.  The parameter BASELINE\n"
    << "# (double, meters) is the physical distance (always positive) between\n"
    << "# the antenna phase centers.  The along-track componenet of the\n"
    << "# baseline is assumed to be zero.  The parameter BASELINEANGLE_DEG\n"
    << "# (double, degrees) is the angle between the antenna phase centers\n"
    << "# with respect to the local horizontal.  Suppose the interferogram is\n"
    << "# s1*conj(s2).  The baseline angle is defined as the angle of antenna2\n"
    << "# above the horizontal line extending from antenna1 towards the side\n"
    << "# of the SAR look direction.  Thus, if the baseline angle minus the\n"
    << "# look angle is less than -pi/2 or greater than pi/2, the topographic\n"
    << "# height increases with increasing elevation.  The units of\n"
    << "# BASELINEANGLE_RAD are radians.\n"
    << "#BASELINE                150.0 (example)\n"
    << "#BASELINEANGLE_DEG       225.0 (example)\n"
    << "#BASELINEANGLE_RAD      3.92699 (example)\n"
    << "BASELINE                " << BASELINE << "\n"
    << "BASELINEANGLE_DEG       " << BASELINEANGLE_DEG << "\n"
    << "#BASELINEANGLE_RAD      " << BASELINEANGLE_RAD << "\n"
    << "\n"
    << "# If the BPERP parameter is given, the baseline angle is taken to be\n"
    << "# equal to the look angle (mod pi) at midswath, and the length of the\n"
    << "# baseline is set accordingly.  Particular attention must be paid to\n"
    << "# the sign of this parameter--it should be negative if increasing\n"
    << "# phase implies increasing topographic height.\n"
    << "#BPERP          -150.0 (example)\n"
    << "#BPERP          " << Bperp << "\n"
    << "\n"
    << "# The transmit mode should be either REPEATPASS or PINGPONG if both\n"
    << "# antennas transmitted and both received (REPEATPASS and PINGPONG have\n"
    << "# the same effect); the transmit mode should be SINGLEANTENNATRANSMIT\n"
    << "# if only one antenna was used to transmit while both antennas\n"
    << "# received.  In single-antenna-transmit mode, the baseline is\n"
    << "# effectively halved.  This parameter is ignored for cost modes other\n"
    << "# than topography.\n"
    << "TRANSMITMODE    REPEATPASS\n"
    << "\n"
    << "# Slant range from platform to first range bin in input data file\n"
    << "# (double, meters).  Be sure to modify this parameter if the input\n"
    << "# file is extracted from a larger scene.  The parameter does not need\n"
    << "# to be modified is snaphu is unwrapping only a subset of the input file.\n"
    << "#NEARRANGE       831000.0 (example)\n"
    << "NEARRANGE       " 
    << setiosflags(ios::fixed) << setw(20) << setprecision(2)
    << NEARRANGE << "\n"
    << setw(20) << setprecision(10)
    << "\n"
    << "# Slant range and azimuth pixel spacings of input interferogram after\n"
    << "# any multilook averaging.  This is not the same as the resolution.\n"
    << "# (double, meters).\n"
    << "#DR              8.0 (example)\n"
    << "#DA              20.0 (example)\n"
    << "DR              " << DR << "\n"
    << "DA              " << DA << "\n"
    << "\n"
    << "# Single-look slant range and azimuth resolutions.  This is not the\n"
    << "# same as the pixel spacing.  (double, meters).\n"
    << "#RANGERES        10.0 (example)\n"
    << "#AZRES           6.0 (example)\n"
    << "RANGERES        " << RANGERES << "\n"
    << "AZRES           " << AZRES    << "\n"
    << "\n"
    << "# Wavelength (double, meters).\n"
    << "#LAMBDA          0.0565647 (example)\n"
    << "LAMBDA          " << master.wavelength << "\n"
    << "\n"
    << "# Number of real (not necessarily independent) looks taken in range and\n"
    << "# azimuth to form the input interferogram (long).\n"
    << "NLOOKSRANGE     " << interferogram.multilookP << "\n"
    << "NLOOKSAZ        " << interferogram.multilookL << "\n"
    << "\n"
    << "\n"
    << "# Equivalent number of independent looks (double, dimensionless) that were\n"
    << "# used to generate correlation file if one is specified.  This parameter\n"
    << "# is ignored if the correlation data are generated by the interferogram\n"
    << "# and amplitude data.\n"
    << "#"
    << "# The equivalent number of independent looks is approximately equal to the\n"
    << "# real number of looks divided by the product of range and azimuth\n"
    << "# resolutions, and multiplied by the product of the single-look range and\n"
    << "# azimuth pixel spacings.  It is about 0.53 times the number of real looks\n"
    << "# for ERS data processed without windowing.\n"
    << "# (BK: taken from Curtis example config file.)\n"
    << "NCORRLOOKS      23.8\n"
    << "\n"
    << "# Number of looks that should be taken in range and azimuth for estimating\n"
    << "# the correlation coefficient from the interferogram and the amplitude\n"
    << "# data.  These numbers must be larger than NLOOKSRANGE and NLOOKSAZ.\n"
    << "# The actual numbers used may be different since we prefer odd integer\n"
    << "# multiples of NLOOKSRANGE and NLOOKSAZ (long).  These numbers are ignored\n"
    << "# if a separate correlation file is given as input.\n"
    << "# (BK: taken from Curtis example config file.)\n"
    << "NCORRLOOKSRANGE 3\n"
    << "NCORRLOOKSAZ    15\n"
        << "\n"
        << "################\n"
        << "# Tile control #\n"
        << "################\n"
        << "\n"
        << "# Parameters in this section describe how the input files will be \n"
        << "# tiled.  This is mainly used for tiling, in which different \n"
        << "# patches of the interferogram are unwrapped separately.\n"
        << "\n"
        << "# Number of rows and columns of tiles into which the data files are\n"
        << "# to be broken up.\n"
        << "NTILEROW           " << unwrapinput.ntilerow << "\n"
        << "NTILECOL           " << unwrapinput.ntilecol << "\n"
        << "\n"
        << "# Overlap, in pixels, between neighboring tiles.\n"
        << "ROWOVRLP           " << unwrapinput.rowovrlp << "\n"
        << "COLOVRLP           " << unwrapinput.colovrlp << "\n"
        << "\n"
        << "# Maximum number of child processes to start for parallel tile\n"
        << "# unwrapping.\n"
        << "NPROC              " << unwrapinput.nproc << "\n"
        << "\n"
        << "# Cost threshold to use for determining boundaries of reliable regions\n"
        << "# (long, dimensionless; scaled according to other cost constants).\n"
        << "# Larger cost threshold implies smaller regions---safer, but\n"
        << "# more expensive computationally.\n"
        << "TILECOSTTHRESH   "  << unwrapinput.tilecostthresh << "\n"
        << "\n"
        << "# Minimum size (long, pixels) of a reliable region in tile mode.  \n"
        << "# MINREGIONSIZE             100\n"
        << "\n"
        << "# Extra weight applied to secondary arcs on tile edges.\n"
        << "# TILEEDGEWEIGHT    2.5\n"
        << "\n"
        << "# Maximum flow magnitude (long) whose cost will be stored in the secondary \n"
        << "# cost lookup table.  Secondary costs larger than this will be approximated\n"
        << "# by a quadratic function.\n"
        << "# SCNDRYARCFLOWMAX  8\n"
        << "\n"
        << "# The program will remove temporary tile files if this is set.\n"
        << "# RMTMPTILE                 FALSE\n"
        << "\n"
        << "# If this is set to anything besides FALSE, the program will skip\n"
        << "# the unwrapping step and only assemble temporary tile files from a previous \n"
        << "# invocation saved in the directory whose name is given here.  The tile size \n"
        << "# parameters and file names must be the same.\n"
        << "# ASSEMBLEONLY              tiledir\n"
    << "\n"
    << "# End of snaphu configuration file\n"
    << "\n";
  snaphuconfigfile.close();
  PROGRESS << "Configuration file for snaphu: " << configfile << " created." << ends;
  PROGRESS.print();

  if (unwrapinput.snaphu_dumponlyconf == false)
    {
      // ______ Call snaphu! ______
      INFO.print("System call for running snaphu follows:");
      INFO.print(basecmdstring);
      system(basecmdstring);
    }
  else
    {
      INFO.print("System call for running snaphu follows:");
      INFO.print(basecmdstring);
      INFO.print("UW_SNAPHU_DUMPONLYCONF was ON. Please run SNAPHU manually.");
    }

  // ====== Write info to resultfiles ======
  ofstream scratchlogfile("scratchlogunwrap", ios::out | ios::trunc);
  scratchlogfile
    << "\n\n*********************************************************"
    << "\n**** UNWRAP                                         *****"
    << "\n*********************************************************"
    << "\nData_output_file: \t\t"
    <<  unwrapinput.fouint
    << "\nData_output_format: \t\t"
    << "hgt or real4"
    << "\nProgram for unwrappng: \t\t"
    <<  prog
    << "\n* End_unwrap:"
    << "\n*********************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchresunwrap", ios::out | ios::trunc);
  bk_assert(scratchresfile,"unwrap: scratchresunwrap",__FILE__,__LINE__);
  scratchresfile 
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_unwrap]
    << "\n*******************************************************************"
    << "\nData_output_file:                     \t"
    <<  unwrapinput.fouint
    << "\nData_output_format:                   \t";
    // bugfix #%// Bert Kampes, 25-May-2005: format 
    if (unwrapinput.oformatflag == FORMATHGT)
      scratchresfile << "hgt";
    else if (unwrapinput.oformatflag == FORMATR4)
      scratchresfile << "real4";
    else
      throw(unhandled_case_error);
  scratchresfile 
    << "\nFirst_line (w.r.t. original_master):  \t"
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  interferogram.win.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  interferogram.win.pixhi
    << "\nMultilookfactor_azimuth_direction:    \t"
    <<  interferogram.multilookL
    << "\nMultilookfactor_range_direction:      \t"
    <<  interferogram.multilookP
    << "\nNumber of lines (multilooked):        \t"
    <<  filelines 
    << "\nNumber of pixels (multilooked):       \t"
    <<  filepixels
    << "\nProgram used for unwrapping:          \t"
    <<  prog
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_unwrap] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();

  // ______ Tidy up ______
  PROGRESS.print("Finished snaphu_unwrap.");
  } // END snaphu_unwrap

