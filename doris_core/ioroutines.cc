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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/ioroutines.cc,v $
 * $Revision: 3.25 $
 * $Date: 2005/10/06 11:09:20 $
 * $Author: kampes $
 *
 * implementation of io routines.
 ****************************************************************/


#include "constants.hh"         // global constants
#include "ioroutines.hh"        // header of this file
#include "utilities.hh"         // ispower2
#include "slcimage.hh"          // my slc image class
#include "productinfo.hh"       // my 'products' class
#include "exceptions.hh"                 // my exceptions class



#include <fstream>              // for file streams
#include <strstream>            // for memory stream
#include <cstring>              // for strcmp etc.
#include <iomanip>              // for setprecision etc.
#include <cstdlib>              // exit
#include <cctype>               // isspace()
#include <cstdarg>              // for ellipses in functions
#include <algorithm>            // generic max
#include <cstdio>               // some compilers, remove function
#include <ctime>                // some compilers, clock functions
#include <sstream>              // int2str [MA]
char *strptime(const char *s, const char  *format,  struct tm *tm);




/****************************************************************
 *    printcpu                                                  *
 *                                                              *
 * Prints cpu and wallclock to screen as DEBUG                  *
 * Printcpu(true) should be called at start main.               *
 *                                                              *
 * input:                                                       *
 *  - bool init: first time true, default false                 *
 * output:                                                      *
 *  - screen: cpu and wallclock time                            *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void printcpu(
        bool init)
  {
  TRACE_FUNCTION("printcpu (BK 11-Dec-1998)")
  time_t nseconds = time(NULL);
  char *tijd = ctime(&nseconds);                      // includes newline
  INFO << "Current time: " << tijd;
  INFO.print();

  //static int32 cput_last   = clock();
  //static real4 cput_total  = 0.;
  //int32        cput_process = clock() - cput_last;
  //static int32 wct_total   = 0;               // wallclock, init at first call
  //static int32 wct_last    = time(NULL);      // wallclock, init at first call
  //int32        wct_process = time(NULL) - wct_last;

  static clock_t cput_last    = clock();// init only first call
  clock_t        cput_process = clock() - cput_last;
  static real4   cput_total   = 0.0;// init only first call

  static time_t wct_last    = time(NULL);       // wallclock, init at first call
  time_t        wct_process = time(NULL) - wct_last;
  static int32  wct_total   = 0;                // wallclock, init at first call

  if (init)                                     // return if only initialization
    return;

  cput_total += (real4(cput_process))/real4(CLOCKS_PER_SEC);
  cput_last   = clock();
  //wct_total  += wct_process;
  wct_total  += int32(wct_process);
  wct_last    = time(NULL);

  const int32 cputmin = int32(floor(cput_total/60));
  const real4 cputsec = cput_total - cputmin*60;
  DEBUG << " cputime used for process: \t"
       <<  setw(6) << real4(cput_process)/real4(CLOCKS_PER_SEC)
       << " sec (total: "
       << setw(4) << cputmin << " min "
       << setw(3) << cputsec << " sec)\n"
       << "\t   wallclock: \t\t\t"
       << setw(6) << wct_process
       << " sec (total: "
       << setw(4) << wct_total/60 << " min "
       << setw(3) << wct_total%60 << " sec)";
  DEBUG.print();
  cerr << "total cpu: \t"
       << setw(4) << cputmin << " min "
       << setw(3) << cputsec << " sec\n";
  } // END printcpu



/****************************************************************
 *    inittest                                                  *
 *                                                              *
 * Some initial tests and info for portability.                 *
 *                                                              *
 * input:                                                       *
 *  - void                                                      *
 * output:                                                      *
 *  - void, exits if trouble                                    *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void inittest()
  {
  TRACE_FUNCTION("inittest (BK 11-Dec-1998)")
  // ______ Some info ______
  #ifdef __DEBUG
    DEBUG.print("\"__DEBUG\"                 defined (extra verbose)");
  #else
    DEBUG.print("\"__DEBUG\"                 not defined (extra verbose)");
  #endif
  #ifdef __DEBUGMAT1
    DEBUG.print("\"__DEBUGMAT1\"             defined (range checking matrix class)");
  #else
    DEBUG.print("\"__DEBUGMAT1\"             not defined (range checking matrix class)");
  #endif
  #ifdef __DEBUGMAT2
    DEBUG.print("\"__DEBUGMAT2\"             defined (verbose matrix class)");
  #else
    DEBUG.print("\"__DEBUGMAT2\"             not defined (verbose matrix class)");
  #endif
  #ifdef __USE_FFTW_LIBRARY__
    DEBUG.print("\"__USE_FFTW_LIBRARY__\"    defined (using fftw routines)");
  #else
    DEBUG.print("\"__USE_FFTW_LIBRARY__\"    not defined (not using fftw routines)");
  #endif
  #ifdef __USE_VECLIB_LIBRARY__
    DEBUG.print("\"__USE_VECLIB_LIBRARY__\"  defined (using library routines)");
  #else
    DEBUG.print("\"__USE_VECLIB_LIBRARY__\"  not defined (using internal routines)");
  #endif
  #ifdef __USE_LAPACK_LIBRARY__
    DEBUG.print("\"__USE_LAPACK_LIBRARY__\"  defined (using library routines)");
  #else
    DEBUG.print("\"__USE_LAPACK_LIBRARY__\"  not defined (using internal routines)");
  #endif
  #ifdef __GplusplusCOMPILER__
    DEBUG.print("\"__GplusplusCOMPILER__\"   defined (GNU c++ compiler)");
    WARNING.print("NOT USED ANYMORE. since v3.8");
  #else
    DEBUG.print("\"__GplusplusCOMPILER__\"   not defined (no GNU c++ compiler)");
  #endif
  #ifdef __GNUC__
    DEBUG << "\"__GNUC__\":               " << __GNUC__;
    DEBUG.print();
  #endif
  #ifdef __GNUC_MINOR__
    DEBUG << "\"__GNUC_MINOR__\":         " << __GNUC_MINOR__;
    DEBUG.print();
  #endif
  // ___ Testing endianness ___
  DEBUG.print("Testing endianness.");
  int i;
  int littleendian=0;// thus big assumed
  char test[64]={0};
  for (i=0;i<32;i++) test[i]=0;
  test[0]=0x01;
  if (0x0001==*(int *)(&test[0])) littleendian=1;
  #ifdef __X86PROCESSOR__
    DEBUG.print("\"__X86PROCESSOR__\"        defined (little Endian machine)");
    if (littleendian == 1)
      INFO.print("Little Endian machine defined and this is correct.");
    else
      WARNING.print("Little Endian machine defined and this is NOT correct.");
  #else
    DEBUG.print("\"__X86PROCESSOR__\"        not defined (big Endian machine)");
    if (littleendian == 0)
      INFO.print("Big Endian machine defined and this is correct.");
    else
      WARNING.print("Big Endian machine defined and this is NOT correct.");
  #endif

  #ifdef __NO_STRPTIME
    DEBUG.print("\"__NO_STRPTIME\"           defined (using internal routine)");
  #else
    DEBUG.print("\"__NO_STRPTIME\"           not defined (using library routine)");
  #endif
  // ___ Testing strptime function ___
  DEBUG.print("Testing strptime function.");
  struct tm tm_ref;
  char utc_ref[100] = "05-JAN-1985 01:02:03.000";
  strptime(utc_ref,"%d-%b-%Y %T",&tm_ref);
  int status = 0;
  if (tm_ref.tm_mday != 5)  status=1;
  if (tm_ref.tm_mon  != 0)  status=1;
  if (tm_ref.tm_year != 85) status=1;
  if (tm_ref.tm_hour != 1)  status=1;
  if (tm_ref.tm_min  != 2)  status=1;
  if (tm_ref.tm_sec  != 3)  status=1;
  if (status == 0)
    INFO.print("strptime function works fine.");
  else
    WARNING.print("strptime function seems NOT TO BE OK.");


  // ______ Some info ______
  if (sizeof(int16) != 2)
    {
    PRINT_ERROR("code: 900: sizeof int16(short) != 2: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(int32) != 4)
    {
    PRINT_ERROR("code: 900: sizeof int32(int) != 4: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(uint) != 4)
    {
    PRINT_ERROR("code: 900: sizeof uint(unsigned int) != 4: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(real8) != 8)
    {
    PRINT_ERROR("code: 900: sizeof real8(double) != 8: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(real4) != 4)
    {
    PRINT_ERROR("code: 900: sizeof real4(float) != 4: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(compli16) != 4)
    {
    PRINT_ERROR("code: 900: sizeof compli16(complex short) != 4: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(compli32) != 8)
    {
    PRINT_ERROR("code: 900: sizeof compli32(complex int) != 8: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(complr4) != 8)
    {
    PRINT_ERROR("code: 900: sizeof complr4(complex float) != 8: see typedefs in constants.h")
    throw(some_error);
    }
  if (sizeof(complr8) != 16)
    {
    PRINT_ERROR("code: 900: sizeof complr16(complex double) != 16: see typedefs in constants.h")
    throw(some_error);
    }
  if (int32(7.5) != 7)
    {
    PRINT_ERROR("code: 900: it is assumed that int(7.5)==7")
    throw(some_error);
    }
  if (complr4(1.1) != complr4(1.1,0))
    {
    PRINT_ERROR("code: ???: it is assumed that complr4(7.5) == 7.5 + 0.0i")
    throw(some_error);
    }
  DEBUG.print("test somewhere also memory,diskspace.");

#ifdef __DEBUG
  // ___ Testing cn class ___
  DEBUG.print("Testing coordinate class:");
  cn P; P.test();
  // ___ Testing window class ___
  DEBUG.print("Testing window class:");
  window W; W.test();
  // ___ Testing ellips class ___
  DEBUG.print("Testing ellips class:");
  input_ell E; E.test();
#endif
  } // END inittest


/****************************************************************
 *    doinitwrite                                               *
 *                                                              *
 * Decides to call initwrite or not                             *
 *  necessary because no call for image (slave or master).      *
 * if something will be processed then true                     *
 * input:                                                       *
 *  - input processes                                           *
 *  - imageid: =2 master, =3 slave                              *
 * output:                                                      *
 *  - bool true or false                                        *
 * future: THIS ROUTINE SHOULD BE UPDATED AS WE GO ALONG        *
 *                                                              *
 *    Bert Kampes, 15-Dec-1998                                  *
 *                                                              *
 * Mahmut: make slc result files equivalent                     *
 ****************************************************************/
bool doinitwrite(
        input_gen &generalinput,
        int16 imageid)
  {
  TRACE_FUNCTION("doinitwrite (BK 15-Dec-1998)")
  // ______Check input and decide (normally true)______
  if (imageid == MASTERID)
    {
    if (!(generalinput.process[pr_m_readfiles]  ||
          generalinput.process[pr_m_crop]       ||
          generalinput.process[pr_m_oversample] ||
          generalinput.process[pr_m_porbits]    ||
	  generalinput.process[pr_m_morbits]    || //[HB]
          generalinput.process[pr_m_simamp]     || //[MA] 2008
          generalinput.process[pr_m_mtiming]    || //[MA] 2008
          generalinput.process[pr_m_resample]   || //[MA] 200903, fake entry to make slc.res files equivalent. Note: no dependency defined in check routine.
          generalinput.process[pr_m_filtazi]    ||
          generalinput.process[pr_m_filtrange]  ||
          generalinput.process[pr_m_EXTRA]        ))
      return false;
    }
  else if (imageid == SLAVEID)
    {
    if (!(generalinput.process[pr_s_readfiles]  ||
          generalinput.process[pr_s_crop]       ||
          generalinput.process[pr_s_oversample] ||
          generalinput.process[pr_s_porbits]    ||
	  generalinput.process[pr_s_morbits]    || //[HB]
          generalinput.process[pr_s_simamp]     || //[MA] 200903. fake entry
          generalinput.process[pr_s_mtiming]    || //[MA] 200903. fake entry
          generalinput.process[pr_s_resample]   ||
          generalinput.process[pr_s_filtazi]    ||
          generalinput.process[pr_s_filtrange]  ||
          generalinput.process[pr_s_EXTRA]        ))
      return false;
    }
  else if (imageid == INTERFID)
    {
    if (!(generalinput.process[pr_i_coarse]      ||
          generalinput.process[pr_i_coarse2]     ||
          generalinput.process[pr_i_fine]        ||
          generalinput.process[pr_i_timing]      || //[FvL]
          generalinput.process[pr_i_demassist]   || //[FvL]
          generalinput.process[pr_i_coregpm]     ||
          generalinput.process[pr_i_comprefpha]  ||
          generalinput.process[pr_i_subtrrefpha] ||
          generalinput.process[pr_i_comprefdem]  ||
          generalinput.process[pr_i_subtrrefdem] ||
          generalinput.process[pr_i_interfero]   ||
          generalinput.process[pr_i_coherence]   ||
          generalinput.process[pr_i_filtphase]   ||
          generalinput.process[pr_i_unwrap]      ||
	  generalinput.process[pr_i_estorbits]   || //[HB]
          generalinput.process[pr_i_slant2h]     ||
          generalinput.process[pr_i_geocoding]   ||
          generalinput.process[pr_i_dinsar]      ||
          generalinput.process[pr_i_EXTRA2]        ))
      return false;
    }
  else
    {
    PRINT_ERROR("code 901: wrong input.")
    throw(some_error);
    }
  return true;
  }  // END doinitwrite


/****************************************************************
 *    initwrite                                                 *
 *                                                              *
 * Writes some general info to logfile and resultfile           *
 *                                                              *
 * input:                                                       *
 *  - file: name of file                                        *
 *  - fileid:                                                   *
 *      1:logfile, 2:masterres, 3:slaveres, 4:interferogram     *
 * output:                                                      *
 *  - void, file is updated                                     *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void initwrite(
        const char* file,
        int16 fileid)
  {
  TRACE_FUNCTION("initwrite (BK 11-Dec-1998)")
  #ifndef SWNAME
  DEBUG.print("SWNAME NOT DEFINED.");
  #endif
  #ifndef SWVERSION
  DEBUG.print("SWVERSION NOT DEFINED.");
  #endif


// ______Initialize time______
// #ifdef __GplusplusCOMPILER__
// g++ 2.7.2.2 wants int, though long is exactly same
// use __GNUC__ = 2, __GNUC_MINOR__ = 7, __GNUG__ = 2
// can be seen with g++ -v dummy.c
//#if __GNUC_MINOR__ == 7
//  int32 nseconds = time(NULL);
//#else // at least version 2.95.2 and HP aCC
//  long nseconds = time(NULL);
//#endif
  time_t nseconds = time(NULL);
  char *tijd = ctime(&nseconds);                      // includes newline

// ______Write logfile if not exists______
  bool ofilenew = true;
  if (existed(file)) ofilenew = false;                // [MA] TODO [hi]: fix a potential bug when ifg.res file is missing, it is regenerated partially.

  ofstream ofile;
  ofile.open(file, ios::out | ios::app);
  bk_assert(ofile,file,__FILE__,__LINE__);


  DEBUG.print("Writing (updating) of file (initwrite).");
  if (ofilenew)
    {
    ofile
          << "=====================================================";
    switch (fileid)
      {
      case LOGID:
        ofile << "\n LOGFILE:                  " << file;
        break;
      case MASTERID:
        ofile << "\n MASTER RESULTFILE:        " << file;
        break;
      case SLAVEID:
        ofile << "\n SLAVE RESULTFILE:         " << file;
        break;
      case INTERFID:
        ofile << "\n INTERFEROGRAM RESULTFILE: " << file;
        break;
      default:
        PRINT_ERROR("panic: code 901:initwrite, internal fileid not ok!")
        throw(unhandled_case_error);
      }
    ofile << "\n\nCreated by:             " // << getlogin() // screen caused trouble with this // [MA] try getUserName() on Windows
          << "\nInSAR Processor:        " << SWNAME
          << "\nVersion:                " << SWVERSION
    #if defined (__DEBUG) || defined (__DEBUGMAT1) || defined (__DEBUGMAT2)
          << " (debug)"
    #else
          << " (optimal)"
    #endif
    #ifdef __USE_FFTW_LIBRARY__
          << "\nFFTW library:           " << "used"
    #else
          << "\nFFTW library:           " << "not used"
    #endif
    #ifdef __USE_VECLIB_LIBRARY__
          << "\nVECLIB library:         " << "used"
    #else
          << "\nVECLIB library:         " << "not used"
    #endif
    #ifdef __USE_LAPACK_LIBRARY__
          << "\nLAPACK library:         " << "used"
    #else
          << "\nLAPACK library:         " << "not used"
    #endif
          << "\nCompiled at:            " << __DATE__ << " " << __TIME__;
    #ifdef __GNUC__
      ofile << "\nBy GNU gcc:             " << __GNUC__;
      #ifdef __GNUC_MINOR__
        ofile << "." << __GNUC_MINOR__;
        #ifdef __GNUG__
          ofile << "." << __GNUG__;
        #endif
      #endif
    #endif
    ofile << "\nFile creation at:       " << tijd
          << "\n -------------------------------------------------------"
          << "\n| Delft Institute of Earth Observation & Space Systems  |"
          << "\n|          Delft University of Technology               |"
          << "\n|              http://doris.tudelft.nl                  |"
          << "\n|                                                       |"
          << "\n| Author: (c) TUDelft - DEOS Radar Group                |"
          << "\n -------------------------------------------------------\n\n";

// ======Write process control in result file======
    switch (fileid)
      {

// ______These may not be changed without changing check process control______
      case MASTERID:
        ofile << "\nStart_process_control\n"
              << processcontrol[pr_m_readfiles]   <<   " \t\t0\n"
              << processcontrol[pr_m_porbits]     <<   " \t0\n"
	      << processcontrol[pr_m_morbits]     <<   " \t\t0\n"    //[HB]
              << processcontrol[pr_m_crop]        <<   " \t\t\t0\n"
              << processcontrol[pr_m_simamp]      <<   " \t\t0\n"    //[MA] 2008
              << processcontrol[pr_m_mtiming]     <<   " \t\t0\n"    //[MA] 2008
              << processcontrol[pr_m_oversample]  <<   " \t\t0\n"
              << processcontrol[pr_m_resample]    <<   " \t\t0\n"    //[MA] 2009. fake entry to make slc.res files equivalent
              << processcontrol[pr_m_filtazi]     <<   " \t\t0\n"
              << processcontrol[pr_m_filtrange]   <<   " \t\t0\n"
              << processcontrol[pr_m_EXTRA]       <<   " \t\t0\n"
              << "End_process_control\n";
        break;

// ______These may not be changed without changing check process control______
      case SLAVEID:
        ofile << "\nStart_process_control\n"
              << processcontrol[pr_s_readfiles]   <<   " \t\t0\n"
              << processcontrol[pr_s_porbits]     <<   " \t0\n"
	      << processcontrol[pr_s_morbits]     <<   " \t\t0\n"    //[HB]
              << processcontrol[pr_s_crop]        <<   " \t\t\t0\n"
              << processcontrol[pr_s_simamp]      <<   " \t\t0\n"    //[MA] 2009. fake entry
              << processcontrol[pr_s_mtiming]     <<   " \t\t0\n"    //[MA] 2009. fake entry no processing is defined
              << processcontrol[pr_s_oversample]  <<   " \t\t0\n"
              << processcontrol[pr_s_resample]    <<   " \t\t0\n"
              << processcontrol[pr_s_filtazi]     <<   " \t\t0\n"
              << processcontrol[pr_s_filtrange]   <<   " \t\t0\n"
              << processcontrol[pr_s_EXTRA]       <<   " \t\t0\n"
              << "End_process_control\n";
        break;

// ______These may not be changed without changing check process control______
      case INTERFID:
        ofile << "\nStart_process_control\n"
              << processcontrol[pr_i_coarse]      <<   " \t\t0\n"
              << processcontrol[pr_i_coarse2]     <<   " \t\t0\n"
              << processcontrol[pr_i_fine]        <<   " \t\t0\n"
              << processcontrol[pr_i_timing]      <<   " \t\t0\n" //[FvL]
              << processcontrol[pr_i_demassist]   <<   " \t\t0\n" //[FvL]
              << processcontrol[pr_i_coregpm]     <<   " \t\t0\n"
              << processcontrol[pr_i_interfero]   <<   " \t\t0\n"
              << processcontrol[pr_i_coherence]   <<   " \t0\n"
              << processcontrol[pr_i_comprefpha]  <<   " \t\t0\n"
              << processcontrol[pr_i_subtrrefpha] <<   " \t\t0\n"
              << processcontrol[pr_i_comprefdem]  <<   " \t\t0\n"
              << processcontrol[pr_i_subtrrefdem] <<   " \t\t0\n"
              << processcontrol[pr_i_filtphase]   <<   " \t\t0\n"
              << processcontrol[pr_i_unwrap]      <<   " \t\t0\n"
	      << processcontrol[pr_i_estorbits]   <<   " \t\t0\n" //[HB]
              << processcontrol[pr_i_slant2h]     <<   " \t\t0\n"
              << processcontrol[pr_i_geocoding]   <<   " \t\t0\n"
              << processcontrol[pr_i_dinsar]      <<   " \t\t\t0\n"
              << processcontrol[pr_i_EXTRA2]      <<   " \t\t0\n"
              << "End_process_control\n";
        break;

      default:
        ; //do nothing
      } // switch fileid
    }

  else                                                          // file already existed
    {
    ofile << endl << endl
          << "   *====================================================================*\n"
          << "   |                                                                    |\n"
          << "       Following part is appended at: " << tijd // \n included in tijd...
          << "                 Using Doris "          << SWVERSION <<       "          \n"
          << "   |                                                                    |\n"
          << "   *--------------------------------------------------------------------*\n\n";
    } //if else file existence

// ______Tidy up______
  ofile.close();
  } // END initwrite


/****************************************************************
 *    updatefile                                                *
 *                                                              *
 * A scratchfile is appended to a logfile.                      *
 * The scratchfile is deleted after copying.                    *
 *                                                              *
 * input:                                                       *
 *  - scratchfile                                               *
 *  - logfile                                                   *
 * output:                                                      *
 *  - void, updated logfile, scratchfile is deleted             *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 * changed to 4*127 after problem.
 #%// BK 07-Mar-2001
 ****************************************************************/
void updatefile(
        const char* sfile,
        const char* lfile)
  {
  TRACE_FUNCTION("updatefile (BK 11-Dec-1998)")

  char                  dummyline[4*ONE27];
  //ifstream scratchfile(sfile, ios::in | ios::nocreate);
  ifstream scratchfile(sfile, ios::in);
  bk_assert(scratchfile,sfile,__FILE__,__LINE__);

  ofstream logfile(lfile, ios::out | ios::app);
  bk_assert(logfile,lfile,__FILE__,__LINE__);

// ______Start writing______
  scratchfile.getline(dummyline,4*ONE27,'\n');            // to prevent twice last line
  while (scratchfile)
    {
    logfile << dummyline << endl;
    scratchfile.getline(dummyline,4*ONE27,'\n');
    }

  // ______ Add current time ______
  time_t nseconds = time(NULL);
  char *tijd = ctime(&nseconds);                      // includes newline
  logfile << "\n   Current time: " << tijd; // \n included in tijd...

// ______Tidy up______
  logfile.close();
  scratchfile.close();
  if (remove(sfile))                                    // remove file
    WARNING.print("Could not remove file.");

  DEBUG << "File: " << lfile
       << " has been updated by file: " << sfile << " (removed.)";
  DEBUG.print();
  } // END updatefile


/****************************************************************
 *    getanswer                                                 *
 *                                                              *
 * Waits for user to press a key.                               *
 * (Mainly) used in interactive mode.                           *
 * should not wait for enter but getch() doesnt work !!         *
 * (dos only?)
 *                                                              *
 * future: no enter required                                    *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void getanswer(
        )
  {
  TRACE_FUNCTION("getanswer (BK 11-Dec-1998)")
  char dummychar;
  cerr << "\n Press <ENTER> to continue.";
  cin.unsetf(ios::skipws);      // ignore ws (just enter)
  cin >> dummychar;
  cerr << " continuing...\n";
  } // END getanswer


/****************************************************************
 *    readres                                                   *
 *                                                              *
 * Pattern is searched in file (1st word),                      *
 * next (skipwords) words are ignored (default=0)               *
 * the following word is returned as returnword,                *
 *  sizeofrw includes '\0'                                      *
 *                                                              *
 * eg: in main c8[9] as readres(c8,sizeof(c8),file,pat)         *
 *                                                              *
 * input:                                                       *
 *  - pointer to returnword                                     *
 *  - size of returnword (including '\0')                       *
 *  - file name to search                                       *
 *  - pattern to search for                                     *
 *  - number of words to skip before desired returnword (0=def  *
 * output:                                                      *
 *  - true if found, returnword is filled                       *
 *                                                              *
 * future: should include identifier as well                    *
 * better?: file.readsome(pattern)                              *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 *    Mahmut Arikan, 17-Apr-2009 skip matching lines            *
 ****************************************************************/
bool readres(
        char* returnword,
        const int16 sizeofrw,
        const char* file,
        const char* pattern,
        const int16 skipwords,      //=0 default
        int16 skiplines)            //=0 default
  {
  TRACE_FUNCTION("readres (BK 11-Dec-1998)")
  char                  dummyline[4*ONE27];
  char                  word[4*ONE27];
  bool                  foundword = false;
  ifstream infile(file, ios::in);
  bk_assert(infile,file,__FILE__,__LINE__);

  // ======Search infile======
  returnword[0] = '\0';                        // garanty word is empty
  char ch       = '\t';
  while (infile)
    {
    infile >> word;
    if (strcmp(pattern,word))                  // no pattern match.
      {
      infile.getline(dummyline,4*ONE27,'\n');    // goto next line.
      }
    else                                       // pattern match.
      {
      if ( skiplines == 0)
        {
        for (register int32 i=0; i<skipwords; i++)
          infile >> word;                        // goto correct position: word
        int32 tin;
        while (isspace(tin=ch))                  // goto correct position: trailing blanks
          {
          infile.get(ch);
          if (ch == '\n')                        // no argument after key
            {
            WARNING.print("Found identifier, but not so desired word!, exiting routine...");
            break;
            }
          }
        }
       else
         {
          skiplines -= 1 ;  // [MA]
          continue;
         }

      // ______ Actual read word ______
      foundword = true;
      //infile.putback(ch);
      //infile.read(returnword,sizeofrw-1);    // read desired word including blanks
      //returnword[sizeofrw]='\0'; ?
      // ______ Test bert ______
      for (int ii=0; ii<sizeofrw; ++ii)
        {
        if (ch == '\n')                        // no argument after key
          {
          returnword[ii] = '\0';
          break;
          }
        else
          {
          returnword[ii] = ch;
          }
        infile.get(ch);
        }
      //returnword[sizeofrw-1]='\0';
      break;                                   // file
      }                                        // else
    }                                          // file

  // ______ Give some info ______
  if (foundword)
    {
    DEBUG << "read: \"" << returnword << "\" as word " << skipwords << " after \""
         << pattern << "\"";
    DEBUG.print();
    }
  else
    {
    ERROR << "readres: could not find pattern: \"" << pattern
         << "\" in file: " << file;
    PRINT_ERROR(ERROR.get_str());
    throw(file_error);
    }

  infile.close();
  return foundword;
  } // END readres


/****************************************************************
 *    updateprocesscontrol                                      *
 *                                                              *
 * Update the process controls flags at begin result file       *
 *                                                              *
 * input:                                                       *
 *  - file: name of file                                        *
 * output:                                                      *
 *  - void, file is updated                                     *
 *                                                              *
 *    Bert Kampes, 14-Dec-1998                                  *
 ****************************************************************/
void updateprocesscontrol(
        const char* file,
        int16 fileid)
  {
  TRACE_FUNCTION("updateprocesscontrol (BK 14-Dec-1998)")
//____RaffaeleNutricato START MODIFICATION SECTION 5
//  int16                 checkprocess[NUMPROCESSES]; // already processed?
  int16                 checkprocess[NUMPROCESSES+1]; // Raffaele Nutricato, Hopefully this fixes bug!
//____RaffaeleNutricato END MODIFICATION SECTION 5

  char                  dummyline[4*ONE27];
  //bool                processed = false;      // process control in resultfile
  //bool                processerror = false;   // to avoid error in middle of copying

  // ______Check existense, input______
  if (!existed(file))
    return;

// ======Update process controls======
  ofstream tmpfile("scratchcopy", ios::out | ios::trunc);          // temporary copy
  bk_assert(tmpfile,"updateprocesscontrol: scratchcopy",__FILE__,__LINE__);
  //ifstream resfile(file, ios::in | ios::nocreate);
  ifstream resfile(file, ios::in);
  bk_assert(resfile,file,__FILE__,__LINE__);

// ______copy resfile to tmpfile and fill checkprocess______
  for (register int32 i=0;i<NUMPROCESSES;i++)
    checkprocess[i]=0;
  resfile.getline(dummyline,4*ONE27,'\n');                // to prevent twice last line
  while(resfile)
    {
    tmpfile << dummyline << endl;                       // copy line
    fillcheckprocess(dummyline,checkprocess,fileid);      // check line to fill ...
    resfile.getline(dummyline,4*ONE27,'\n');              // get next line
    }
  resfile.close();
  tmpfile.close();


// ______Update process_controls______
  //ifstream tmpfile2("scratchcopy", ios::in | ios::nocreate);  // temporary copy
  ifstream tmpfile2("scratchcopy", ios::in);    // temporary copy
  bk_assert(tmpfile2,"updateprocesscontrols: scratchcopy",__FILE__,__LINE__);
  ofstream resfile2(file, ios::out | ios::trunc);                  // do replace !
  bk_assert(resfile2,file,__FILE__,__LINE__);


// ______Copy back, up to process control______
  while (strcmp(dummyline,"Start_process_control"))
    {
    tmpfile2.getline(dummyline,4*ONE27,'\n');
    resfile2 << dummyline << endl;
    }

// ______Update process controls______
  char word[EIGHTY]=" ";
  for (;;)                                      // forever
    {
    tmpfile2 >> word;
    resfile2 << word;
    tmpfile2.getline(dummyline,4*ONE27,'\n');     // go to next line
    if (!strcmp(word,"End_process_control"))
      {
      resfile2 << endl;
      break;
      }

    switch (fileid)
      {
    case MASTERID:
      if      (!strcmp(word,processcontrol[pr_m_readfiles]))
        (checkprocess[pr_m_readfiles]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_porbits]))
        (checkprocess[pr_m_porbits])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_morbits]))
        (checkprocess[pr_m_morbits])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n"; //[HB]
      else if (!strcmp(word,processcontrol[pr_m_crop]))
        (checkprocess[pr_m_crop])   ? resfile2 << "\t\t\t1\n" : resfile2 << "\t\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_simamp]))
        (checkprocess[pr_m_simamp])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";      //[MA] 2008
      else if (!strcmp(word,processcontrol[pr_m_mtiming]))
        (checkprocess[pr_m_mtiming])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";      //[MA] 2008
      else if (!strcmp(word,processcontrol[pr_m_oversample]))
        (checkprocess[pr_m_oversample])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_resample]))                                 //[MA] 2009, fake entry no processing is defined
        (checkprocess[pr_m_resample])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_filtazi]))
        (checkprocess[pr_m_filtazi])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_filtrange]))
        (checkprocess[pr_m_filtrange])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_m_EXTRA]))
        (checkprocess[pr_m_EXTRA])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else
        {
        ERROR << "PANIC: forgotten to update routine? " << word
             << " not recognized in master resultfile.";
        PRINT_ERROR(ERROR.get_str())
        throw(unhandled_case_error);
        }
      break; // switch

    case SLAVEID:
      if      (!strcmp(word,processcontrol[pr_s_readfiles]))
        (checkprocess[pr_s_readfiles]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_porbits]))
        (checkprocess[pr_s_porbits])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_morbits]))
        (checkprocess[pr_s_morbits])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n"; //[HB]
      else if (!strcmp(word,processcontrol[pr_s_crop]))
        (checkprocess[pr_s_crop])   ? resfile2 << "\t\t\t1\n" : resfile2 << "\t\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_simamp]))
        (checkprocess[pr_s_simamp])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";      //[MA] 2009, fake entry no processing is defined
      else if (!strcmp(word,processcontrol[pr_s_mtiming]))
        (checkprocess[pr_s_mtiming])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";      //[MA] 2009, fake entry no processing is defined
      else if (!strcmp(word,processcontrol[pr_s_oversample]))
        (checkprocess[pr_s_oversample])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_resample]))
        (checkprocess[pr_s_resample])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_filtazi]))
        (checkprocess[pr_s_filtazi])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_filtrange]))
        (checkprocess[pr_s_filtrange])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_s_EXTRA]))
        (checkprocess[pr_s_EXTRA])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else
        {
        ERROR << "PANIC: forgotten to update routine? " << word
             << " not recognized in slave resultfile.";
        PRINT_ERROR(ERROR.get_str())
        throw(unhandled_case_error);
        }
      break; // switch

    case INTERFID:
      if      (!strcmp(word,processcontrol[pr_i_coarse]))
        (checkprocess[pr_i_coarse])    ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_coarse2]))
        (checkprocess[pr_i_coarse2])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_fine]))
        (checkprocess[pr_i_fine])      ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_timing]))
        (checkprocess[pr_i_timing])      ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n"; //[FvL]
      else if (!strcmp(word,processcontrol[pr_i_demassist]))
        (checkprocess[pr_i_demassist])      ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n"; //[FvL]
      else if (!strcmp(word,processcontrol[pr_i_coregpm]))
        (checkprocess[pr_i_coregpm])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_interfero]))
        (checkprocess[pr_i_interfero]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_coherence]))
        (checkprocess[pr_i_coherence]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_comprefpha]))
        (checkprocess[pr_i_comprefpha]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_subtrrefpha]))
        (checkprocess[pr_i_subtrrefpha]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_comprefdem]))
        (checkprocess[pr_i_comprefdem]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_subtrrefdem]))
        (checkprocess[pr_i_subtrrefdem]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_filtphase]))
        (checkprocess[pr_i_filtphase]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_unwrap]))
        (checkprocess[pr_i_unwrap])    ? resfile2 << "\t\t\t1\n" : resfile2 << "\t\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_estorbits]))
        (checkprocess[pr_i_estorbits])  ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n"; //[HB]
      else if (!strcmp(word,processcontrol[pr_i_slant2h]))
        (checkprocess[pr_i_slant2h])   ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_geocoding]))
        (checkprocess[pr_i_geocoding]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_dinsar]))
        (checkprocess[pr_i_dinsar]) ? resfile2 << "\t\t\t1\n" : resfile2 << "\t\t\t0\n";
      else if (!strcmp(word,processcontrol[pr_i_EXTRA2]))
        (checkprocess[pr_i_EXTRA2]) ? resfile2 << "\t\t1\n" : resfile2 << "\t\t0\n";
      else
        {
        ERROR << "PANIC: forgotten to update routine? " << word
             << " not recognized in interferogram resultfile.";
        PRINT_ERROR(ERROR.get_str())
        throw(unhandled_case_error);
        }
      break; // switch

    default:
      PRINT_ERROR("PANIC, impossible wrong.")
      throw(unhandled_case_error);
      } // switch
    } // forever


// ======Copy rest of file======
  tmpfile2.getline(dummyline,4*ONE27,'\n');       // get line
  while(tmpfile2)
    {
    resfile2 << dummyline << endl;
    tmpfile2.getline(dummyline,4*ONE27,'\n');
    }
  resfile2.close();
  tmpfile2.close();

// ______Tidy up______
//  if (processerror)
//    ERROR(__FILE__,__LINE__,"updateprocesscontrols: resfile some problem.");
  if (remove("scratchcopy"))                            // remove file
    WARNING.print("Could not remove file: scratchcopy.");
  } // END updateprocesscontrols


/****************************************************************
 *    checkprocessing                                           *
 *                                                              *
 * Checks if requested processing is possible.                  *
 * checkprocess[i] contains already processed                    *
 *  (master index is mis-used to indicate)                      *
 * inputgeneral.process[i] contains requests                     *
 * requirements for processing can be found in software folder  *
 * general if step is requested:                                *
 *  0) check already                                            *
 *  1) check requirements (other steps/files)                   *
 *  2) set checkprocessed true                                   *
 *                                                              *
 * input:                                                       *
 * output:                                                      *
 *  - void (no exit hopefully)                                  *
 *  - checkprocess is filled with what is already processed      *
 *                                                              *
 *BUG: filled with what will be processed as well...            *
 *BUG??: not correctly passed back ???                          *
 *                                                              *
 *    Bert Kampes, 14-Dec-1998                                  *
 ****************************************************************/
void checkprocessing(
        const input_gen &generalinput,
        int16 checkprocess[NUMPROCESSES])
  {
  TRACE_FUNCTION("checkprocessing (BK 14-Dec-1998)")
  char          dummyline[4*ONE27] = " ";
  //char          word[EIGHTY]     = " ";
  const char*   file;
  int16         checkprocess2[NUMPROCESSES+1];   // check with filecontent

  checkprocess2[NUMPROCESSES]=0;
  register int32 i;
  for (i=0;i<NUMPROCESSES;i++)
    {
    checkprocess[i]=0;
    checkprocess2[i]=0;
    }

// ______Fill checkprocess______
  for (int32 fileid=MASTERID; fileid<=INTERFID; fileid++)
    {
    if (fileid == MASTERID)
      file=generalinput.m_resfile;
    else if (fileid == SLAVEID)
      file=generalinput.s_resfile;
    else if (fileid == INTERFID)
      file=generalinput.i_resfile;
    else
      {
      PRINT_ERROR("checkprocessing, PANIC, wrong input")
      throw(input_error);
      }

    if (existed(file))
      {
      int32 breakit=0;
      while(breakit<3)                          // try 2 recover in case of error
        {
        ifstream resfile;// required to put in loop for g++v3.2// BK 08-Apr-2003
        resfile.open(file, ios::in);
        bk_assert(resfile,file,__FILE__,__LINE__);

        int32 linecnt=0;
        while (strcmp(dummyline,"Start_process_control"))
          {
          resfile.getline(dummyline,4*ONE27,'\n');
          DEBUG << "read line: " << dummyline << ends;
          DEBUG.print();
          linecnt++;
          if (linecnt==100)
            {
            WARNING << "Checked first 100 lines, did not find: \"Start_process_control\" in file: " << file ;
            WARNING.print();
            break;
            }
          }

        // ______ Read processcontrols into array checkprocess ______
        linecnt=0;
        while (strcmp(dummyline,"End_process_control"))
          {
          resfile.getline(dummyline,4*ONE27,'\n');
          DEBUG << "read line: " << dummyline;
          DEBUG.print();
          fillprocessed(dummyline,checkprocess,fileid);
          linecnt++;
          if (linecnt==100)
            {
            WARNING << "Checked first 100 lines, did not find: \"End_process_control\"   in file: " << file ;
            WARNING.print();
            break;
            }
          }

        // ______ Read resultsections in file for array checkprocess2 ______
        resfile.getline(dummyline,4*ONE27,'\n');          // read line
        while(resfile)
          {
          fillcheckprocess(dummyline,checkprocess2,fileid);       // check line
          resfile.getline(dummyline,4*ONE27,'\n');                // read line
          DEBUG << "read line: " << dummyline;
          DEBUG.print();
          }
        resfile.close();

        // ______ Check resultsections with process control flags ______
        bool dofixprocesscontrol=false;
        for (i=0;i<NUMPROCESSES;i++)
          {
          if (checkprocess[i] != checkprocess2[i])
            {
            dofixprocesscontrol=true;
            WARNING << "Step: " << i << " (" << processcontrol[i]
                 << ") ";
            if (checkprocess[i]==1)
              WARNING << "in process control flag, but result is not in \""
                   << file << "\".";
            else
              WARNING << "not in process control flag, but result is in \""
                   << file << "\".";
            WARNING.print();
            if (generalinput.interactive) getanswer();
            } // something is wrong
          } // for all steps

        // ______ Check if repairs have to be made ______
        if (dofixprocesscontrol)
          {
          if (breakit == 1)
            {
            cerr << "\nAlready tried to fix process controls. Should I try again?\n";
            getanswer();
            }
          updateprocesscontrol(file, fileid);                   // repair routine
          breakit += 1;                                         // only one try
          }
        else // nothing strange
          breakit = 10;                                         // > stop condition
      } // try to repair
    } // existed(file)
  } // for fileid


// ______ Keep checkprocess like this (only what is actually in resultfile) ______
// ______ and update a new array ______
  int16 checkprocesstmp[NUMPROCESSES];
  for (i=0; i<NUMPROCESSES; i++)
    checkprocesstmp[i]=checkprocess[i];

// ====== Check with requested processes ======
// ______ The order is important (BK 26mar2001) ______
// ______ This means that to avoid warnings, filtrange should be after coarse corr.
// ______ Master ______
  if (generalinput.process[pr_m_readfiles])                       // requested
    checkrequest(pr_m_readfiles,checkprocesstmp,0);               // no requirements
  if (generalinput.process[pr_m_crop])                            // requested
    checkrequest(pr_m_crop,checkprocesstmp,
                 1,pr_m_readfiles);                               // required (RN)?
    //checkrequest(pr_m_crop,checkprocesstmp,0);                  // no requirements
  if (generalinput.process[pr_m_oversample])                      // requested
    checkrequest(pr_m_oversample,checkprocesstmp,1,pr_m_crop);    // oversample requires a cropped image
  if (generalinput.process[pr_m_porbits])                         // requested
    checkrequest(pr_m_porbits,checkprocesstmp, 1,pr_m_readfiles); // required for time info
    if (checkprocess2[NUMPROCESSES])                              // extra check
      DEBUG.print("orbits from leader file will be deleted.");
  if (generalinput.process[pr_m_morbits])                         // requested [HB]
    checkrequest(pr_m_morbits,checkprocesstmp,                    // required
		 1,pr_m_readfiles);
  if (generalinput.process[pr_m_simamp])                          // requested [MA]
    checkrequest(pr_m_simamp,checkprocesstmp, 1,pr_m_crop);       // amplitude simulation requires crop step
  if (generalinput.process[pr_m_mtiming])                         // requested [MA]
    checkrequest(pr_m_mtiming,checkprocesstmp, 1,pr_m_simamp);       // correlation with simulated amplitude requires a simulated amplitude image
  if (generalinput.process[pr_m_filtazi])                         // requested
    checkrequest(pr_m_filtazi,checkprocesstmp,
                 1,pr_m_crop);                                    // required
  if (generalinput.process[pr_m_EXTRA])                           // requested
    checkrequest(pr_m_EXTRA,checkprocesstmp,0);                   // no requirements

// ______ Slave ______
  if (generalinput.process[pr_s_readfiles])                     // requested
    checkrequest(pr_s_readfiles,checkprocesstmp,0);             // no requirements
  if (generalinput.process[pr_s_crop])                          // requested
    checkrequest(pr_s_crop,checkprocesstmp,
                 1,pr_s_readfiles);                             // required for check
  if (generalinput.process[pr_s_oversample])                    // requested
    checkrequest(pr_s_oversample,checkprocesstmp,1,pr_s_crop);  // oversample requires a cropped image
  if (generalinput.process[pr_s_porbits])                       // requested
    checkrequest(pr_s_porbits,checkprocesstmp,
                 1,pr_s_readfiles);                             // required for time info
  if (generalinput.process[pr_s_morbits])                         // requested  [HB]
    checkrequest(pr_s_morbits,checkprocesstmp,                    // required
		 1,pr_s_readfiles);
    if (checkprocess2[NUMPROCESSES])                            // extra check
      DEBUG.print("orbits from leader file will be deleted.");
  if (generalinput.process[pr_s_filtazi])                       // requested
    checkrequest(pr_s_filtazi,checkprocesstmp,
                 1,pr_s_crop);                                  // required
  if (generalinput.process[pr_s_EXTRA])                         // requested
    checkrequest(pr_s_EXTRA,checkprocesstmp,0);                 // no requirements

// ______ Interferogram ______
  if (generalinput.process[pr_i_coarse])                        // requested coarse orbits
    checkrequest(pr_i_coarse,checkprocesstmp,
    2,pr_m_readfiles,pr_s_readfiles);                           // required [MA] crop --> readfiles
  if (generalinput.process[pr_i_coarse2])                       // requested coarse correlation
    checkrequest(pr_i_coarse2,checkprocesstmp,
    2,pr_m_crop,pr_s_crop);                                     // required
  // orbits...
  if (generalinput.process[pr_m_filtrange])                     // requested
    checkrequest(pr_m_filtrange,checkprocesstmp,
                 2,pr_m_crop,pr_i_coarse2);                  // adviced
  // adaptive...
  //if (generalinput.process[pr_m_filtrange])                     // requested
  //  checkrequest(pr_m_filtrange,checkprocesstmp,
  //               2,pr_m_crop,pr_s_resample);                    // required
  //if (generalinput.process[pr_s_filtrange])                   // requested
  //  checkrequest(pr_s_filtrange,checkprocesstmp,
  //               1,pr_s_crop);                                // required
  if (generalinput.process[pr_i_fine])                          // requested
    checkrequest(pr_i_fine,checkprocesstmp,
    2,pr_m_crop,pr_s_crop);
  if (generalinput.process[pr_i_timing])                    // requested
    checkrequest(pr_i_timing,checkprocesstmp,
    2,pr_i_coarse,pr_i_fine); //[FvL]
  if (generalinput.process[pr_i_demassist])                    // requested
    checkrequest(pr_i_demassist,checkprocesstmp,
    3,pr_m_crop,pr_s_crop,pr_i_fine); //[FvL]
  if (generalinput.process[pr_i_coregpm])                       // requested
    checkrequest(pr_i_coregpm,checkprocesstmp,
                 1,pr_i_fine);

  // this should go here...
  // BK 24-Aug-2000
  if (generalinput.process[pr_s_resample])                      // requested
    checkrequest(pr_s_resample,checkprocesstmp,
                 1,pr_i_coregpm);                               // required

  if (generalinput.process[pr_i_interfero])                     // requested
    checkrequest(pr_i_interfero,checkprocesstmp,
    1,pr_s_resample);
  if (generalinput.process[pr_i_coherence])                     // requested
    checkrequest(pr_i_coherence,checkprocesstmp,
    1,pr_s_resample);
  if (generalinput.process[pr_i_comprefpha])                    // requested
    checkrequest(pr_i_comprefpha,checkprocesstmp,
    2,pr_m_readfiles,pr_s_readfiles);                           // orbits required
  if (generalinput.process[pr_i_subtrrefpha])                   // requested
    checkrequest(pr_i_subtrrefpha,checkprocesstmp,
    1,pr_i_interfero);                                          // required
    // 2,pr_i_comprefpha,pr_i_interfero);                       // required
  if (generalinput.process[pr_i_comprefdem])                    // requested
    checkrequest(pr_i_comprefdem,checkprocesstmp,
    2,pr_m_readfiles,pr_s_readfiles);                           // orbits required
  if (generalinput.process[pr_i_subtrrefdem])                   // requested
    checkrequest(pr_i_subtrrefdem,checkprocesstmp,
    2,pr_i_comprefdem,pr_i_interfero);                          // required
  if (generalinput.process[pr_i_filtphase])                     // requested
    checkrequest(pr_i_filtphase,checkprocesstmp,
    1,pr_i_interfero);                                          // required
  if (generalinput.process[pr_i_unwrap])                        // requested
    checkrequest(pr_i_unwrap,checkprocesstmp,
    1,pr_i_interfero);                                          // required
  if (generalinput.process[pr_i_estorbits])                     // requested [HB]
    checkrequest(pr_i_estorbits,checkprocesstmp,                // required
		 1,pr_i_subtrrefpha);
  if (generalinput.process[pr_i_slant2h])                       // requested
    checkrequest(pr_i_slant2h,checkprocesstmp,
    1,pr_i_unwrap);                                             // required
  if (generalinput.process[pr_i_geocoding])                     // requested
    checkrequest(pr_i_geocoding,checkprocesstmp,
    1,pr_i_slant2h);                                            // required
  if (generalinput.process[pr_i_dinsar])                        // requested
    checkrequest(pr_i_dinsar,checkprocesstmp,
    1,pr_i_interfero);                                          // required
  if (generalinput.process[pr_i_EXTRA2])                        // requested
    checkrequest(pr_i_EXTRA2,checkprocesstmp,0);                // no requirements

// ______Tidy up______
  } // END checkprocessing


/****************************************************************
 *    checkrequest                                              *
 *                                                              *
 * Checks if requested processing step can be processed.        *
 *                                                              *
 * input:                                                       *
 *  - step to be processed                                      *
 *  - array of already processed steps                          *
 *  - optional: Number of req, and,                             *
 *                   steps that are required for this step      *
 * output:                                                      *
 *  - (updated) alreadyprocess[]                                *
 *                                                              *
 *                                                              *
 *    Bert Kampes, 16-Dec-1998                                  *
 ****************************************************************/
void checkrequest(
        int16 step,
        int16 alreadyprocess[NUMPROCESSES],
        ...)
  {
  TRACE_FUNCTION("checkrequest (BK 16-Dec-1998)")

// ______Check if step is already processed and update already______
  if (alreadyprocess[step])
    {
    ERROR << "Results of step: "
         << step << " (" << processcontrol[step]
         << ") already in result file.";
    ERROR << "\n TIP    : use \'doris.rmstep.sh\' to cleanup " <<  processcontrol[step]
          << " entries in result file." ;  // [MA] TODO report result filename file  or general.?
    PRINT_ERROR(ERROR.get_str())
    throw(input_error);
    }
  alreadyprocess[step]=1;                        // set to processed

  // ______Check with requirements______
  // ______ It seemed this has to be int instea of int16 on redhat 7.0 linux? ______
  // but why? BK 21-Nov-2000
  va_list arglist;                              // use ellipses
  va_start(arglist,alreadyprocess);
  /* *** SOME compiler required the second form, though it seems wrong,
     *** in order to compile doris comment out the second, and put a comment before
     *** the first form. */
  // seems that while passing '...' type is converted to int, so use that here...
  //int16 N = va_arg(arglist, int16);             // number of arguments=first ellipses
  int16 N = va_arg(arglist, int);

  int16 requiredstep;
  for (register int32 i=0; i<N; i++)
    {
    /* *** SOME compiler required the second form, though it seems wrong,
       *** in order to compile doris comment out the second, and put a comment before
       *** the first form. */
    //requiredstep = va_arg(arglist, int16);
    requiredstep = va_arg(arglist, int);

    if (!(alreadyprocess[requiredstep]))  // [MA] switched from warning to error to eliminate bugs.
      {
      ERROR << "Requested step: "
            << step << " (" << processcontrol[step]
            << ") seems impossible, because step "
            << requiredstep << " (" << processcontrol[requiredstep]
            << ") is not in resultfile.";
      PRINT_ERROR(ERROR.get_str())
      throw(input_error);// exit
      }
    }
  va_end(arglist);
  } // END checkrequest


/****************************************************************
 *    fillcheckprocess                                          *
 *                                                              *
 * Fills processing array for checking.                         *
 * Line comes from file, if processed is finished,              *
 *  checkprocess is filled                                      *
 * Looks for End_*:_NORMAL                                      *
 *                                                              *
 * input:                                                       *
 *  - line from result file                                     *
 *  - (empty) process control array                             *
 * output:                                                      *
 *  - (filled) process control array                            *
 *     normally m_pr* is filled if there is a choice            *
 *     except for pr_orbits: m_* =getorb ;s_* =readleader       *
 *     if there are more methods, the first one is flagged      *
 *                                                              *
 *    Bert Kampes, 16-Dec-1998                                  *
 ****************************************************************/
void fillcheckprocess(
        const char *line,
        int16 checkprocess[NUMPROCESSES+1],
        int16 fileid)
  {
  TRACE_FUNCTION("fillcheckprocess (BK 16-Dec-1998)")
  if (line[0] != '*')                                   // always start with *
    return;

  // ______ Set up compare strings ______
  // could be done smarter...
  char *endnormal[NUMPROCESSES];
  char mread[4*ONE27];
  strcpy(mread,"* End_");
  strcat(mread,processcontrol[pr_m_readfiles]);
  strcat(mread,"_NORMAL");
  endnormal[pr_m_readfiles]=&mread[0];

  char mcrop[4*ONE27];
  strcpy(mcrop,"* End_");
  strcat(mcrop,processcontrol[pr_m_crop]);
  strcat(mcrop,"_NORMAL");
  endnormal[pr_m_crop]=&mcrop[0];

//____RaffaeleNutricato START MODIFICATION SECTION 10
  char moversample[4*ONE27];
  strcpy(moversample,"* End_");
  strcat(moversample,processcontrol[pr_m_oversample]);
  strcat(moversample,"_NORMAL");
  endnormal[pr_m_oversample]=&moversample[0];
//____RaffaeleNutricato END MODIFICATION SECTION 10

  char mporbits[4*ONE27];
  strcpy(mporbits,"* End_");
  strcat(mporbits,processcontrol[pr_m_porbits]);
  strcat(mporbits,"_NORMAL");
  endnormal[pr_m_porbits]=&mporbits[0];

  char mmorbits[4*ONE27];              // [HB]
  strcpy(mmorbits,"* End_");
  strcat(mmorbits,processcontrol[pr_m_morbits]);
  strcat(mmorbits,"_NORMAL");
  endnormal[pr_m_morbits]=&mmorbits[0];

  char msimamp[4*ONE27];             // [MA]
  strcpy(msimamp,"* End_");
  strcat(msimamp,processcontrol[pr_m_simamp]);
  strcat(msimamp,"_NORMAL");
  endnormal[pr_m_simamp]=&msimamp[0];

  char mtiming[4*ONE27];             // [MA]
  strcpy(mtiming,"* End_");
  strcat(mtiming,processcontrol[pr_m_mtiming]);
  strcat(mtiming,"_NORMAL");
  endnormal[pr_m_mtiming]=&mtiming[0];

  char mfiltazi[4*ONE27];
  strcpy(mfiltazi,"* End_");
  strcat(mfiltazi,processcontrol[pr_m_filtazi]);
  strcat(mfiltazi,"_NORMAL");
  endnormal[pr_m_filtazi]=&mfiltazi[0];

  char mfiltrange[4*ONE27];
  strcpy(mfiltrange,"* End_");
  strcat(mfiltrange,processcontrol[pr_m_filtrange]);
  strcat(mfiltrange,"_NORMAL");
  endnormal[pr_m_filtrange]=&mfiltrange[0];

  char mEXTRA[4*ONE27];
  strcpy(mEXTRA,"* End_");
  strcat(mEXTRA,processcontrol[pr_m_EXTRA]);
  strcat(mEXTRA,"_NORMAL");
  endnormal[pr_m_EXTRA]=&mEXTRA[0];

  char sread[4*ONE27];
  strcpy(sread,"* End_");
  strcat(sread,processcontrol[pr_s_readfiles]);
  strcat(sread,"_NORMAL");
  endnormal[pr_s_readfiles]=&sread[0];

  char scrop[4*ONE27];
  strcpy(scrop,"* End_");
  strcat(scrop,processcontrol[pr_s_crop]);
  strcat(scrop,"_NORMAL");
  endnormal[pr_s_crop]=&scrop[0];

//____RaffaeleNutricato START MODIFICATION SECTION 11
  char soversample[4*ONE27];
  strcpy(soversample,"* End_");
  strcat(soversample,processcontrol[pr_s_oversample]);
  strcat(soversample,"_NORMAL");
  endnormal[pr_s_oversample]=&soversample[0];
//____RaffaeleNutricato END MODIFICATION SECTION 11

  char sporbits[4*ONE27];
  strcpy(sporbits,"* End_");
  strcat(sporbits,processcontrol[pr_s_porbits]);
  strcat(sporbits,"_NORMAL");
  endnormal[pr_s_porbits]=&sporbits[0];

  char smorbits[4*ONE27];              // [HB]
  strcpy(smorbits,"* End_");
  strcat(smorbits,processcontrol[pr_s_morbits]);
  strcat(smorbits,"_NORMAL");
  endnormal[pr_s_morbits]=&smorbits[0];

  char sfiltazi[4*ONE27];
  strcpy(sfiltazi,"* End_");
  strcat(sfiltazi,processcontrol[pr_s_filtazi]);
  strcat(sfiltazi,"_NORMAL");
  endnormal[pr_s_filtazi]=&sfiltazi[0];

  char sfiltrange[4*ONE27];
  strcpy(sfiltrange,"* End_");
  strcat(sfiltrange,processcontrol[pr_s_filtrange]);
  strcat(sfiltrange,"_NORMAL");
  endnormal[pr_s_filtrange]=&sfiltrange[0];

  char sresample[4*ONE27];
  strcpy(sresample,"* End_");
  strcat(sresample,processcontrol[pr_s_resample]);
  strcat(sresample,"_NORMAL");
  endnormal[pr_s_resample]=&sresample[0];

  char sEXTRA[4*ONE27];
  strcpy(sEXTRA,"* End_");
  strcat(sEXTRA,processcontrol[pr_s_EXTRA]);
  strcat(sEXTRA,"_NORMAL");
  endnormal[pr_s_EXTRA]=&sEXTRA[0];

  char icoarse[4*ONE27];
  strcpy(icoarse,"* End_");
  strcat(icoarse,processcontrol[pr_i_coarse]);
  strcat(icoarse,"_NORMAL");
  endnormal[pr_i_coarse]=&icoarse[0];

  char icoarse2[4*ONE27];
  strcpy(icoarse2,"* End_");
  strcat(icoarse2,processcontrol[pr_i_coarse2]);
  strcat(icoarse2,"_NORMAL");
  endnormal[pr_i_coarse2]=&icoarse2[0];

  char ifine[4*ONE27];
  strcpy(ifine,"* End_");
  strcat(ifine,processcontrol[pr_i_fine]);
  strcat(ifine,"_NORMAL");
  endnormal[pr_i_fine]=&ifine[0];

  char itiming[4*ONE27]; //[FvL]
  strcpy(itiming,"* End_");
  strcat(itiming,processcontrol[pr_i_timing]);
  strcat(itiming,"_NORMAL");
  endnormal[pr_i_timing]=&itiming[0];

  char idemassist[4*ONE27]; //[FvL]
  strcpy(idemassist,"* End_");
  strcat(idemassist,processcontrol[pr_i_demassist]);
  strcat(idemassist,"_NORMAL");
  endnormal[pr_i_demassist]=&idemassist[0];

  char icoregpm[4*ONE27];
  strcpy(icoregpm,"* End_");
  strcat(icoregpm,processcontrol[pr_i_coregpm]);
  strcat(icoregpm,"_NORMAL");
  endnormal[pr_i_coregpm]=&icoregpm[0];

  char iinterfero[4*ONE27];
  strcpy(iinterfero,"* End_");
  strcat(iinterfero,processcontrol[pr_i_interfero]);
  strcat(iinterfero,"_NORMAL");
  endnormal[pr_i_interfero]=&iinterfero[0];

  char icoherence[4*ONE27];
  strcpy(icoherence,"* End_");
  strcat(icoherence,processcontrol[pr_i_coherence]);
  strcat(icoherence,"_NORMAL");
  endnormal[pr_i_coherence]=&icoherence[0];

  char icomprefpha[4*ONE27];
  strcpy(icomprefpha,"* End_");
  strcat(icomprefpha,processcontrol[pr_i_comprefpha]);
  strcat(icomprefpha,"_NORMAL");
  endnormal[pr_i_comprefpha]=&icomprefpha[0];

  char isubtrrefpha[4*ONE27];
  strcpy(isubtrrefpha,"* End_");
  strcat(isubtrrefpha,processcontrol[pr_i_subtrrefpha]);
  strcat(isubtrrefpha,"_NORMAL");
  endnormal[pr_i_subtrrefpha]=&isubtrrefpha[0];

  char icomprefdem[4*ONE27];
  strcpy(icomprefdem,"* End_");
  strcat(icomprefdem,processcontrol[pr_i_comprefdem]);
  strcat(icomprefdem,"_NORMAL");
  endnormal[pr_i_comprefdem]=&icomprefdem[0];

  char isubtrrefdem[4*ONE27];
  strcpy(isubtrrefdem,"* End_");
  strcat(isubtrrefdem,processcontrol[pr_i_subtrrefdem]);
  strcat(isubtrrefdem,"_NORMAL");
  endnormal[pr_i_subtrrefdem]=&isubtrrefdem[0];

  char ifiltphase[4*ONE27];
  strcpy(ifiltphase,"* End_");
  strcat(ifiltphase,processcontrol[pr_i_filtphase]);
  strcat(ifiltphase,"_NORMAL");
  endnormal[pr_i_filtphase]=&ifiltphase[0];

  char iunwrap[4*ONE27];
  strcpy(iunwrap,"* End_");
  strcat(iunwrap,processcontrol[pr_i_unwrap]);
  strcat(iunwrap,"_NORMAL");
  endnormal[pr_i_unwrap]=&iunwrap[0];

  char iestorbits[4*ONE27];                     // [HB]
  strcpy(iestorbits,"* End_");
  strcat(iestorbits,processcontrol[pr_i_estorbits]);
  strcat(iestorbits,"_NORMAL");
  endnormal[pr_i_estorbits]=&iestorbits[0];

  char islant2h[4*ONE27];
  strcpy(islant2h,"* End_");
  strcat(islant2h,processcontrol[pr_i_slant2h]);
  strcat(islant2h,"_NORMAL");
  endnormal[pr_i_slant2h]=&islant2h[0];

  char igeocoding[4*ONE27];
  strcpy(igeocoding,"* End_");
  strcat(igeocoding,processcontrol[pr_i_geocoding]);
  strcat(igeocoding,"_NORMAL");
  endnormal[pr_i_geocoding]=&igeocoding[0];

  char idinsar[4*ONE27];
  strcpy(idinsar,"* End_");
  strcat(idinsar,processcontrol[pr_i_dinsar]);
  strcat(idinsar,"_NORMAL");
  endnormal[pr_i_dinsar]=&idinsar[0];

  char iEXTRA2[4*ONE27];
  strcpy(iEXTRA2,"* End_");
  strcat(iEXTRA2,processcontrol[pr_i_EXTRA2]);
  strcat(iEXTRA2,"_NORMAL");
  endnormal[pr_i_EXTRA2]=&iEXTRA2[0];

// ======Fill checkprocessarray======
  if (fileid == MASTERID)
    {
    //<< "\n* End_" << processcontrol[pr_i_coarse] << "_NORMAL"
    //if (!strcmp(line,"* End_readfiles:_NORMAL"))
    if (!strcmp(line,endnormal[pr_m_readfiles]))
      checkprocess[pr_m_readfiles]=1;
    else if (!strcmp(line,"* End_leader_datapoints:_NORMAL"))   // fixed string...
      checkprocess[NUMPROCESSES]=1;
    else if (!strcmp(line,endnormal[pr_m_porbits]))
      checkprocess[pr_m_porbits]=1;
    else if (!strcmp(line,endnormal[pr_m_morbits]))            // [HB]
      checkprocess[pr_m_morbits]=1;
    else if (!strcmp(line,endnormal[pr_m_crop]))
      checkprocess[pr_m_crop]=1;
    else if (!strcmp(line,endnormal[pr_m_simamp]))             // [MA]
      checkprocess[pr_m_simamp]=1;
    else if (!strcmp(line,endnormal[pr_m_mtiming]))            // [MA]
      checkprocess[pr_m_mtiming]=1;
    else if (!strcmp(line,endnormal[pr_m_oversample]))
      checkprocess[pr_m_oversample]=1;
    else if (!strcmp(line,endnormal[pr_m_filtazi]))
      checkprocess[pr_m_filtazi]=1;
    else if (!strcmp(line,endnormal[pr_m_filtrange]))
      checkprocess[pr_m_filtrange]=1;
    else if (!strcmp(line,endnormal[pr_m_EXTRA]))
      checkprocess[pr_m_EXTRA]=1;
    }

  else if (fileid == SLAVEID)
    {
    if (!strcmp(line,endnormal[pr_s_readfiles]))
      checkprocess[pr_s_readfiles]=1;
    else if (!strcmp(line,"* End_leader_datapoints:_NORMAL"))   // fixed string...
      checkprocess[NUMPROCESSES]=1;
    else if (!strcmp(line,endnormal[pr_s_porbits]))
      checkprocess[pr_s_porbits]=1;
    else if (!strcmp(line,endnormal[pr_s_morbits]))            // [HB]
      checkprocess[pr_s_morbits]=1;
    else if (!strcmp(line,endnormal[pr_s_crop]))
      checkprocess[pr_s_crop]=1;
    else if (!strcmp(line,endnormal[pr_s_oversample]))
      checkprocess[pr_s_oversample]=1;
    else if (!strcmp(line,endnormal[pr_s_filtazi]))
      checkprocess[pr_s_filtazi]=1;
    else if (!strcmp(line,endnormal[pr_s_filtrange]))
      checkprocess[pr_s_filtrange]=1;
    else if (!strcmp(line,endnormal[pr_s_resample]))
      checkprocess[pr_s_resample]=1;
    else if (!strcmp(line,endnormal[pr_s_EXTRA]))
      checkprocess[pr_s_EXTRA]=1;
    }

  else if (fileid == INTERFID)
    {
    // ______Interferogram______
    if (!strcmp(line,endnormal[pr_i_coarse]))
      checkprocess[pr_i_coarse]=1;
    else if (!strcmp(line,endnormal[pr_i_coarse2]))
      checkprocess[pr_i_coarse2]=1;
    else if (!strcmp(line,endnormal[pr_i_fine]))
      checkprocess[pr_i_fine]=1;
    else if (!strcmp(line,endnormal[pr_i_timing])) //[FvL]
      checkprocess[pr_i_timing]=1;
    else if (!strcmp(line,endnormal[pr_i_demassist])) //[FvL]
      checkprocess[pr_i_demassist]=1;
    else if (!strcmp(line,endnormal[pr_i_coregpm]))
      checkprocess[pr_i_coregpm]=1;
    else if (!strcmp(line,endnormal[pr_i_interfero]))
      checkprocess[pr_i_interfero]=1;
    else if (!strcmp(line,endnormal[pr_i_coherence]))
      checkprocess[pr_i_coherence]=1;
    else if (!strcmp(line,endnormal[pr_i_comprefpha]))
      checkprocess[pr_i_comprefpha]=1;
    else if (!strcmp(line,endnormal[pr_i_subtrrefpha]))
      checkprocess[pr_i_subtrrefpha]=1;
    else if (!strcmp(line,endnormal[pr_i_comprefdem]))
      checkprocess[pr_i_comprefdem]=1;
    else if (!strcmp(line,endnormal[pr_i_subtrrefdem]))
      checkprocess[pr_i_subtrrefdem]=1;
    else if (!strcmp(line,endnormal[pr_i_filtphase]))
      checkprocess[pr_i_filtphase]=1;
    else if (!strcmp(line,endnormal[pr_i_unwrap]))
      checkprocess[pr_i_unwrap]=1;
    else if (!strcmp(line,endnormal[pr_i_estorbits]))  // [HB]
      checkprocess[pr_i_estorbits]=1;
    else if (!strcmp(line,endnormal[pr_i_slant2h]))
      checkprocess[pr_i_slant2h]=1;
    else if (!strcmp(line,endnormal[pr_i_geocoding]))
      checkprocess[pr_i_geocoding]=1;
    else if (!strcmp(line,endnormal[pr_i_dinsar]))
      checkprocess[pr_i_dinsar]=1;
    else if (!strcmp(line,endnormal[pr_i_EXTRA2]))
      checkprocess[pr_i_EXTRA2]=1;
    }
  else
    {
    PRINT_ERROR("wrong input.")
    throw(unhandled_case_error);
    }
  } // END fillcheckprocess


/****************************************************************
 *    fillprocessed                                             *
 *                                                              *
 * Fills processing array for checking.                         *
 * Line comes from file, if processed is finished,              *
 *  checkprocess is filled                                      *
 * Reads process controls:       1/0                            *
 *                                                              *
 * input:                                                       *
 *  - line from result file                                     *
 *  - (empty) process control array                             *
 *                                                              *
 * output:                                                      *
 *  - (filled) process control array                            *
 *     normally m_pr* is filled if there is a choice            *
 *     except for pr_orbits: m_* =getorb ;s_* =readleader       *
 *     if there are more methods, the first one is flagged      *
 *                                                              *
 * See documentation page 14                                    *
 *    Bert Kampes, 16-Dec-1998                                  *
 ****************************************************************/
void fillprocessed(
        const char *line,
        int16 checkprocess[NUMPROCESSES+1],
        int16 fileid)
  {
  TRACE_FUNCTION("fillprocessed (BK 16-Dec-1998)")
  bool  space=false;                            // line must contain 2 words
  char  ch='\t';                                // isspace
  int32         tin;

  int32 linesz = strlen(line);                  // w/o \0
  if (linesz > 30)
    return;
  char  word[4*ONE27];                            // should be enough
// ______ Disect line ______
  for (register int32 i=0;i<linesz;i++)
    {
    ch = line[i];
    word[i]=ch;
    if (isspace(tin=ch))
      {
      space = true;
      if (line[i-1]==':')                               // process_control: 1/0
        {
        word[i]='\0';                                   // replace space by \0
        break;                                          // for
        }
      else
        {
        return;
        }
      }
    }
  if (!space) return;                           // must be a space in line

  int16 processflag;
  if (line[linesz-1]=='0')
    processflag=0;
  else if (line[linesz-1]=='1')
    processflag=1;
  else
    return;

  // ====== Fill process control ======
  switch (fileid)
    {
    case MASTERID:
      if      (!strcmp(word,processcontrol[pr_m_readfiles]))
        checkprocess[pr_m_readfiles]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_porbits]))
        checkprocess[pr_m_porbits]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_morbits]))       // [HB]
        checkprocess[pr_m_morbits]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_crop]))
        checkprocess[pr_m_crop]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_simamp]))        // [MA]
        checkprocess[pr_m_simamp]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_mtiming]))       // [MA]
        checkprocess[pr_m_mtiming]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_oversample]))
        checkprocess[pr_m_oversample]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_filtazi]))
        checkprocess[pr_m_filtazi]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_filtrange]))
        checkprocess[pr_m_filtrange]=processflag;
      else if (!strcmp(word,processcontrol[pr_m_EXTRA]))
        checkprocess[pr_m_EXTRA]=processflag;
      else
        {
         DEBUG << "Not a correct entry: \"" << word << "\" in master process control (ignored)"; // [MA]
         DEBUG.print();
        }
      break;

    case SLAVEID:
      if      (!strcmp(word,processcontrol[pr_s_readfiles]))
        checkprocess[pr_s_readfiles]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_porbits]))
        checkprocess[pr_s_porbits]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_morbits]))   // [HB]
        checkprocess[pr_s_morbits]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_crop]))
        checkprocess[pr_s_crop]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_oversample]))
        checkprocess[pr_s_oversample]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_filtazi]))
        checkprocess[pr_s_filtazi]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_filtrange]))
        checkprocess[pr_s_filtrange]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_resample]))
        checkprocess[pr_s_resample]=processflag;
      else if (!strcmp(word,processcontrol[pr_s_EXTRA]))
        checkprocess[pr_s_EXTRA]=processflag;
      else
        {
         DEBUG << "Not a correct entry: \"" << word << "\" in slave process control (ignored)"; // [MA]
         DEBUG.print();
        }
      break;

    // ______ Interferogram ______
    case INTERFID:
      if      (!strcmp(word,processcontrol[pr_i_coarse]))
        checkprocess[pr_i_coarse]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_coarse2]))
        checkprocess[pr_i_coarse2]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_fine]))
        checkprocess[pr_i_fine]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_timing])) //[FvL]
        checkprocess[pr_i_timing]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_demassist])) //[FvL]
        checkprocess[pr_i_demassist]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_coregpm]))
        checkprocess[pr_i_coregpm]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_interfero]))
        checkprocess[pr_i_interfero]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_coherence]))
        checkprocess[pr_i_coherence]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_comprefpha]))
        checkprocess[pr_i_comprefpha]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_subtrrefpha]))
        checkprocess[pr_i_subtrrefpha]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_comprefdem]))
        checkprocess[pr_i_comprefdem]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_subtrrefdem]))
        checkprocess[pr_i_subtrrefdem]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_filtphase]))
        checkprocess[pr_i_filtphase]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_unwrap]))
        checkprocess[pr_i_unwrap]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_estorbits])) // [HB]
        checkprocess[pr_i_estorbits]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_slant2h]))
        checkprocess[pr_i_slant2h]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_geocoding]))
        checkprocess[pr_i_geocoding]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_dinsar]))
        checkprocess[pr_i_dinsar]=processflag;
      else if (!strcmp(word,processcontrol[pr_i_EXTRA2]))
        checkprocess[pr_i_EXTRA2]=processflag;
      else
        {
         DEBUG << "Not a correct entry: \"" << word << "\" in interferrogram process control (ignored)"; // [MA]
         DEBUG.print();
        }
      break;

    default:
      PRINT_ERROR("wrong input.")
      throw(unhandled_case_error);
    } // switch fileid
  } // END fillprocessed



/****************************************************************
 *    filelines                                                 *
 *                                                              *
 * Returns number of \n (enters).                               *
 * input:                                                       *
 *  - filename                                                  *
 * output:                                                      *
 *  - number of eols.                                           *
 *                                                              *
 *    Bert Kampes, 10-Jun-1999                                  *
 ****************************************************************/
int32 filelines(
        const char *file)
  {
  TRACE_FUNCTION("filelines (BK 10-Jun-1999)")
  char dummyline[4*ONE27];
  register int32 nlines = 0;
  //ifstream ifile(file, ios::in | ios::nocreate);
  ifstream ifile(file, ios::in);
  bk_assert(ifile,file,__FILE__,__LINE__);
  while (ifile)
    {
    ifile.getline(dummyline,4*ONE27,'\n');
    nlines++;
    }
  ifile.clear();// g++ v3.2 may have a problem else
  ifile.close();
  return nlines-1;
  }  // END filelines



/****************************************************************
 *    existed                                                   *
 *                                                              *
 * Returns true if file exists.                                 *
 * input:                                                       *
 *  - filename                                                  *
 * output:                                                      *
 *  - true if existi, false if not.                             *
 *                                                              *
 *    Bert Kampes, 24-Dec-1998                                  *
 ****************************************************************/
bool existed(
        const char *file)
  {
  bool exists = false;
  ifstream ifile(file, ios::in);
  if (ifile)
    {
    exists = true;
    DEBUG << "file: " << file << " does exist.";
    }
  else
    {
    DEBUG << "file: " << file << " does not exist.";
    }
  DEBUG.print();

  ifile.close();
  return exists;
  } // END exist



/****************************************************************
 *    removedatleader                                           *
 *                                                              *
 * Part containing datapoints from leader is removed            *
 *   from resultfile. (because precise orbits will be there)    *
 *                                                              *
 * input:                                                       *
 *  - file                                                      *
 * output:                                                      *
 *  - void, updated logfile,                                    *
 *                                                              *
 *    Bert Kampes, 07-Jan-1999                                  *
 ****************************************************************/
void removedatleader(
        const char* file)
  {
  TRACE_FUNCTION("removedatleader (BK 07-Jan-1999)")
  // ______Copy relevant part to tmp file______
  // ______ then rename file to old name ______
  //ifstream ifile(file, ios::in | ios::nocreate);
  ifstream ifile(file, ios::in);
  bk_assert(ifile,file,__FILE__,__LINE__);
  ofstream tmpfile("scratchtmp", ios::out | ios::trunc);
  bk_assert(tmpfile,"removedatleader: scratchtmp",__FILE__,__LINE__);

  // ______Copy upto data section______
  char                  dummyline[4*ONE27];
  ifile.getline(dummyline,4*ONE27,'\n');          // to prevent twice last line
  while (ifile && strcmp(dummyline,"*_Start_leader_datapoints")) // fixed string...
    {
    tmpfile << dummyline << endl;
    ifile.getline(dummyline,4*ONE27,'\n');
    }
  if (!ifile) return;                                   // file eof?

  // ______Dont copy data section______
  while (ifile && strcmp(dummyline,"* End_leader_datapoints:_NORMAL"))  // fixed string..
    {
    ifile.getline(dummyline,4*ONE27,'\n');
    }
  ifile.getline(dummyline,4*ONE27,'\n');

  // ______Copy rest of file______
  while (ifile)
    {
    tmpfile << dummyline << endl;
    ifile.getline(dummyline,4*ONE27,'\n');
    }

  // ______Tidy up______
  ifile.close();
  tmpfile.close();
  if (rename("scratchtmp",file))                        // rename file
    {
    ERROR << "Could not rename file: scratchtmp to: " << file;
    PRINT_ERROR(ERROR.get_str())
    throw(file_error);
    }
  } // END removedatleader



/****************************************************************
 *    filesize                                                  *
 *                                                              *
 * Returns uint filesize of (closed) file in bytes.             *
 * input:                                                       *
 *  - filename                                                  *
 * output:                                                      *
 *  - size (B)                                                  *
 *                                                              *
 *    Bert Kampes, 08-Jan-1999                                  *
 ****************************************************************/
// inline uint filesize(
uint filesize(
        const char* file)
  {
  TRACE_FUNCTION("filesize (BK 08-Jan-1999)")
  //ifstream ifile(file, ios::in | ios::ate | ios::nocreate);
  // problem of g++ with ate, app does work ok, but
  //ifstream ifile(file, ios::in | ios::nocreate);
  ifstream ifile(file, ios::in);
#ifdef __DEBUG
  bk_assert(ifile,file,__FILE__,__LINE__);
#endif
  ifile.seekg(0,ios::end);                      // pointer to end...
  // closed automatically
  return ifile.tellg();
  } // END filesize



/****************************************************************
 * getoverlap                                                   *
 *                                                              *
 * overlap of 2 windows in same coord. system                   *
 #%// BK 21-Sep-2000
 ****************************************************************/
window getoverlap(
        const window &master,
        const window &slave)
  {
  TRACE_FUNCTION("getoverlap (BK 10-Mar-1999)")
  window overlap = slave;
  if (master.linelo>overlap.linelo) overlap.linelo=master.linelo;
  if (master.linehi<overlap.linehi) overlap.linehi=master.linehi;
  if (master.pixlo >overlap.pixlo)  overlap.pixlo =master.pixlo;
  if (master.pixhi <overlap.pixhi)  overlap.pixhi =master.pixhi;
  return overlap;
  } // getoverlap



/****************************************************************
 * getoverlap                                                   *
 *                                                              *
 * compute approx. rectangular overlap (master coord.)          *
 * between master/slave with help of transformation polynomial  *
 * Bert Kampes 10-Mar-99                                        *
 ****************************************************************/
window getoverlap(
        const slcimage      &master,
        const slcimage      &slave,
        const matrix<real8> &cpmL,
        const matrix<real8> &cpmP)
  {
  TRACE_FUNCTION("getoverlap (BK 10-Mar-1999)")

// ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;

  INFO << "getoverlap: polynomial normalized by factors: "
       << minL << " " << maxL << " " << minP << " " << maxP
       << " to [-2,2]";
  INFO.print();

// ______offset = A(slave system) - A(master system)______
// ______Corners of slave in master system______
// ______Offsets for slave corners (approx.)______
// ______ approx: defined as offset = f(l,p)_M in master system not slave.
  real8 approxoffL = cpmL(0,0);                         // zero order term;
  real8 approxoffP = cpmP(0,0);                         // zero order term;

//  real8 sL00 = slave.currentwindow.linelo -
//                polyval(slave.currentwindow.linelo - approxoffL,
//                        slave.currentwindow.pixlo  - approxoffP, cpmL);
// ______ Use normalized polynomial ______
  const real8 sL00 = slave.currentwindow.linelo -
       polyval(normalize(real8(slave.currentwindow.linelo)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixlo) -approxoffP,minP,maxP),
               cpmL);
  const real8 sP00 = slave.currentwindow.pixlo -
       polyval(normalize(real8(slave.currentwindow.linelo)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixlo) -approxoffP,minP,maxP),
               cpmP);
  const real8 sL0N = slave.currentwindow.linelo -
       polyval(normalize(real8(slave.currentwindow.linelo)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixhi) -approxoffP,minP,maxP),
               cpmL);
  const real8 sP0N = slave.currentwindow.pixhi -
       polyval(normalize(real8(slave.currentwindow.linelo)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixhi) -approxoffP,minP,maxP),
               cpmP);
  const real8 sLN0 = slave.currentwindow.linehi -
       polyval(normalize(real8(slave.currentwindow.linehi)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixlo) -approxoffP,minP,maxP),
               cpmL);
  const real8 sPN0 = slave.currentwindow.pixlo -
       polyval(normalize(real8(slave.currentwindow.linehi)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixlo) -approxoffP,minP,maxP),
               cpmP);
  const real8 sLNN = slave.currentwindow.linehi -
       polyval(normalize(real8(slave.currentwindow.linehi)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixhi) -approxoffP,minP,maxP),
               cpmL);
  const real8 sPNN = slave.currentwindow.pixhi -
       polyval(normalize(real8(slave.currentwindow.linehi)-approxoffL,minL,maxL),
               normalize(real8(slave.currentwindow.pixhi) -approxoffP,minP,maxP),
               cpmP);


//  // ______make window rectangular + int______
//  // ______window is uint type______
//    if (sL00 < 0.) sL00 = 0.;
//    if (sP00 < 0.) sP00 = 0.;
//    window slaveinmaster = {ceil(max(sL00,sL0N)),
//                            floor(min(sLN0,sLNN)),
//                            ceil(max(sP00,sPN0)),
//                            floor(min(sP0N,sPNN))};


  // ______Corners of overlap master,slave in master system______
  window win;
  win.linelo = max(int32(master.currentwindow.linelo),
                   int32(ceil(max(sL00,sL0N))));
  win.linehi = min(int32(master.currentwindow.linehi),
                   int32(floor(min(sLN0,sLNN))));
  win.pixlo  = max(int32(master.currentwindow.pixlo),
                   int32(ceil(max(sP00,sPN0))));
  win.pixhi  = min(int32(master.currentwindow.pixhi),
                   int32(floor(min(sP0N,sPNN))));

#ifdef __DEBUG
  DEBUG.print("Finished getoverlap");
  DEBUG.print("Approximate overlap master/slave window:");
  win.disp();
#endif
  return win;
  } // END getoverlap


/****************************************************************
 * getoverlap                                                   *
 *                                                              *
 * compute rectangular overlap between master and slave         *
 * (master coordinate system)                                   *
 * Freek van Leijen, 22-SEP-2007                                *
 ****************************************************************/
window getoverlap(
        const slcimage      &master,
        const slcimage      &slave,
        const real8         &Npointsd2,
        const real8         &timing_L,
        const real8         &timing_P)
  {
    TRACE_FUNCTION("getoverlap (FvL 22-SEP-07)")

    real8 ml0 = master.currentwindow.linelo;
    real8 mlN = master.currentwindow.linehi;
    real8 mp0 = master.currentwindow.pixlo;
    real8 mpN = master.currentwindow.pixhi;

    real8 sl00 = slave.currentwindow.linelo+slave.slavemasteroffsets.l00+Npointsd2-timing_L;
    real8 sp00 = slave.currentwindow.pixlo+slave.slavemasteroffsets.p00+Npointsd2-timing_P;
    real8 sl0N = slave.currentwindow.linelo+slave.slavemasteroffsets.l0N+Npointsd2-timing_L;
    real8 sp0N = slave.currentwindow.pixhi+slave.slavemasteroffsets.p0N-Npointsd2-timing_P;
    real8 slN0 = slave.currentwindow.linehi+slave.slavemasteroffsets.lN0-Npointsd2-timing_L;
    real8 spN0 = slave.currentwindow.pixlo+slave.slavemasteroffsets.pN0+Npointsd2-timing_P;
    real8 slNN = slave.currentwindow.linehi+slave.slavemasteroffsets.lNN-Npointsd2-timing_L;
    real8 spNN = slave.currentwindow.pixhi+slave.slavemasteroffsets.pNN-Npointsd2-timing_P;

    matrix<real8> mh1sv1(2,1), mh1sv2(2,1), mh2sv1(2,1), mh2sv2(2,1),
      mv1sh1(2,1), mv1sh2(2,1), mv2sh1(2,1), mv2sh2(2,1);
    lineintersect(ml0,mp0,ml0,mpN,sl00,sp00,slN0,spN0,mh1sv1);
    lineintersect(ml0,mp0,ml0,mpN,sl0N,sp0N,slNN,spNN,mh1sv2);
    lineintersect(mlN,mp0,mlN,mpN,sl00,sp00,slN0,spN0,mh2sv1);
    lineintersect(mlN,mp0,mlN,mpN,sl0N,sp0N,slNN,spNN,mh2sv2);
    lineintersect(ml0,mp0,mlN,mp0,sl00,sp00,sl0N,sp0N,mv1sh1);
    lineintersect(ml0,mp0,mlN,mp0,slN0,spN0,slNN,spNN,mv1sh2);
    lineintersect(ml0,mpN,mlN,mpN,sl00,sp00,sl0N,sp0N,mv2sh1);
    lineintersect(ml0,mpN,mlN,mpN,slN0,spN0,slNN,spNN,mv2sh2);

    real8 overlap_l0 = max(max(max(max(max(max(ml0,sl00),sl0N),mh1sv1(0,0)),mh1sv2(0,0)),mv1sh1(0,0)),mv2sh1(0,0));
    real8 overlap_p0 = max(max(max(max(max(max(mp0,sp00),spN0),mh1sv1(1,0)),mh2sv1(1,0)),mv1sh1(1,0)),mv1sh2(1,0));
    real8 overlap_lN = min(min(min(min(min(min(mlN,slN0),slNN),mh2sv1(0,0)),mh2sv2(0,0)),mv1sh2(0,0)),mv2sh2(0,0));
    real8 overlap_pN = min(min(min(min(min(min(mpN,sp0N),spNN),mh1sv2(1,0)),mh2sv2(1,0)),mv2sh1(1,0)),mv2sh2(1,0));

    // ______Corners of overlap master,slave in master system______
    window overlap;
    overlap.linelo = int32(ceil(overlap_l0));
    overlap.linehi = int32(floor(overlap_lN));
    overlap.pixlo = int32(ceil(overlap_p0));
    overlap.pixhi = int32(floor(overlap_pN));

    return overlap;
  } // END getoverlap


/****************************************************************
 * lineintersect                                                *
 *                                                              *
 * compute intersection point of two line segments              *
 * (master coordinate system)                                   *
 * Freek van Leijen, 22-SEP-2007                                *
 ****************************************************************/
void lineintersect(
                   const real8   &ax,
                   const real8   &ay,
                   const real8   &bx,
                   const real8   &by,
                   const real8   &cx,
                   const real8   &cy,
                   const real8   &dx,
                   const real8   &dy,
                   matrix<real8> &exy)
      {
        TRACE_FUNCTION("lineintersect (FvL 22-SEP-2007)")

        real8 u1 = bx-ax;
        real8 u2 = by-ay;
        real8 v1 = dx-cx;
        real8 v2 = dy-cy;
        real8 w1 = ax-cx;
        real8 w2 = ay-cy;

        real8 s = (v2*w1-v1*w2)/(v1*u2-v2*u1);
        exy(0,0) = ax+s*u1;
        exy(1,0) = ay+s*u2;
      } // END lineintersect


/****************************************************************
 *    readcoeff                                                 *
 *                                                              *
 * Pattern is searched in file (1st word),                      *
 * After Pattern the Ncoeffs must follow                        *
 * next (Ncoeff) lines are assumed to contain the coefficients  *
 *                                                              *
 * e.g.: readcoeff(resfile,"Degree_flat:",9)                    *
 *                                                              *
 * input:                                                       *
 *  - file name to search                                       *
 *  - pattern to search for                                     *
 *                                                              *
 * output:                                                      *
 *  - matrix<real8> coefficients(Nc , 1)                        *
 *                                                              *
 *    Bert Kampes, 12-Mar-1999                                  *
 ****************************************************************/
matrix<real8> readcoeff(
        const char* file,
        const char* pattern,
        const int16 Ncoeff)
  {
  TRACE_FUNCTION("readcoeff (BK 12-Mar-1999)")
  char                  dummyline[ONE27];
  char                  word[EIGHTY];
  bool                  foundword = false;
  matrix<real8>         coeffs(Ncoeff,1);               // store coefficients

  ifstream infile(file, ios::in);
  bk_assert(infile,file,__FILE__,__LINE__);

  // ====== Search infile ======
  while (infile)
    {
    infile >> word;
    if (strcmp(pattern,word))                           // no pattern match.
      {
      infile.getline(dummyline,4*ONE27,'\n');             // goto next line.
      }
    else                                                // pattern match.
      {
      infile.getline(dummyline,ONE27,'\n');             // start data
      foundword = true;
      for (register int32 i=0;i<Ncoeff;i++)
        {
        infile >> coeffs(i,0);
        infile.getline(dummyline,ONE27,'\n');           // goto next line.
        }
      break;                                            // file
      }                                                 // else
    }                                                   // file
  infile.close();

  if (!foundword)
    {
    ERROR << "readcoeff: file: " << file
         << ": could not find string \"" << pattern << "\".";
    PRINT_ERROR(ERROR.get_str());
    throw(file_error);
    }
  else
    {
    INFO << "read: " << Ncoeff << " coefficients after: \""
         << pattern << "\"";
    INFO.print();
    }
  return coeffs;
  } // END readcoeff


/****************************************************************
 *    readcoeff                                                 *
 *                                                              *
 * Pattern is searched in file (1st word),                      *
 * After Pattern the Ncoeffs must follow                        *
 * next (Ncoeff) lines are assumed to contain the coefficients  *
 *                                                              *
 * e.g.: readcoeff(resfile,"Degree_flat:",9)                    *
 *                                                              *
 * input:                                                       *
 *  - file name to search                                       *
 *  - pattern to search for                                     *
 *                                                              *
 * output:                                                      *
 *  - matrix<real8> coefficients(Nc , 1)                        *
 *                                                              *
 *    Bert Kampes, 12-Mar-1999                                  *
 ****************************************************************/
matrix<real8> readnormcoeff(
        const char* file,
        const char* pattern)
  {
  TRACE_FUNCTION("readcoeff (BK 12-Mar-1999)")
  char                  dummyline[ONE27];
  char                  word[EIGHTY];
  bool                  foundword = false;
  matrix<real8>         coeffs(2,1);               // store coefficients

  ifstream infile(file, ios::in);
  bk_assert(infile,file,__FILE__,__LINE__);

  // ====== Search infile ======
  while (infile)
    {
    infile >> word;
    if (strcmp(pattern,word))                           // no pattern match.
      {
      infile.getline(dummyline,4*ONE27,'\n');             // goto next line.
      }
    else                                                // pattern match.
      {
     
        foundword = true;
     
        infile >> coeffs(0,0)>>coeffs(1,0);
     
      break;    
      break;                                            // file
      }                                                 // else
    }                                                   // file
  infile.close();

  if (!foundword)
    {
    ERROR << "readcoeff: file: " << file
         << ": could not find string \"" << pattern << "\".";
    PRINT_ERROR(ERROR.get_str());
    throw(file_error);
    }
  else
    {
    INFO << "read: " << 2 << " coefficients after: \""
         << pattern << "\"";
    INFO.print();
    }
  return coeffs;
  } // END readcoeff


/****************************************************************
 *    openfstream                                               *
 *                                                              *
 *  Open an input file stream by name                           *
 * input:                                                       *
 *  - fstream                                                   *
 *  - filename                                                  *
 * output:                                                      *
 *  - opened file, check yourself with assert                   *
 *                                                              *
 *    Bert Kampes, 11-Sep-2004                                  *
 ****************************************************************/
void openfstream(
        ifstream &stream,
        const char* ifilename)
  {
  TRACE_FUNCTION("openfstream (BK 11-Sep-2004)")
  INFO << "Opening input file: " << ifilename;
  INFO.print();
  DEBUG << "Opening input file: " << ifilename;
  DEBUG.print();
  #ifdef __NO_IOS_BINARY__
  stream.open(ifilename, ios::in);
  #else
  stream.open(ifilename, ios::in | ios::binary);
  #endif
  } // END openfstream



/****************************************************************
 *    openfstream                                               *
 *                                                              *
 *  Open an output file stream by name                          *
 * input:                                                       *
 *  - ofstream                                                  *
 *  - filename                                                  *
 *  - overwrite bool                                            *
 * output:                                                      *
 *  - opened file, check yourself with assert                   *
 *                                                              *
 *    Bert Kampes, 11-Sep-2004                                  *
 ****************************************************************/
void openfstream(
        ofstream &stream,
        const char* ofilename,
        const bool overwrite)
  {
  TRACE_FUNCTION("openfstream (BK 11-Sep-2004)")
  DEBUG << "Opening output file: " << ofilename;
  DEBUG.print();
  if (overwrite == true)
    {
    #ifdef __NO_IOS_BINARY__
    stream.open(ofilename, ios::out | ios::trunc);
    #else
    stream.open(ofilename, ios::out | ios::binary | ios::trunc);
    #endif
    }
  else
    {
    if (existed(ofilename) == true)
      {
      ERROR << "output file \": "
            << ofilename << "\" exists, use OVERWRITE ON";
      PRINT_ERROR(ERROR.get_str())
      throw(file_error);
      }
    else
      {
      #ifdef __NO_IOS_BINARY__
      stream.open(ofilename, ios::out);
      #else
      stream.open(ofilename, ios::out | ios::binary);
      #endif
      }
    }
  } // END openfstream



/****************************************************************
 *    bk_assert                                                 *
 *                                                              *
 *  Check opened input fstream                                  *
 * input:                                                       *
 *  - fstream                                                   *
 *  - filedescription                                           *
 *  - __FILE__                                                  *
 *  - __LINE__                                                  *
 * output:                                                      *
 *  - exit or ok                                                *
 *                                                              *
 *    Bert Kampes, 04-Jun-1999                                  *
 ****************************************************************/
void bk_assert(
        const ifstream &stream,
        const char* ofilename,
        const char* callingfile,
        int32 linenumber)
  {
  TRACE_FUNCTION("bk_assert (BK 04-Jun-1999)")
  DEBUG << "Opened input file:  " << ofilename;
  DEBUG.print();
  if (!stream)
    {
    ERROR << "Unable to open input file: "      << ofilename
         << " in source file: "                 << callingfile
         << " at line: "                        << linenumber;
    PRINT_ERROR(ERROR.get_str());
    ERROR << "Missing a required file: "      << ofilename;  // [MA] 201004
    PRINT_ERROR(ERROR.get_str());
    throw(file_error);
    }
  } // END bk_assert



/****************************************************************
 *    bk_assert                                                 *
 *                                                              *
 *  Check opened output fstream                                 *
 * input:                                                       *
 *  - fstream                                                   *
 *  - filedescription                                           *
 *  - __FILE__                                                  *
 *  - __LINE__                                                  *
 * output:                                                      *
 *  - exit or ok                                                *
 *                                                              *
 *    Bert Kampes, 04-Jun-1999                                  *
 ****************************************************************/
void bk_assert(
        const ofstream &stream,
        const char* ofilename,
        const char* callingfile,
        int32 linenumber)
  {
  TRACE_FUNCTION("bk_assert (BK 04-Jun-1999)")
  DEBUG << "Opened output file: " << ofilename;
  DEBUG.print();
  if (!stream)
    {
    ERROR << "Unable to open output file: "     << ofilename
         << " in source file: "                 << callingfile
         << " at line: "                        << linenumber
         << " (OVERWRIT OFF/non existing directory?)";
    PRINT_ERROR(ERROR.get_str());
    throw(file_error);
    }
  } // END bk_assert



/****************************************************************
 *    tolower                                                   *
 *                                                              *
 * Convert string to lower case                                 *
 * input:                                                       *
 *  - __LINE__                                                  *
 * output:                                                      *
 *  - converted lowercase string                                *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void tolower(char *s)
  {
  TRACE_FUNCTION("tolower (BK 06-Sep-1999)")
  #ifdef WIN32
    s = _strlwr(s);// Jia
  #else
  while (*s != '\0')
    *s++ = tolower(*s);// cctype
  #endif
  }



/****************************************************************
 *    toupper                                                   *
 *                                                              *
 * Convert string to upper case                                 *
 * input:                                                       *
 *  - __LINE__                                                  *
 * output:                                                      *
 *  - converted uppercase string                                *
 *                                                              *
 *    Bert Kampes, 06-Sep-1999                                  *
 ****************************************************************/
void toupper(char *s)
  {
  TRACE_FUNCTION("toupper (BK 06-Sep-1999)")
  #ifdef WIN32
    s = _strupr(s);// Jia
  #else
  while (*s != '\0') {
    //*s++ = toupper(*s);  // [AV] do not works correctly with g++ 4.8 
    *s = std::toupper(*s);                 // cctype
    ++s;
  }
  #endif
  }


/****************************************************************
 *    int2str                                                   *
 *                                                              *
 * Convert integer to string                                    *
 * input:                                                       *
 *  - integer                                                   *
 * output:                                                      *
 *  - converted to string                                       *
 *                                                              *
 *    Mahmut Arikan, 23-Oct-2008                                *
 ****************************************************************/
string int2str(const int &integer)
        {
        TRACE_FUNCTION("int2str (MA 23-Oct-2008)")
        ostringstream datastream;
        datastream << integer << flush;
        return (datastream.str());
        } // End of int2str

