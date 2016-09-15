/*
 @file   excpetions.cc exception handling for Doris InSAR processor
 @brief  exception handling for Doris InSAR processor
*/
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

#include <csignal>                      // function signal()
#include <iostream>                     // cerr
#include <cstdlib>                      // exit [MA]
#include "exceptions.hh"                // class



// ====== Globals to throw everywhere, e.g., throw(some_error) ======
SOME_ERROR           some_error;// can be thrown from all programs
INPUT_ERROR          input_error;// can be thrown from all programs
FILE_ERROR           file_error;// can be thrown from all programs
MEMORY_ERROR         memory_error;// can be thrown from all programs
UNHANDLED_CASE_ERROR unhandled_case_error;// can be thrown from all programs
ARGUMENT_ERROR       argument_error;// can be thrown from all programs
KEYWORD_ERROR        keyword_error;// can be thrown from all programs
USAGE_ERROR          usage_error;// can be thrown from all programs
// TODO: SLiu above not required see end of exceptions.hh



/*********************************************************************
 * @brief exception handler for floating point exception
 *********************************************************************/
/* function: CatchSignals()
 * ------------------------
 * Traps common signals that by default cause the program to abort.  
 * Sets (pointer to function) Handler as the signal handler for all.
 * Note that SIGKILL usually cannot be caught.  No return value.
 */   // following is from snaphu:
void CatchSignals(void (*SigHandler)(int)){
#ifndef WIN32
  signal(SIGHUP,SigHandler);// Hang-up signal
  signal(SIGQUIT,SigHandler);
  signal(SIGPIPE,SigHandler);
  signal(SIGALRM,SigHandler);
  signal(SIGBUS,SigHandler);
#endif
  signal(SIGINT,SigHandler);
  signal(SIGILL,SigHandler);
  signal(SIGABRT,SigHandler);
  signal(SIGFPE,SigHandler);// floating point exception
  signal(SIGSEGV,SigHandler);// segmentation fault: introduces when compiled with -O in gcc4?
  signal(SIGTERM,SigHandler);
}



// following is based on snaphu code: but needs some work.
/* function: SetDump()
 * -------------------
 * Set the global variable dumpresults_global to TRUE if SIGINT or SIGHUP
 * signals recieved.  Also sets requestedstop_global if SIGINT signal 
 * received.  This function should only be called via signal() when 
 * a signal is caught.
 */
void handle_signal(int signum)
  {
  switch (signum)
    {
#ifndef WIN32
    case SIGHUP:
      cout << "Caught SIGHUP: Hang-up signal." << endl;
      cerr << "Caught SIGHUP: Hang-up signal." << endl;
      break;
    case SIGQUIT:
      cout << "Caught SIGQUIT: Quit signal." << endl;
      cerr << "Caught SIGQUIT: Quit signal." << endl;
      break;
    case SIGPIPE:
      cout << "Caught SIGPIPE: ? signal." << endl;
      cerr << "Caught SIGPIPE: ? signal." << endl;
      break;
    case SIGALRM:
      cout << "Caught SIGALRM: Alarm signal." << endl;
      cerr << "Caught SIGALRM: Alarm signal." << endl;
      break;
    case SIGBUS:
      cout << "Caught SIGBUS: Bus error (accessing memory incorrectly)?" << endl;
      cerr << "Caught SIGBUS: Bus error (accessing memory incorrectly)?" << endl;
      break;
#endif
    case SIGINT:
      cout << "Caught SIGINT: User interupt signal." << endl;
      cerr << "Caught SIGINT: User interupt signal." << endl;
      exit(1);
      break;
    case SIGFPE:
      cout << "Caught SIGFPE: floating point exception, zero division, etc." << endl;
      cerr << "Caught SIGFPE: floating point exception, zero division, etc." << endl;
      break;
    case SIGILL:
      cout << "Caught SIGILL: ? signal." << endl;
      cerr << "Caught SIGILL: ? signal." << endl;
      break;
    case SIGABRT:
      cout << "Caught SIGABRT: Abort signal." << endl;
      cerr << "Caught SIGABRT: Abort signal." << endl;
      break;
    case SIGSEGV:
      cout << "Caught SIGSEGV: Segmentation fault." << endl;
      cerr << "Caught SIGSEGV: Segmentation fault." << endl;
      exit(1);
      break;
    case SIGTERM:
      cout << "Caught SIGTERM: ? signal." << endl;
      cerr << "Caught SIGTERM: ? signal." << endl;
      break;
    default:
      cout << "Caught an unknown signal.  Signum = " << signum << endl;
      cerr << "Caught an unknown signal.  Signum = " << signum << endl;
    }
  }




