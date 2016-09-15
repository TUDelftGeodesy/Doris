/*
 @file   excpetions.hh exception handling for Doris InSAR processor
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

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

using namespace std;                    // BK 29-Mar-2003, new compiler?
                                        // TODO SLiu see constants.hh

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <csignal>                      // function signal()
#include <cstring>                      // strcpy



/*********************************************************************
 * @brief exception handler for floating point exception
 *********************************************************************/
void CatchSignals(void (*SigHandler)(int));
void handle_signal(int signum);


/*********************************************************************
 * @brief exception handler class
 *********************************************************************/
// #include "exceptions.hh"  // this file..  
// int main(
//        int argc,
//        char* argv[])
// {
// // ___ catch math errors ___
// CatchSignals(handle_signal);
//
// try // --- start trying -----------------
//   {
//   code...
//   if (error) throw(some_error);// object
//   if (error) throw("sasas");// string
//   }
// }// --- end of try block; now catch thrown exceptions -------------
// catch(EXCEPTION& error)// catch errors of EXCEPTION class
//   {
//   cerr << "i caught an error!" << endl;
//   cerr << "i think: " << (const char*)error << endl;
//   exit(1);
//   }
// catch(const char* error_string)// catch handled errors
//   {
//   cerr << "i caught an error_string!" << endl;
//   cerr << "it is: " << error_string << endl;
//   exit(2);
//   }
// catch(...) // catches other errors
//   {
//   cerr << "i caught an unhandled error!" << endl;
//   exit(3);
//   }
//
// cout << "\n\nNormal termination.\nThank you for using Doris.\n\n";
// return int(0);
// } // END main


// ______ Base class for all error exceptions that can be caught ______
class EXCEPTION
  {
  public:
    EXCEPTION()          {};
    virtual ~EXCEPTION() {};
    operator const char*() const {return(get_error_string());};
    // virtual operator const char*() const {return(get_error_string());}; // suggested by SLiu
    virtual const char* get_error_string() const {return("generic error");};
  };
// ______ Now all errors follow ______
// ______ some error ______
class SOME_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];
  public:
    SOME_ERROR() {strcpy(err_str,"specific error");};
    virtual ~SOME_ERROR() {};
    operator const char*() const {return(get_error_string());};   // overloading, why?
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ some input error ______
class INPUT_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    INPUT_ERROR() {strcpy(err_str,"input error");};
    virtual ~INPUT_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ some file error ______
class FILE_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    FILE_ERROR() {strcpy(err_str,"file error");};
    virtual ~FILE_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ memory error ______
class MEMORY_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    MEMORY_ERROR() {strcpy(err_str,"memory error");};
    virtual ~MEMORY_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ unhandled case error ______
class UNHANDLED_CASE_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    UNHANDLED_CASE_ERROR() {strcpy(err_str,"unhandled case error");};
    virtual ~UNHANDLED_CASE_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ unhandled case error ______
class ARGUMENT_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    ARGUMENT_ERROR() {strcpy(err_str,"wrong input argument(s) to function");};
    virtual ~ARGUMENT_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ keyword error ______
class KEYWORD_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    KEYWORD_ERROR() {strcpy(err_str,"incorrect keyword");};
    virtual ~KEYWORD_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };
// ______ usage request error ______
class USAGE_ERROR : public EXCEPTION
  {
  private:
    char err_str[64];// make part of base class?
  public:
    USAGE_ERROR() {strcpy(err_str,"done");};
    virtual ~USAGE_ERROR() {};
    operator const char*() const {return(get_error_string());};
    virtual const char* get_error_string() const {return(err_str);};
  };


// ====== Globals to throw everywhere, e.g., throw(some_error) ======
extern SOME_ERROR       some_error;// can be thrown from all programs
extern INPUT_ERROR      input_error;// can be thrown from all programs
extern FILE_ERROR       file_error;// can be thrown from all programs
extern MEMORY_ERROR     memory_error;// can be thrown from all programs
extern UNHANDLED_CASE_ERROR unhandled_case_error;// can be thrown from all programs
extern ARGUMENT_ERROR   argument_error;// can be thrown from all programs
extern KEYWORD_ERROR    keyword_error;// can be thrown from all programs
extern USAGE_ERROR      usage_error;// can be thrown from all programs


#endif // EXCEPTIONS_H


