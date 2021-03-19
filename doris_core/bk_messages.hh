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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/bk_messages.hh,v $ *
 * $Revision: 3.12 $ *
 * $Date: 2005/08/24 10:03:18 $ *
 * $Author: kampes $ *
 *
 * definition of bk_message class *
 #%// BK 10-Apr-2003
 ****************************************************************/

#ifndef BK_MESSAGES_H
#define BK_MESSAGES_H
using namespace std;

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>    // cout
#include <iomanip>     // setw()
#include <strstream>   // mem.buf
#include <cstring>     // strcpy()
#include <cstdlib>     // exit()

// #if majorversion>3 and g++
//#include <stringstream>   // new mem.buf
//#include <string>



// ---------------------------------------------------------------- 
// --- Choice was either to make something like:
// --- message ERROR;
// --- message DEBUG;
// --- message INFO;
// --- with functions like: ERROR.doprint(1); ERROR.setidentifyer("ERROR"); 
// ---                      ERROR << "dsaada" << 4.5 << setw(4) << nw;
// ---                      ERROR.print();
// --- with separate bufs, or to combine them in a single object
// --- like:
// --- messages logging;
// --- logging.setlevel("INFO");
// --- logging << "dsaada" << 4.5 << setw(4) << nw;
// --- logging.info();
// --- The first is nice, different separated bufs, and easy to change 
// --- existing code.  However,
// --- The latter is better with Doris, since it allows things like
// --- if (good) then logging.info() else logging.warning();
// --- But by using construction like 
// --- WARNING.print(INFO.get_str()); this is solved
// --- #%// BK 10-Apr-2003
//
// most important member functions:
//   get_str():       leave buffer unaffected, pointer to a copy
//   print("sss"):    leave buffer unaffected, print a copy
//   print():         leave buffer unaffected, terminate buffer, rewind buffer;
// this means that you have to rrewind the buffer yourself if you misuse it
// as storage, like TRACE.rewind();
// INFO.precision(10); /* permanent */  INFO << setp(13)<<q;// once
// INFO.width(11);
// INFO.reset();
// ---------------------------------------------------------------- 
// EXAMPLE USEAGE:
// ---------------------------------------------------------------- 
// in code simply change with a define calls to ERROR to give a WARNING summary, etc.
// #define printERROR cout<<"this is defined...\n"; warning.summary(); error<<" ("<<__FILE__<<": "<<__LINE__<<")"; error.print();
//
// ---------------------------------------------------------------- 
// in globals, put them as extern to access them everywhere
// #include "globals.h"
// --- declaration of global info message object ---
// bk_messages error;
//bk_messages INFO;
//bk_messages DEBUG;
//int main( int argc, char* argv[])
//  {
//  DEBUG.setidentifyer("DEBUG");
//  INFO.setidentifyer("INFO");
//  INFO.doprint(1);
//  INFO.dostderr(0);
//  INFO.bellrings(2);
//  DEBUG << 421;
//  DEBUG << " dsadsaa: " << 433 ;
//  DEBUG << 421.1;
//  DEBUG << setw(20) << double(421.1) << ";;;";
//  DEBUG.print();
//  INFO << "this is a info message";
//  INFO.print();
//  INFO.print("this is a another info message");
//  INFO.summary();
//  return 0;
//  }



// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
class bk_messages
  {
  // ______ Private data of bk_message class ______
  private:
    bool print_cout;//    switch to print or not to cout
    bool print_cerr;//    switch to print or not to cerr
    bool do_exit_;//      switch to exit after print or not
    int  ringbell_;//     number of times to ring bell.
    int  strwidth;//      original str.width();
    int  strprecision;//  original str.precision();
    char name_[10];//     identifyer, <9 chars (PROGRESS)
    ostrstream buf;//     message buffer stream
    long long nummessages;//    counter
    char first_messages[6][256];// summary of first 6 warnings displayed at end or at exit;

  // ______ Public function of bk_message class ______
  public:
    // ___Constructor/Destructor___
    bk_messages()
      {
      //cerr << "YOU ARE IN THE CONSTRUCTOR\n";
      print_cout  = true;//  default do print;
      print_cerr  = false;// default no print to cerr;
      do_exit_    = false;// default no exit after print
      ringbell_   = 0;//     default no bell;
      strwidth    = buf.width();
      strprecision= buf.precision();
      strcpy(name_,"M??????");// default identifyer
      nummessages = 0;
      // --- allocate enough space, then freeze it?, and keep this
      // --- stream as is ---
      // --- this is required, at g++, and can't really hurt ---
      buf.seekp(0);
      buf << ends
          << "This is the constructor...This is the constructor..."
          << "This is the constructor...This is the constructor..."
          << "This is the constructor...This is the constructor..."
          << "This is the constructor...This is the constructor..."
          << "This is the constructor...This is the constructor..."
          << ends;
      buf.seekp(0);// reset pointer
      //      //buf.rdbuf()->freeze(0); // Unfreeze
      //      //buf.rdbuf()->freeze(); // freeze string, keep allocated mem
      //      buf.freeze();// freeze string, keep allocated memory ???
      }

    // ___Constructor/Destructor___
    ~bk_messages()                 //{;}// nothing to destruct
      {
      delete []buf.str();
      }// dealloc

    // ___Helper functions___
    void terminate()            {buf<<ends;}
    void rewind()               {buf.seekp(0);}
    char* get_id()              {return name_;}
    void dostderr(bool s)       {print_cerr = s;}
    void doprint(bool s)        {print_cout = s;}
    void doexit(bool s)         {do_exit_   = s;}
    void bellrings(int s)       {ringbell_  = s;}
    int  width()          const {return int(buf.width());}
    void width(int s)           {buf.width(s);}
    int  precision()      const {return int(buf.precision());}
    void precision(int s)       {buf.precision(s);}
    //void setiosflags(long s)    {buf.setiosflags(s);}
    //void resetiosflags(long s)  {buf.resetiosflags(s);}
    //void unsetf(long s)         {buf.unsetf(s);} // is this standard?
    void reset()                {// what does resetiosflags do?
                                 // add all flags?
                                 buf << resetiosflags(
                                        ios::fixed       |
                                        ios::showpoint   | 
                                        ios::skipws      |
                                        ios::showpos     |
                                        ios::scientific  |
                                        ios::left        |
                                        ios::right       |
                                        ios::floatfield  |
                                        ios::adjustfield);
                                 width(strwidth);
                                 precision(strprecision);
                                 terminate();
                                 rewind();
                                }

    // ___Set a name to prepend to message string___
    void setidentifyer(const char *id)
      {
      // ___ only use first 9 of ID string ___
      if (strlen(id)>=10)
        {
        strncpy(name_,id, 9);
        name_[9] = '\0';// terminate id
        }
      else
        {
        strcpy(name_,id);// identifyer
        }
      // ___ default message behavior for some identifyers ___
      if (!strcmp(id,"ERROR"))
        {
        do_exit_   = true;// default message behavior
        print_cerr = true;
        ringbell_  = 3;
        }
      else if (!strcmp(id,"WARNING"))
        {
        ringbell_ = 2;
        print_cerr = true;
        }
      else if (!strcmp(id,"PROGRESS"))
        {
        ringbell_ = 1;
        print_cerr = true;
        }
      }

    // ___ problem that buf.str() empties string buffer, thus do cp ___
    // ___ append a null, since common use will be to get a copy after written fully ___
    // ___ this returns a pointer to buf, thus, that may be changed afterwards ___
    // ___ use immediately? ___
    //char* get_str() {return buf.str();}
    char* get_str()
      {
      char *s = buf.str();
      buf.freeze(0);// unfreeze, since str() freezes it (does not seem to work?) 
      return s;
      }

    // ___What it is all about, print the message___
    // ___the buffer pointer is set to (0), but the message is not deleted___
    // ___so people can access the buffer multiple times___
    // ___but you cannot add to it later...____
    void print()
      {
      // --- Allways rewind string to prevent unlimited growth ---
      // --- Only terminate it when new info is there ------------
      //cerr << "pcount: " << buf.pcount() << endl;
      if (buf.pcount()>0) {terminate(); rewind();}
      print(buf.str());
      buf.freeze();
      }

    // ___What it is all about, print the message___
    // ___the buffer pointer is set to (0), but the message is not deleted___
    // ___so people can access the buffer multiple times___
    // ___but you cannot add to it later...____
    void print(const char *s)
      {
      if (nummessages<6) strcpy(first_messages[nummessages],s);
      nummessages++;
      // --- only print if level has been set ----------------------
      if(print_cout==true)
        {
        if (do_exit_==true) 
          cout << endl << "!!! ABNORMAL TERMINATION !!!" << endl;
        cout << setw(8) << setiosflags(ios::left) << name_ << ": " << s << endl;
        #ifndef WIN32
        if (print_cerr==true)
          cerr << setw(8) << setiosflags(ios::left) << name_ << ": " << s << endl;
        #endif
        // --- Ring bell -----------------------------------------
        for (int i=0;i<ringbell_;++i) cerr << "\a";
        // --- Exit ----------------------------------------------
        if (do_exit_==true) 
          {
          /* with a define you can print other stuff at exit...
           * #define printERROR cout<<"this is defined...\n"; \
           * warning.summary(); error<<" ("<<__FILE__<<": "
           *  <<__LINE__<<")"; error.print();
           */
          cerr << endl << "!!! ABNORMAL TERMINATION !!!" << endl;
          exit(9999);
          // Bert Kampes, 31-Mar-2005: throw relevant exception
          //cerr << "throwing string exception" << endl;
          //throw("!!! ABNORMAL TERMINATION !!!");// handle exception as string
          }
        }
      }


    // --- SUMMARY OF FIRST 6 MESSAGES ----------------------------
    // --- on stderr and stdout -----------------------------------
    void summary()
      {
      for (int i=0;i<ringbell_;++i) cerr << "\a";
      cout << endl << " --- " << setiosflags(ios::left) << name_ << " SUMMARY ---" << endl;
      #ifndef WIN32
      //if (print_cerr==true)
      cerr << endl << " --- " << setiosflags(ios::left) << name_ << " SUMMARY ---" << endl;
      #endif
      switch (nummessages)
        {
        case 0:
          cout << "There were no messages." << endl;
          cerr << "There were no messages." << endl;
        break;
        case 1:
          cout << "There was 1 message:" << endl
               << " 1: " << first_messages[0] << endl;
          cerr << "There was 1 message:" << endl
               << " 1: " << first_messages[0] << endl;
        break;
        default:
          cout << "There were " << nummessages << " messages:" << endl;
          cerr << "There were " << nummessages << " messages:" << endl;
          for (int i=0; i<nummessages; ++i)
            {
            if (i==6) 
              {
              cout << " " << i+1 << ": [...]" << endl;
              cerr << " " << i+1 << ": [...]" << endl;
              break; // only first six
              }
            cout << " " << i+1 << ": " << first_messages[i] << endl;
            cerr << " " << i+1 << ": " << first_messages[i] << endl;
            }
        }// switch nummessages
     }

    // ___Operator to build string from different types___
    // ostream& operator << (const int    X) {return buf << X;} 
    template <class Type> 
      ostream& operator << (const Type X) {return buf << X;}

    // ___ e.g., cout  << info; ___
    // ___ e.g., ofile << info; ___
    // friend ostream& operator << (ostream& s)
    // ostream& operator << (ostream& s, const bk_message &X)
    //    {
    //      if(print_cout==true)
    //        {
    //       X.terminate();
    //  s << X.getid() << ": " << X.getstr() << endl;
    //  X.rewind();
    //        buf << ends;
    //        s << name_ << ": " << buf.str() << endl;
    //        buf.seekp(0);
    //        }
    //      return s;
    //      }

  };// eof class


#endif // BK_MESSAGES_H

