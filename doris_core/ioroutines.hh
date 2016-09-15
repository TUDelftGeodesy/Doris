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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/ioroutines.hh,v $
 * $Revision: 3.9 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Definitions of general input/output routines.
 ****************************************************************/
/****************************************************************
 * Note that the strstream class using a buffer and stream output memory
 * can be replaced by the following construction (BK, april 2003):
 * -----------------------------------------------------
 * static ostringstream     message;
 * message << "Hi Bert. " << 1.1 << ends;
 * cout << "message.str(): " << message.str() << "\n";
 * message.str(string());// empty buffer
 * -----------------------------------------------------
 ****************************************************************/


#ifndef IOROUTINES_H
#define IOROUTINES_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                 // global ONE27, inputstructs etc.
#include "readinput.hh"                 // input structs
#include "matrixbk.hh"                  // my matrix class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class

#include <cstdlib>                      // exit



//  // ====== Prototypes and inlined functions ======
// ______ Screenoutput ______
void printcpu(bool init=false);         // prototype


// ====== Prototypes ======
// ______ Initialization ______
void inittest();                        // Some initial tests


// ______ Find out if initwrite is necessary ______
bool doinitwrite(
        input_gen &generalinput,
        int16 imageid);


// ______ Write initial data to files ______
// ______ id=1: logfile, 2:master result, 3: slaveresult ______
void initwrite(
        const char* file,
        int16 fileid);


// ______ Append filein to fileout and remove in ______
void updatefile(
        const char* filein,
        const char* fileout);


// ______ Search file for result and return that ______
bool readres(
        char* returnword,
        const int16 sizeofreturnword, 
        const char* file, 
        const char* patern, 
        const int16 napatern=0,
              int16 skiplines=0);  


// ______ Read coefficients from resultfile ______
matrix<real8> readcoeff(
        const char* file,
        const char* pattern,
        const int16 Ncoefficients);

// ______ Read normalization coefficients from ifgs result file________
matrix<real8> readnormcoeff(
        const char* file,
        const char* pattern);

// ______ Updates process_control in resultfiles ______
void updateprocesscontrol(
        const char* filein,
        int16 fileid);


// ______ Check if only ones results per step ______
void checkprocessing(
        const input_gen &generalinput,
        int16 alreadyprocessed[NUMPROCESSES]);


// ______ Check if requested processing is logical ______
void checkrequest(
        int16 step,
        int16 alreadyprocess[NUMPROCESSES], ...);


// ______ Check if process control flags are written correctly _____
void fillcheckprocess(
        const char *line,
        int16 checkprocess[NUMPROCESSES],
        int16 fileid);


// ______ Fill array what is already processed ______
void fillprocessed(
        const char *line,
        int16 checkprocess[NUMPROCESSES],
        int16 fileid);


// ______ Removes data section from resultfile ______
void removedatleader(
        const char *file);


// ______ True if file existed ______
bool existed(
        const char *file);


// ______ Returns filesize ______
uint filesize(
        const char *file);


// ______ Count number of eols ______
int32 filelines(
        const char *file);


// ______ Pause for interactive processing ______
void getanswer(
        );


// ______ Open an input file stream ______
void openfstream(
        ifstream &str,
        const char* ifilename);


// ______ Open an output file stream ______
void openfstream(
        ofstream &str,
        const char* ofilename,
        const bool overwrit);


// ______ Check if file is opened correctly ______
void bk_assert(
        const ifstream &str,
        const char* ifilename,
        const char* callingfilename = "?",
        int32 linenumber = 0);


// ______ Check if file is opened correctly ______
void bk_assert(
        const ofstream &str,
        const char* ofilename,
        const char* callingfilename = "?",
        int32 linenumber = 0);


// _____ Convert a string to lower/upper case ______
void tolower(char *s);
void toupper(char *s);

// _____ Convert an integer to string ______
string int2str(const int &integer); // [MA]




//  // ______ (partially) Fills struct image/inter-info ______
//  void fillslcimage(
//      slcimage &image,
//      const char *file);


//  // ______ Struct slcimage: update it after certain step ______
//  void updateslcimage(
//      slcimage &image,
//      const char *file,
//      const char *iden);


// ______ Fill struct interferogram info ______
//void fillproductinfo(
//      productinfo         &interferogram, 

// ______ compute overlap in same system ______
window getoverlap(
        const window        &master,
        const window        &slave);

// compute overlap approx with transf. polynomial between slave2master coord.system
window getoverlap(
        const slcimage      &master,
        const slcimage      &slave,
        const matrix<real8> &cpmL,
        const matrix<real8> &cpmP);

window getoverlap(
        const slcimage      &master,
        const slcimage      &slave,
        const real8         &Npointsd2,
        const real8         &timing_L,
        const real8         &timing_P); //[FvL]

void lineintersect(
        const real8   &ax,
        const real8   &ay,
        const real8   &bx,
        const real8   &by,
        const real8   &cx,
        const real8   &cy,
        const real8   &dx,
        const real8   &dy,
        matrix<real8> &exy); //[FvL]

//  // ______ Struct productinfo: fill info after certain step ______
//  void fillproductinfo(
//      productinfo &interferogram,
//      const char  *file,
//      const char  *iden);




#endif // IOROUTINES_H



