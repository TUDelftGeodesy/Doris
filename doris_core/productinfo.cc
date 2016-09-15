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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/productinfo.cc,v $        *
 * $Revision: 3.13 $                                            *
 * $Date: 2005/08/01 15:51:19 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of product info class.                        *
 * - data filling/updating.                                     * 
 * - reading in matrices.                                       *
 * - etc.                                                       *
 #%// BK 01-Sep-2000
 ****************************************************************/


#include "constants.hh"                 // DEBUG
#include "productinfo.hh"               // declarations, constants.h
#include "ioroutines.hh"                // printcpu
#include "exceptions.hh"                 // my exceptions class



/****************************************************************
 *    productinfo::fillproductinfo                              *
 * Fills data of productinfo object, read from resultfile after *
 * steps that produced a product, identified by 'iden'.         *
 * 13 lines are checked for info, first occurence counts.       *
 * exit if 8 are found, or 15 lines are checked.                *
 *  and read information upto string ":_NORMAL" (i.e., only read  *
 *  appropriate section                                         *
 *                                                              *
 *  "Data_output_file:"                                         *
 *  "Data_output_format:"                                       *
 *  "Multilookfactor_azimuth_direction:"                        *
 *  "Multilookfactor_range_direction:"                          *
 *  "First_line"  (w.r.t. original):                            *
 *  "Last_line"  (w.r.t. original):                             *
 *  "First_pixel"  (w.r.t. original):                           *
 *  "Last_pixel"  (w.r.t. original):                            *
 *                                                              *
 * input:                                                       *
 *  - resultfilename                                            *
 *  - identifier                                                *
 * output:                                                      *
 *  - (updated) productinfo                                     *
 *                                                              *
 *    Bert Kampes, 22-Dec-1998                                  *
 *  Dont know why found1 is used more then ones?                *
 *  BK 13-april-2000                                            *
 #%// BK 01-Sep-2000 class implementation                       *
 * added check for filesize                                     *
 #%// BK 08-Mar-2001                                            *
 * added stop after ":_NORMAL" found                            *
 #%// Bert Kampes, 01-Aug-2005                                  *
 #%// Mahmut Arikan, 02/2010 ONE27 --> 2*ONE27                  *
 ****************************************************************/
void productinfo::fillproductinfo(
        const char* resultfile,
        const char* iden)
  {
  TRACE_FUNCTION("fillproductinfo (BK 01-Sep-2000)")
  char                  word[4*ONE27]=" ";      // MA 4*ONE27
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
      break;
      }
    }

// ______Check if section has been found______
  if (!foundiden)
    {
    ERROR << "(fillproductinfo). identifier: \"" << iden
         << "\" not found in file: "            << resultfile;
    PRINT_ERROR(ERROR.get_str())
    throw(file_error);
    }

// ======Extract info (file name/format and window size)======
  int32 numlineschecked = 0;                    // number of lines checked after iden
  bool foundfilename    = false;
  bool foundfileformat  = false;
  bool foundfirstline   = false;
  bool foundlastline    = false;
  bool foundfirstpixel  = false;
  bool foundlastpixel   = false;
  bool foundmlL         = false;
  bool foundmlP         = false;

  while(resfile)                                 // i.e. rest of file
    {
    resfile.getline(dummyline,4*ONE27,'\n');       // next line
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
    if (foundfilename  && foundfileformat &&
        foundfirstline && foundlastline && foundfirstpixel && foundlastpixel &&
        foundmlL       && foundmlP) break; 
    //if (numlineschecked == 13) break;
    // report that 13 was too little for some output sections
    // so increased it to 15, BK 07-05-2002
    if (numlineschecked == 15) break;
    numlineschecked++;

    // ___ Check for strings with parameters ___
    resfile >> word;                             // read word
    if (!strcmp(word,"Data_output_file:"))
      {
      if (!foundfilename)
        {
        foundfilename = true;
        resfile >> file;
        DEBUG << "String: \"Data_output_file:\", \t\t\tvalue: " << file;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Data_output_file:\" found again (ignored).");
      }
    else if (!strcmp(word,"Data_output_format:"))
      {
      if (!foundfileformat)
        {
        foundfileformat = true;
        resfile >> word;
        if (!strcmp(word,"complex_short"))
          formatflag = FORMATCI2;
        else if (!strcmp(word,"complex_real4"))
          formatflag = FORMATCR4;
        else if (!strcmp(word,"real4"))
          formatflag = FORMATR4;
        else if (!strcmp(word,"hgt"))
          formatflag = FORMATHGT;
        else
          {
          PRINT_ERROR("wrong format specifier (impossible?)")
          throw(unhandled_case_error);
          }
        DEBUG << "String: \"Data_output_format:\", \t\tvalue: " << formatflag;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Data_output_format:\" found again (ignored).");
      }

    else if (!strcmp(word,"Multilookfactor_azimuth_direction:"))
      {
      if (!foundmlL)
        {
        foundmlL = true;
        resfile >> multilookL;
        DEBUG << "String: \"Multilookfactor_azimuth_direction:\", value: "
             << multilookL;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Multilookfactor_azimuth_direction:\" found again (ignored).");
      }
    else if (!strcmp(word,"Multilookfactor_range_direction:"))
      {
      if (!foundmlP)
        {
        foundmlP = true;
        resfile >> multilookP;
        DEBUG << "String: \"Multilookfactor_range_direction:\", \tvalue: "
             << multilookP;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Multilookfactor_range_direction:\" found again (ignored).");
      }

    else if (!strcmp(word,"First_line"))                // (w.r.t. original):
      {
      if (!foundfirstline)
        {
        foundfirstline = true;
        resfile >> word >> word >> win.linelo ;
        DEBUG << "String: \"First_line:\", \t\t\tvalue: " << win.linelo;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"First_line:\" found again (ignored).");
      }

    else if (!strcmp(word,"Last_line"))                 // (w.r.t. original):
      {
      if  (!foundlastline)
        {
        foundlastline = true;
        resfile >> word >> word >> win.linehi ;
        DEBUG << "String: \"Last_line:\", \t\t\tvalue: " << win.linehi;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Last_line:\" found again (ignored).");
      }

    else if (!strcmp(word,"First_pixel"))               // (w.r.t. original):
      {
      if (!foundfirstpixel)
        {
        foundfirstpixel = true;
        resfile >> word >> word >> win.pixlo ;
        DEBUG << "String: \"First_pixel:\", \t\t\tvalue: " << win.pixlo;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"First_pixel:\" found again (ignored).");
      }

    else if (!strcmp(word,"Last_pixel"))                // (w.r.t. original):
      {
      if (!foundlastpixel)
        {
        foundlastpixel = true;
        resfile >> word >> word >> win.pixhi ;
        DEBUG << "String: \"Last_pixel:\", \t\t\tvalue: " << win.pixhi;
        DEBUG.print();
        }
      else
        WARNING.print("String: \"Last_pixel:\" found again (ignored).");
      }
    } // while resultfile


  // ______ Check filesize with format/dimensions (BK 08-Mar-2001) ______
  // BK 08-Mar-2001
  //ifstream tmpfile(file, ios::in | ios::nocreate);
  ifstream tmpfile(file, ios::in);
  if (tmpfile)
    {
    tmpfile.seekg(0,ios::end); // internal filesize, normal one exists if not exists
    // uint filesizetrue  = tmpfile.tellg();
    const streamoff &filesizetrue  = tmpfile.tellg();   // [MA] file > 4GB support, this fix eliminates wrong warning
    int16 bytesperelem = 4;
    if (formatflag==FORMATCI2) bytesperelem=4;
    if (formatflag==FORMATCR4) bytesperelem=8;
    if (formatflag==FORMATR4) bytesperelem=4;
    if (formatflag==FORMATI2) bytesperelem=2;
    if (formatflag==FORMATI2_BIGENDIAN) bytesperelem=2;
    if (formatflag==FORMATR8) bytesperelem=8;
    if (formatflag==FORMATHGT) bytesperelem=8;
    // int32 filesizecomp = int32(win.lines()/multilookL) *
    uint64 filesizecomp = (uint64)(win.lines()/multilookL) *
                                 (win.pixels()/multilookP) *
                                               bytesperelem;
    DEBUG << "Checking format/dimensions file=" << file;
    DEBUG.print();
    if (filesizecomp != filesizetrue)
      {
      WARNING << "File: \'" << file << "\' has wrong fileformat or dimensions"   // [MA] reorganized + add file
           << ": bytesperpix="    << bytesperelem
           << ", #l="             << win.lines()/multilookL      // [MA]
           << ", #p="             << win.pixels()/multilookP
           << "; size_computed=" << filesizecomp << "B" << " v."
           << " size_ondisk="    << filesizetrue << "B";
      WARNING.print();
      }
    else
      {
      DEBUG.print("Fileformat and dimensions are checked and ok.");
      }
    tmpfile.close();
    } // stream ok
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
  DEBUG.print("finished fillproductinfo");
  DEBUG.print("content of struct:");
  showdata();
#endif
  } // END fillproductinfo


/****************************************************************
 * productinfo::readphase                                       *
 *  read data from file in a real4 matrix                       *
 *  phase is extracted from hgt, not unwrapped set to NaN       *
 *  real4 is read 'as is'.                                      *
 *  complex real4: angle is taken.                              *
 *  hgt: phase is read, set to NaN if amplitude==0              *
 * cutout window in own system starting at line 1,              *
 * (no offset win)                                              *
 * This function is useful for reading products: complex interf.*
 * unwrapped interf. and to get phase for s2h etc.              *
 #%// BK 22-Sep-2000                                            *
 ****************************************************************/
matrix<real4> productinfo::readphase(
        window cutout) const               // window to be read
  {
  TRACE_FUNCTION("productinfo::readphase (BK 22-Sep-2000)")
  // ______ Check if cutout < #lines/pixels ______
  const int32 linesondisk  = int32(win.lines()/multilookL);
  const int32 pixelsondisk = int32(win.pixels()/multilookP);

  #ifdef __DEBUGMAT2
    if (cutout.linelo < 1)
      {
      PRINT_ERROR("readphase: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.linehi>linesondisk)
      {
      PRINT_ERROR("readphase: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixlo  < 1)
      {
      PRINT_ERROR("readphase: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixhi>pixelsondisk)
      {
      PRINT_ERROR("readphase: cutout window larger than what's on disk!")
      throw(input_error);
      }
  #endif

  // ______ Open file ______
//#ifdef __NO_IOS_BINARY__
//  ifstream ifile(file, ios::in);
//#else
//  ifstream ifile(file, ios::in | ios::binary);
//#endif
  ifstream ifile;
  openfstream(ifile,file);
  bk_assert(ifile,file,__FILE__,__LINE__);

  matrix<real4> Result(cutout.lines(),cutout.pixels());

  // ====== Actual read data ======
  switch (formatflag)
    {
    case FORMATR4:
      {
      // why not use readfile function?
      DEBUG.print("reading real4 from file");
      matrix<real4> LINE(1,cutout.pixels());            // phase
      for (int32 line=cutout.linelo; line<=cutout.linehi; ++line)
        {
        uint64 start = (uint64)(cutout.pixlo-1 + pixelsondisk*(line-1))
                      * sizeof(real4);
        ifile.seekg(start,ios::beg);
        ifile >> LINE;
        Result.setrow(line-cutout.linelo,LINE);         // not the fastest way...
                                                        // but cannot write directly
                                                        // to private data ... (bk)
        }
      break;
      }

    case FORMATCR4:
      {
      DEBUG.print("reading complex real4 from file (get phase)");
      const window dummyoffset(1,99999,1,99999);        // only 1's used -> no offset
      matrix<complr4> TMP(cutout.lines(),cutout.pixels());
      readfile(TMP,file,linesondisk,cutout,dummyoffset);
      Result = angle(TMP);
      break;
      }

    // ______ hgt: first amplitude band, then phase, then next line ______
    case FORMATHGT:
      {
      DEBUG.print("reading hgt (band interleaved) from file (get phase)");
      // ______ Read in one line, set to NaN, fill result ______
      matrix<real4> LINE(1,win.pixels()*2);     // disk phase/ampl band interleaved
      for (int32 line=cutout.linelo; line<=cutout.linehi; ++line)
        {
        const int32 start = (line-1) * pixelsondisk * 2 * sizeof(real4);
        ifile.seekg(start,ios::beg);            // file pointer to 1st pix of line
        ifile >> LINE;
        for (int32 pix=cutout.pixlo; pix<=cutout.pixhi; ++pix)
          {
          // ______ Check amplitude defined as 0 if not unwrapped ok. ______
          if (LINE(0,pix-1) == 0.)              // ampl. band
            Result(line-cutout.linelo,pix-cutout.pixlo) = NaN;
          else                                                  // phase band
            Result(line-cutout.linelo,pix-cutout.pixlo) = LINE(0,pix-1+pixelsondisk);
          }
        }
      break;
      }

    default:
      PRINT_ERROR("readphase::not correct format on file.")
      throw(unhandled_case_error);
    }

  ifile.close();
  return Result;
  } // END readphase


/****************************************************************
 * productinfo::readdata                                        *
 *  read data from file in a complex real4 matrix               *
 *  complex real4                                               *
 * cutout window in own system starting at line 1,              *
 * (no offset win)                                              *
 * This function is useful for reading products: complex interf.*
 * unwrapped interf. and to get phase for s2h etc.              *
 #%// BK 22-Sep-2000 (readphase)                                *
 #%// Batu 30-JUL-2007                                          *
 ****************************************************************/
matrix<complr4> productinfo::readdata(
        window cutout) const               // window to be read
  {
  TRACE_FUNCTION("productinfo::readdata (Batu 30-JUL-2007)")
  // ______ Check if cutout < #lines/pixels ______
  const int32 linesondisk  = int32(win.lines()/multilookL);
  const int32 pixelsondisk = int32(win.pixels()/multilookP);

  #ifdef __DEBUGMAT2
    if (cutout.linelo < 1)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.linehi>linesondisk)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixlo  < 1)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixhi>pixelsondisk)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
  #endif

  // ______ Open file ______
//#ifdef __NO_IOS_BINARY__
//  ifstream ifile(file, ios::in);
//#else
//  ifstream ifile(file, ios::in | ios::binary);
//#endif
  ifstream ifile;
  openfstream(ifile,file);
  bk_assert(ifile,file,__FILE__,__LINE__);

  matrix<complr4> Result(cutout.lines(),cutout.pixels());

  // ====== Actual read data ======
      DEBUG.print("reading complex real4 from file (readdata)");
      const window dummyoffset(1,99999,1,99999);        // only 1's used -> no offset
      readfile(Result,file,linesondisk,cutout,dummyoffset);

  ifile.close();
  return Result;
  } // END readdata


/****************************************************************
 * productinfo::readdatar4                                      *
 *  read data from file in a complex real4 matrix               *
 *  complex real4                                               *
 * cutout window in own system starting at line 1,              *
 * (no offset win)                                              *
 * This function is useful for reading products: complex interf.*
 * unwrapped interf. and to get phase for s2h etc.              *
 #%// BK 22-Sep-2000 (readphase)                                *
 #%// Batu 30-JUL-2007                                          *
 #%// MA   19-NOV-2008 [TODO] template                          *
 ****************************************************************/

matrix<real4> productinfo::readdatar4(
        window cutout) const               // window to be read
  {
  TRACE_FUNCTION("productinfo::readdatar4 (MA & BO 19-NOV-2008)")
  // --- Log debug info ---
  DEBUG << "Reading file:     " << file;
  DEBUG.print();
  DEBUG << "Formatflag:       " << formatflag;
  DEBUG.print();
  DEBUG << "Currentwindow:    ";  win.disp();     // appends to DEBUG and prints
                                                  // [MA] tricked since productinfo lacks currentwindow()
  DEBUG << "Window from file: ";  cutout.disp();  // appends to DEBUG and prints
  //
  
  // ______ Check if cutout < #lines/pixels ______
  const int32 linesondisk  = int32(win.lines()/multilookL);
  const int32 pixelsondisk = int32(win.pixels()/multilookP);

  #ifdef __DEBUGMAT2
    if (cutout.linelo < 1)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.linehi>linesondisk)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixlo  < 1)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
    if (cutout.pixhi>pixelsondisk)
      {
      PRINT_ERROR("readdata: cutout window larger than what's on disk!")
      throw(input_error);
      }
  #endif

  ifstream ifile;
  openfstream(ifile,file);
  bk_assert(ifile,file,__FILE__,__LINE__);
  ifile.close();

  //matrix<real4> Result(cutout.lines(),cutout.pixels());

  // ====== Actual read data ======
  switch (formatflag)
    {
    case FORMATR4:                  // [MA] bigendian support
      {
      matrix<real4> Result(cutout.lines(),cutout.pixels());
      DEBUG.print("reading real4 from file (readdata)");
      //const window dummyoffset(1,99999,1,99999);      // only 1's used -> no offset
      //readfile(Result,file,linesondisk,cutout,dummyoffset);
      readfile(Result,file,linesondisk,cutout,win);
      return Result;
      break;
      }
//   case FORMATR8:
//     {
//     matrix<real8> Result(cutout.lines(),cutout.pixels());
//     DEBUG.print("reading real8 from file (readdata)");
//     readfile(Result,file,linesondisk,cutout,win);
//     return Result;
//     break;
//     }
//   case FORMATCR4:
//     {
//     matrix<complr4> Result(cutout.lines(),cutout.pixels());
//     DEBUG.print("reading complex real4 from file (readdata)");
//     readfile(Result,file,linesondisk,cutout,win);
//     return Result;
//     break;
//     }
    default:
      PRINT_ERROR("readdata::not correct format on file.")
      throw(file_error);
   }
  //return Result;
  } // END readdatar4
