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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readdata.cc,v $
 * $Revision: 3.37 $
 * $Date: 2005/10/06 11:09:20 $
 * $Author: kampes $
 *
 * routines for initial work, reading SLC,
 * writing to internal format.
 * see e.g: http://earth.esa.int:81/sarslc
 ****************************************************************/


#include "matrixbk.hh"
#include "constants.hh"
#include "slcimage.hh"                  // my slc image class
#include "ioroutines.hh"                // ?
#include "exceptions.hh"                // my exceptions class
#include "coregistration.hh"            // for kernel raised cosine

#include <fstream>                      // for streams
#include <iomanip>                      // setw
#include <cstdlib>                      // atoi, system
#include <cstring>                      // strcmp
#include <cstdio>                       // some compilers, remove function
#include <ctime>                        // some compilers, strptime function
#ifdef WIN32
  #include "winsock2.h"                 // Jia changed this
#else
  #include <netinet/in.h>                 // ntohl byteorder x86-HP unix
#endif
#include <algorithm>        // max,min
#include "utilities.hh"     // sinc() cr4toci2()
char *strptime(const char *s, const char  *format,  struct tm *tm);



/****************************************************************
 *    julday()                                                  *
 * mm JAN=1, FEB=2, etc; id: day=1,2,3.., iyyy=1996             *
 *                                                              *
 #%// Bert Kampes, 07-Apr-2005
 ****************************************************************/
#define IGREG (15+31L*(10+12L*1582))
int32 julday(int32 id, int32 mm, int32 iyyy)
  {
  int32 jul;
  int32 ja,jy,jm;
  if (iyyy==0)
    {
    //PRINT_ERROR("julday: error")
    //throw(some_error);
    WARNING.print("julday: year=0; impossible; continuing (only for Btemp computation)");
    return 0;
    }
  if (iyyy<0) ++iyyy;
  if (mm>2)
    {
    jy=iyyy;
    jm=mm+1;
    }
  else
    {
    jy=iyyy-1;
    jm=mm+13;
    }
  jul = int32(floor(365.25*jy)+floor(30.6001*jm)+id+1720995);
  if (id+31L*(mm+12L*iyyy) >= IGREG)
    {
    ja   = int32(0.01*jy);
    jul += 2-ja+int32(0.25*ja);
    }
  return jul;
  }
#undef IGREG



/****************************************************************
 *    readvolume                                                *
 *                                                              *
 * reads volumefile                                             *
 *  and writes to scratchlogfile and resfile                    *
 *  checks if master is same as slave by id                     *
 * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
 *      doc:er-is-epo-gs-5902.3                                 *
 http://earth.esa.int/rootcollection/sysutil/01008.html
 *                                                              *
 * input:                                                       *
 *  - struct: input options for readfiles                       *
 *  - 3 checks with other (slave or master) volumefile          *
 * output:                                                      *
 *  - void, scratchfiles: ?                                     *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 * if radarsat, skip record number 4 of 360 bytes.              *
 * according to specs                                           *
 #%// Bert Kampes, 02-Aug-2004                                  *
 #%// Davide Nitti (Don), 11-Nov-2008 Reader update for ALOS    *
 ****************************************************************/
void readvolume(
        input_readfiles &readfiles_arg,
        const char* checkvol1,
        const char* checkvol2,
        const char* checkvol3)
  {
  TRACE_FUNCTION("readvolume (BK 11-Dec-1998)")
  const int16           sizeb1 = 1,             // some constants for reading
                        sizeb4 = 4,             //   binary fields.
                        sizei4 = 4,             //   binary fields.
                        sizei8 = 8,
                        sizea8 = 8,
                        sizea12 = 12,
                        sizea16 = 16,
                        sizea28 = 28,
                        sizea40 = 40,
                        sizea60 = 60;
  uint                  lenrec1,                // length of record1
                        lenrec2,                // length of record2
                        lenrec3,                // length of record3
                        lenrec4,                // length of record4
                        numpointrec,            // volfile...
                        numrec,                 // volfile...
                        numvol;                 // volfile...
  char                  c4dummy[5],             // correctly 5 for \0
                        c8date[9],
                        c8time[9],
                        c8agency[9],            // correctly 9 for \0
                        c8nlins[9],             // correctly 9 for \0
                        c12logvol[13],
                        c12country[13],
                        c12facility[13],        // correctly 13 for \0
                        c16physid[17],
                        c16logvolid[17],
                        c16setid[17],
                        c16checkfilename[17],
                        c16dataref[17],         // correctly 17 for \0
                        c28leaderrefclass[29],
                        c28datarefclass[29],
                        c40typespec[41],
                        c40physvolid[41],
                        c40sceneid[41],
                        c40sceneloc[41],
                        c60product[61];
  // --- Check for Atlantis processor (RSAT?) #%// Bert Kampes, 02-Aug-2004 ---
  uint rec_seq;// type B4
  unsigned char rec_sub1, rec_type, rec_sub2, rec_sub3;// type B1

// ======Open files======
  // ___ check if opened correctly, if not, try to use uppercase
  // ___ from SCENE1/lea_01.001 #%// BK 27-Nov-2003
  ifstream volumefile;
  openfstream(volumefile,readfiles_arg.volfile);
  bk_assert(volumefile,readfiles_arg.volfile,__FILE__,__LINE__);

// ______First check control arguments, allow check to be wrong ______
  DEBUG.print("readvol: check on filename consistency.");
  for (register int32 i=0; i<999; i++)                          // try some times
    {
    volumefile.seekg(44,ios::beg);
    volumefile.read((char*)&c16physid,sizea16);         // physical logical volume ID
    c16physid[16]='\0';
    volumefile.read((char*)&c16logvolid,sizea16);       // logical volume ID
    c16logvolid[16]='\0';
    volumefile.read((char*)&c16setid,sizea16);          // volume set ID
    c16setid[16]='\0';

    if (readfiles_arg.sensor_id==SLC_ALOS)   // added by don
      {
      char c8logvoltime[9];
      volumefile.seekg(120,ios::beg);
      volumefile.read((char*)&c8logvoltime,sizea8);             // Logical volume creation time
      c8logvoltime[8]='\0';
      for (register int32 idx=0; idx<8; idx++)
        c16physid[idx+8]=c8logvoltime[idx];
      }

    if ((!strcmp(checkvol1,c16physid)) &&                       // check failed
        (!strcmp(checkvol2,c16logvolid)) &&
        (!strcmp(checkvol3,c16setid)))
      {
      if (i==0) // first time
        {
        WARNING.print("Volume file of master and slave seem to be the same, please change cd.");
        getanswer();                                    // wait
        }

      else if (i==1)    // second time
        {
        WARNING << "ID of volume file 1: "
             << c16physid << ": " << c16logvolid << ": " << c16setid;
        WARNING.print();
        WARNING << "ID of volume file 2: "
             << checkvol1 << ":" << checkvol2 << ":" << checkvol3;
        WARNING.print();
        WARNING.print("Next time I will assume they are different.");
        getanswer();                                    // wait
        }

      else // if i == 2, enough tries
        {
        break;
        }
      }
    else                                                // check passed
      {
      break;
      }
    }


  // ======Read volumefile======
  // --- RECORD 1 ---
  DEBUG.print("record 1 of volume file.");
  volumefile.seekg(0,ios::beg);// place pointer at beginning of file
  volumefile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  volumefile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  volumefile.read((char*)&rec_type,sizeb1);// record type code
  volumefile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  volumefile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("Expecting record 1 with code {192,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // --- Read important fields ---
  volumefile.seekg(8,ios::beg);
  volumefile.read((char*)&lenrec1,sizeb4);              // length of record1
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  if (lenrec1 != 360 )
    {
    WARNING.print("probably something wrong here, byte order x86 (?).");
    WARNING << "readvolume: length of record 1 = \""
         <<  lenrec1 << "\"; expected \"360\" for ESA SLC (full scene).";
    WARNING.print();
    }
  volumefile.seekg(32,ios::beg);
  volumefile.read((char*)&c12logvol,sizea12);           // logical volume etc.
  c12logvol[12]='\0';
  volumefile.seekg(112,ios::beg);
  volumefile.read((char*)&c8date,sizea8);               // generating date
  c8date[8]='\0';
  volumefile.read((char*)&c8time,sizea8);               // generating time
  c8time[8]='\0';
  volumefile.read((char*)&c12country,sizea12);          // generating country
  c12country[12]='\0';
  volumefile.read((char*)&c8agency,sizea8);             // generating agency
  c8agency[8]='\0';
  volumefile.read((char*)&c12facility,sizea12);         // generating facility
  c12facility[12]='\0';
  volumefile.read((char*)&c4dummy,sizei4);              // #pointer records in vol.
  c4dummy[4]='\0';
  numpointrec = atoi(c4dummy);

  // if (readfiles_arg.sensor_id!=SLC_RSAT)
  // Modified by LG for reading ALOS Fine
  if ((readfiles_arg.sensor_id!=SLC_RSAT) && (readfiles_arg.sensor_id!=SLC_ALOS) )
    {
    if (numpointrec!=2)
      {
      ERROR << "readvolume: number of pointer records = \""
           <<  c4dummy << "\"; expected \"2\" for ESA SLC (full scene)." << ends;
      WARNING.print(ERROR.get_str());// Jan Kianicka report atlantis processor
      WARNING.print("but for Atlantis processor this may be correct?");
      ERROR.reset();
      }
    }
  else
    if (numpointrec!=3)
      WARNING.print("expected 3 pointer records for Radarsat?");

  volumefile.read((char*)&c4dummy,sizei4);              // #records in vol.
  c4dummy[4]='\0';
  numrec = atoi(c4dummy);
  // Modified by LG for reading ALOS Fine
  if (readfiles_arg.sensor_id == SLC_ALOS)
  {
          if (numrec!=1)
      {
      ERROR << "readvolume: number of records = \""
           <<  c4dummy << "\"; expected \"1\" for JAXA SLC (full scene)." << ends;
      WARNING.print(ERROR.get_str());
      ERROR.reset();
      }
  }
  else if (readfiles_arg.sensor_id!=SLC_RSAT)
    {
    if (numrec!=4)
      {
      ERROR << "readvolume: number of records = \""
           <<  c4dummy << "\"; expected \"4\" for ESA SLC (full scene)." << ends;
      WARNING.print(ERROR.get_str());
      ERROR.reset();
      }
    }
  else
    if (numrec!=5)
      WARNING.print("expected 5 records for Radarsat?");
  volumefile.read((char*)&c4dummy,sizei4);              // #vol.
  c4dummy[4]='\0';
  numvol = atoi(c4dummy);
  if (readfiles_arg.sensor_id!=SLC_RSAT)
    if (numvol!=1)
      {
      ERROR << "readvolume: number of volumes = \""
           <<  c4dummy << "\"; expected \"1\" for ESA SLC (full scene)." << ends;
      WARNING.print(ERROR.get_str());
      ERROR.reset();
      }


  // --- RECORD 2 ---
  const uint startrec2 = lenrec1;
  volumefile.seekg(startrec2,ios::beg);
  volumefile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  volumefile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  volumefile.read((char*)&rec_type,sizeb1);// record type code
  volumefile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  volumefile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("Expecting record 2 with code {219,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // --- Read important fields ---
  volumefile.seekg(startrec2+8,ios::beg);               //  leader file pointer rec.
  DEBUG.print("record 2 of volume file.");
  volumefile.read((char*)&lenrec2,sizeb4);              // length of record2
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.
  if (lenrec2 != 360 )
    {
    WARNING.print("probably something wrong here, byte order x86. (?)");
    WARNING << "readvolume: length of record 2 = \""
         <<  lenrec2 << "\"; expected \"360\" for ESA SLC (full scene).";
    WARNING.print();
    }
  volumefile.seekg(startrec2+20,ios::beg);
  volumefile.read((char*)&c16checkfilename,sizea16);    // referenced file name
  c16checkfilename[16]='\0';
  volumefile.read((char*)&c28leaderrefclass,sizea28);   // referenced file class
  c28leaderrefclass[28]='\0';


  // --- RECORD 3 ---
  DEBUG.print("record 3 of volume file.");
  const uint startrec3 = lenrec1 + lenrec2;
  volumefile.seekg(startrec3,ios::beg);
  volumefile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  volumefile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  volumefile.read((char*)&rec_type,sizeb1);// record type code
  volumefile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  volumefile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("Expecting record 3 with code {219,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // --- Read important fields ---
  volumefile.seekg(startrec3+8,ios::beg);               //  data file pointer rec.
  volumefile.read((char*)&lenrec3,sizeb4);              // length of record3
  lenrec3 = ntohl(lenrec3);     // bk 6 jul 2000, byteorder x86 machines.
  if (lenrec3 != 360 )
    {
    WARNING.print("probably something wrong here, byte order x86. (?)");
    WARNING << "readvolume: length of record 3 = \""
         <<  lenrec3 << "\"; expected \"360\" for ESA SLC (full scene).";
    WARNING.print();
    }
  volumefile.seekg(startrec3+20,ios::beg);
  volumefile.read((char*)&c16dataref,sizea16);          // referenced file name
  c16dataref[16]='\0';
  volumefile.read((char*)&c28datarefclass,sizea28);     // referenced file class
  c28datarefclass[28]='\0';
  volumefile.seekg(startrec3+100,ios::beg);             // numlines for checking
  volumefile.read((char*)&c8nlins,sizei8);
  c8nlins[8]='\0';


  // --- RECORD 4 ---
  uint startrec4 = lenrec1 + lenrec2 + lenrec3;
  //if (readfiles_arg.sensor_id==SLC_RSAT)
  //   startrec4=startrec4+360;// skip trailer record for RSAT
  DEBUG.print("record 4 of volume file.");
  DEBUG << "readvolume::rec4: start at byte: " << startrec4;
  DEBUG.print();
  // --- Check this record, RSAT has for us useless trailer info as rec.4 ---
  // --- while next record of 360 is the text record of ERS ---
  // ---start test for ATLANTIS (RSAT) ------------------------------------------
  DEBUG.print("VMP: Expecting record 4 with code {18,63,18,18}");
  volumefile.seekg(startrec4,ios::beg);
  volumefile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  volumefile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  volumefile.read((char*)&rec_type,sizeb1);// record type code
  volumefile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  volumefile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==63 && int(rec_sub2)==18 && int(rec_sub3)==18)
    {
    DEBUG.print("This is the expected text record with code {18,63,18,18}");
    readfiles_arg.sar_processor = SARPR_VMP;// set determined sar processor/format
    }
  else
    {
    //WARNING.print("This is not the expected text record, trying next one.")
    readfiles_arg.sar_processor = SARPR_ATL;// set determined sar processor/format

        // Modified by LG for reading ALOS Fine
        if ( (readfiles_arg.sensor_id == SLC_ALOS) && (int(rec_type)==192) )
                         readfiles_arg.sar_processor = SARPR_JAX;
        else if (readfiles_arg.sensor_id!=SLC_RSAT)
      {
      DEBUG.print("This is NOT the expected text record with code {18,63,18,18}");
      DEBUG.print("I assume in following Atlantis processed data (ERS or RSAT).");
      DEBUG.print("If this is not the case, please contact doris_users@tudelft.nl");
      }
    volumefile.seekg(startrec4+8,ios::beg);             //  text record
    volumefile.read((char*)&lenrec4,sizeb4);            // length of record4
    lenrec4   = ntohl(lenrec4); // bk 6 jul 2000, byteorder x86 machines.
    startrec4 = startrec4 + lenrec4;// hopefully there is a next record
    }
  // ---end test for ATLANTIS (RSAT) -------------------------------------------

  // --- Read important fields ---
  volumefile.seekg(startrec4+8,ios::beg);               //  text record
  volumefile.read((char*)&lenrec4,sizeb4);              // length of record4
  lenrec4 = ntohl(lenrec4);     // bk 6 jul 2000, byteorder x86 machines.
  if (lenrec4 != 360 )
    {
    WARNING.print("probably something wrong here, byte order x86. (?)");
    WARNING << "readvolume: length of record 3 = \""
         <<  lenrec4 << "\"; expected \"360\" for ESA SLC (full scene).";
    WARNING.print();
    }
  volumefile.seekg(startrec4+16,ios::beg);
  volumefile.read((char*)&c40typespec,sizea40);         // product type specifier
  c40typespec[40]='\0';
  volumefile.read((char*)&c60product,sizea60);          // loc&date product gen.
  c60product[60]='\0';
  volumefile.read((char*)&c40physvolid,sizea40);        // physical vol id
  c40physvolid[40]='\0';
  volumefile.read((char*)&c40sceneid,sizea40);          // scene id
  c40sceneid[40]='\0';
  volumefile.read((char*)&c40sceneloc,sizea40);         // scene loc
  c40sceneloc[40]='\0';

  volumefile.close();


// ====== Write results to scratch files ======
  ofstream scratchlogfile("scratchlogvol", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"scratchlogvol",__FILE__,__LINE__);

// ______Write information to scratchfiles______
  scratchlogfile << "\n*******************************************************************"
                 << "\n* EXTRACTED DATA FROM VOLUME FILE: "
                 <<  readfiles_arg.volfile << " *"
                 << "\n*******************************************************************"

                 << "\n\nVolume descriptor record"
                 << "\n------------------------"
                 << "\nLogical volume generating facility software "
                 << "\n +release and revision level: \t\t\t"
                 <<  c12logvol
                 << "\nID of physical volume containing "
                 << "\n +this volume descriptor: \t\t\t"
                 <<  c16physid
                 << "\nLogical volume identifier: \t\t\t"
                 <<  c16logvolid
                 << "\nVolume set identifier: \t\t\t\t"
                 <<  c16setid
                 << "\nLogical volume creation date (YYYYMMDD): \t"
                 <<  c8date
                 << "\nLogical volume creation time (HHMMSSDD): \t"     // DD=deci-secs
                 <<  c8time
                 << "\nLogical volume generation country: \t\t"
                 <<  c12country
                 << "\nLogical volume generating agency: \t\t"
                 <<  c8agency
                 << "\nLogical volume generation facility: \t\t"
                 <<  c12facility

                 << "\n\nLeader file pointer record"
                 << "\n--------------------------"
                 << "\nReferenced file name: \t\t\t\t"
                 <<  c16checkfilename
                 << "\nReferenced file class: \t\t\t\t"
                 <<  c28leaderrefclass

                 << "\n\nData file pointer record"
                 << "\n------------------------"
                 << "\nReferenced file name: \t\t\t\t"
                 <<  c16dataref
                 << "\nReferenced file class: \t\t\t\t"
                 <<  c28datarefclass
                 << "\nNumber of records in referenced file: \t\t"
                 <<  c8nlins

                 << "\n\nText record"
                 << "\n-----------"
                 << "\nProduct type specifier: \t\t\t"
                 <<  c40typespec
                 << "\nLocation and date/time of product creation: \t"
                 <<  c60product
                 << "\nPhysical volume identification: \t\t"
                 <<  c40physvolid
                 << "\nScene identification: \t\t\t\t"
                 <<  c40sceneid
                 << "\nScene location: \t\t\t\t"
                 <<  c40sceneloc
                 << "\nEND VOLUME FILE "
                 <<  readfiles_arg.volfile
                 <<  endl << endl;
  scratchlogfile.close();

  ofstream scratchresfile("scratchresvol", ios::out | ios::trunc);
  bk_assert(scratchresfile,"scratchresvol",__FILE__,__LINE__);
  scratchresfile
    << "\n*******************************************************************";
  if (readfiles_arg.fileid==MASTERID)
    scratchresfile << "\n*_Start_" << processcontrol[pr_m_readfiles];
  if (readfiles_arg.fileid==SLAVEID)
    scratchresfile << "\n*_Start_" << processcontrol[pr_m_readfiles];
  scratchresfile
    << "\n*******************************************************************"
    << "\nVolume file: \t\t\t\t\t"
    <<  readfiles_arg.volfile
    << "\nVolume_ID: \t\t\t\t\t"
    <<  c16physid
    << "\nVolume_identifier: \t\t\t\t"
    <<  c16logvolid
    << "\nVolume_set_identifier: \t\t\t\t"
    <<  c16setid
    << "\n(Check)Number of records in ref. file: \t\t"
    <<  atoi(c8nlins);
  // --- write new string to determine sar processor ---
  // --- this is used in cropping the data ---
  if (readfiles_arg.sar_processor==SARPR_VMP)
    scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tVMP";
  if (readfiles_arg.sar_processor==SARPR_ATL)
    scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tATLANTIS";
  if (readfiles_arg.sar_processor==SARPR_TUD)
    scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tTUD";
  // start_added_by_don
  if (readfiles_arg.sar_processor==SARPR_JAX)
    scratchresfile << "\nSAR_PROCESSOR:  \t\t\t\tJAX";
  // end_added_by_don
  if (readfiles_arg.sensor_id==SLC_RSAT)
    scratchresfile
      << "\nProduct type specifier: \t\t\t"
      <<  "RSAT";// else it is just "PRODUCT:", but it is used in CROP step
  // start_added_by_don
  else if (readfiles_arg.sensor_id==SLC_ALOS)
    scratchresfile
      << "\nProduct type specifier: \t\t\t"
      <<  "PRODUCT ALOS";// else it is just "PRODUCT:", but it is used in CROP step
  // end_added_by_don
  else
    scratchresfile
      << "\nProduct type specifier: \t\t\t"
      <<  c40typespec;
  scratchresfile
    << "\nLogical volume generating facility: \t\t\t"
    <<  c12facility
    << "\nLogical volume creation date: \t\t\t"
    <<  c8date
    << "\nLocation and date/time of product creation: \t"
    <<  c60product
    << "\nScene identification: \t\t\t\t"
    <<  c40sceneid
    << "\nScene location: \t\t\t\t"
    <<  c40sceneloc
    << endl;
  scratchresfile.close();

// ______Tidy up______
  PROGRESS.print("readvolume finished.");
  } // END READVOLUME



/****************************************************************
 *    readleader                                                *
 *                                                              *
 * reads leaderfile                                             *
 *  and writes to scratchlogfile, scratchresfile                *
 * checks with volumefile #lines                                *
 *                                                              *
 * input:                                                       *
 *  - struct with arguments for step0                           *
 *  - check with volumefile                                     *
 * output:                                                      *
 *  - scratchlogfile                                            *
 *  - scratchresfile                                            *
 *  - scratchdatapoints                                         *
 * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
 *      doc:er-is-epo-gs-5902.3                                 *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 * Included RSAT format based on document of ASF                *
 #%// Bert Kampes, 03-Aug-2004                                  *
 #%// Davide Nitti (Don), 11-Nov-2008  fixes for doppler        *
 #     coefficient unit for Radarsat1 and ALOS                  *
 ****************************************************************/
void readleader(
        input_readfiles &readfiles_arg,
        const int32 checklines)
  {
  TRACE_FUNCTION("readleader (BK 11-Dec-1998)")
  const int16           sizea2 = 2,
                        sizeb1 = 1,             // some constants for reading
                        sizeb4 = 4,             // some constants for reading
                        sizei4 = 4,             //   binary fields.
                        sizea4 = 4,             //   binary fields.
                        sizei8 = 8,
                        sizea8 = 8,
                        sizea12 = 12,
                        sizea16 = 16,
                        sizef16 = 16,
                        sizei16 = 16,
                        sized22 = 22,
                        sizea24 = 24,
                        sizea32 = 32,
                        sizea64 = 64;
  uint                  lenrec1,                // length of record1
                        lenrec2,                // length of record2
                        lenrec3,                // length of record3
                        lenrec4,                // bc length of record4
                        lenrec5,                // bc/gk length of record5
                        lenrec6,                // gk length of record6
                        lenrec7;                // bk rsat record
 char                   c2motioncomp[3],
                        c4dummy[5],             // correctly 5 for \0
                        c4year[5],
                        c4month[5],
                        c4day[5],
                        c4dayofyear[5],
                        c4conversion[5],
                        c4compression[5],
                        c4clutterlock[5],
                        c4autofocus[5],
                        c8dummy[9],             // dummy
                        c8orbitnr[9],
                        c8platformlat[9],
                        c8platformlon[9],
                        c8platformheading[9],
                        c8clockangle[9],
                        c8incidence[9],
                        c8freq[9],
                        c8systemid[9],
                        c8versionid[9],
                        c8satclockstep[9],
                        c8extindex[9],
                        c8qperch[9],
                        c8timeline[9],
                        c8timepix[9],
                        c8linecontent[9],
                        c12qdesc[13],
                        c16dummy[17],
                        c16lat11[17], c16lon11[17],
                        c16lat1N[17], c16lon1N[17],
                        c16latNN[17], c16lonNN[17],
                        c16latN1[17], c16lonN1[17],
                        c16leafilename[17],
                        c16centerlat[17],
                        c16centerlon[17],
                        c16centerheading[17],
                        c16ellipsoid[17],
                        c16semimajor[17],
                        c16semiminor[17],
                        c16GM[17],
                        c16J2[17],
                        c16J3[17],
                        c16J4[17],
                        c16scenelength[17],
                        c16scenewidth[17],
                        c16platformid[17],
                        c16wavelength[17],
                        c16facilityid[17],
                        c16numpix[17],
                        c16numlin[17],
                        c16interpix[17],
                        c16interlin[17],
                        c16pulse[17],
                        c16ampconst[17],
                        c16amplinear[17],
                        c16ampquadratic[17],
                        c16ampcubic[17],
                        c16ampquartic[17],
                        c16phaseconst[17],
                        c16phaselinear[17],
                        c16phasequadratic[17],
                        c16phasecubic[17],
                        c16phasequartic[17],
                        c16sattimecode[17],
                        c16samplingrate[17],
                        c16rangedelay[17],
                        c16ranpulselen[17],
                        c16dci[17],
                        c16dcq[17],
                        c16boresight[17],
                        c16imbalance[17],
                        c16prf[17],
                        c16looksazi[17],
                        c16looksrange[17],
                        c16bandazi[17],
                        c16bandrange[17],
                        c16bandazitot[17],
                        c16bandrangetot[17],
                        c16inputsource[17],
                        c16resrange[17],
                        c16resazi[17],
                        c16linespace[17],
                        c16pixspace[17],
                        c16atdoppcconst[17],
                        c16atdoppclinear[17],
                        c16atdoppcquadratic[17],
                        c16xtdoppcconst[17],
                        c16xtdoppclinear[17],
                        c16xtdoppcquadratic[17],
                        c16atdopprconst[17],
                        c16atdopprlinear[17],
                        c16atdopprquadratic[17],
                        c16xtdopprconst[17],
                        c16xtdopprlinear[17],
                        c16xtdopprquadratic[17],
                        c16rcompdes[17],
                        c16zd1strange[17],
                        c16zdcenrange[17],
                        c16zdlstrange[17],
                        c16orien[17],
                        c16platincl[17],
                        c16platascn[17],
                        c16geocenter[17],
                        c16platalt[17],
                        c16plathead[17],
                        c16platgs[17],
                        c16refmajor[17],
                        c16refminor[17],
                        c16ltposerr[17],
                        c16ctposerr[17],
                        c16rposerr[17],
                        c22dummy[23],
                        c22seconds[23],
                        c22interval[23],
                        c22gmha[23],
                        c24zd1stazitime[25],
                        c24zdcenazitime[25],
                        c24zdlstazitime[25],
                        c32sceneref[33],
                        c32scenetime[33],
                        c32sensorid[33],
                        c32typespec[33],
                        c32algid[33],
                        c32projection[33],
                        c32sattime[33],
                        c32weightrange[33],
                        c32weightazi[33],
                        c32refellips[33],
                        c64rcs[65];

  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
  uint rec_seq;// type B4
  unsigned char rec_sub1, rec_type, rec_sub2, rec_sub3;// type B1



// ======Open files======
  ifstream leaderfile;
  openfstream(leaderfile,readfiles_arg.leaderfile);
  bk_assert(leaderfile,readfiles_arg.leaderfile,__FILE__,__LINE__);

  // ======Read leaderfile======
  // --- RECORD 1 ---
  leaderfile.seekg(0,ios::beg);// place pointer at beginning of file
  leaderfile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("Expecting record 1 with code {63,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // --- Read important data ---
  leaderfile.seekg(8,ios::beg);                         // file descriptor record
  leaderfile.read((char*)&lenrec1,sizeb4);              // length of record1
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 1 of leader file: length: " << lenrec1;
  DEBUG.print();
  if (lenrec1 != 720 )
    {
    WARNING.print("probably something wrong here, byte order x86.");
    WARNING << "readleader: length of record 1 = \""
         <<  lenrec1 << "\"; expected \"720\" for ESA SLC (full scene).";
    WARNING.print();
    }
  leaderfile.seekg(48,ios::beg);
  leaderfile.read((char*)&c16leafilename,sizea16);      // file name
  c16leafilename[16]='\0';


  // --- RECORD 2 ---
  DEBUG.print("readleader::reading record 2 of leader file.");
  uint startrec2 = lenrec1;
  leaderfile.seekg(startrec2,ios::beg);// place pointer at beginning of record
  leaderfile.read((char*)&rec_seq,sizeb4);// record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 2 with code {10,10,31,20}");
  DEBUG.print("RSAT: Expecting record 2 with code {18,10,18,20}");
  DEBUG.print("RSAT record length is 4096, ERS 1886, but");
  DEBUG.print("ERS contains more info on zero doppler times, etc.");
  DEBUG.print("RSAT seems to have that info in the data file.");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // ___ Read important parameters ___
  leaderfile.seekg(startrec2+8,ios::beg);               // slc data set summary record
  leaderfile.read((char*)&lenrec2,sizeb4);              // length of record2
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 2: start at: " << lenrec1 << "; length: " << lenrec2;
  DEBUG.print();
  //if (readfiles_arg.sensor_id==SLC_RSAT)
  if (readfiles_arg.sar_processor==SARPR_ATL)
    {
    if (lenrec2 != 4096)
      {
      WARNING.print("SARPR_ATL (RSAT) has 4096 record length, but other value found?");
      }
    }
  // start_added_by_don
  else if (readfiles_arg.sar_processor==SARPR_JAX)
    {
    if (lenrec2 != 4096)
      {
      WARNING.print("SARPR_JAX (ALOS) has 4096 record length, but other value found?");
      }
    }
  // end_added_by_don
  else
    {
    if (lenrec2 != 1886)
      {
      WARNING.print("SARPR_ATL (RSAT) has 4096 record length");
      WARNING.print("probably something wrong here, byte order on x86?");
      WARNING << "readleader: length of record 2 = \""
           <<  lenrec2 << "\"; expected \"1886\" for ESA SLC (full scene).  continuing";
      WARNING.print();
      }
    }


  // ______Scene parameters______
 // Modified by LG for reading ALOS Fine
 if (readfiles_arg.sensor_id == SLC_ALOS)
 {
         leaderfile.seekg(startrec2+20,ios::beg);
         leaderfile.read((char*)&c32sceneref,sizea32);
         leaderfile.seekg(startrec2+20+sizea32+1,ios::beg);
 }
 else
 {
         leaderfile.seekg(startrec2+36,ios::beg);
         leaderfile.read((char*)&c32sceneref,sizea32);          // scene ref. number
 }
  c32sceneref[32]='\0';
  leaderfile.read((char*)&c32scenetime,sizea32);        // scene center time
  c32scenetime[32]='\0';
  leaderfile.seekg(startrec2+116,ios::beg);
  leaderfile.read((char*)&c16centerlat,sizef16);        // centre latitude
  c16centerlat[16]='\0';
  leaderfile.read((char*)&c16centerlon,sizef16);        // centre longitude
  c16centerlon[16]='\0';
  leaderfile.read((char*)&c16centerheading,sizef16);    // center true heading
  c16centerheading[16]='\0';
  leaderfile.read((char*)&c16ellipsoid,sizea16);        // ell. designated
  c16ellipsoid[16]='\0';
  leaderfile.read((char*)&c16semimajor,sizef16);        // ell. semi major
  c16semimajor[16]='\0';
  leaderfile.read((char*)&c16semiminor,sizef16);        // ell. semi minor
  c16semiminor[16]='\0';
  leaderfile.read((char*)&c16GM,sizef16);               // GM
  c16GM[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // dummy
  leaderfile.read((char*)&c16J2,sizef16);               // J2
  c16J2[16]='\0';
  leaderfile.read((char*)&c16J3,sizef16);               // J3
  c16J3[16]='\0';
  leaderfile.read((char*)&c16J4,sizef16);               // J4
  c16J4[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // dummy
  leaderfile.read((char*)&c16dummy,sizef16);            // dummy
  leaderfile.read((char*)&c8dummy,sizei8);              // center line#
  c8dummy[8]='\0';
  uint scenecenterline = atoi(c8dummy);
  leaderfile.read((char*)&c8dummy,sizei8);              // center pixel#
  c8dummy[8]='\0';
  uint scenecenterpixel = atoi(c8dummy);
  leaderfile.read((char*)&c16scenelength,sizef16);      // scene length
  c16scenelength[16]='\0';
  leaderfile.read((char*)&c16scenewidth,sizef16);       // scene width
  c16scenewidth[16]='\0';

// ______General mission / sensor parameters______
  leaderfile.seekg(startrec2+396,ios::beg);
  leaderfile.read((char*)&c16platformid,sizea16);       // platform mission id
  c16platformid[16]='\0';
  leaderfile.read((char*)&c32sensorid,sizea32);         // sensor id
  c32sensorid[32]='\0';
  leaderfile.read((char*)&c8orbitnr,sizea8);            // orbit number
  c8orbitnr[8]='\0';
  leaderfile.read((char*)&c8platformlat,sizea8);        // platform latitude
  c8platformlat[8]='\0';
  leaderfile.read((char*)&c8platformlon,sizea8);        // platform longitude
  c8platformlon[8]='\0';
  leaderfile.read((char*)&c8platformheading,sizea8);    // platform heading
  c8platformheading[8]='\0';
  leaderfile.read((char*)&c8clockangle,sizea8);         // sensor clock angle
  c8clockangle[8]='\0';
  leaderfile.read((char*)&c8incidence,sizea8);          // incidence angle
  c8incidence[8]='\0';
  leaderfile.read((char*)&c8freq,sizea8);               // radar frequency
  c8freq[8]='\0';
  leaderfile.read((char*)&c16wavelength,sizea16);       // radar wavelength
  c16wavelength[16]='\0';
  leaderfile.read((char*)&c2motioncomp,sizea2);         // indicator for compensation
  c2motioncomp[2]='\0';
  leaderfile.read((char*)&c16pulse,sizea16);            // range pulse code specifier
  c16pulse[16]='\0';
  leaderfile.read((char*)&c16ampconst,sizef16);         // amplitude constant term
  c16ampconst[16]='\0';
  leaderfile.read((char*)&c16amplinear,sizef16);        // amplitude linear term
  c16amplinear[16]='\0';
  leaderfile.read((char*)&c16ampquadratic,sizef16);     // amplitude quadrati term
  c16ampquadratic[16]='\0';
  leaderfile.read((char*)&c16ampcubic,sizef16);         // amplitude cubic term
  c16ampcubic[16]='\0';
  leaderfile.read((char*)&c16ampquartic,sizef16);       // amplitude quartic term
  c16ampquartic[16]='\0';
  leaderfile.read((char*)&c16phaseconst,sizef16);       // phase constant term
  c16phaseconst[16]='\0';
  leaderfile.read((char*)&c16phaselinear,sizef16);      // phase linear term
  c16phaselinear[16]='\0';
  leaderfile.read((char*)&c16phasequadratic,sizef16);   // phase quadratic term
  c16phasequadratic[16]='\0';
  leaderfile.read((char*)&c16phasecubic,sizef16);       // phase cubicterm
  c16phasecubic[16]='\0';
  leaderfile.read((char*)&c16phasequartic,sizef16);     // phase quartic term
  c16phasequartic[16]='\0';
  leaderfile.read((char*)&c8extindex,sizei8);           // chirp extraction
  c8extindex[8]='\0';
  leaderfile.read((char*)&c8dummy,sizei8);              // spare
  leaderfile.read((char*)&c16samplingrate,sizef16);     // range sampling rate
  c16samplingrate[16]='\0';
  leaderfile.read((char*)&c16rangedelay,sizef16);       // delay
  c16rangedelay[16]='\0';
  leaderfile.read((char*)&c16ranpulselen,sizef16);      // range pulselength
  c16ranpulselen[16]='\0';
  leaderfile.read((char*)&c4conversion,sizea4);         // flag
  c4conversion[4]='\0';
  leaderfile.read((char*)&c4compression,sizea4);        // flag
  c4compression[4]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // reserved
  leaderfile.read((char*)&c16dummy,sizef16);            // reserved
  leaderfile.read((char*)&c8qperch,sizei8);             // quantization
  c8qperch[8]='\0';
  leaderfile.read((char*)&c12qdesc,sizea12);            // quantization description
  c12qdesc[12]='\0';
  leaderfile.read((char*)&c16dci,sizef16);              // bias for i comp.
  c16dci[16]='\0';
  leaderfile.read((char*)&c16dcq,sizef16);              // bias for q comp.
  c16dcq[16]='\0';
  leaderfile.read((char*)&c16imbalance,sizef16);        // gain imbalance i&q
  c16imbalance[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // spare
  leaderfile.read((char*)&c16dummy,sizef16);            // spare
  leaderfile.read((char*)&c16dummy,sizef16);            // reserved
  leaderfile.read((char*)&c16boresight,sizef16);        // antenna
  c16boresight[16]='\0';
  leaderfile.read((char*)&c4dummy,sizea4);              // reserved
  leaderfile.read((char*)&c16prf,sizef16);              // pulse repetition frequency
  c16prf[16]='\0';

// ______Sensor specific parameters______
  leaderfile.seekg(startrec2+982,ios::beg);
  leaderfile.read((char*)&c16sattimecode,sizei16);      // sat time code
  c16sattimecode[16]='\0';
  leaderfile.read((char*)&c32sattime,sizea32);          // sat time
  c32sattime[32]='\0';
  leaderfile.read((char*)&c8satclockstep,sizei8);       // sat clock step length
  c8satclockstep[8]='\0';

// ______General processing parameters______
  leaderfile.seekg(startrec2+1046,ios::beg);
  leaderfile.read((char*)&c16facilityid,sizea16);       // proc. facility id
  c16facilityid[16]='\0';
  leaderfile.read((char*)&c8systemid,sizea8);           // proc. system id
  c8systemid[8]='\0';
  leaderfile.read((char*)&c8versionid,sizea8);          // proc. version id
  c8versionid[8]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // dummy
  leaderfile.read((char*)&c16dummy,sizef16);            // dummy
  leaderfile.read((char*)&c32typespec,sizea32);         // produkt type spec.
  c32typespec[32]='\0';
  leaderfile.read((char*)&c32algid,sizea32);            // proc. alg. id
  c32algid[32]='\0';
  leaderfile.read((char*)&c16looksazi,sizef16);         // number of looks
  c16looksazi[16]='\0';
  leaderfile.read((char*)&c16looksrange,sizef16);       // number of looks
  c16looksrange[16]='\0';
  leaderfile.read((char*)&c16bandazi,sizef16);          // bandwidth
  c16bandazi[16]='\0';
  leaderfile.read((char*)&c16bandrange,sizef16);        // bandwidth
  c16bandrange[16]='\0';
  leaderfile.read((char*)&c16bandazitot,sizef16);       // bandwidth
  c16bandazitot[16]='\0';
  leaderfile.read((char*)&c16bandrangetot,sizef16);     // bandwidth
  c16bandrangetot[16]='\0';
  leaderfile.read((char*)&c32weightazi,sizea32);        // weighting function
  c32weightazi[32]='\0';
  leaderfile.read((char*)&c32weightrange,sizea32);      // weighting function
  c32weightrange[32]='\0';
  leaderfile.read((char*)&c16inputsource,sizea16);      // data input
  c16inputsource[16]='\0';
  leaderfile.read((char*)&c16resrange,sizef16);         // resolution
  c16resrange[16]='\0';
  leaderfile.read((char*)&c16resazi,sizef16);           // resolution
  c16resazi[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // reserved
  leaderfile.read((char*)&c16dummy,sizef16);            // reserved
  leaderfile.read((char*)&c16atdoppcconst,sizef16);     // along track centroid
  c16atdoppcconst[16]='\0';
  leaderfile.read((char*)&c16atdoppclinear,sizef16);    // along track centroid
  c16atdoppclinear[16]='\0';
  leaderfile.read((char*)&c16atdoppcquadratic,sizef16); // along track centroid
  c16atdoppcquadratic[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // spare
  leaderfile.read((char*)&c16xtdoppcconst,sizef16);     // cross track centroid
  c16xtdoppcconst[16]='\0';
  leaderfile.read((char*)&c16xtdoppclinear,sizef16);    // cross track centroid
  c16xtdoppclinear[16]='\0';
  leaderfile.read((char*)&c16xtdoppcquadratic,sizef16); // cross track centroid
  c16xtdoppcquadratic[16]='\0';
  leaderfile.read((char*)&c8timepix,sizea8);            // time direction
  c8timepix[8]='\0';
  leaderfile.read((char*)&c8timeline,sizea8);           // time direction
  c8timeline[8]='\0';
  leaderfile.read((char*)&c16atdopprconst,sizef16);     // along track rate
  c16atdopprconst[16]='\0';
  leaderfile.read((char*)&c16atdopprlinear,sizef16);    // along track rate
  c16atdopprlinear[16]='\0';
  leaderfile.read((char*)&c16atdopprquadratic,sizef16); // along track rate
  c16atdopprquadratic[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // spare
  leaderfile.read((char*)&c16xtdopprconst,sizef16);     // cross track rate
  c16xtdopprconst[16]='\0';
  leaderfile.read((char*)&c16xtdopprlinear,sizef16);    // cross track rate
  c16xtdopprlinear[16]='\0';
  leaderfile.read((char*)&c16xtdopprquadratic,sizef16); // cross track rate
  c16xtdopprquadratic[16]='\0';
  leaderfile.read((char*)&c16dummy,sizef16);            // spare
  leaderfile.read((char*)&c8linecontent,sizea8);        // indicator
  c8linecontent[8]='\0';
  leaderfile.read((char*)&c4clutterlock,sizea4);        // flag
  c4clutterlock[4]='\0';
  leaderfile.read((char*)&c4autofocus,sizea4);          // flag
  c4autofocus[4]='\0';
  leaderfile.read((char*)&c16linespace,sizef16);        //
  c16linespace[16]='\0';
  leaderfile.read((char*)&c16pixspace,sizef16);         //
  c16pixspace[16]='\0';
  leaderfile.read((char*)&c16rcompdes,sizea16);         // range compression designator
  c16rcompdes[16]='\0';



  // --- RSAT does not fill this spare part of the record (blanks) ---
  DEBUG.print("Assuming user has specified method RSAT if it is RSAT.");
  DEBUG.print("although we could also easily detect it before");
  bool skipmapprojrecord = false;       // some problem s with old IPAF?
  uint numdatapoints=99999;
  char  c16incangle1strange[17],                        //bc
        c16incanglecenrange[17],                        //bc
        c16incanglelstrange[17],                        //bc
        calK[17],                                       //gk
        repplspwr[17];                                  //gk
  char  c4numvalid[5];                                  //bk
  char  c4numinvalid[5];                                //bk
  strcpy(c16incangle1strange,     "skipped");
  strcpy(c16incanglecenrange,     "skipped");
  strcpy(c16incanglelstrange,     "skipped");
  strcpy(calK,                    "skipped");
  strcpy(repplspwr,               "skipped");
  strcpy(c4numvalid,              "999");
  strcpy(c4numinvalid,            "999");
  matrix<real8> STATE;// has to be declared at large scope, used later.
  matrix<real8> STATE_INERTIAL;// has to be declared at large scope, used later.


//if (readfiles_arg.sar_processor!=SARPR_ATL)
// seems ERS data are same, also with Atlantis, so only check for RSAT here.
if (readfiles_arg.sensor_id!=SLC_RSAT)
  {
  DEBUG.print("Reading rest of this record, which for RSAT is empty.");

  // ______Sensor specific local use segment______
  leaderfile.seekg(startrec2+1766,ios::beg);
  leaderfile.read((char*)&c16zd1strange,sizef16);       // zero doppler 1st pixel
  c16zd1strange[16]='\0';
  leaderfile.read((char*)&c16zdcenrange,sizef16);       // zero doppler centre pixel
  c16zdcenrange[16]='\0';
  leaderfile.read((char*)&c16zdlstrange,sizef16);       // zero doppler last pixel 2way
  c16zdlstrange[16]='\0';
  leaderfile.read((char*)&c24zd1stazitime,sizea24);     // zero doppler 1st pixel
  c24zd1stazitime[24]='\0';
  leaderfile.read((char*)&c24zdcenazitime,sizea24);     // zero doppler 1st pixel
  c24zdcenazitime[24]='\0';
  leaderfile.read((char*)&c24zdlstazitime,sizea24);     // zero doppler 1st pixel
  c24zdlstazitime[24]='\0';


  // --- RECORD 3 ---
  DEBUG.print("record 3 of leader file (ERS).");
  uint startrec3 = lenrec1 + lenrec2;
  leaderfile.seekg(startrec3,ios::beg);//     map projection data record
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 3 with code {10,20,31,20}");
  DEBUG.print("RSAT: Expecting record 3 with code {18,60,18,20}");
  DEBUG.print("RSAT record length is 4096, ERS 1886, but");
  DEBUG.print("ERS contains more info on zero doppler times, etc.");
  DEBUG.print("RSAT seems to have that in data file");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();

  // ___ Read important information ___
  leaderfile.seekg(startrec3+8,ios::beg);               //  map projection data record
  leaderfile.read((char*)&lenrec3,sizeb4);              // length of record3
  lenrec3 = ntohl(lenrec3);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 3: start at: " << startrec3 << "; length: " << lenrec3;
  DEBUG.print();

  //if (lenrec3 != 1620 ) // quarter scene 1046 (?)
  if (lenrec3 < 1620 ) // quarter scene 1046 (?)
    {
    WARNING.print("Probably something wrong here, byte order on x86?");
    WARNING.print("you may also be trying to read leader of ESA.RAW not ESA.SLC.");
    WARNING.print("format (or quarter scene).");
    WARNING.print("We try to do something useful still, but be careful.");
    WARNING.print("Skipping map projection record, not in ESA.RAW leaderfile format.");
    WARNING << "readleader: length of record 3 = \""
         <<  lenrec3 << "\"; expected \"1620\" for ESA SLC (full scene).";
    WARNING.print();
    lenrec3 = 0;
    skipmapprojrecord = true;
    }
  // --- Trying to read anyway if it is longer ---
  if (lenrec3 > 1620 ) // atlantis? rec2 is 4096 not 1886 (?)
    {
    //skipmapprojrecord = true;// no: want to read it sometimes...
    WARNING.print("record length longer than expected.  Atlantis?.  reading map.");
    WARNING.print("seems not ok. assuming mapprojection record not present");
    WARNING.print("trying to read further with next record");
    startrec3 = lenrec1+lenrec2+lenrec3;// try to continue with rec4??? rec3 seems platform data
    lenrec3 = 0;// so next time this record is read, but with platform data format...
    }

  // skip map projection record if not present ...
  // BK 18-Jul-2000
  if (skipmapprojrecord) // not present, do not read.
    {
    WARNING.print("Skipping map projection record. Not found.");
    WARNING.print("Continuing with platform position data.");
    }
  else // ESA.SLC format, do read map proj. record.
    {
    // ______ Map projection general information ______
    leaderfile.seekg(startrec3+28,ios::beg);
    leaderfile.read((char*)&c32projection,sizea32);     // map proj. descr.
    c32projection[32]='\0';
    leaderfile.read((char*)&c16numpix,sizei16);         // numpixels
    c16numpix[16]='\0';
    leaderfile.read((char*)&c16numlin,sizei16);         // numlines
    c16numlin[16]='\0';


    // ______ Continue ______
    leaderfile.read((char*)&c16interpix,sizef16);       // dist inter-pixel
    c16interpix[16]='\0';
    leaderfile.read((char*)&c16interlin,sizef16);       // dist inter-lines
    c16interlin[16]='\0';
    leaderfile.read((char*)&c16orien,sizef16);          // orientation at output
    c16orien[16]='\0';
    leaderfile.read((char*)&c16platincl,sizef16);       // actual platform inclination
    c16platincl[16]='\0';
    leaderfile.read((char*)&c16platascn,sizef16);       // actual ascending node
    c16platascn[16]='\0';
    leaderfile.read((char*)&c16geocenter,sizef16);      //
    c16geocenter[16]='\0';
    leaderfile.read((char*)&c16platalt,sizef16);        // altitude
    c16platalt[16]='\0';
    leaderfile.read((char*)&c16platgs,sizef16);         // ground speed
    c16platgs[16]='\0';
    leaderfile.read((char*)&c16plathead,sizef16);       // heading
    c16plathead[16]='\0';
    leaderfile.read((char*)&c32refellips,sizea32);      // ellipsoid
    c32refellips[32]='\0';
    leaderfile.read((char*)&c16refmajor,sizef16);       // semi major
    c16refmajor[16]='\0';
    leaderfile.read((char*)&c16refminor,sizef16);       // semi minor
    c16refminor[16]='\0';

    // ______ Coordinates of four corner points ______
    leaderfile.seekg(startrec3+1072,ios::beg);
    leaderfile.read((char*)&c16lat11,sizef16);          // lat. 1st line 1st pix.
    c16lat11[16]='\0';
    leaderfile.read((char*)&c16lon11,sizef16);          // lon. 1st line 1st pix.
    c16lon11[16]='\0';
    leaderfile.read((char*)&c16lat1N,sizef16);          // lat. 1st line last pix.
    c16lat1N[16]='\0';
    leaderfile.read((char*)&c16lon1N,sizef16);          // lon. 1st line last pix.
    c16lon1N[16]='\0';
    leaderfile.read((char*)&c16latNN,sizef16);          // lat. last line last pix.
    c16latNN[16]='\0';
    leaderfile.read((char*)&c16lonNN,sizef16);          // lon. last line last pix.
    c16lonNN[16]='\0';
    leaderfile.read((char*)&c16latN1,sizef16);          // lat. last line 1st pix.
    c16latN1[16]='\0';
    leaderfile.read((char*)&c16lonN1,sizef16);          // lon. last line 1st pix.
    c16lonN1[16]='\0';
  } // end skip map proj.


  // --- RECORD 4 ---
  DEBUG.print("record 4 of leader file.");
  uint startrec4=lenrec1+lenrec2+lenrec3;
  leaderfile.seekg(startrec4+8,ios::beg);               //  slc platform position data record
  leaderfile.read((char*)&lenrec4,sizeb4);              // length of record4
  lenrec4 = ntohl(lenrec4);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 4: start at: " << startrec4
        << "; length (variable): " << lenrec4;
  DEBUG.print();


// ______ Positional data points ______
  leaderfile.seekg(startrec4+140,ios::beg);
  leaderfile.read((char*)&c4dummy,sizei4);              // number of data points
  c4dummy[4]='\0';
  numdatapoints = atoi(c4dummy);
  leaderfile.read((char*)&c4year,sizei4);               // year
  c4year[4]='\0';
  leaderfile.read((char*)&c4month,sizei4);              // month
  c4month[4]='\0';
  leaderfile.read((char*)&c4day,sizei4);                // day
  c4day[4]='\0';
  leaderfile.read((char*)&c4dayofyear,sizei4);          // day of year
  c4dayofyear[4]='\0';
  leaderfile.read((char*)&c22seconds,sized22);          // sec
  c22seconds[22]='\0';
  leaderfile.read((char*)&c22interval,sized22);         // interval time
  c22interval[22]='\0';
  leaderfile.read((char*)&c64rcs,sizea64);              // ref. coord. system
  c64rcs[64]='\0';
  leaderfile.read((char*)&c22gmha,sized22);             // greenwich mean hour angle
  c22gmha[22]='\0';
  leaderfile.read((char*)&c16ltposerr,sizef16);         // along track pos. error
  c16ltposerr[16]='\0';
  leaderfile.read((char*)&c16ctposerr,sizef16);         // across track pos. error
  c16ctposerr[16]='\0';
  leaderfile.read((char*)&c16rposerr,sizef16);          // radial pos. error
  c16rposerr[16]='\0';

// ______ Read statevector of data points ______
  if (numdatapoints==0) // test, sometimes on linux problems
    {
    WARNING.print("numdatapoints=0: probalbly something wrong with filepointer");
    WARNING.print("But continuing, be very careful!");
    numdatapoints=1;                                    // arbitrary...
    STATE.resize(6,numdatapoints);                      // storage statevector
    }
  else
    {
    STATE.resize(6,numdatapoints);                      // storage statevector
    for (register int32 i=0;i<numdatapoints;i++)        // number of data points
      {
      leaderfile.seekg(startrec4+386+i*132,ios::beg);   // start ith datarecord
      for (register int32 j=0;j<6;j++)                  // read x,y,z,xdot,ydot,zdot
        {
        leaderfile.read((char*)&c22dummy,sized22);
        c22dummy[22]='\0';
        STATE(j,i)=atof(c22dummy);
        }
      }
    }

  // --- RECORD 5 (Bianca Cassee) slc facility related data record [general type]) ---
  DEBUG.print("record 5 of leader file.");              //bc 18 dec 2003
  uint startrec5=lenrec1+lenrec2+lenrec3+lenrec4;       //bc
  leaderfile.seekg(startrec5+8,ios::beg);       //slc facility related
  leaderfile.read((char*)&lenrec5,sizeb4);      //bc length of record4
  lenrec5 = ntohl(lenrec5);                     //byteorder x86 machines.
  DEBUG << "readleader::record 5: start at: " << startrec5
        << "; length: " << lenrec5;
  DEBUG.print();


  // ______ Calibration information  ______              //bc
  leaderfile.seekg(startrec5+582,ios::beg);             //bc
  leaderfile.read((char*)&c16incangle1strange,sizef16); //bc
  c16incangle1strange[16]='\0';                         //bc
  leaderfile.read((char*)&c16incanglecenrange,sizef16); //bc
  c16incanglecenrange[16]='\0';                         // gk bc
  leaderfile.read((char*)&c16incanglelstrange,sizef16); //bc
  c16incanglelstrange[16]='\0';                         //bc
  leaderfile.seekg(startrec5+662,ios::beg);             //gk
  leaderfile.read((char*)&calK,sizef16);                //gk
  calK[16]='\0';                                        //gk
  leaderfile.seekg(startrec5+566,ios::beg);             //gk
  leaderfile.read((char*)&repplspwr,sizef16);           //gk
  repplspwr[16]='\0';                                   //gk
  // ___ BK: number of invalid samples at end (DPAF/UKPAF problem) ___
  char  c4SWSTflag[5];                                  //bk
  char  c4SWSTchange[5];                                //bk
  char  c4missingrawlines[5];                           //bk
  char  c4validperline[5];                              //bk
  leaderfile.seekg(startrec5+98,ios::beg);
  leaderfile.read((char*)&c4SWSTflag,sizei4);           // numsamples
  c4SWSTflag[4]='\0';
  if (atoi(c4SWSTflag) != 0)
    WARNING.print("SWST was not constant");
  leaderfile.seekg(startrec5+138,ios::beg);
  leaderfile.read((char*)&c4SWSTchange,sizei4);         // numsamples
  c4SWSTchange[4]='\0';
  if (atoi(c4SWSTchange) != 0)
    {
    WARNING << "c4SWSTchange: " << c4SWSTchange;
    WARNING.print();
    }
  leaderfile.seekg(startrec5+142,ios::beg);
  leaderfile.read((char*)&c4missingrawlines,sizei4);            // numsamples
  c4missingrawlines[4]='\0';
  if (atoi(c4missingrawlines) != 0)
    {
    WARNING << "c4missingrawlines: " << c4missingrawlines;
    WARNING.print();
    }
  leaderfile.seekg(startrec5+1722,ios::beg);
  leaderfile.read((char*)&c4validperline,sizei4);               // numsamples
  c4validperline[4]='\0';
  DEBUG << "c4validperline: " << c4validperline;
  DEBUG.print();

  /*
  leaderfile.seekg(startrec5+XXX,ios::beg);
  leaderfile.read((char*)&c4numvalid,sizei4);           // numsamples
  c4numvalid[4]='\0';
  leaderfile.read((char*)&c4numinvalid,sizei4);         // zero padded at end
  c4numinvalid[4]='\0';
  if (atoi(c4numinvalid) != 0)
    {
    WARNING << "Number of invalid samples " << c4numinvalid
            << " not 0: rsr may be wrongly computed.";
    WARNING.print();
    }
  */

  //  // --- RECORD 6 (slc facility related data record [pcs type]) ---
  //  DEBUG.print("record 6 of leader file.");                 //gk 28 jan 2004
  //  uint startrec6=lenrec1+lenrec2+lenrec3+lenrec4+lenrec5;  //gk
  //  leaderfile.seekg(startrec6+8,ios::beg);       //slc facility related
  //  leaderfile.read((char*)&lenrec6,sizeb4);      //gk length of record5
  //  lenrec6 = ntohl(lenrec6);                     //byteorder x86 machines.
  //  DEBUG << "readleader::record 6: start at: " << startrec6
  //        << "; length: " << lenrec6;
  //  DEBUG.print();

  }// endif ERS/RSAT
else//RSAT method specified
  {
  skipmapprojrecord = true;// RSAT does not have this here
  // --- RECORD 3 (skip) ---
  DEBUG.print("record 3 of RSAT leader file (data quality).");
  uint startrec3 = lenrec1 + lenrec2;
  leaderfile.seekg(startrec3,ios::beg);//
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 3 with code {10,20,31,20}");
  DEBUG.print("RSAT: Expecting record 3 with code {18,60,18,20}");
  DEBUG.print("RSAT record length should be 1620, ERS 1620");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  leaderfile.seekg(startrec3+8,ios::beg);//
  leaderfile.read((char*)&lenrec3,sizeb4);// length of record3
  lenrec3 = ntohl(lenrec3);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 3: start at: " << startrec3
        << "; length: " << lenrec3;
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==60 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected RSAT record with code {18,60,18,20}");
  else
    WARNING.print("This is NOT the expected 3rd RSAT record with code {18,60,18,20}");


  // --- RECORD 4 (skip) ---
  DEBUG.print("record 4 of RSAT leader file (data hist).");
  uint startrec4 = lenrec1 + lenrec2 + lenrec3;
  leaderfile.seekg(startrec4,ios::beg);//
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 4 with code {10,30,31,20}");
  DEBUG.print("RSAT: Expecting record 4 with code {18,70,18,20}");
  DEBUG.print("RSAT record length should be 16920");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  leaderfile.seekg(startrec4+8,ios::beg);//
  leaderfile.read((char*)&lenrec4,sizeb4);// length of record4
  lenrec4 = ntohl(lenrec4);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 4: start at: " << startrec4
        << "; length: " << lenrec4;
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==70 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected RSAT record with code {18,70,18,20}");
  else
    WARNING.print("This is NOT the expected 4th RSAT record with code {18,70,18,20}");


  // --- RECORD 5 (skip) ---
  DEBUG.print("record 5 of RSAT leader file (data hist).");
  uint startrec5 = lenrec1 + lenrec2 + lenrec3 + lenrec4;
  leaderfile.seekg(startrec5,ios::beg);//
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 5 with code {10,200,31,50}");
  DEBUG.print("RSAT: Expecting record 5 with code {18,70,18,20}");
  DEBUG.print("RSAT record length should be 16920");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  leaderfile.seekg(startrec5+8,ios::beg);//
  leaderfile.read((char*)&lenrec5,sizeb4);// length of record5
  lenrec5 = ntohl(lenrec5);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 5: start at: " << startrec5
        << "; length: " << lenrec5;
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==70 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected RSAT record with code {18,70,18,20}");
  else
    WARNING.print("This is NOT the expected 5th RSAT record with code {18,70,18,20}");


  // --- RECORD 6 (read a bit) ---
  DEBUG.print("record 6 of RSAT leader file (detailed proc. param.).");
  uint startrec6 = lenrec1 + lenrec2 + lenrec3 + lenrec4 + lenrec5;
  leaderfile.seekg(startrec6,ios::beg);//
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS:  Expecting record 6 with code {10,200,31,50}");
  DEBUG.print("RSAT: Expecting record 6 with code {18,120,18,20}");
  DEBUG.print("RSAT record length should be 7726");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  leaderfile.seekg(startrec6+8,ios::beg);//
  leaderfile.read((char*)&lenrec6,sizeb4);// length of record5
  lenrec6 = ntohl(lenrec6);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 6: start at: " << startrec6
        << "; length: " << lenrec6;
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==120 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected RSAT record with code {18,120,18,20}");
  else
    WARNING.print("This is NOT the expected 6th RSAT record with code {18,120,18,20}");

  DEBUG.print("todo: read asc/desc, BEAM, etc.");
  // TODO: Annex B-11: LEA record 6
  // if center lat not known, it seems here too, ignore them


  // --- RECORD 7 (read platform data) ---
  DEBUG.print("record 7 of RSAT leader file (platform pos. data).");
  DEBUG.print("this is exact same as for ERS, but in later record.");
  uint startrec7 = lenrec1 + lenrec2 + lenrec3 + lenrec4 + lenrec5 + lenrec6;
  leaderfile.seekg(startrec7,ios::beg);//
  leaderfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  leaderfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  leaderfile.read((char*)&rec_type,sizeb1);// record type code
  leaderfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  leaderfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("RSAT: Expecting record 7 with code {18,30,18,20}");
  DEBUG.print("RSAT record length should be 8960");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  leaderfile.seekg(startrec7+8,ios::beg);//
  leaderfile.read((char*)&lenrec7,sizeb4);// length of record5
  lenrec7 = ntohl(lenrec7);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG << "readleader::record 7: start at: " << startrec7
        << "; length: " << lenrec7;
  DEBUG.print();
  if (int(rec_sub1)==18 && int(rec_type)==30 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected RSAT record with code {18,30,18,20}");
  else
    WARNING.print("This is NOT the expected 7th RSAT record with code {18,30,18,20}");
  // ______ Keppler elements may be useful to get some sort of orbit data points ______
  // ______ by evaluating them myself, in earth fixed system ______
  char  c16semimajor[17];
  char  c16inclination[17];
  char  c16eccentricity[17];
  char  c16argofperi[17];
  char  c16lonofnode[17];
  char  c16meananomaly[17];
  // Modified by LG for reading JERS Fine
  if ( readfiles_arg.sensor_id == SLC_ALOS )
                startrec7 =     startrec3;
  leaderfile.seekg(startrec7+44,ios::beg);
  leaderfile.read((char*)&c16semimajor,sizef16);
  c16semimajor[16]='\0';
  leaderfile.read((char*)&c16inclination,sizef16);
  c16inclination[16]='\0';
  leaderfile.read((char*)&c16eccentricity,sizef16);
  c16eccentricity[16]='\0';
  leaderfile.read((char*)&c16argofperi,sizef16);
  c16argofperi[16]='\0';
  leaderfile.read((char*)&c16lonofnode,sizef16);
  c16lonofnode[16]='\0';
  leaderfile.read((char*)&c16meananomaly,sizef16);
  c16meananomaly[16]='\0';
  // Modified by LG for reading JERS Fine
  if ( readfiles_arg.sensor_id != SLC_ALOS )
  {
    INFO.print("Orbit Kepplerian elements follow:");
    INFO << "c16semimajor [km]:     " << c16semimajor;
    INFO.print();
    INFO << "c16inclination [rad]:  " << c16inclination;
    INFO.print();
    INFO << "c16eccentricity [-]:   " << c16eccentricity;
    INFO.print();
    INFO << "c16argofperi [rad]:    " << c16argofperi;
    INFO.print();
    INFO << "c16lonofnode [rad]:    " << c16lonofnode;
    INFO.print();
    INFO << "c16meananomaly [rad]:  " << c16meananomaly;
    INFO.print();
  }
  //real8 semimajor = atof(c16semimajor);

  // ______ Positional data points ______
  leaderfile.seekg(startrec7+140,ios::beg);
  leaderfile.read((char*)&c4dummy,sizei4);              // number of data points
  c4dummy[4]='\0';
  numdatapoints = atoi(c4dummy);
  leaderfile.read((char*)&c4year,sizei4);               // year of first data point
  c4year[4]='\0';
  leaderfile.read((char*)&c4month,sizei4);              // month of first data point
  c4month[4]='\0';
  leaderfile.read((char*)&c4day,sizei4);                // day of first data point
  c4day[4]='\0';
  leaderfile.read((char*)&c4dayofyear,sizei4);          // daynumber of year
  c4dayofyear[4]='\0';
  leaderfile.read((char*)&c22seconds,sized22);          // sec of day of first point
  c22seconds[22]='\0';
  leaderfile.read((char*)&c22interval,sized22);         // interval time
  c22interval[22]='\0';
  leaderfile.read((char*)&c64rcs,sizea64);              // ref. coord. system
  c64rcs[64]='\0';
  leaderfile.read((char*)&c22gmha,sized22);             // greenwich mean hour angle
  c22gmha[22]='\0';
  leaderfile.read((char*)&c16ltposerr,sizef16);         // along track pos. error
  c16ltposerr[16]='\0';
  leaderfile.read((char*)&c16ctposerr,sizef16);         // across track pos. error
  c16ctposerr[16]='\0';
  leaderfile.read((char*)&c16rposerr,sizef16);          // radial pos. error
  c16rposerr[16]='\0';

  // --- Maybe Julian Day is useful for conversion ---
  INFO << "Seconds of day of first data point: " << c22seconds;
  INFO.print();
  INFO << "Day of data:        " << c4day;
  INFO.print();
  INFO << "Month of data:      " << c4month;
  INFO.print();
  INFO << "Year of data:       " << c4year;
  INFO.print();
  int32 jd_statevector = julday(atoi(c4day),atoi(c4month),atoi(c4year));
  INFO << "Julian Day of data: " << jd_statevector;
  INFO.print();


  // ______ Read statevector of data points ______
  STATE.resize(6,numdatapoints);                        // storage statevector
  for (register int32 i=0;i<numdatapoints;i++)          // number of data points
    {
    leaderfile.seekg(startrec7+386+i*132,ios::beg);     // start ith datarecord
    for (register int32 j=0;j<6;j++)                    // read x,y,z,xdot,ydot,zdot
      {
      leaderfile.read((char*)&c22dummy,sized22);
      c22dummy[22]='\0';
      STATE(j,i)=atof(c22dummy);
      }
    }
  // --- if this is an inertial system (rsat) we need to convert to earth fixed ---
  // --- Moreover: we want to reduce the number of data points to local arc -------
  WARNING.print("Convert orbit data to earth fixed (please check this).");
  if (!(strcmp(c64rcs,"INERTIAL")))
    INFO.print("Inertial system for orbit: transforming to Earth fixed.");
  else
    WARNING.print("system for orbit: transforming to Earth fixed (?).");
  // This angle refers to ascending node of equator crossing (time), ie. where Z=0
  // not to first state vectors time!  so following is wrong... ?
  // note that we have also the Keppler elements, etc.
  real8 csi2cse  = -1.0;// 1 for initial to earth fixed, -1 for back transform
  real8 GMHA     = deg2rad(atof(c22gmha));// Greenwich Mean Hour Angle in [deg]
  real8 earthrot = 7.292115856E-5;// [rad/s] 2pi rad in 23h56'04.091" (siderial day)
  // can we simply use the given Greenwich mean angle to rotate around Z?
  // or do we need to take precesion.nutation polarwobble into account?
  // to what (time) does the annotated GMHA refer?
  DEBUG << "GMHA [rad]: " << GMHA;
  DEBUG.print();
  DEBUG << "Convertion from inertial to earth fixed [1/-1]: " << csi2cse;
  DEBUG.print();
  DEBUG << "earthrot [rad/s]: " << earthrot;
  DEBUG.print();
  // --- Create a new state vector matrix ---
  // --- these computation could be checked by using TIEPOINT, the center point -----
  // --- Limit to only 2 points before, 2 after scene, and use polyfit(3) interp. ---
  // --- Still this seems to have a large error... ---
  STATE_INERTIAL = STATE;// copy the Z values
  for (register int32 i=0;i<numdatapoints;i++)  // number of data points
    {
    real8 dt       = real8(i)*atof(c22interval);// since annotated GMHA???
    real8 angle    = csi2cse*(GMHA+earthrot*dt);// current angle of Greenwich
    DEBUG << "current angle for this data point [deg]: " << rad2deg(angle);
    DEBUG.print();
    STATE(0,i)     = cos(angle)*STATE_INERTIAL(0,i)-sin(angle)*STATE_INERTIAL(1,i);// x
    STATE(1,i)     = sin(angle)*STATE_INERTIAL(0,i)+cos(angle)*STATE_INERTIAL(1,i);// y
    //STATE(2,i)     = STATE_INERTIAL(2,i);// Z is same
    STATE(3,i)     = -1.0*csi2cse*earthrot*
                     (sin(angle)*STATE_INERTIAL(0,i)+cos(angle)*STATE_INERTIAL(1,i))+
                      cos(angle)*STATE_INERTIAL(3,i)-sin(angle)*STATE_INERTIAL(4,i);// xdot
    STATE(4,i)     = csi2cse*earthrot*
                     (cos(angle)*STATE_INERTIAL(0,i)-sin(angle)*STATE_INERTIAL(1,i))+
                      sin(angle)*STATE_INERTIAL(3,i)+cos(angle)*STATE_INERTIAL(4,i);// ydot
    //STATE(5,i)     = STATE_INERTIAL(5,i);// Zdot is same ?
    }
  INFO.print("info on fast/slow time is obtained in readdat, not here.");
  WARNING.print("todo: reduce the number of orbital data to arc around area only");
  // --- Limit the data points to 2 points before/after: 32 minutes ---
  // --- better yet, interpolate using keppler to 5 points 5 seconds dt over take ---
}

  leaderfile.close();



// ====== Compute prf and rsr and use these (based on data is better) ======
  // Modified by LG for reading ALOS Fine
  if ( (readfiles_arg.sar_processor==SARPR_ATL) || (readfiles_arg.sar_processor == SARPR_JAX ) )
    skipmapprojrecord = true;
  real8 prfcomputed;            // [Hz]  written to result file
  real8 rsrcomputed;            // [MHz] written to result file
  if (skipmapprojrecord==true)  // not present in old? ESA.SLC versions?
    {
    strcpy(c32projection,"skipped");
    strcpy(c16numpix,    "skipped");
    strcpy(c16numlin,    "skipped");
    //strcpy(c16numpix,     c16scenewidth);
    //strcpy(c16numlin,     c16scenelength);
    strcpy(c16interpix,  "skipped");
    strcpy(c16interlin,  "skipped");
    strcpy(c16orien,     "skipped");
    strcpy(c16platincl,  "skipped");
    strcpy(c16platascn,  "skipped");
    strcpy(c16geocenter, "skipped");
    strcpy(c16platalt,   "skipped");
    strcpy(c16platgs,    "skipped");
    strcpy(c16plathead,  "skipped");
    strcpy(c32refellips, "skipped");
    strcpy(c16refmajor,  "skipped");
    strcpy(c16refminor,  "skipped");
    strcpy(c16lat11,     "skipped");
    strcpy(c16lon11,     "skipped");
    strcpy(c16lat1N,     "skipped");
    strcpy(c16lon1N,     "skipped");
    strcpy(c16latNN,     "skipped");
    strcpy(c16lonNN,     "skipped");
    strcpy(c16latN1,     "skipped");
    strcpy(c16lonN1,     "skipped");
    WARNING.print("Using nominal values for PRF, RSR, not computed.");
    WARNING.print(" + because map projection record not present(?).");
    prfcomputed = atof(c16prf);                 // Hz
    // start_added_by_don
    if (readfiles_arg.sensor_id==SLC_ALOS)
      prfcomputed /= 1e3;
    // end_added_by_don
    rsrcomputed = atof(c16samplingrate);        // MHz
    }
  else
    {
    // ______ Check prf ______
    struct tm tijdstart;
    struct tm tijdend;
    strptime(c24zd1stazitime,"%d-%b-%Y %T",&tijdstart);
    strptime(c24zdlstazitime,"%d-%b-%Y %T",&tijdend);
    int32 index=0;
    while (c24zd1stazitime[index] != '.')
      ++index;
    char fracsec1st[25];
    strcpy(fracsec1st,&c24zd1stazitime[index+1]);
    index=0;
    while (c24zdlstazitime[index] != '.')
      ++index;
    char fracseclst[25];
    strcpy(fracseclst,&c24zdlstazitime[index+1]);
    const real8 ta1 = atof(fracsec1st)/1000. +
      real8(tijdstart.tm_sec + 60*tijdstart.tm_min + 3600*tijdstart.tm_hour);
    const real8 taN = atof(fracseclst)/1000. +
      real8(tijdend.tm_sec + 60*tijdend.tm_min + 3600*tijdend.tm_hour);
    // BK 28-Sep-2000: numlin -1 !
    prfcomputed = (atof(c16numlin) - 1.) / (taN-ta1);

    // ______ compute rsr ______
    //  const real8 tr1     = atof(c16zd1strange);              // 2way ms
    //  const real8 trN     = atof(c16zdlstrange);              // 2way ms
    // BK 28-Sep-2000: numpix -1 !
    rsrcomputed = 0.001 * (atof(c16numpix)-1.0) /
                  (atof(c16zdlstrange)-atof(c16zd1strange));    // MHz
    } // else skipped map projection record

  // BK 28-Oct-2003, for Atlantis processor
  // ___ Check if rsr is in MHz, assume about 20 MHz ---
  //if (rsrcomputed>10000000.0 && rsrcomputed<30000000.0)
  if (rsrcomputed>10000.0)
    {
    WARNING.print("Adapting RSR from Hz to MHz.");
    WARNING << "Old value rsrcomputed: " << rsrcomputed
            << "; new value rsrcomputed: " << rsrcomputed/1e6;
    WARNING.print();
    rsrcomputed/=1e6;
    }
  real8 rangebandwidth = atof(c16bandrangetot);
  // start_added_by_don
  if (readfiles_arg.sensor_id==SLC_ALOS)
    rangebandwidth /= 1e3;
  // end_added_by_don
  if (rangebandwidth < 1.0)// e.g. JERS not annotated.
    {
    WARNING.print("computed rangebandwidth < 1.0 ignoring this.");
    rangebandwidth = rsrcomputed;// set this to no oversampling
    }
  //if (rangebandwidth>10000000.0 && rangebandwidth<30000000.0)
  if (rangebandwidth>10000.0)
    {
    WARNING.print("Adapting RBW from Hz to MHz.");
    WARNING << "Old value rangebandwidth: " << rangebandwidth
            << "; new value rangebandwidth: " << rangebandwidth/1e6;
    WARNING.print();
    rangebandwidth/=1e6;
    }


// ====== Write information to scratchfiles ======
  ofstream scratchlogfile("scratchloglea", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"readleader: scratchloglea",__FILE__,__LINE__);

  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* EXTRACTED DATA FROM LEADER FILE: "
    <<  readfiles_arg.leaderfile << " *"
    << "\n*******************************************************************"

    << "\n\nFile descriptor record"
    << "\n----------------------"
    << "\nFile name: "
    <<  c16leafilename

    << "\n\nSLC data set summary record: scene parameters"
    << "\n---------------------------------------------"
    << "\nScene reference number: \t\t\t\t" // (e.g. orbit-frame number): "
    <<  c32sceneref
    << "\nScene centre time (UTC) <YYYYMMDDhhmmssttt>: \t\t"
    <<  c32scenetime
    << "\nProcessed scene centre geodetic latitude"
    << "\n +(positive for North latitude) (degrees): \t\t"
    <<  c16centerlat
    << "\nProcessed scene centre longitude"
    << "\n +(negative for West longitude) (degrees): \t\t"
    <<  c16centerlon
    << "\nProcessed scene centre true heading (degrees): \t\t"
    <<  c16centerheading
    << "\nEllipsoid designator: \t\t\t\t\t"
    <<  c16ellipsoid
    << "\nEllipsoid semimajor axis (km): \t\t\t\t"
    <<  c16semimajor
    << "\nEllipsoid semiminor axis (km): \t\t\t\t"
    <<  c16semiminor
    << "\nEarth mass times grav. const. (M.G) (kg.m/s^2): \t"
    <<  c16GM
    << "\nEllipsoid J2 parameter: \t\t\t\t"
    <<  c16J2
    << "\nEllipsoid J3 parameter: \t\t\t\t"
    <<  c16J3
    << "\nEllipsoid J4 parameter: \t\t\t\t"
    <<  c16J4
    << "\nScene center line number (the line number at "
    << "\n +the scene centre including zero fill): \t\t"
    <<  scenecenterline
    << "\nScene center pixel number (the pixel number at"
    << "\n +the scene centre including zero fill): \t\t"
    <<  scenecenterpixel
    << "\nProcessed scene length incl. zero fill (km): \t\t"
    <<  c16scenelength
    << "\nProcessed scene width incl. zero fill (km): \t\t"
    <<  c16scenewidth

    << "\n\nSLC data set summary record: general mission/sensor parameters"
    << "\n--------------------------------------------------------------"
    << "\nSensor platform mission identifier: \t\t\t"
    <<  c16platformid
    << "\nSensor ID, mode of op. for this channel: \t\t"
    <<  c32sensorid
    << "\nOrbit number: \t\t\t\t\t\t\t"
    <<  c8orbitnr
    << "\nSensor platform geodetic latitude at nadir"
    << "\n +corresponding to scene centre"
    << "\n +(positive for North latitude) (degrees): \t\t\t"
    <<  c8platformlat
    << "\nSensor platform geodetic longitude at nadir"
    << "\n +corresponding to scene centre"
    << "\n +(negative for West longitude) (degrees): \t\t\t"
    <<  c8platformlon
    << "\nSensor platform heading at nadir"
    << "\n +corresponding to scene centre"
    << "\n +(clockwise positive from North) (degrees): \t\t\t"
    <<  c8platformheading
    << "\nSensor clock angle as measured relative to"
    << "\n +sensor platform flight direction (degrees): \t"
    <<  c8clockangle
    << "\nIncidence angle at scene centre as derived"
    << "\n +from sensor platform orientation (degrees): \t"
    <<  c8incidence
    << "\nRadar frequency (GHz): \t\t\t\t\t\t"
    <<  c8freq
    << "\nRadar wavelength (meters): \t\t\t"
    <<  c16wavelength
    << "\nMotion compensator identifier"
    << "\n +(00=no, 01=on board, 10=in processor, 11=both): \t"
    <<  c2motioncomp
    << "\nRange pulse specifier: \t\t\t\t\t"
    <<  c16pulse
    << "\nNominal range pulse (chirp) amplitude coefficient,"
    << "\n +constant term: \t\t\t\t\t"
    <<  c16ampconst
    << "\n +linear term (sec-1): \t\t\t\t\t"
    <<  c16amplinear
    << "\n +quadratic term (sec-2): \t\t\t\t"
    <<  c16ampquadratic
    << "\n +cubic term (sec-3): \t\t\t\t\t"
    <<  c16ampcubic
    << "\n +quartic term (sec-4): \t\t\t\t"
    <<  c16ampquartic
    << "\nNominal range pulse (chirp) phase coefficient,"
    << "\n +constant term (cycles): \t\t\t\t"
    <<  c16phaseconst
    << "\n +linear term (Hz): \t\t\t\t\t"
    <<  c16phaselinear
    << "\n +quadratic term (Hz/sec): \t\t\t\t"
    <<  c16phasequadratic
    << "\n +cubic term (Hz/sec2): \t\t\t\t"
    <<  c16phasecubic
    << "\n +quartic term (Hz/sec3): \t\t\t\t"
    <<  c16phasequartic
    << "\nDown linked chirp extraction index (samples): \t\t\t"
    <<  c8extindex
    << "\nRange sampling rate (MHz): \t\t\t\t"
    <<  c16samplingrate
    << "\nRange gate delay at early ege (in time)"
    << "\n +at the start of the image (microsec): \t"
    <<  c16rangedelay
    << "\nRange pulse length (microsec): \t\t\t\t"
    <<  c16ranpulselen
    << "\nBase band conversion flag: \t\t\t"
    <<  c4conversion
    << "\nRange compressed flag (YES=range compressed data): \t"
    <<  c4compression
    << "\nQuantization per channel I & Q"
    << "\n +(5I 5Q/6I 6Q for OGRB/OBRC) (bits): \t\t\t\t"
    <<  c8qperch
    << "\nQuantizer descriptor: \t\t\t\t\t"
    <<  c12qdesc
    << "\nDC Bias for I-component (actual value): \t\t"
    <<  c16dci
    << "\nDC Bias for Q-component (actual value): \t\t"
    <<  c16dcq
    << "\nGain imbalance for I & Q (actual value): \t\t"
    <<  c16imbalance
    << "\nAntenna mechanical boresight angle"
    << "\n +relative to platform vertical axis: \t\t\t"
    <<  c16boresight
    << "\nPulse Repetition Frequency (PRF) (actual value): \t"
    <<  c16prf

    << "\n\nSLC data set summary record: sensor specific parameters"
    << "\n-------------------------------------------------------"

    << "\nSatellite encoded binary time code: \t\t\t"
    <<  c16sattimecode
    << "\nSatellite clock time (UTC) (YYYYMMDDhhmmssttt): \t"
    <<  c32sattime
    << "\nSatellite clock step length (nanosec): \t\t\t\t"
    <<  c8satclockstep

    << "\n\nSLC data set summary record: general processing parameters"
    << "\n----------------------------------------------------------"
    << "\nProcessing facility identifier: \t\t\t"
    <<  c16facilityid
    << "\nProcessing system identifier: \t\t\t\t"
    <<  c8systemid
    << "\nProcessing version identifier: \t\t\t\t"
    <<  c8versionid
    << "\nProduct type specifier: \t\t\t\t"
    <<  c32typespec
    << "\nProcessing algorithm identifier: \t\t\t"
    <<  c32algid
    << "\nNominal number of looks processed in azimuth (looks): \t"
    <<  c16looksazi
    << "\nNominal number of looks processed in range (looks): \t"
    <<  c16looksrange
    << "\nBandwidth per look in azimuth (null-to-null) (Hz): \t"
    <<  c16bandazi
    << "\nBandwidth per look in range (MHz): \t\t\t"
    <<  c16bandrange
    << "\nTotal processor bandwidth in azimuth (Hz): \t\t"
    <<  c16bandazitot
    << "\nTotal processor bandwidth in range (MHz): \t\t"
    <<  c16bandrangetot
    << "\nWeighting function designator in azimuth: \t\t"
    <<  c32weightazi
    << "\nWeighting function designator in range: \t\t"
    <<  c32weightrange
    << "\nData input source: \t\t\t\t\t"
    <<  c16inputsource
    << "\nNominal resolution in range (3-dB width) (m): \t\t"
    <<  c16resrange
    << "\nNominal resolution in azimuth (3-dB width) (m): \t"
    <<  c16resazi
    << "\nAlong track Doppler frequency centroid"
    << "\n +at early edge of image"
    << "\n +constant term (Hz): \t\t\t\t\t"
    <<  c16atdoppcconst
    << "\n +linear term (Hz/sec): \t\t\t\t"
    <<  c16atdoppclinear
    << "\n +quadratic term (Hz/sec/sec): \t\t\t\t"
    <<  c16atdoppcquadratic
    << "\nCross track Doppler frequency centroid"
    << "\n +at early edge of image"
    << "\n +constant term (Doppler centroid) (Hz): \t\t"
    <<  c16xtdoppcconst
    << "\n +linear term (Slope of Doppler centroid) (Hz/sec): \t"
    <<  c16xtdoppclinear
    << "\n +quadratic term (Hz/sec/sec): \t\t\t\t"
    <<  c16xtdoppcquadratic
    << "\nTime direction indicator along pixel direction: \t"
    <<  c8timepix
    << "\nTime direction indicator along line direction: \t\t"
    <<  c8timeline
    << "\nAlong track Doppler frequency rate"
    << "\n +at early edge of image"
    << "\n +constant term (Hz/sec): \t\t\t\t"
    <<  c16atdopprconst
    << "\n +linear term (Hz/sec/sec): \t\t\t\t"
    <<  c16atdopprlinear
    << "\n +quadratic term (Hz/sec/sec/sec): \t\t\t"
    <<  c16atdopprquadratic
    << "\nCross track Doppler frequency rate"
    << "\n +at early edge of image"
    << "\n +constant term (Azimuth FM rate) (Hz/sec): \t\t"
    <<  c16xtdopprconst
    << "\n +linear term (Slope of Azimuth FM rate) (Hz/sec/sec): \t"
    <<  c16xtdopprlinear
    << "\n +quadratic term (Hz/sec/sec/sec): \t\t\t"
    <<  c16xtdopprquadratic
    << "\nLine content indicator: \t\t\t\t"
    <<  c8linecontent
    << "\nClutterlock applied flag: \t\t\t\t"
    <<  c4clutterlock
    << "\nAutofocussing applied flag: \t\t\t\t"
    <<  c4autofocus
    << "\nLine spacing (m): \t\t\t\t\t"
    <<  c16linespace
    << "\nPixel spacing (in range) (m): \t\t\t\t"
    <<  c16pixspace
    << "\nProcessor range compression designator: \t\t"
    <<  c16rcompdes

    << "\n\nSLC data set summary record: sensor specific local use segment"
    << "\n--------------------------------------------------------------"
    << "\nZero-doppler range time (two-way)"
    << "\n +of first range pixel (millisec): \t\t\t"
    <<  c16zd1strange
    << "\n +of centre range pixel (millisec): \t\t\t"
    <<  c16zdcenrange
    << "\n +of last range pixel (millisec): \t\t\t"
    <<  c16zdlstrange
    << "\nZero-doppler azimuth time"
    << "\n +of first azimuth pixel (UTC): \t\t"
    <<  c24zd1stazitime
    << "\n +of centre azimuth pixel (UTC): \t\t"
    <<  c24zdcenazitime
    << "\n +of last azimuth pixel (UTC): \t\t\t"
    <<  c24zdlstazitime

    << "\n\nMap projection data record: general information"
    << "\n-----------------------------------------------"
    << "\nMap projection descriptor: \t\t\t\t"
    <<  c32projection
    << "\nNumber of pixels per line of image: \t\t\t"
    <<  c16numpix
    << "\nNumber of lines: \t\t\t\t\t"
    <<  c16numlin
    << "\nNominal inter-pixel distance in output scene (m): \t"
    <<  c16interpix
    << "\nNominal inter-line distance in output scene (m): \t"
    <<  c16interlin
    << "\nOrientation at output scene centre"
    << "\n +[for geodetic products this is simply the"
    << "\n +convergence of the meridians, ie: the angle"
    << "\n +between geographic north and map grid north"
    << "\n +(angle of projection axis from true North)] (degrees): "
    <<  c16orien
    << "\nActual platform orbital inclination (degrees): \t\t"
    <<  c16platincl
    << "\nActual ascending node (longitude at equator) (degrees): "
    <<  c16platascn
    << "\nGeocentre to platform distance at input scene centre: \t"
    <<  c16geocenter
    << "\nPlatform geodetic altitude over the ellipsoid: \t\t"
    <<  c16platalt
    << "\nGround speed at nadir at input scene centre time: \t"
    <<  c16platgs
    << "\nPlatform heading at nadir"
    << "\n +corresponding to scene centre (degrees): \t\t"
    <<  c16plathead
    << "\nName of reference ellipsoid: \t\t\t\t"
    <<  c32refellips
    << "\nSemimajor axis of ref.ellipsoid (km): \t\t\t"
    <<  c16refmajor
    << "\nSemiminor axis of ref.ellipsoid (km): \t\t\t"
    <<  c16refminor

    << "\n\nMap projection data record: coordinates of four corner points"
    << "\n-------------------------------------------------------------"
    << "\n1st line 1st pixel geodetic latitude"
    << "\n +(positive for North latitude) (degrees): \t\t"
    <<  c16lat11
    << "\n1st line 1st pixel geodetic longitude"
    << "\n +(negative for West longitude) (degrees): \t\t"
    <<  c16lon11
    << "\n1st line last pixel geodetic latitude (degrees): \t"
    <<  c16lat1N
    << "\n1st line last pixel geodetic longitude (degrees): \t"
    <<  c16lon1N
    << "\nlast line last pixel geodetic latitude (degrees): \t"
    <<  c16latNN
    << "\nlast line last pixel geodetic longitude (degrees): \t"
    <<  c16lonNN
    << "\nlast line 1st pixel geodetic latitude (degrees): \t"
    <<  c16latN1
    << "\nlast line 1st pixel geodetic longitude (degrees): \t"
    <<  c16lonN1

    << "\n\nSLC platform position data record: positional data points"
    << "\n---------------------------------------------------------"
    << "\nNumber of data points: \t\t\t\t"
    <<  numdatapoints
    << "\nYear of data point <YYYY>: \t\t\t"
    <<  c4year
    << "\nMonth of data point <$$MM>: \t\t\t"
    <<  c4month
    << "\nDay of data point <$$DD>: \t\t\t"
    <<  c4day
    << "\nDay in the year <GMT> (jan 1st=day 1): \t\t"
    <<  c4dayofyear
    << "\nSeconds of day of data: \t\t\t"
    <<  c22seconds
    << "\nTime interval between data points: \t\t"
    <<  c22interval
    << "\nReference coordinate system: \t\t\t"
    <<  c64rcs
    << "\nGreenwich mean hour angle (degrees): \t\t"
    <<  c22gmha
    << "\nAlong track position error (meters): \t\t"
    <<  c16ltposerr
    << "\nAcross track position error (meters): \t\t"
    <<  c16ctposerr
    << "\nRadial position error (meters): \t\t"
    <<  c16rposerr

    << "\n\nSLC facility related data record [general type]: calibration information"
    << "\n------------------------------------------------------------------------"
    << "\nIncidence angle at first range pixel (at mid-azimuth): \t"
    <<  c16incangle1strange //bc
    << "\nIncidence angle at centre range pixel (at mid-azimuth): "
    <<  c16incanglecenrange //bc
    << "\nIncidence angle at last range pixel (at mid-azimuth): \t"
    <<  c16incanglelstrange //bc
    << "\nAbsolute calibration constant K: \t"
    <<  calK     //gk
    << "\nLenrec6: \t"
    <<  lenrec6  //gk
    << "\nReplica Pulse Power: \t"
    <<  repplspwr//gk
    << "\nNumber of valid samples:   \t\t\t" << c4numvalid
    << "\nNumber of invalid samples: \t\t\t" << c4numinvalid
    <<  endl;
  int32 point;

  // --- STATE VECTORS in INERTIAL SYSTEM ---
  // only convert these if RSAT and ATLANTIS?
//if (readfiles_arg.sensor_id==SLC_RSAT)
if (readfiles_arg.sensor_id==SLC_RSAT && readfiles_arg.sar_processor==SARPR_ATL)
  {
  scratchlogfile << "\nRSAT: INERTIAL system orbital data as in leaderfile";
  for (register int32 k=0;k<numdatapoints;k++)          // number of data points
    {
    point = k+1;
    real8 secofday = atof(c22seconds) + real8(k)*atof(c22interval);
  scratchlogfile << "\nSLC platform position data record: data point: " << point
                 << "\n------------------------------------------------\n"
                 << point << " data point - Seconds of day (s):    \t\t"
                 << setprecision(13)
                 << secofday << endl
                 << point << " data point - Position vector X (m): \t\t"
                 << setprecision(13)
                 << STATE_INERTIAL(0,k) << endl
                 << point << " data point - Position vector Y (m): \t\t"
                 << setprecision(13)
                 << STATE_INERTIAL(1,k) << endl
                 << point << " data point - Position vector Z (m): \t\t"
                 << setprecision(13)
                 << STATE_INERTIAL(2,k) << endl
                 << point << " data point - Velocity vector X (mm/s): \t"
                 << setprecision(13)
                 << STATE_INERTIAL(3,k) << endl
                 << point << " data point - Velocity vector Y (mm/s): \t"
                 << setprecision(13)
                 << STATE_INERTIAL(4,k) << endl
                 << point << " data point - Velocity vector Z (mm/s): \t"
                 << setprecision(13)
                 << STATE_INERTIAL(5,k) << endl;
    }
  scratchlogfile << "\nRSAT: and converted to earth fixed using GMHA:";
  }

  // --- STATE VECTORS in EARTH FIXED SYSTEM ---
  for (register int32 k=0;k<numdatapoints;k++)          // number of data points
    {
    point = k+1;
    real8 secofday = atof(c22seconds) + real8(k)*atof(c22interval);
  scratchlogfile << "\nSLC platform position data record: data point: " << point
                 << "\n------------------------------------------------\n"
                 << point << " data point - Seconds of day (s):    \t\t"
                 << setprecision(13)
                 << secofday << endl
                 << point << " data point - Position vector X (m): \t\t"
                 << setprecision(13)
                 << STATE(0,k) << endl
                 << point << " data point - Position vector Y (m): \t\t"
                 << setprecision(13)
                 << STATE(1,k) << endl
                 << point << " data point - Position vector Z (m): \t\t"
                 << setprecision(13)
                 << STATE(2,k) << endl
                 << point << " data point - Velocity vector X (mm/s): \t"
                 << setprecision(13)
                 << STATE(3,k) << endl
                 << point << " data point - Velocity vector Y (mm/s): \t"
                 << setprecision(13)
                 << STATE(4,k) << endl
                 << point << " data point - Velocity vector Z (mm/s): \t"
                 << setprecision(13)
                 << STATE(5,k) << endl;
    }

  scratchlogfile << "\nEND LEADER FILE "
                 <<  readfiles_arg.leaderfile
                 <<  endl << endl;
  scratchlogfile.close();



  // ___ "repair" bandwidth if not given, e.g., in JERS? ___
  if (atof(c16bandazitot) < 1.0)
    {
    WARNING.print("could not find azimuth band width, using prf...");
    strcpy(c16bandazitot,c16prf);
    }
  if (atof(c16bandrangetot) < 1.0)
    {
    WARNING.print("could not find range band width, using rsr...");
    strcpy(c16bandrangetot,c16samplingrate);
    }

  real8 wavelength_computed = atof(c16wavelength);// seems freq. not ok for RSAT
  // --- seems not all information there for Atlantis CEOS ---
  //if (readfiles_arg.sensor_id!=SLC_RSAT)
  if (readfiles_arg.sar_processor!=SARPR_ATL)
    {
    INFO << "Computing wavelength from radar freq. " << c8freq << "GHz";// GHz
    INFO.print();
    wavelength_computed = (0.000000001*SOL/atof(c8freq));// seems more reliable, BK 03/04
    INFO << "Annotated wavelength: " << c16wavelength << " re-computed: " << wavelength_computed;
    INFO.print();
    }

  // Email KKMohanty: time indicator is different then line counter, i.e, image
  // can be flipped up/down left/right.  The utilities converting line number to
  // az_time can thus be confused.  Repair this.
  //  Time direction indicator along pixel direction:         INCREASE or DECREASE
  //  Time direction indicator along line direction:          INCREASE or DECREASE
  // BK 20-Apr-2004
  INFO.print("Checking time indicator along line direction.");
  // ___ check indicator in line direction ___
  if (!strcmp(c8timeline,"INCREASE"))
    {
    INFO.print("Time indicator along line direction is INCREASE, this is OK.");
    }
  else if (!strcmp(c8timeline,"DECREASE"))
    {
    WARNING.print("Time indicator along line direction is DECREASE: should adapt times.");
    WARNING.print("TODO: done for RSAT, not for ERS.");
    // ___ Adapt the time of first line/last line, such that orbit is correctly ___
    // ___ interpolated and line2ta works OK ___
    }
  else
    {
    WARNING.print("unexpected: time indicator along line direction not INCREASE nor DECREASE.");
    }
  // ___ check indicator in pixel direction ___
  if (!strcmp(c8timepix,"INCREASE"))
    {
    INFO.print("Time indicator along pixel direction is INCREASE, this is OK.");
    }
  else if (!strcmp(c8timepix,"DECREASE"))
    {
    WARNING.print("Time indicator along pixel direction is DECREASE: should adapt times.");
    WARNING.print("TODO: done for RSAT, not for ERS.");
    }
  else
    {
    WARNING.print("unexpected: time indicator along pixel direction not INCREASE nor DECREASE.");
    }



  ofstream scratchresfile("scratchreslea", ios::out | ios::trunc);
  bk_assert(scratchresfile,"readleader: scratchreslea",__FILE__,__LINE__);
  // Modified by LG for reading ALOS Fine
  if(!strcmp(c16centerlat,"                "))
          strcpy(c16centerlat,"0");
  if(!strcmp(c16centerlon,"                "))
          strcpy(c16centerlon,"0");

  scratchresfile
    << "Leader file:                                 \t"
    <<  readfiles_arg.leaderfile
    << "\nSensor platform mission identifer:         \t"
    <<  c32sensorid
    << "\nScene_centre_latitude:                     \t"
    <<  c16centerlat
    << "\nScene_centre_longitude:                    \t"
    <<  c16centerlon
    //<< "\nScene_centre_heading:                    \t";   // MA
    << "\nScene_centre_heading";   // MA
    if (readfiles_arg.sar_processor==SARPR_ATL)
      {
      scratchresfile
      << "_true:                    \t"
      << c16centerheading;                                 // MA for VMP true heading  + radarsat
      }
    else
      {
      scratchresfile
       << ":                    \t"   // MA
      << c16plathead;                                      // heading more digits no true, PGS doesn't have true heading but VMP does. TODO put switch for VMP and PGS
      }
    //<< c8platformheading;                               // heading
    //<< c16centerheading;                                 // MA for VMP true heading
  // start_added_by_don
  if (readfiles_arg.sensor_id == SLC_ALOS)
    {
    scratchresfile
      << "\nRadar_wavelength (m):                      \t"
      <<  c16wavelength;
    }
  else
  // end_added_by_don
    {
    scratchresfile
    << "\nRadar_wavelength (m):                      \t"
    << setprecision(16) // important to have enough digits
    <<  wavelength_computed;
    // <<  c16wavelength;
    }

    // ______ Azimuth ______
  //if (readfiles_arg.sar_processor==SARPR_ATL)
  // start_added_by_don
  if ( (readfiles_arg.sensor_id == SLC_ALOS) || (readfiles_arg.sar_processor==SARPR_ATL))
  // end_added_by_don
    {
    WARNING.print("First_pixel_azimuth_time not in leader, but done in readdat.");
    scratchresfile
      << "\nRSAT_First_pixel_azimuth_time (UTC):            \t"
      <<  "skipped_see_datfile_below";
      //<<  c24zd1stazitime;
    }
  else
    {
    scratchresfile
      << "\nFirst_pixel_azimuth_time (UTC):            \t"
      <<  c24zd1stazitime;
    }
  scratchresfile
    // email to ESA, they said this value is most accurate and should be used.
    //<< "\nPulse_Repetition_Frequency (actual, Hz):   \t"
    //<<  c16prf
    // (but that seems strange to me...)
    << "\nPulse_Repetition_Frequency (computed, Hz): \t"
    << setprecision(16)
    <<  prfcomputed
    << "\nTotal_azimuth_band_width (Hz):             \t"
    <<  c16bandazitot
    << "\nWeighting_azimuth:                         \t"
    <<  c32weightazi
    << "\nXtrack_f_DC_constant (Hz, early edge):     \t"
    <<  c16xtdoppcconst;

  // start_added_by_don
  if ( (readfiles_arg.sar_processor==SARPR_JAX) || (readfiles_arg.sar_processor==SARPR_ATL))
    {
    real8 xtdoppclinear = atof(c16xtdoppclinear) * rsrcomputed * 1e6;
    real8 xtdoppcquadratic = atof(c16xtdoppcquadratic) * sqr(rsrcomputed*1e6);
    scratchresfile
      << "\nXtrack_f_DC_linear (Hz/s, early edge):     \t"
      << setprecision(16)
      <<  xtdoppclinear
      << "\nXtrack_f_DC_quadratic (Hz/s/s, early edge): \t"
      <<  xtdoppcquadratic;
    }
  else
  // start_added_by_don
    {
    scratchresfile
    << "\nXtrack_f_DC_linear (Hz/s, early edge):     \t"
    <<  c16xtdoppclinear
    << "\nXtrack_f_DC_quadratic (Hz/s/s, early edge): \t"
    <<  c16xtdoppcquadratic;
  }

  // ______ Range ______
  //if (readfiles_arg.sar_processor==SARPR_ATL)
  // Modified by LG for reading ALOS Fine
  if ( (readfiles_arg.sensor_id == SLC_ALOS) || (readfiles_arg.sar_processor==SARPR_ATL))
    {
    WARNING.print("Range_time_to_first_pixel not in leader, but done in readdat.");
    scratchresfile
      << "\nRSAT_Range_time_to_first_pixel (2way) (ms):     \t"
      <<  "skipped_see_datfile_below";
      //<<  c16zd1strange;
    }
  else
    {
    real8 range_t1 = atof(c16zd1strange);
    // Bert Kampes, 20-Sep-2005: analysis of Corner Reflectors showed about 4
    // samples offset: systematic error, correct for ERS1 and ERS2.
    // this means that the radarcoding is now better.
    // If you don;t like this, you can use card [MS]_T_RG_ERROR to correct
    // but it is not that important.
    //
    // No, we did not like this and turned it off. FvL, 13-6-2006
    //
    //real8 range_t1 = atof(c16zd1strange);
    //if (readfiles_arg.sar_processor==SARPR_VMP)
    //  {
    //  INFO.print("ERS range time correction: -0.00019235 msec (two-way); ");
    //  range_t1 += -0.00019235;// ERS CR computed by BK ~75m
    //  }
    scratchresfile
      << "\nRange_time_to_first_pixel (2way) (ms):     \t" //<<  c16zd1strange;
      << setprecision(16)
      <<  range_t1;
    }
  scratchresfile
    << "\nRange_sampling_rate (computed, MHz):       \t"
    << setprecision(16)
    <<  rsrcomputed
    << "\nTotal_range_band_width (MHz):               \t"
    <<  rangebandwidth
    << "\nWeighting_range:                            \t"
    <<  c32weightrange
    << endl;

  scratchresfile.setf(ios::fixed | ios::floatfield);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_leader_datapoints"
    << "\n*******************************************************************"
    << "\nt(s)\t\tX(m)\t\tY(m)\t\tZ(m)"
    << "\nNUMBER_OF_DATAPOINTS: \t\t\t"
    <<  numdatapoints
    << endl;

  for (register int32 l=0;l<numdatapoints;l++)          // number of data points
    {
    scratchresfile
      << setprecision(11)                                // [MA] 10 --> 11 since we have iiiii.dddddd decimal digits for day of seconds
      << atof(c22seconds)+real4(l)*atof(c22interval);
    for (register int32 m=0;m<3;m++)                    // no velocities
      {
      scratchresfile
        << " \t"
        << setprecision(13)
        << STATE(m,l);
      }
    scratchresfile << endl;
    }
  scratchresfile
    << "\n*******************************************************************"
    << "\n* End_leader_datapoints:_NORMAL"      // fixed string...
    << "\n*******************************************************************\n";
  scratchresfile.close();


// ______ Check with volumefile (but first write what is read) ______
  if (atoi(c16numlin) != checklines)
    {
    WARNING << "Number of lines of leader file seems "
         << "not to be consistent with volume file: "
         <<  c16numlin << " != " << checklines;
    WARNING.print();
    WARNING.print("This may be caused by wrong? format SLC of this PAF, check parameterpool.");
    }
  if (100*((atof(c16prf)-prfcomputed)/prfcomputed) > 2.)
    WARNING.print("deviation PRF computed from PRF read in leaderfile > 2%");
  if (100*((atof(c16samplingrate)-rsrcomputed)/rsrcomputed) > 2.)
    WARNING.print("deviation range sampling rate computed from RSR read in leaderfile > 2%");


// ______ Write some info ______
  INFO << "Number of lines, pixels:               "
       << c16numlin << " " << c16numpix;
  INFO.print();
  INFO << "Pulse repetition frequency (Hz):       "
       << c16prf << ends;
  INFO.print();
  INFO << "Pulse repetition frequency (computed): "
       << setprecision(16) << prfcomputed;
  INFO.print();
  INFO << "Range sampling rate (Mhz):             "
       << c16samplingrate << ends;
  INFO.print();
  INFO << "Range sampling rate (computed Mhz):    "
       << setprecision(16) << rsrcomputed;
  INFO.print();
  INFO << "UTC of first azimuth line:             "
       << c24zd1stazitime;
  INFO.print();
  INFO << "UTC of last azimuth line:              "
       << c24zdlstazitime;
  INFO.print();
  INFO << "Range time to first pixel (2way ms):   "
       << c16zd1strange;
  INFO.print();
  INFO << "Range time to last pixel (2way ms):    "
       << c16zdlstrange;
  INFO.print();
  INFO << "First corner of image (lat,lon):       "
       << c16lat11 << " " << c16lon11;
  INFO.print();
  INFO << "Second corner of image (lat,lon):      "
       << c16lat1N << " " << c16lon1N;
  INFO.print();
  INFO << "Third corner of image (lat,lon):       "
       << c16latNN << " " << c16lonNN;
  INFO.print();
  INFO << "Fourth corner of image (lat,lon):      "
       << c16latN1 << " " << c16lonN1;
  INFO.print();

  INFO << "Weighting function designator azimuth: "
       << c32weightazi;
  INFO.print();
  toupper(c32weightazi);
  if (strncmp(c32weightazi,"HAMMING",7))
    {
    WARNING << INFO.get_str() << " not HAMMING.";
    WARNING.print();
    }
  INFO << "Weighting function designator range:   "
       << c32weightrange;
  INFO.print();
  toupper(c32weightrange);
  if (strncmp(c32weightrange,"HAMMING",7))
    {
    WARNING << INFO.get_str() << " not HAMMING.";
    WARNING.print();
    }

  PROGRESS.print("readleader finished.");
  } // END READLEADER



/****************************************************************
 *    readnull                                                  *
 *                                                              *
 * reads nullfile                                               *
 *  and writes to scratchlogfile, scratchresfile                *
 * see: annex C ERS SAR.SLC CCTand EXABYTE                      *
 *      doc:er-is-epo-gs-5902.3                                 *
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void readnull(
        const input_readfiles &readfiles_arg)
  {
// ======Open files======
  ifstream nullfile;
  openfstream(nullfile,readfiles_arg.nullfile);
  bk_assert(nullfile,readfiles_arg.nullfile,__FILE__,__LINE__);

  ofstream scratchlogfile("scratchlognull", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"readnull: scratchlognull",__FILE__,__LINE__);

  ofstream scratchresfile("scratchresnull", ios::out | ios::trunc);
  bk_assert(scratchresfile,"readnull: scratchresnull",__FILE__,__LINE__);

// ======Read nullfile======
  DEBUG.print("nullfile should be read?");

// ______Tidy up______
  nullfile.close();
  scratchlogfile.close();
  scratchresfile.close();
  PROGRESS.print("readnull finished.");
  } // END READNUL



/****************************************************************
 *    readdat                                                   *
 *                                                              *
 * Extract info from data file,                                 *
 *  see also appendix C of CD-R distribution (ESA).             *
 * checks with volumefile #lines                                *
 * checks some things like numchannels                          *
 *                                                              *
 * info is written to scratchresfile?                           *
 *                                                              *
 *    Bert Kampes, 21-Dec-1998                                  *
 ****************************************************************/
void readdat(
        input_readfiles &readfiles_arg,
        const int32 checklines)
  {
  TRACE_FUNCTION("readdat (BK 21-Dec-1998)")
  const int16           sizeb4 = 4,             // some constants for reading
                        sizeb1 = 1,
                        sizei4 = 4,             //   binary fields.
                        sizei6 = 6,
                        sizei8 = 8;
  uint                  numchannels,            // number of channels (=1?)
                        leftborder,             // offset (=0?)
                        rightborder,            // offset (=0?)
                        bottomborder,           // offset (=0?)
                        topborder;              // offset (=0?)
  uint                  numdatarec,             // number of SAR DATA records
                        numlines,               // number of lines
                        numpixels,              // number of pixels
                        numbytesdata,           // bytes per line
                        lendatarec2;            // SAR DATA records length
  char                  c4[5],                  // correctly 5 for \0
                        c6[7],                  // correctly 7 for \0
                        c8[9];                  // correctly 9 for \0


  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
  uint rec_seq;// type B4
  unsigned char rec_sub1, rec_type, rec_sub2, rec_sub3;// type B1



// ______Open files______
  ifstream datfile;
  openfstream(datfile,readfiles_arg.datfile);
  bk_assert(datfile,readfiles_arg.datfile,__FILE__,__LINE__);


  // --- RECORD 1 ---
  DEBUG.print("record 1 of data file (ERS and RSAT).");
  datfile.seekg(0,ios::beg);//
  datfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  datfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  datfile.read((char*)&rec_type,sizeb1);// record type code
  datfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  datfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("ERS/RSAT: Expecting record 1 with code {63,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  if (int(rec_sub1)==63 && int(rec_type)==192 && int(rec_sub2)==18 && int(rec_sub3)==18)
    DEBUG.print("This is the expected record with code {63,192,18,18}");
  else
    WARNING.print("This is NOT the expected record with code {63,192,18,18}");
  uint                  lenrec1;                // length of record1
  datfile.seekg(8,ios::beg);//
  datfile.read((char*)&lenrec1,sizeb4);// length of record
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("ERS:  Expecting record 1 with length 10012");
  DEBUG.print("RSAT: Expecting record 1 with length 16252");
  DEBUG << "readdat::record 1: start at: " << 0
        << "; length: " << lenrec1;
  DEBUG.print();


  // ====== Get general info ======
  datfile.seekg(180,ios::beg);
  datfile.read((char*)&c6,sizei6);              // number of SAR DATA records (lines)
  c6[6]='\0';
    numdatarec = atoi(c6);
  datfile.read((char*)&c6,sizei6);              // SAR DATA record length
  c6[6]='\0';
    lendatarec2 = atoi(c6);
  datfile.seekg(232,ios::beg);                  // SAR Related data
  datfile.read((char*)&c4,4);
  c4[4]='\0';
  numchannels = atoi(c4);
  datfile.read((char*)&c8,sizei8);
  c8[8]='\0';
  numlines = atoi(c8);

  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    leftborder = atoi(c4);
  datfile.read((char*)&c8,sizei8);              // number of pixels
    c8[8]='\0';
    numpixels = atoi(c8);
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    rightborder = atoi(c4);
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    topborder = atoi(c4);
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    bottomborder = atoi(c4);

  datfile.seekg(280,ios::beg);                  // Record data
  datfile.read((char*)&c8,sizei8);
  c8[8]='\0';
  numbytesdata = atoi(c8);


// ======Write to scratchfiles======
  ofstream scratchlogdat("scratchlogdat", ios::out | ios::trunc);
  bk_assert(scratchlogdat,"readdat: scratchlogdat",__FILE__,__LINE__);

  scratchlogdat << "\n\n*******************************************************************"
                << "\n* EXTRACTED DATA FROM DATA FILE: "
                <<  readfiles_arg.datfile << " *"
                << "\n*******************************************************************"
                << "\nNumber of SAR channels in file:         \t"
                <<  numchannels
                << "\nNumber of lines:                        \t"
                <<  numlines
                << "\nNumber of pixels:                       \t"
                <<  numpixels
                << "\nNumber of left border pixels per line:  \t"
                <<  leftborder
                << "\nNumber of right border pixels per line: \t"
                <<  rightborder
                << "\nNumber of topborder lines:              \t"
                <<  topborder
                << "\nNumber of bottom border lines:          \t"
                <<  bottomborder
                << "\nNumber of bytes per data group:         \t"
                <<  numbytesdata
                <<  endl;
  scratchlogdat.close();


  // --- write RESULTFILE, add info for RSAT which was not in LEADER ---
  ofstream scratchresdat("scratchresdat", ios::out | ios::trunc);
  bk_assert(scratchresdat,"readdat: scratchresdat",__FILE__,__LINE__);
  scratchresdat << "Datafile: \t\t\t\t\t"
                <<  readfiles_arg.datfile
                << "\nNumber_of_lines_original: \t\t\t"
                <<  numlines
                << "\nNumber_of_pixels_original: \t\t\t"
                <<  numpixels;

  // --- if rsat, read datfile further, write az/rg time ----------------
//if (readfiles_arg.sensor_id==SLC_RSAT)
//if (readfiles_arg.sar_processor==SARPR_ATL)
  if ((readfiles_arg.sar_processor==SARPR_ATL) || (readfiles_arg.sar_processor==SARPR_JAX)) //LG
  {
  // --- RECORD 2 ---
  uint startrec2 = lenrec1;
  DEBUG.print("record 2 of data file (RSAT).");
  datfile.seekg(startrec2,ios::beg);//
  datfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);// Bert Kampes, 07-Apr-2005
  datfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  datfile.read((char*)&rec_type,sizeb1);// record type code
  datfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  datfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("RSAT: Expecting record 2 with code {50,11,18,20}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  if (int(rec_sub1)==50 && int(rec_type)==11 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected record with code {50,11,18,20}");
  else
    WARNING.print("This is NOT the expected record with code {50,11,18,20}");
  uint                  lenrec2;                // length of record2
  datfile.seekg(startrec2+8,ios::beg);//
  datfile.read((char*)&lenrec2,sizeb4);// length of record
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("RSAT: Expecting record 2 with length variable");
  DEBUG << "readdat::record 2: start at: " << startrec2
        << "; length: " << lenrec2;
  DEBUG.print();
  uint startrec3 = lenrec1+lenrec2;
  uint startrecN = lenrec1+(numlines-1)*lenrec2;// start of last record
  // --- azimuth time to first line (depends on decrease/increase): ---
  uint zdmsecofday1 = 99999;// B4
  uint zdmsecofday2 = 99999;// B4
  uint zdmsecofdayN = 99999;// B4
  datfile.seekg(startrec2+44,ios::beg);//
  datfile.read((char*)&zdmsecofday1,sizeb4);//  range to first pix
  zdmsecofday1 = ntohl(zdmsecofday1);// Bert Kampes, 07-Apr-2005
  datfile.seekg(startrec3+44,ios::beg);//
  datfile.read((char*)&zdmsecofday2,sizeb4);//  range to first pix
  zdmsecofday2 = ntohl(zdmsecofday2);// Bert Kampes, 07-Apr-2005
  datfile.seekg(startrecN+44,ios::beg);//
  datfile.read((char*)&zdmsecofdayN,sizeb4);//  range to first pix
  zdmsecofdayN = ntohl(zdmsecofdayN);// Bert Kampes, 07-Apr-2005
  INFO << "zdmsecofday1: " << zdmsecofday1;
  INFO.print();
  DEBUG << "zdmsecofday2: " << zdmsecofday2;
  DEBUG.print();
  INFO << "zdmsecofdayN: " << zdmsecofdayN;
  INFO.print();
  // --- Check PRF/RSR ---
  real8 prf_check = real8(numlines-1)/(abs(real8(zdmsecofday1)-real8(zdmsecofdayN))/1000.0);
  INFO << "PRF check (computed [Hz]): " << prf_check;
  INFO.print();


  // format should be: "22-AUG-1997 18:22:10.246"
  if (zdmsecofday1 < zdmsecofdayN)// increase, use ZD time of first line
    {
    INFO.print("INCREASING azimuth time detected RSAT");
    }
  else // decrease: use first-numlines*prf? or read last record.
    {
    INFO.print("DECREASING azimuth time detected (flip image up-down) RSAT");
    zdmsecofday1 = zdmsecofdayN;// use last record (line) as first line (flip)
    startrec2    = startrecN;// use last record (line) as first line (flip)
    }
  // --- Format the string UTC ---
  uint acq_year = 1999;
  uint acq_day  = 19;
  datfile.seekg(startrec2+36,ios::beg);//
  datfile.read((char*)&acq_year,sizeb4);//  range to first pix
  acq_year = ntohl(acq_year);// Bert Kampes, 07-Apr-2005
  datfile.read((char*)&acq_day,sizeb4);//  range to first pix
  acq_day = ntohl(acq_day);// Bert Kampes, 07-Apr-2005
  INFO << "acq_year: " << acq_year;
  INFO.print();
  INFO << "acq_day: " << acq_day;
  INFO.print();
  // --- fill tm struct using daynumber with strptime, re-format it using strftime ---
  char datestring[13];// e.g., "25-Jan-1999"
  struct tm tm_tmp;
  char buf[9]; // e.g., "1999 191";
  sprintf(buf,"%4d %03d", acq_year,acq_day);
  strptime(buf,"%Y %j",&tm_tmp);
  int q = strftime(datestring,12,"%d-%b-%Y",&tm_tmp);// q: number of bytes written
  INFO << "numbytes in datestring: (12?) " << q;
  INFO.print();

  // --- Fill important part: (sec of day) ---
  char c24zd1stazitime[25];
  int32 hour = int32(real8(zdmsecofday1)/1000.0/60.0/60.0);// floor
  int32 min  = int32((real8(zdmsecofday1)-hour*1000.0*60.0*60.0)/1000.0/60.0);// floor
  int32 sec  = int32((real8(zdmsecofday1)-hour*1000.0*60.0*60.0-min*1000.0*60.0)/1000.0);// floor
  int32 msec = int32((real8(zdmsecofday1)-hour*1000.0*60.0*60.0-min*1000.0*60.0-1000.0*sec));// floor
  //sprintf(c24zd1stazitime,"01-JAN-1990 %02d:%02d:%02d.%03d", hour,min,sec,msec);
  sprintf(c24zd1stazitime,"%11s %02d:%02d:%02d.%03d", datestring, hour,min,sec,msec);
  INFO << "c24zd1stazitime: " << c24zd1stazitime;
  INFO.print();

  // --- range time to first pixel is distance -------------------------
  uint range1st = 99999;//
  uint rangelst = 99999;
  // start_added_by_don
  real8 zd1strange;
  if (readfiles_arg.sar_processor==SARPR_JAX)
    {
    datfile.seekg(startrec2+116,ios::beg);//
    datfile.read((char*)&range1st,sizeb4);//  range to first pix
    range1st = ntohl(range1st);// Bert Kampes, 07-Apr-2005
    zd1strange = 2000.0*range1st/SOL;
    INFO << "range1st: " << range1st;
    INFO.print();
     }
  if (readfiles_arg.sar_processor==SARPR_ATL)
  // end_added_by_don
    {
  datfile.seekg(startrec2+64,ios::beg);//
  datfile.read((char*)&range1st,sizeb4);//  range to first pix
  range1st = ntohl(range1st);// Bert Kampes, 07-Apr-2005
  datfile.seekg(startrec2+72,ios::beg);//
  datfile.read((char*)&rangelst,sizeb4);//  range to last pix
  rangelst = ntohl(rangelst);// Bert Kampes, 07-Apr-2005
  zd1strange = (range1st<rangelst) ? 2000.0*range1st/SOL : 2000.0*rangelst/SOL;
  INFO << "range1st: " << range1st;
  INFO.print();
  INFO << "rangelst: " << rangelst;
  INFO.print();
  // --- Check PRF/RSR ---
  real8 rsr_check = real8(numpixels-1)/abs(2000000.0*(real8(range1st)-real8(rangelst))/SOL);
  INFO << "RSR check (computed [MHz]): " << rsr_check;
  INFO.print();
  }
  //char dummydate[] = "01-JAN-1990 ";// not used except maybe getorb later
  //WARNING.print("RSAT: using a dummy date for orbit, only secofday important.");
  scratchresdat
    << "\nFirst_pixel_azimuth_time (UTC):            \t"
    //    << "22-AUG-1997 18:22:10.246"
    << c24zd1stazitime
    << "\nRange_time_to_first_pixel (2way) (ms):     \t"
    << setprecision(16)
    <<  zd1strange;
  }

  scratchresdat << "\n*******************************************************************";
  if (readfiles_arg.fileid == MASTERID)
    scratchresdat <<  "\n* End_" << processcontrol[pr_m_readfiles] << "_NORMAL";
  if (readfiles_arg.fileid == SLAVEID)
    scratchresdat <<  "\n* End_" << processcontrol[pr_s_readfiles] << "_NORMAL";
  scratchresdat << "\n*******************************************************************"
                << endl;
  scratchresdat.close();
  datfile.close();

// ______Tidy up______
  if (numchannels != 1)                         // ??
    {
    WARNING << "code 904: Number of channels in file: "
         << readfiles_arg.datfile << " = "
         << numchannels << " != 1 ";
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }

// ______ Check with volumefile ______
  if (numlines != checklines)
    {
    WARNING << "code 902: data file: "
         << readfiles_arg.datfile
         << " numlin=" << numlines
         << " vs.  volume file: "
         << readfiles_arg.volfile
         << " numlin=" << checklines;
    WARNING.print();
    WARNING.print(" +this means data and volume file seem not to correspond.");
    }

// ______ Check with previous section ______
  if (numlines != numdatarec)
    {
    WARNING << "code 904: Number of lines seems not to be consistent in file: "
         << readfiles_arg.datfile << " : " << numlines << " != " << numdatarec;
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if (bottomborder != topborder != leftborder != rightborder != 0)
    {
    WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: "
         << leftborder << "," << rightborder << "," << bottomborder << ","
         << topborder << " in file: " << readfiles_arg.datfile;
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }

  // start_added_by_don
  if (readfiles_arg.sar_processor==SARPR_JAX)
     {
     if ((numbytesdata / 8) != numpixels)
       {
       WARNING << "code 904AAA: Number of pixels seems to be inconsistent in file: "
            << readfiles_arg.datfile << ": "
     << numpixels << " != " << (numbytesdata / 8);
       WARNING.print();
       WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
       }
     }
  else
     {
  if ((numbytesdata / 4) != numpixels)
    {
    WARNING << "code 904: Number of pixels seems to be inconsistent in file: "
         << readfiles_arg.datfile << ": "
         << numpixels << " != " << (numbytesdata / 4);
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  }
  PROGRESS.print("readdat (header info) finished.");
  } // END readdat



/****************************************************************
 *    writeslc                                                  *
 *                                                              *
 * Inputfile in ceos slc format is converted to                 *
 *  raw format outputfile.                                      *
 * Some info is extracted from header,                          *
 *  see also appendix C of CD-R distribution (ESA).             *
 * checks with volumefile #lines                                *
 *                                                              *
 * A buffer (BUFFERMEMSIZE) is used to increase speed.          *
 * Input is the inputfile and outputfile,                       *
 *  an identifier for this run,                                 *
 *  the area which should be converted (default is total image) *
 * info is written to scratchresfile                            *
 http://earth.esa.int/rootcollection/sysutil/01008.html
 *                                                              *
 *    Bert Kampes, 11-Dec-1998                                  *
 ****************************************************************/
void writeslc(
        const input_gen &generalinput,
        const input_crop &crop_arg,
        const int       checklines)
  {
  const int16           sizeb4 = 4,             // some constants for reading
                        sizei4 = 4,             //   binary fields.
                        sizei6 = 6,
                        sizei8 = 8;
  uint                  lenrec1;                // length of general record1
  uint                  lenrec2;                // (nominal) length of data records
  char                  c4[5],                  // correctly 5 for \0
                        c6[7],                  // correctly 7 for \0
                        c8[9];                  // correctly 9 for \0

  // ______ Write some info ______
  TRACE_FUNCTION("writeslc (BK 11-Dec-1998)")
  PROGRESS.print("Start cropping slc data.");
  #ifdef __X86PROCESSOR__
  INFO.print("Swapping Big Endian (CEOS input) to Little Endian (your platform).");
  #else
  INFO.print("NO byte swapping performed, you must be on Big Endian platform.");
  #endif

  // ______ Open files ______
  ifstream datfile;
  openfstream(datfile,crop_arg.filein1);
  bk_assert(datfile,crop_arg.filein1,__FILE__,__LINE__);

  // ====== Get data such as recordlength ======
  datfile.seekg(8,ios::beg);
  datfile.read((char*)&lenrec1,sizeb4);         // length of record1
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("record 1 of data file.");
  if (lenrec1 != 19612 )  // quarter scene 10012? JERS 22196?
    {
    WARNING << "writeslc: length of record 1 = \""
         <<  lenrec1 << "\"; expected \"19612\" for ERS SLC (CEOS, full scene).";
    WARNING.print();
    DEBUG.print("10012 seems ERS quarter scene?");
    DEBUG.print("22196 seems JERS scene?");
    }

  datfile.seekg(180,ios::beg);
  datfile.read((char*)&c6,sizei6);              // number of SAR DATA records (lines)
    c6[6]='\0';
    const uint numdatarec = atoi(c6);
  DEBUG << "numdatarec: " << numdatarec;
  DEBUG.print();
  //datfile.seekg(186,ios::beg);
  datfile.read((char*)&c6,sizei6);              // SAR DATA record length
    c6[6]='\0';
    const uint lendatarec2 = atoi(c6);
  DEBUG << "lendatarec2: " << lendatarec2;
  DEBUG.print();
  datfile.seekg(232,ios::beg);                  // SAR Related data
  datfile.read((char*)&c4,4);
    c4[4]='\0';
    const uint numchannels = atoi(c4);
  DEBUG << "numchannels: " << numchannels;
  DEBUG.print();

  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    uint numlines = atoi(c8);
  DEBUG << "numlines: " << numlines;
  DEBUG.print();

  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint leftborder = atoi(c4);
  DEBUG << "leftborder: " << leftborder;
  DEBUG.print();
  datfile.read((char*)&c8,sizei8);              // number of pixels
    c8[8]='\0';
    uint numpixels = atoi(c8);
  DEBUG << "numpixels: " << numpixels;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint rightborder = atoi(c4);
  DEBUG << "rightborder: " << rightborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint topborder = atoi(c4);
  DEBUG << "topborder: " << topborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint bottomborder = atoi(c4);
  DEBUG << "bottomborder: " << bottomborder;
  DEBUG.print();

  datfile.seekg(280,ios::beg);                  // Record data
  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    const uint numbytesdata = atoi(c8);
  DEBUG << "numbytesdata: " << numbytesdata;
  DEBUG.print();


// ====== Check with volumefile / internal ======
  if (numlines != checklines)
    {
    WARNING << "code 902: data file: "
         << crop_arg.filein1
         << " numlin=" << numlines
         << " vs. volume file: "
         << crop_arg.filein1
         << " numlin=" << checklines;
    WARNING.print();
    WARNING.print(" +this means data and volume file seem not to correspond.");
    }

// ______ Check with previous section ______
  if (numlines != numdatarec)
    {
    WARNING << "code 904: Number of lines seems not to be consistent in file: "
         << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if ((numbytesdata / 4) != numpixels)
    {
    WARNING << "code 904: Number of pixels seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << numpixels << " != " << (numbytesdata / 4);
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }


// ====== Start copy input to output (raw) format with buffer======
// ______ Check and process optional offset parameters______
// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
  uint linestart  = 1;                                  // counters for loops
  uint lineend    = numlines;
  uint pixelstart = 1;
  uint pixelend   = numpixels;                          // only for resultfile

  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    {
    window tempdbow(crop_arg.dbow.linelo, crop_arg.dbow.linehi,
                    crop_arg.dbow.pixlo,  crop_arg.dbow.pixhi);
    if (crop_arg.dbow.linehi>numlines)
      {
      WARNING << "Specified input DBOW linehi > numlines: "
           << crop_arg.dbow.linehi << " > " << numlines
           << ". I set linehi = " << numlines;
      WARNING.print();
      tempdbow.linehi=numlines;
      }
    if (crop_arg.dbow.pixhi>numpixels)
      {
      WARNING << "Specified input DBOW pixhi > numpixels: "
           << crop_arg.dbow.pixhi << " > " << numpixels
           << ". I set pixhi = " << numpixels;
      WARNING.print();
      tempdbow.pixhi=numpixels;
      }
// ______ Only hi values are possibly adapted, low is a constant ______
    numlines   = tempdbow.linehi - crop_arg.dbow.linelo + 1;
    numpixels  = tempdbow.pixhi  - crop_arg.dbow.pixlo  + 1;

    linestart  = crop_arg.dbow.linelo;
    lineend    = tempdbow.linehi;
    pixelstart = crop_arg.dbow.pixlo;
    pixelend   = tempdbow.pixhi;                        // only for resultfile
    }

// ______ Note complex<short> not in ANSI c ______
  matrix <int16>        LINE(1,2*numpixels);            // size of short

// ====== Process requested lines ======
  ofstream datoutfile;
  openfstream(datoutfile,crop_arg.fileout1,generalinput.overwrit);
  bk_assert(datoutfile,crop_arg.fileout1,__FILE__,__LINE__);

  // ______ info on data, to avoid X86 problems ______
  // ______ according to CEOS specs, byte 13-16 is first complex pixel, etc. ______
  datfile.seekg(lenrec1+12,ios::beg);
  matrix <int16> TMPSHORT(1,2);
  datfile >> TMPSHORT;          // read in first complex pixel for test
  real8 tmpmag = sqrt(
    real8(int16(ntohs(TMPSHORT(0,0)))*int16(ntohs(TMPSHORT(0,0)))) +
    real8(int16(ntohs(TMPSHORT(0,1)))*int16(ntohs(TMPSHORT(0,1)))));
  DEBUG << "First complex element in datafile: ("
       << int16(ntohs(TMPSHORT(0,0))) << ","
       << int16(ntohs(TMPSHORT(0,1)))
       << "); mag = " << tmpmag;
  DEBUG.print();
  if (tmpmag > 10000.)
    {
    WARNING.print(DEBUG.get_str());
    WARNING.print("this is a byteorder problem on X86? (use ntohs)");
    }
  DEBUG << "TEST: (realpart): " << TMPSHORT(0,0)
       << ", (imagpart): " << TMPSHORT(0,1);
  DEBUG.print();
  DEBUG << "TEST: htons(realpart): " << htons(TMPSHORT(0,0))
       << ", htons(imagpart): " << htons(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: ntohs(realpart): " << ntohs(TMPSHORT(0,0))
       << ", ntohs(imagpart): " << ntohs(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: short int(ntohs(realpart)): " << int16(ntohs(TMPSHORT(0,0)))
       << ", (imagpart): " << int16(ntohs(TMPSHORT(0,1)));
  DEBUG.print();


  // ====== perline is faster than perbuffer, less memory etc. BK1998 ======
  datfile.seekg(lenrec1+(linestart-1)*lendatarec2+8,ios::beg);
  datfile.read((char*)&lenrec2,sizeb4);         // length of first record
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.

  if (lenrec2 != lendatarec2)
    {
    ERROR << "code 904: Length of datarecords seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << lenrec2 << " != " << lendatarec2;
    WARNING.print(ERROR.get_str());
    ERROR.reset();
    }

  const int32 TEN        = 10;
  const int32 TENPERCENT = int32((.5*TEN+numlines)/TEN);        // number of lines
  int32 percentage       = 0;                                   // initialization
  const int32 tmpstart   = lenrec1+12-lendatarec2+(pixelstart-1)*4;     // sizeof=4
  for (register int32 linecnt=linestart; linecnt<=lineend; linecnt++)
    {
    if (!((linecnt-linestart)%TENPERCENT))
      {
      PROGRESS << "WRITESLC: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += TEN;
      }
    datfile.seekg(tmpstart+linecnt*lendatarec2,ios::beg);
    datfile    >> LINE;
    // ______ BK 13 July 2000: swapbytes for X86 (intel) linux cpus ______
    #ifdef __X86PROCESSOR__
    for (int ii=0; ii<LINE.pixels(); ++ii)
      LINE(0,ii) = int16(ntohs(LINE(0,ii)));    // changed from htons 171100 BK
    #endif
    datoutfile << LINE;
    }
  datfile.close();                                      // close files
  datoutfile.close();



// ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
    //<< "\n*_Start_crop:\t\t\t"
    //<<  crop_arg.idcrop
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"

// ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  linestart
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  lineend
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  pixelstart
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  pixelend
    << "\nNumber of lines (non-multilooked): \t\t" <<  lineend-linestart+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  pixelend-pixelstart+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();

// ______ Tidy up do checks here ______
  if (numchannels != 1)                         // ??
    {
    WARNING << "code 904: Number of channels in file: "
         << crop_arg.filein1 << " = "
         << numchannels << " != 1 ";
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if (bottomborder != topborder != leftborder != rightborder != 0)
    {
    WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: "
         << leftborder << "," << rightborder << "," << bottomborder << ","
         << topborder << " in file: " << crop_arg.filein1;
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  PROGRESS.print("WRITESLC: 100%");
  } // END writeslc

/****************************************************************
 *    gammaprocessor_crop                                       *
 *    crop function for PROC_GAM                                *
 *    SLC has to be in pixel interleaved Complex Short Int frmt *
 *    meaning: real part = 2bytes, imag part= 2 bytes 		*
 *                                                              *
 *    Batuhan Osmanoglu, 15-DEC-2009, after BKampes of course ;)*
 ****************************************************************/
void gammaprocessor_crop(
        const input_gen &generalinput,
        const slcimage &master,
        const input_crop &crop_arg)
  {
  const int16           sizeb4 = 4,             // some constants for reading
                        sizei4 = 4,             //   binary fields.
                        sizei6 = 6,
                        sizei8 = 8;
  uint                  lenOneLine;                // (nominal) length of data records
  char                  c4[5],                  // correctly 5 for \0
                        c6[7],                  // correctly 7 for \0
                        c8[9];                  // correctly 9 for \0

  // ______ Write some info ______
  TRACE_FUNCTION("gammaprocessor_crop (Batu 15-Dec-2009)")
  PROGRESS.print("Start cropping Gamma focused slc data.");
  INFO.print("Data assumed in the host platform format. NO byte swapping performed.");


  uint numlines=master.originalwindow.linehi-master.originalwindow.linelo+1;
  uint numpixels=master.originalwindow.pixhi-master.originalwindow.pixlo+1;
  // ====== Check data such as recordlength ======
  if (numlines <= 0)
    {
    WARNING << "Number of lines (less than zero?) ="
         <<  numlines
         <<  "Check input result file.";
    WARNING.print();
    }
  if (numpixels <= 0)
    {
    WARNING << "Number of pixels (less than zero?) ="
         <<  numpixels
         <<  "Check input result file.";
    WARNING.print();
    }


// ====== Start copy input to output (raw) format with buffer======
// ______ Check and process optional offset parameters______
// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
  uint linestart  = 1;                                  // counters for loops
  uint lineend    = numlines;
  uint pixelstart = 1;
  uint pixelend   = numpixels;                          // only for resultfile

  lenOneLine      = numpixels*4;			// each complex pixel is 4 bytes.BO.

  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    {
    window tempdbow(crop_arg.dbow.linelo, crop_arg.dbow.linehi,
                    crop_arg.dbow.pixlo,  crop_arg.dbow.pixhi);
    if (crop_arg.dbow.linehi>numlines)
      {
      WARNING << "Specified input DBOW linehi > numlines: "
           << crop_arg.dbow.linehi << " > " << numlines
           << ". I set linehi = " << numlines;
      WARNING.print();
      tempdbow.linehi=numlines;
      }
    if (crop_arg.dbow.pixhi>numpixels)
      {
      WARNING << "Specified input DBOW pixhi > numpixels: "
           << crop_arg.dbow.pixhi << " > " << numpixels
           << ". I set pixhi = " << numpixels;
      WARNING.print();
      tempdbow.pixhi=numpixels;
      }
// ______ Only hi values are possibly adapted, low is a constant ______
    numlines   = tempdbow.linehi - crop_arg.dbow.linelo + 1;
    numpixels  = tempdbow.pixhi  - crop_arg.dbow.pixlo  + 1;

    linestart  = crop_arg.dbow.linelo;
    lineend    = tempdbow.linehi;
    pixelstart = crop_arg.dbow.pixlo;
    pixelend   = tempdbow.pixhi;                        // only for resultfile
    }

// ______ Note complex<short> not in ANSI c ______
  matrix <int16>        LINE(1,2*numpixels);            // size of short

  ifstream datfile;
  openfstream(datfile,crop_arg.filein1);
  bk_assert(datfile,crop_arg.filein1,__FILE__,__LINE__);

// ====== Process requested lines ======
  ofstream datoutfile;
  openfstream(datoutfile,crop_arg.fileout1,generalinput.overwrit);
  bk_assert(datoutfile,crop_arg.fileout1,__FILE__,__LINE__);

  // ______ info on data, to avoid X86 problems ______
  // ______ go to the beginning of file to first complex pixel, etc. ______
  datfile.seekg(0, ios::beg);
  matrix <int16> TMPSHORT(1,2);
  datfile >> TMPSHORT;          // read in first complex pixel for test
  real8 tmpmag = sqrt(
    real8(int16(ntohs(TMPSHORT(0,0)))*int16(ntohs(TMPSHORT(0,0)))) +
    real8(int16(ntohs(TMPSHORT(0,1)))*int16(ntohs(TMPSHORT(0,1)))));
  DEBUG << "First complex element in datafile: ("
       << int16(ntohs(TMPSHORT(0,0))) << ","
       << int16(ntohs(TMPSHORT(0,1)))
       << "); mag = " << tmpmag;
  DEBUG.print();
  if (tmpmag > 10000.)
    {
    WARNING.print(DEBUG.get_str());
    WARNING.print("this is a byteorder problem on X86? (use ntohs)");
    }
  DEBUG << "TEST: (realpart): " << TMPSHORT(0,0)
       << ", (imagpart): " << TMPSHORT(0,1);
  DEBUG.print();
  DEBUG << "TEST: htons(realpart): " << htons(TMPSHORT(0,0))
       << ", htons(imagpart): " << htons(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: ntohs(realpart): " << ntohs(TMPSHORT(0,0))
       << ", ntohs(imagpart): " << ntohs(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: short int(ntohs(realpart)): " << int16(ntohs(TMPSHORT(0,0)))
       << ", (imagpart): " << int16(ntohs(TMPSHORT(0,1)));
  DEBUG.print();


  // ====== perline is faster than perbuffer, less memory etc. BK1998 ======
  const int32 TEN        = 10;
  const int32 TENPERCENT = int32((.5*TEN+numlines)/TEN);        // number of lines
  int32 percentage       = 0;                                   // initialization
  const int32 offset     = (pixelstart-1)*4;     // sizeof=4
  for (register int32 linecnt=linestart; linecnt<=lineend; linecnt++)
    {
    if (!((linecnt-linestart)%TENPERCENT))
      {
      PROGRESS << "WRITESLC: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += TEN;
      }
    datfile.seekg((linecnt-1)*lenOneLine+offset, ios::beg);
    datfile    >> LINE;
    datoutfile << LINE;
    }
  datfile.close();                                      // close files
  datoutfile.close();

// ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
    //<< "\n*_Start_crop:\t\t\t"
    //<<  crop_arg.idcrop
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"

// ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  linestart
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  lineend
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  pixelstart
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  pixelend
    << "\nNumber of lines (non-multilooked): \t\t" <<  lineend-linestart+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  pixelend-pixelstart+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();

  PROGRESS.print("WRITESLC: 100%");
  } // END gammaprocessor_crop


/****************************************************************
 *    envisatdump_data                                          *
 *                                                              *
 * Via a system call to the c-program envisat_dump_data, the    *
 * SLC data is wrtten to file.  The resfile is here created.    *
 * envisatdumpdata writes SLC data out in host order.           *
 * it is important that crop_arg.dbow is correctly filled.      *
 *    Bert Kampes, 16-JUN-2003                                  *
 ****************************************************************/
void envisat_dump_data(
        const input_crop &crop_arg)
  {
  // ______ Write some info ______
  TRACE_FUNCTION("envisat_dump_data (BK 16-Jun-2003)")
  // ______ Build command ______
  // ______ make sure l0 etc. are correctly defined ______
  // ____ assume these are filled correctly ___
  INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "envisat_dump_data " << crop_arg.filein1
         << " " << crop_arg.fileout1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "envisat_dump_data " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command the envisat data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  system(cmd);// this does the work
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
  } // END envisat_dump_data

// AP_VV functions, later migrate to above one with a switch

void envisat_dump_VV(
        const input_crop &crop_arg)
  {
  // ______ Write some info ______
  TRACE_FUNCTION("envisat_dump_VV (MCC 16-Jun-2003)")
  // ______ Build command ______
  // ______ make sure l0 etc. are correctly defined ______
  // ____ assume these are filled correctly ___
  INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "envisat_dump_VV " << crop_arg.filein1
         << " " << crop_arg.fileout1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "envisat_dump_VV " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command the envisat data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  system(cmd);// this does the work
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
  } // END envisat_dump_VV


void envisat_dump_HH(
        const input_crop &crop_arg)
  {
  // ______ Write some info ______
  TRACE_FUNCTION("envisat_dump_HH (BK 16-Jun-2003)")
  // ______ Build command ______
  // ______ make sure l0 etc. are correctly defined ______
  // ____ assume these are filled correctly ___
  INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "envisat_dump_HH " << crop_arg.filein1
         << " " << crop_arg.fileout1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "envisat_dump_HH " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command the envisat data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  system(cmd);// this does the work
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
  } // END envisat_dump_HH


/****************************************************************
 *    tsxdump_data                                              *
 *                                                              *
 * Via a system call to the python tsx_dump_data, the           *
 * SLC data is written to file.  The resfile is here created.   *
 * tsxdumpdata writes SLC data out in host order.               *
 * it is important that crop_arg.dbow is correctly filled.      *
 *                                                              *
 * Dependencies: GDAL                                           *
 ****************************************************************/
void tsx_dump_data(
       const input_crop &crop_arg)
{
  // ______ Write some info ______
  TRACE_FUNCTION("tsx_dump_data (PM 06-Apr-2009)")
    // ______ Build command ______
    // ______ make sure l0 etc. are correctly defined ______
    // ____ assume these are filled correctly ___
    int16 status = 0;    // [MA] check exit status of system calls for proper error handling
    INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "tsx_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1
         //<< " " << crop_arg.dbow.linelo - 1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         //<< " " << crop_arg.dbow.pixlo - 1
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "tsx_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command TSX cosar data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  status=system(cmd);// this does the work
  if (status != 0)                                                          // [MA] TODO make it a function
    {
    ERROR << "tsx_dump_data.py: failed with exit code: " << status;
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
} // END tsx_dump_data


/****************************************************************
 *    rs2dump_data                                              *
 *                                                              *
 * Dump Radarsat-2 TIFF file to raw                             *
 * via a system call to the python rs2_dump_data, the           *
 * SLC data is written to file.  The resfile is here created.   *
 * rs2dumpdata writes SLC data out in host order.               *
 * it is important that crop_arg.dbow is correctly filled.      *
 *                                                              *
 * Dependencies: GDAL, PYTHON                                   *
 *                                                              *
 *                                                              *
 ****************************************************************/
void rs2_dump_data(
       const input_crop &crop_arg)
{
  // ______ Write some info ______
  TRACE_FUNCTION("rs2_dump_data (MA,PM 04-Oct-2009)")
    // ______ Build command ______
    // ______ make sure l0 etc. are correctly defined ______
    // ____ assume these are filled correctly ___
    int16 status = 0;                                                       // [MA] check exit status of system calls for proper error handling
    INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "rs2_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1
         //<< " " << crop_arg.dbow.linelo - 1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         //<< " " << crop_arg.dbow.pixlo - 1
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "rs2_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command RS2 data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  status=system(cmd);// this does the work
  if (status != 0)                                                          // [MA] TODO make it a function
    {
    ERROR << "rs2_dump_data.py: failed with exit code: " << status;
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
} // END rs2_dump_data


/****************************************************************
 *    cskdump_data                                              *
 *                                                              *
 * Dump Cosmo-Skymed HD5 file to raw                            *
 * via a system call to the python csk_dump_data, the           *
 * SLC data is written to file.  The resfile is here created.   *
 * cskdumpdata writes SLC data out in host order.               *
 * it is important that crop_arg.dbow is correctly filled.      *
 *                                                              *
 * Dependencies: HD5, PYTHON                                    *
 *                                                              *
 *                                                              *
 ****************************************************************/
void csk_dump_data(
       const input_crop &crop_arg)
{
  // ______ Write some info ______
  TRACE_FUNCTION("csk_dump_data (MA,PD 23-Jun-2010)")
    // ______ Build command ______
    // ______ make sure l0 etc. are correctly defined ______
    // ____ assume these are filled correctly ___
    int16 status = 0;                                                       // [MA] check exit status of system calls for proper error handling
    INFO.reset();
  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    INFO << "csk_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1
         //<< " " << crop_arg.dbow.linelo - 1
         << " " << crop_arg.dbow.linelo
         << " " << crop_arg.dbow.linehi
         //<< " " << crop_arg.dbow.pixlo - 1
         << " " << crop_arg.dbow.pixlo
         << " " << crop_arg.dbow.pixhi << ends;
  else
    INFO << "csk_dump_data.py " << crop_arg.filein1
         << " " << crop_arg.fileout1 << ends;
  char cmd[512];// command string
  strcpy(cmd, INFO.get_str());
  INFO.print("With following command CSK data was cropped.");
  INFO.print(cmd);
  PROGRESS.print("system call may take some time...");
  status=system(cmd);// this does the work
  if (status != 0)                                                          // [MA] TODO make it a function
    {
    ERROR << "csk_dump_data.py: failed with exit code: " << status;
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  INFO.reset();
  INFO.print();

  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"
    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linelo
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.linehi
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixlo
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  crop_arg.dbow.pixhi
    << "\nNumber of lines (non-multilooked): \t\t" <<  crop_arg.dbow.linehi-crop_arg.dbow.linelo+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  crop_arg.dbow.pixhi-crop_arg.dbow.pixlo+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();
} // END csk_dump_data


/****************************************************************
 *    radarsat_dump_data                                        *
 *                                                              *
 * Inputfile in ceos slc format is converted to                 *
 *  raw format outputfile.                                      *
 * image is flipped if time direction is not INCREASE,INCREASE  *
 * crop is thus done either geometric correct, or if specified  *
 #%// Bert Kampes, 04-Aug-2004
 * DBOW a bit strange ... ?
 #%// Fix for DBOW for cropping radarsat images - Andy  Mar,2009
 ****************************************************************/
void radarsat_dump_data(
        const input_gen &generalinput,
        const input_crop &crop_arg)
  {
  const int16           sizeb4 = 4,             // some constants for reading
                        sizeb1 = 1,             //   binary fields.
                        sizei4 = 4,             //   binary fields.
                        sizei6 = 6,
                        sizei8 = 8;
  uint                  lenrec1;                // length of general record1
  uint                  lenrec2;                // (nominal) length of data records
  char                  c4[5],                  // correctly 5 for \0
                        c6[7],                  // correctly 7 for \0
                        c8[9];                  // correctly 9 for \0
  // --- Check for RSAT #%// Bert Kampes, 02-Aug-2004 ---
  uint rec_seq;// type B4
  unsigned char rec_sub1, rec_type, rec_sub2, rec_sub3;// type B1

  // ______ Write some info ______
  TRACE_FUNCTION("radarsat_dump_data (Bert Kampes 04-Aug-2004)")
  PROGRESS.print("Start cropping slc data for RADARSAT.");
  #ifdef __X86PROCESSOR__
  INFO.print("Swapping Big Endian (CEOS input) to Little Endian (your platform).");
  #else
  INFO.print("NO byte swapping performed, you must be on Big Endian platform.");
  #endif

  // ______ Open files ______
  ifstream datfile;
  openfstream(datfile,crop_arg.filein1);
  bk_assert(datfile,crop_arg.filein1,__FILE__,__LINE__);

  // ====== Get data such as recordlength ======
  // --- RECORD 1 ---
  DEBUG.print("record 1 of data file (ERS and RSAT).");
  datfile.seekg(0,ios::beg);//
  datfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);     // bk 6 jul 2000, byteorder x86 machines.
  datfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  datfile.read((char*)&rec_type,sizeb1);// record type code
  datfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  datfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("RSAT: Expecting record 1 with code {63,192,18,18}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  if (int(rec_sub1)==63 && int(rec_type)==192 && int(rec_sub2)==18 && int(rec_sub3)==18)
    DEBUG.print("This is the expected record with code {63,192,18,18}");
  else
    WARNING.print("This is NOT the expected record with code {63,192,18,18}");
  datfile.seekg(8,ios::beg);//
  datfile.read((char*)&lenrec1,sizeb4);// length of record
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("RSAT: Expecting record 1 with length 16252");
  DEBUG << "radarsat_dump_data::record 1: start at: " << 0
        << "; length: " << lenrec1;
  DEBUG.print();

  // --- Get some info in RECORD 1 ---
  datfile.seekg(180,ios::beg);
  datfile.read((char*)&c6,sizei6);              // number of SAR DATA records (lines)
    c6[6]='\0';
    const uint numdatarec = atoi(c6);
  DEBUG << "numdatarec: " << numdatarec;
  DEBUG.print();
  //datfile.seekg(186,ios::beg);
  datfile.read((char*)&c6,sizei6);              // SAR DATA record length
    c6[6]='\0';
    const uint lendatarec2 = atoi(c6);
  DEBUG << "lendatarec2: " << lendatarec2;
  DEBUG.print();
  datfile.seekg(232,ios::beg);                  // SAR Related data
  datfile.read((char*)&c4,4);
    c4[4]='\0';
    const uint numchannels = atoi(c4);
  DEBUG << "numchannels: " << numchannels;
  DEBUG.print();

  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    uint numlines = atoi(c8);
  DEBUG << "numlines: " << numlines;
  DEBUG.print();

  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint leftborder = atoi(c4);
  DEBUG << "leftborder: " << leftborder;
  DEBUG.print();
  datfile.read((char*)&c8,sizei8);              // number of pixels
    c8[8]='\0';
    uint numpixels = atoi(c8);
  DEBUG << "numpixels: " << numpixels;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint rightborder = atoi(c4);
  DEBUG << "rightborder: " << rightborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint topborder = atoi(c4);
  DEBUG << "topborder: " << topborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint bottomborder = atoi(c4);
  DEBUG << "bottomborder: " << bottomborder;
  DEBUG.print();

  datfile.seekg(280,ios::beg);                  // Record data
  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    const uint numbytesdata = atoi(c8);
  DEBUG << "numbytesdata: " << numbytesdata;
  DEBUG.print();

// ______ Check with previous section ______
  if (numlines != numdatarec)
    {
    WARNING << "code 904: Number of lines seems not to be consistent in file: "
         << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if ((numbytesdata / 4) != numpixels)
    {
    WARNING << "code 904: Number of pixels seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << numpixels << " != " << (numbytesdata / 4);
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }


// ====== Start copy input to output (raw) format with buffer======
// ______ Check and process optional offset parameters______
// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
  uint linestart  = 1;                  // counters for loops, first=1;
  uint lineend    = numlines;
  uint pixelstart = 1;                  // first pix is 1
  uint pixelend   = numpixels;          // only for resultfile
  uint orig_numlines = numlines;
  uint orig_numpixels = numpixels;

  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    {
    WARNING.print("cropping data may be difficult, due to INC/DECREASE storage");
    window tempdbow(crop_arg.dbow.linelo, crop_arg.dbow.linehi,
                    crop_arg.dbow.pixlo,  crop_arg.dbow.pixhi);
    if (crop_arg.dbow.linehi>numlines)
      {
      WARNING << "Specified input DBOW linehi > numlines: "
           << crop_arg.dbow.linehi << " > " << numlines
           << ". I set linehi = " << numlines;
      WARNING.print();
      tempdbow.linehi=numlines;
      }
    if (crop_arg.dbow.pixhi>numpixels)
      {
      WARNING << "Specified input DBOW pixhi > numpixels: "
           << crop_arg.dbow.pixhi << " > " << numpixels
           << ". I set pixhi = " << numpixels;
      WARNING.print();
      tempdbow.pixhi=numpixels;
      }
    // ______ Only hi values are possibly adapted, low is a constant ______
    numlines   = tempdbow.linehi - crop_arg.dbow.linelo + 1;
    numpixels  = tempdbow.pixhi  - crop_arg.dbow.pixlo  + 1;

    linestart  = crop_arg.dbow.linelo;
    lineend    = tempdbow.linehi;
    pixelstart = crop_arg.dbow.pixlo;
    pixelend   = tempdbow.pixhi;                        // only for resultfile
    }


  // --- Find out if data is stored increasing in line/pix ---
  // --- RECORD 2 ---
  uint startrec2 = lenrec1;
  DEBUG.print("record 2 of data file (RSAT).");
  datfile.seekg(startrec2,ios::beg);//
  datfile.read((char*)&rec_seq,sizeb4);//  record number
  rec_seq = ntohl(rec_seq);     // bk 6 jul 2000, byteorder x86 machines.
  datfile.read((char*)&rec_sub1,sizeb1);// first record sub type code
  datfile.read((char*)&rec_type,sizeb1);// record type code
  datfile.read((char*)&rec_sub2,sizeb1);// second record sub type code
  datfile.read((char*)&rec_sub3,sizeb1);// third record sub type code
  DEBUG.print("RSAT: Expecting record 2 with code {50,11,18,20}");
  DEBUG << "rec_seq: " << rec_seq
        << "; rec_sub1: " << int(rec_sub1)
        << "; rec_type: " << int(rec_type)
        << "; rec_sub2: " << int(rec_sub2)
        << "; rec_sub3: " << int(rec_sub3);
  DEBUG.print();
  if (int(rec_sub1)==50 && int(rec_type)==11 && int(rec_sub2)==18 && int(rec_sub3)==20)
    DEBUG.print("This is the expected record with code {50,11,18,20}");
  else
    WARNING.print("This is NOT the expected record with code {50,11,18,20}");
  datfile.seekg(startrec2+8,ios::beg);//
  datfile.read((char*)&lenrec2,sizeb4);// length of record
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("RSAT: Expecting record 2 with length variable");
  DEBUG << "radarsat_dump_data::record 2: start at: " << startrec2
        << "; length: " << lenrec2;
  DEBUG.print();
  uint startrec3 = lenrec1+lenrec2;
  uint startrecN = lenrec1+(numlines-1)*lenrec2;// start of last record
  // --- azimuth time to first line (depends on decrease/increase): ---
  uint zdmsecofday1 = 99999;// B4
  uint zdmsecofday2 = 99999;// B4
  uint zdmsecofdayN = 99999;// B4
  datfile.seekg(startrec2+44,ios::beg);//
  datfile.read((char*)&zdmsecofday1,sizeb4);//  range to first pix
  zdmsecofday1 = ntohl(zdmsecofday1);   // bk 6 jul 2000, byteorder x86 machines.
  datfile.seekg(startrec3+44,ios::beg);//
  datfile.read((char*)&zdmsecofday2,sizeb4);//  range to first pix
  zdmsecofday2 = ntohl(zdmsecofday2);   // bk 6 jul 2000, byteorder x86 machines.
  datfile.seekg(startrecN+44,ios::beg);//
  datfile.read((char*)&zdmsecofdayN,sizeb4);//  range to first pix
  zdmsecofdayN = ntohl(zdmsecofdayN);   // bk 6 jul 2000, byteorder x86 machines.
  INFO << "zdmsecofday1: " << zdmsecofday1;
  INFO.print();
  DEBUG << "zdmsecofday2: " << zdmsecofday2;
  DEBUG.print();
  INFO << "zdmsecofdayN: " << zdmsecofdayN;
  INFO.print();
  bool increasing_line = true;// assume this
  // --- I assume linestart is smaller than end, so swap them if required ---
  if (zdmsecofday1 < zdmsecofdayN)// increase, use ZD time of first line
    {
    INFO.print("INCREASE line direction detected (OK).");
    }
  else // decreasing lines: flip up-down
    {
    WARNING.print("DECREASE line direction detected: I will flip up-down the RSAT data");
    increasing_line = false;
    }

  // --- range time to first pixel is distance -------------------------
  uint range1st = 99999;//
  uint rangelst = 99999;
  datfile.seekg(startrec2+64,ios::beg);//
  datfile.read((char*)&range1st,sizeb4);//  range to first pix
  range1st = ntohl(range1st);   // bk 6 jul 2000, byteorder x86 machines.
  datfile.seekg(startrec2+72,ios::beg);//
  datfile.read((char*)&rangelst,sizeb4);//  range to last pix
  rangelst = ntohl(rangelst);   // bk 6 jul 2000, byteorder x86 machines.
  INFO << "range1st: " << range1st;
  INFO.print();
  INFO << "rangelst: " << rangelst;
  INFO.print();
  bool increasing_pix  = true;// assume this
  if (range1st < rangelst)// increase
    {
    INFO.print("INCREASE pixel direction detected (OK).");
    }
  else // decrease: flip data left-right
    {
    WARNING.print("DECREASE pixel direction detected: I will flip left-right the RSAT data");
    increasing_pix = false;
    }



  // ====== Process requested lines ======
  ofstream datoutfile;
  openfstream(datoutfile,crop_arg.fileout1,generalinput.overwrit);
  bk_assert(datoutfile,crop_arg.fileout1,__FILE__,__LINE__);

  // ______ info on data, to avoid X86 problems ______
  // ______ according to RSAT CEOS specs, byte 192-195 is first complex pixel, etc. ______
  datfile.seekg(lenrec1+192,ios::beg);
  matrix <int16> TMPSHORT(1,2);
  datfile >> TMPSHORT;          // read in first complex pixel for test
  real8 tmpmag = sqrt(
    real8(int16(ntohs(TMPSHORT(0,0)))*int16(ntohs(TMPSHORT(0,0)))) +
    real8(int16(ntohs(TMPSHORT(0,1)))*int16(ntohs(TMPSHORT(0,1)))));
  DEBUG << "First complex element in datafile: ("
       << int16(ntohs(TMPSHORT(0,0))) << ","
       << int16(ntohs(TMPSHORT(0,1)))
       << "); mag = " << tmpmag;
  DEBUG.print();
  if (tmpmag > 10000.)
    {
    WARNING.print(DEBUG.get_str());
    WARNING.print("this is a byteorder problem on X86? (use ntohs)");
    }
  DEBUG << "TEST: (realpart): " << TMPSHORT(0,0)
       << ", (imagpart): " << TMPSHORT(0,1);
  DEBUG.print();
  DEBUG << "TEST: htons(realpart): " << htons(TMPSHORT(0,0))
       << ", htons(imagpart): " << htons(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: ntohs(realpart): " << ntohs(TMPSHORT(0,0))
       << ", ntohs(imagpart): " << ntohs(TMPSHORT(0,1));
  DEBUG.print();
  DEBUG << "TEST: short int(ntohs(realpart)): " << int16(ntohs(TMPSHORT(0,0)))
       << ", (imagpart): " << int16(ntohs(TMPSHORT(0,1)));
  DEBUG.print();


  // --- Simple way pix by pix reading so we can easily flip if req. ---
  datfile.seekg(lenrec1+(linestart-1)*lendatarec2+8,ios::beg);
  datfile.read((char*)&lenrec2,sizeb4);         // length of first record
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.
  if (lenrec2 != lendatarec2)
    {
    ERROR << "code 904: Length of datarecords seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << lenrec2 << " != " << lendatarec2;
    WARNING.print(ERROR.get_str());
    ERROR.reset();
    }


  // ______ Note complex<short> not in ANSI c ______
  int32 percentage       = 0;                           // initialization
  const int32 TENPERCENT = int32((5+numlines)/10);      // number of lines
  for (register int32 linecnt=linestart; linecnt<=lineend; linecnt++)
    {
    // 03/2009 AH
    //int32 line2read = (increasing_line==true) ? linecnt : lineend-(linecnt-linestart);
    int32 line2read = (increasing_line==true) ? linecnt : orig_numlines-linecnt+1;
    for (register int32 pixcnt=pixelstart; pixcnt<=pixelend; pixcnt++)
      {
      // 03/2009 AH
      //int32 pix2read = (increasing_pix==true) ? pixcnt : pixelend-(pixcnt-pixelstart);
      int32 pix2read = (increasing_pix==true) ? pixcnt : orig_numpixels-pixcnt+1;
      // --- set pointer before complex pixel to read ---
      uint tmpstart = lenrec1+192+(line2read-1)*lendatarec2+(pix2read-1)*4;
      datfile.seekg(tmpstart,ios::beg);
      datfile    >> TMPSHORT;
      // ______ BK 13 July 2000: swapbytes for X86 (intel) linux cpus ______
      #ifdef __X86PROCESSOR__
      TMPSHORT(0,0) = int16(ntohs(TMPSHORT(0,0)));// byteswap
      TMPSHORT(0,1) = int16(ntohs(TMPSHORT(0,1)));// byteswap
      #endif
      datoutfile << TMPSHORT;
      }
    if (!((linecnt-linestart)%TENPERCENT))
      {
      PROGRESS << "radarsat_dump_data: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += 10;
      }
    }
  datfile.close();                                      // close files
  datoutfile.close();



// ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_short"

    // ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  linestart
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  lineend
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  pixelstart
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  pixelend
    << "\nNumber of lines (non-multilooked): \t\t" <<  lineend-linestart+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  pixelend-pixelstart+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();

// ______ Tidy up do checks here ______
  if (numchannels != 1)                         // ??
    {
    WARNING << "code 904: Number of channels in file: "
         << crop_arg.filein1 << " = "
         << numchannels << " != 1 ";
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if (bottomborder != topborder != leftborder != rightborder != 0)
    {
    WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: "
         << leftborder << "," << rightborder << "," << bottomborder << ","
         << topborder << " in file: " << crop_arg.filein1;
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  PROGRESS.print("radarsat_dump_data: 100%");
  } // END radarsat_dump_data




//____RaffaeleNutricato START MODIFICATION SECTION 2
/****************************************************************
 *    OversampleSLC                                             *
 *                                                              *
 * Oversamples the SLC by an integer factor.                    *
 * For now only range oversampling is performed.                *
 *    Raffaele Nutricato, 12-Jan-2004                           *
 * Azimuth oversampling with factor 2 added.                    *
 *    Bert Kampes, 30-Jul-2005                                  *
 ***************************************************************/
void  OversampleSLC(
        const input_gen        &generalinput,
        const slcimage         &imageinfo,
        const input_oversample &oversampleinput,
        const int16            fileid)
  {
  TRACE_FUNCTION("OversampleSLC (Raffaele Nutricato 12-Jan-2004)")
  //char infile[EIGHTY];  // Input file which is master/slave.raw renamed as .old
  //char outfile[EIGHTY]; // Output file which is the oversampled version.
  char infile[2*ONE27];  // Input file which is master/slave.raw renamed as .old  // MA
  char outfile[2*ONE27]; // Output file which is the oversampled version.
  strcpy(infile,imageinfo.file);
  strcpy(outfile,oversampleinput.fileoutovs);
  const int32 OsrRange   = oversampleinput.OsrRange;   // Range oversampling ratio.
  const int32 OsrAzimuth = oversampleinput.OsrAzimuth; // Azimuth oversampling ratio.
  const int32 FilterSize = oversampleinput.FilterSize; // Length of the kernel for the oversampling in range.
  if (OsrAzimuth!=1 && OsrAzimuth!=2)
    {
    ERROR.print("oversampling in azimuth: only factor 2");
    throw(some_error);
    }

  // ______Open input file_____
  ifstream ifile;
  openfstream(ifile,infile);
  bk_assert(ifile,infile,__FILE__,__LINE__);

  // ______Open output file_____
  ofstream ofile;
  openfstream(ofile,outfile,generalinput.overwrit);
  bk_assert(ofile,outfile,__FILE__,__LINE__);

  // ______ Compute the size of original cropped image ______
  //const int32 numlines  = imageinfo.currentwindow.lines();
  //const int32 numpixels = imageinfo.currentwindow.pixels();
  const uint numlines  = imageinfo.currentwindow.lines();  // MA
  const uint numpixels = imageinfo.currentwindow.pixels();

  // ______ Define the accuracy of the digital signal processings ______
  // ______ in range (in azimuth I use float, BK ______
  #define RN_DSP_ACCURACY    double
  #define RN_DSP_CPXACCURACY complr8

  // ______ Interpolation kernel section ______
  matrix <RN_DSP_CPXACCURACY> LINE_IN(1,OsrRange*(numpixels-1)+1);// zero alternating
  matrix <RN_DSP_CPXACCURACY> LINE_OUT(1,OsrRange*numpixels);// ovs line
  const int32 interp_size         = FilterSize  * OsrRange - 1;
  matrix <RN_DSP_CPXACCURACY> INTERP_KERNEL(1,interp_size);

  // ______ Generate the range interpolator impulse response ______
  const RN_DSP_ACCURACY invosr    = 1.0/RN_DSP_ACCURACY(OsrRange);
  RN_DSP_ACCURACY interpsamplepos = invosr - RN_DSP_ACCURACY(FilterSize)/2.0;
  for (int32 i=0; i<interp_size; i++)
    {
    INTERP_KERNEL(0,i) = RN_DSP_CPXACCURACY(sinc(interpsamplepos),0);
    interpsamplepos   += invosr;
    }

  // ______ Normalize kernel (BK) ______
  // ______ (i would use 6pn real4 rc kernel, but ok.)
  RN_DSP_ACCURACY qsum = 0.0;
  for (int32 i=0; i<interp_size; ++i)
    qsum += real(INTERP_KERNEL(0,i));
  INTERP_KERNEL /= qsum;// complex divided by float

  // ====== Variables for oversampling in azimuth [BK] ======
  // ______ Factor two oversampling is supported in azimuth ______
  // ______ We use a buffer of 6 lines (oversampled in range) ______
  // ______ these samples are interpolated in between using a ______
  // ______ raised cosine kernel of 6 points (claimed 0.9999) ______
  // ______ A big matrix is used that contains all shifted kernels ______
  const int32 NP_kernel_az = 6;
  matrix<complr4> BUFFER_AZ(NP_kernel_az,LINE_OUT.pixels());// rotating buffer
  matrix<complr4> KERNEL_AZ(NP_kernel_az,LINE_OUT.pixels());// raised cosine
  matrix<complr4> LINE_OUT2(1,LINE_OUT.pixels());// interpolated in azi.
  if (OsrAzimuth!=1)
    {
    INFO.print("Initializing azimuth kernel");
    const real4 CHI = imageinfo.prf/imageinfo.abw;// oversampling factor az
    matrix<real4> x_axis(NP_kernel_az,1);
    for (int32 i=0; i<NP_kernel_az; ++i)
      x_axis(i,0) = -NP_kernel_az/2 + 0.5 + i;// [-2.5 -1.5 -0.5 0.5 1.5 2.5]
    matrix<complr4> tmp_kernel = mat2cr4(rc_kernel(x_axis, CHI, NP_kernel_az));
    DEBUG.print("Normalizing kernel");
    real4 qsum = 0.0;
    for (int32 i=0; i<NP_kernel_az; ++i)
      qsum += real(tmp_kernel(i,0));
    tmp_kernel /= qsum;// complr4
    DEBUG.print("Shifting kernels with Doppler");
    // ___ Doppler centroid is function of range only ____
    // ___ to shift spectrum of convolution kernel to fDC of data, multiply
    // ___ in the space domain with a phase trend of -2pi*t*fdc/prf
    // ___ (to shift back (no need) you would use +fdc), see manual;
    for (int32 x=0; x<LINE_OUT.pixels(); ++x)
      {
      const real4 pix   = real4(imageinfo.currentwindow.pixlo)+real4(x)/2.0;
      const real4 slope = 2.0*PI*imageinfo.pix2fdc(pix)/imageinfo.prf;
      for (int32 i=0; i<NP_kernel_az; ++i)
        {
        // ___ Modify kernel, shift spectrum to fDC ___
        const real4 t  = x_axis(i,0)*slope;// the phase ramp
        KERNEL_AZ(i,x) = tmp_kernel(i,0)*complr4(cos(t),-sin(t));// note '-' (see manual)
        }
      }
    }//if azimuth ovs


  // ====== Loop on the lines to oversample ======
  const int32 TEN        = 10;
  const int32 TENPERCENT = int32(floor((0.5*TEN+numlines)/TEN));  // number of lines
  int32 percentage       = 0;                                     // initialization
  for (register int32 linecnt=0; linecnt<numlines; linecnt++)
    {
    if (!(linecnt%TENPERCENT))
      {
      PROGRESS << "OVERSAMPLESLC: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += TEN;
      }
    // ______ Read input data in larger line ______
    switch (imageinfo.formatflag)
      {
      case FORMATCR4:
        {
        matrix<real4> bufferrreal4(1,1);
        matrix<real4> bufferrimag4(1,1);
        for (int32 ii=0; ii<numpixels; ++ii)
          {
          ifile >> bufferrreal4;
          ifile >> bufferrimag4;
          // ______Generate a zero filled copy of LINE______
          // RN LINE_IN must be cleaned!!!
          LINE_IN(0,OsrRange*ii) = RN_DSP_CPXACCURACY(bufferrreal4(0,0),bufferrimag4(0,0));
          }
        break;
        }
      // ______ Convert first to ci2 before writing to file ______
      case FORMATCI2:
        {
        matrix<int16> bufferrealint16(1,1);
        matrix<int16> bufferimagint16(1,1);
        for (int32 ii=0; ii<numpixels; ++ii)
          {
          ifile >> bufferrealint16;
          ifile >> bufferimagint16;
          // ______Generate a zero filled copy of LINE______
          // RN LINE_IN must be cleaned!!!
          LINE_IN(0,OsrRange*ii) = RN_DSP_CPXACCURACY(bufferrealint16(0,0),bufferimagint16(0,0));
          }
        break;
        }
      default:
        PRINT_ERROR("Unknown input format for the cropped image.");
        throw(unhandled_case_error);
      }// end switch reading input line

    // ______Spatial convolution between LINE_IN and INTERP_KERNEL______
    int jmin, jmax, RN_k, minpos, maxpos;
    RN_k   = 0;
    minpos = (interp_size-1)/2;
    maxpos = (interp_size-1)/2 + (LINE_IN.pixels() - 1) + OsrRange - 1;
    for (int ii=minpos; ii<=maxpos; ii++)
      {
      LINE_OUT(0,RN_k) = 0;
      jmin = max(int32(0), int32(ii-interp_size+1));
      jmax = min(int32(ii),int32(LINE_IN.pixels()-1));
      for (int j=jmin; j<=jmax; j++)
        LINE_OUT(0,RN_k) += LINE_IN(0,j) * INTERP_KERNEL(0,ii-j);
      RN_k++;
      }

    // ______ Oversample in azimuth ______
    // ______ e.g., kernel is 6 points, buffer is 6 lines ______
    // ______ if linecnt==5, then buffer is filled with ______
    // ______ first 6 lines of file.  LINE_OUT2 is the ______
    // ______ interpolated line at BUFFER[2.5,*], i.e., ______
    // ______ to write in correct order each time we write ______
    // ______ BUFFER[2,*] and LINE_OUT2[*];  not 5 and 5.5 ______
    // ______ To correct for this offset, we do: ______
    // ______ linecnt==0: do not write (add at end)
    // ______ linecnt==1: do not write (add at end)
    // ______ linecnt==2: do not write (add at end)
    // ______ linecnt==3: write 0 and 0.5
    // ______ linecnt==4: write 1 and 1.5
    // ______ linecnt==5: write 2 and 2.5
    // ______ linecnt==6: write etc.
    // ______ linecnt==last: write last + 6 extra lines.
    if (OsrAzimuth!=1)//i.e., 2
      {
      // ______ Rotating BUFFER_AZ ______
      // ______ I dont trust shifting pointers so I copy lines ______
      // ______ this is slower, but guarantees continuous in mem. ______
      for (int32 x=0; x<LINE_OUT.pixels(); ++x)
        for (int32 i=0; i<NP_kernel_az-1; ++i)
          BUFFER_AZ(i,x) = BUFFER_AZ(i+1,x);//rotate buffer by copy
      for (int32 x=0; x<LINE_OUT.pixels(); ++x)
        BUFFER_AZ(NP_kernel_az-1,x) = complr4(LINE_OUT(0,x));// add new line
      // ______ Oversample in azimuth (interpolate) ______
      LINE_OUT2 = sum(dotmult(BUFFER_AZ,KERNEL_AZ),1);// at half+0.5
      }

    // ______ Write LINE_OUT in the output file ______
    switch (oversampleinput.oformatflag)
      {
      // ______ Convert first to cr4 before writing to file ______
      case FORMATCR4:
        {
        complr4 buffercr4;
        if (OsrAzimuth!=1)
          {
          if (linecnt>=NP_kernel_az/2)
            {
            // _____ write line 2.0 ______
            for (int ii=0; ii<LINE_OUT.pixels(); ii++)
              {
              buffercr4 = complr4(BUFFER_AZ(NP_kernel_az/2-1,ii));
              ofile.write((char*)&buffercr4,sizeof(complr4));
              }
            // _____ write line 2.5 ______
            for (int ii=0; ii<LINE_OUT.pixels(); ii++)
              {
              buffercr4 = complr4(LINE_OUT2(0,ii));
              ofile.write((char*)&buffercr4,sizeof(complr4));
              }
            // _____ write additional zero lines at end ______
            if (linecnt==numlines-1)//write extra lines at end
              {
              buffercr4 = complr4(0.0,0.0);
              for (int i=0; i<NP_kernel_az; ++i)
                for (int ii=0; ii<LINE_OUT.pixels(); ii++)
                  ofile.write((char*)&buffercr4,sizeof(complr4));
              }
            }
          }
        else // no azimuth oversampling
          {
          for (int ii=0; ii<LINE_OUT.pixels(); ii++)
            {
            buffercr4 = complr4(LINE_OUT(0,ii));
            ofile.write((char*)&buffercr4,sizeof(complr4));
            }
          }
        break;//switch
        }
      // ______ Convert first to ci2 before writing to file ______
      case FORMATCI2:
        {
        compli16 bufferci16;
        if (OsrAzimuth!=1)
          {
          if (linecnt>=NP_kernel_az/2)
            {
            // _____ write line 2.0 ______
            for (int ii=0; ii<LINE_OUT.pixels(); ii++)
              {
              bufferci16 = cr4toci2(complr4(BUFFER_AZ(NP_kernel_az/2-1,ii)));
              ofile.write((char*)&bufferci16,sizeof(compli16));
              }
            // _____ write line 2.5 ______
            for (int ii=0; ii<LINE_OUT.pixels(); ii++)
              {
              bufferci16 = cr4toci2(complr4(LINE_OUT2(0,ii)));
              ofile.write((char*)&bufferci16,sizeof(compli16));
              }
            // _____ write extra line at end ______
            if (linecnt==numlines-1)// write extra lines at end
              {
              bufferci16 = compli16(0,0);
              for (int i=0; i<NP_kernel_az; ++i)
                for (int ii=0; ii<LINE_OUT.pixels(); ii++)
                  ofile.write((char*)&bufferci16,sizeof(compli16));
              }
            }
          }
        else // no azimuth oversampling
          {
          for (int ii=0; ii<LINE_OUT.pixels(); ii++)
            {
            bufferci16 = cr4toci2(complr4(LINE_OUT(0,ii)));
            ofile.write((char*)&bufferci16,sizeof(compli16));
            }
          }
        break;
        }
        default:
          PRINT_ERROR("Unknown output format for the oversampled image.");
          throw(unhandled_case_error);
      } // end switch writing output
    LINE_IN.clean();
    } // end for loop over lines


  // ______ Close files ______
  #undef RN_DSP_ACCURACY
  #undef RN_DSP_CPXACCURACY
  ifile.close();
  ofile.close();


  // ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchoversample", ios::out | ios::trunc);
  bk_assert(scratchresfile,"OversampleSLC: scratchoversample",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
  if (fileid == MASTERID)
  {
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_oversample];
    scratchresfile << "\t\t\t" <<  "master";
  }
  if (fileid == SLAVEID)
    {
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_oversample];
    scratchresfile << "\t\t\t" <<  "slave";
    }
  // ______ updateslcimage greps these ______
  // ______ [BK] write new file size, etc. here ______
  scratchresfile
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  outfile
    << "\nData_output_format: \t\t\t\t";
    if (oversampleinput.oformatflag==FORMATCR4)
      scratchresfile << "complex_real4";
    if (oversampleinput.oformatflag==FORMATCI2)
      scratchresfile << "complex_short";
  scratchresfile
    << "\nFirst_line (w.r.t. ovs_image):       \t\t"
    <<  (imageinfo.currentwindow.linelo-1)*OsrAzimuth+1
    << "\nLast_line (w.r.t. ovs_image):        \t\t"
    <<  imageinfo.currentwindow.linehi*OsrAzimuth
    << "\nFirst_pixel (w.r.t. ovs_image):      \t\t"
    <<  (imageinfo.currentwindow.pixlo-1)*OsrRange+1
    << "\nLast_pixel (w.r.t. ovs_image):       \t\t"
    <<  imageinfo.currentwindow.pixhi*OsrRange
    << "\nMultilookfactor_azimuth_direction:   \t\t"
    << 1.0/OsrAzimuth
    << "\nMultilookfactor_range_direction:     \t\t"
    << 1.0/OsrRange
  // ______ updateslcimage does not grep these ______
    << "\nNumber of lines (oversampled):       \t\t"
    << OsrAzimuth*numlines
    << "\nNumber of pixels (oversampled):      \t\t"
    << OsrRange*numpixels
    << "\n#First_line (w.r.t. original_image): \t\t"
    <<  imageinfo.currentwindow.linelo
    << "\n#Last_line (w.r.t. original_image):  \t\t"
    <<  imageinfo.currentwindow.linehi
    << "\n#First_pixel (w.r.t. original_image):\t\t"
    <<  imageinfo.currentwindow.pixlo
    << "\n#Last_pixel (w.r.t. original_image): \t\t"
    <<  imageinfo.currentwindow.pixhi;
  scratchresfile
    << "\n*******************************************************************";
  if (fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_oversample] << "_NORMAL";
  if (fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_oversample] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();

  PROGRESS.print("OVERSAMPLESLC: 100%");
  } // END OversampleSLC
//____RaffaeleNutricato END MODIFICATION SECTION 2

// Copy code from writeslc
// do some modification
// Modified by LG for reading ALOS Fine
// To read complex real4 data
void palsar_fine_dump_data(
        const input_gen &generalinput,
        const input_crop &crop_arg,
        const int       checklines)
{

  const int16           sizeb4 = 4,             // some constants for reading
                        sizei4 = 4,             //   binary fields.
                        sizei6 = 6,
                        sizei8 = 8;
  uint                  lenrec1;                // length of general record1
  uint                  lenrec2;                // (nominal) length of data records
  char                  c4[5],                  // correctly 5 for \0
                        c6[7],                  // correctly 7 for \0
                        c8[9];                  // correctly 9 for \0

  // ______ Write some info ______
  TRACE_FUNCTION("palsar_fine_dump_data (LG 28-Dec-2005)")
  PROGRESS.print("Start cropping slc data.");

  // ______ Open files ______
  ifstream datfile;
  openfstream(datfile,crop_arg.filein1);
  bk_assert(datfile,crop_arg.filein1,__FILE__,__LINE__);

  // ====== Get data such as recordlength ======
  datfile.seekg(8,ios::beg);
  datfile.read((char*)&lenrec1,sizeb4);         // length of record1
  lenrec1 = ntohl(lenrec1);     // bk 6 jul 2000, byteorder x86 machines.
  DEBUG.print("record 1 of data file.");
  if (lenrec1 != 720 )  // PALSAR data = 720 + Rec_bytes*lines
    {
    WARNING << "palsar_fine_dump_data : length of record 1 = \""
         <<  lenrec1 << "\"; expected \"720\" for PALSAR FINE SLC (CEOS, full scene).";
    WARNING.print();

    }

  datfile.seekg(180,ios::beg);
  datfile.read((char*)&c6,sizei6);              // number of SAR DATA records (lines)
    c6[6]='\0';
    const uint numdatarec = atoi(c6);
  DEBUG << "numdatarec: " << numdatarec;
  DEBUG.print();
  //datfile.seekg(186,ios::beg);
  datfile.read((char*)&c6,sizei6);              // SAR DATA record length
    c6[6]='\0';
    const uint lendatarec2 = atoi(c6);
  DEBUG << "lendatarec2: " << lendatarec2;
  DEBUG.print();
  datfile.seekg(232,ios::beg);                  // SAR Related data
  datfile.read((char*)&c4,4);
    c4[4]='\0';
    const uint numchannels = atoi(c4);
  DEBUG << "numchannels: " << numchannels;
  DEBUG.print();

  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    uint numlines = atoi(c8);
  DEBUG << "numlines: " << numlines;
  DEBUG.print();

  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint leftborder = atoi(c4);
  DEBUG << "leftborder: " << leftborder;
  DEBUG.print();
  datfile.read((char*)&c8,sizei8);              // number of pixels
    c8[8]='\0';
    uint numpixels = atoi(c8);
  DEBUG << "numpixels: " << numpixels;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint rightborder = atoi(c4);
  DEBUG << "rightborder: " << rightborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint topborder = atoi(c4);
  DEBUG << "topborder: " << topborder;
  DEBUG.print();
  datfile.read((char*)&c4,sizei4);
    c4[4]='\0';
    const uint bottomborder = atoi(c4);
  DEBUG << "bottomborder: " << bottomborder;
  DEBUG.print();

  datfile.seekg(280,ios::beg);                  // Record data
  datfile.read((char*)&c8,sizei8);
    c8[8]='\0';
    const uint numbytesdata = atoi(c8);
  DEBUG << "numbytesdata: " << numbytesdata;
  DEBUG.print();


// ====== Check with volumefile / internal ======
// It seems that the lines (N+1) get from volume file is wrong,
  // it seems to be the pixle number in one line
  if (numlines != checklines)
    {
    WARNING << "code 902: data file: "
         << crop_arg.filein1
         << " numlin=" << numlines
         << " vs. volume file: "
         << crop_arg.filein1
         << " numlin=" << checklines;
    WARNING.print();
    WARNING.print(" +this means data and volume file seem not to correspond.");
    }

// ______ Check with previous section ______
  if (numlines != numdatarec)
    {
    WARNING << "code 904: Number of lines seems not to be consistent in file: "
         << crop_arg.filein1 << " : " << numlines << " != " << numdatarec;
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if ((numbytesdata / 8) != numpixels)
    {
    WARNING << "code 904: Number of pixels seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << numpixels << " != " << (numbytesdata / 8);
    WARNING.print();
    WARNING.print(" +this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }


// ====== Start copy input to output (raw) format with buffer======
// ______ Check and process optional offset parameters______
// ______ Lcnlow is corner line, lcnhi is other corner, pcnlow, pixel coord. low etc.
  uint linestart  = 1;                                  // counters for loops
  uint lineend    = numlines;
  uint pixelstart = 1;
  uint pixelend   = numpixels;                          // only for resultfile

  if (crop_arg.dbow.linehi!=0 && crop_arg.dbow.linelo!=0 &&
      crop_arg.dbow.pixhi!=0 && crop_arg.dbow.pixlo!=0)
    {
    window tempdbow(crop_arg.dbow.linelo, crop_arg.dbow.linehi,
                    crop_arg.dbow.pixlo,  crop_arg.dbow.pixhi);
    if (crop_arg.dbow.linehi>numlines)
      {
      WARNING << "Specified input DBOW linehi > numlines: "
           << crop_arg.dbow.linehi << " > " << numlines
           << ". I set linehi = " << numlines;
      WARNING.print();
      tempdbow.linehi=numlines;
      }
    if (crop_arg.dbow.pixhi>numpixels)
      {
      WARNING << "Specified input DBOW pixhi > numpixels: "
           << crop_arg.dbow.pixhi << " > " << numpixels
           << ". I set pixhi = " << numpixels;
      WARNING.print();
      tempdbow.pixhi=numpixels;
      }
// ______ Only hi values are possibly adapted, low is a constant ______
    numlines   = tempdbow.linehi - crop_arg.dbow.linelo + 1;
    numpixels  = tempdbow.pixhi  - crop_arg.dbow.pixlo  + 1;

    linestart  = crop_arg.dbow.linelo;
    lineend    = tempdbow.linehi;
    pixelstart = crop_arg.dbow.pixlo;
    pixelend   = tempdbow.pixhi;                        // only for resultfile
    }

// ______ Note complex<short> not in ANSI c ______
  //
      matrix <real4>    LINE(1,2*numpixels);            // size of real4

// ====== Process requested lines ======
  ofstream datoutfile;
  openfstream(datoutfile,crop_arg.fileout1,generalinput.overwrit);
  bk_assert(datoutfile,crop_arg.fileout1,__FILE__,__LINE__);

  // ______ info on data, to avoid X86 problems ______
  // ______ according to CEOS specs, byte 413 is first complex pixel, etc. ______
  // 720 + 84124*linenum + 412
  datfile.seekg(lenrec1 + 412,ios::beg);

  matrix <real4> TMPREAL4(1,2);
  datfile >> TMPREAL4;          // read in first complex pixel for test

  real8 tmpmag = sqrt(
    real8(real4(ntohl(TMPREAL4(0,0)))*real4(ntohl(TMPREAL4(0,0)))) +
    real8(real4(ntohl(TMPREAL4(0,1)))*real4(ntohl(TMPREAL4(0,1)))));

  DEBUG << "First complex element in datafile: ("
       << real4(ntohl(TMPREAL4(0,0))) << ","
       << real4(ntohl(TMPREAL4(0,1)))
       << "); mag = " << tmpmag;
  DEBUG.print();
  if (tmpmag > 10000.)
    {
    WARNING.print(DEBUG.get_str());
    WARNING.print("this is a byteorder problem on X86? (use ntohs)");
    }
  DEBUG << "TEST: (realpart): " << TMPREAL4(0,0)
       << ", (imagpart): " << TMPREAL4(0,1);
  DEBUG.print();
  DEBUG << "TEST: htons(realpart): " << ntohl(TMPREAL4(0,0))
       << ", htons(imagpart): " << ntohl(TMPREAL4(0,1));
  DEBUG.print();
  DEBUG << "TEST: ntohs(realpart): " << ntohl(TMPREAL4(0,0))
       << ", ntohs(imagpart): " << ntohl(TMPREAL4(0,1));
  DEBUG.print();
  DEBUG << "TEST: short int(ntohs(realpart)): " << real4(ntohl(TMPREAL4(0,0)))
       << ", (imagpart): " << real4(ntohl(TMPREAL4(0,1)));
  DEBUG.print();


  // ====== perline is faster than perbuffer, less memory etc. BK1998 ======

  datfile.seekg(lenrec1+(linestart-1)*lendatarec2 + 8,ios::beg);
  datfile.read((char*)&lenrec2,sizeb4);         // length of first record
  lenrec2 = ntohl(lenrec2);     // bk 6 jul 2000, byteorder x86 machines.

  if (lenrec2 != lendatarec2)
    {
    ERROR << "code 904: Length of datarecords seems to be inconsistent in file: "
         << crop_arg.filein1 << ": "
         << lenrec2 << " != " << lendatarec2;
    WARNING.print(ERROR.get_str());
    ERROR.reset();
    }

  const int32 TEN        = 10;
  const int32 TENPERCENT = int32((.5*TEN+numlines)/TEN);        // number of lines
  int32 percentage       = 0;                                   // initialization
  const int32 tmpstart   = lenrec1+412-lendatarec2+(pixelstart-1)*8;    // sizeof=8
  char  pD,*pc;
  for (register int32 linecnt=linestart; linecnt<=lineend; linecnt++)
    {
    if (!((linecnt-linestart)%TENPERCENT))
      {
      PROGRESS << "WRITESLC: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += TEN;
      }
    datfile.seekg(tmpstart+linecnt*lendatarec2,ios::beg);
    datfile    >> LINE;
    // ______ LG 28 DEC 2005: swapbytes for X86 (intel) linux cpus ______
    #ifdef __X86PROCESSOR__
    for (int ii=0; ii<LINE.pixels(); ++ii)
        {
                pc = (char*)&LINE(0,ii);
                pD = *pc; *pc = *(pc+3); *(pc+3) = pD;
                pD = *(pc+1);  *(pc+1) = *(pc+2); *(pc+2) = pD;
    //  LINE(0,ii) = LINE(0,ii)));      // changed from htons 171100 BK
        }
    #endif
    datoutfile << LINE;
    }
  datfile.close();                                      // close files
  datoutfile.close();



// ====== Write results to scratchfile ======
  ofstream scratchresfile("scratchres2raw", ios::out | ios::trunc);
  bk_assert(scratchresfile,"writeslc: scratchres2raw",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************\n";
    //<< "\n*_Start_crop:\t\t\t"
    //<<  crop_arg.idcrop
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_m_crop];
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "*_Start_" << processcontrol[pr_s_crop];
  scratchresfile
    << "\t\t\t" <<  crop_arg.idcrop
    << "\n*******************************************************************"
    << "\nData_output_file: \t\t\t\t"
    <<  crop_arg.fileout1
    << "\nData_output_format: \t\t\t\t"
    << "complex_real4"

// ______ updateslcimage greps these ______
    << "\nFirst_line (w.r.t. original_image): \t\t"
    <<  linestart
    << "\nLast_line (w.r.t. original_image): \t\t"
    <<  lineend
    << "\nFirst_pixel (w.r.t. original_image): \t\t"
    <<  pixelstart
    << "\nLast_pixel (w.r.t. original_image): \t\t"
    <<  pixelend
    << "\nNumber of lines (non-multilooked): \t\t" <<  lineend-linestart+1
    << "\nNumber of pixels (non-multilooked): \t\t" <<  pixelend-pixelstart+1
    << "\n*******************************************************************";
  if (crop_arg.fileid == MASTERID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_m_crop] << "_NORMAL";
  if (crop_arg.fileid == SLAVEID)
    scratchresfile <<  "\n* End_" << processcontrol[pr_s_crop] << "_NORMAL";
  scratchresfile
    << "\n*******************************************************************"
    <<  endl;
  scratchresfile.close();

// ______ Tidy up do checks here ______
  if (numchannels != 1)                         // ??
    {
    WARNING << "code 904: Number of channels in file: "
         << crop_arg.filein1 << " = "
         << numchannels << " != 1 ";
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  if (bottomborder != topborder != leftborder != rightborder != 0)
    {
    WARNING << "code 904: Not implemented: offset border: left,right,bottom,top: "
         << leftborder << "," << rightborder << "," << bottomborder << ","
         << topborder << " in file: " << crop_arg.filein1;
    WARNING.print();
    WARNING.print("this means SLC FORMAT IS DIFFERENT THEN EXPECTED.");
    }
  PROGRESS.print("WRITESLC: 100%");
} // end palsar_fine_dump_data

