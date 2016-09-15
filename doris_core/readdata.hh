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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/readdata.hh,v $
 * $Revision: 3.13 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * routines for reading and logging volumefile,
 *  leaderfile, nullfile and datafile.
 * writing of data to raw format outputfile.
 ****************************************************************/


#ifndef READDATA_H
#define READDATA_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "constants.hh"                                 // typedefs
#include "readinput.hh"                                 // input structs

// ______prototypes______
void readvolume(
        input_readfiles &readfiles_arg,
        const char      *chk1,
        const char      *chk2,
        const char      *chk3);

void readleader(
        input_readfiles &readfiles_arg,
        const int32      check);                        // process leader file

void readnull(
        const input_readfiles &readfiles_arg);          // process null file

void readdat(
        input_readfiles &readfiles_arg,
        const int32      check);

// Modified by LG for reading ALOS Fine
void palsar_fine_dump_data(
                           const input_gen &generalinput,               
                           const input_crop &writeslc_arg,
                           const int32   check);                        

void writeslc(
        const input_gen &generalinput,                  // mem/ overwrite
        const input_crop &writeslc_arg,
        const int32      check);                        // process data file

void envisat_dump_data(
        const input_crop &writeslc_arg);

void envisat_dump_VV(
        const input_crop &writeslc_arg);

void envisat_dump_HH(
        const input_crop &writeslc_arg);

void tsx_dump_data(
                   const input_crop &writeslc_arg);

void rs2_dump_data(
                   const input_crop &writeslc_arg);

void csk_dump_data(
                   const input_crop &writeslc_arg);

void radarsat_dump_data(
        const input_gen &generalinput,                  // mem/ overwrite
        const input_crop &writeslc_arg);

// BO.20100917
void gammaprocessor_crop(
        const input_gen &generalinput,                  // mem/ overwrite
	const slcimage 	 &master,
        const input_crop &writeslc_arg);

void  OversampleSLC(
       const input_gen        &generalinput,
       const slcimage         &imageinfo,
       const input_oversample &oversampleinput,
       const int16            fileid);

#endif // READDATA_H
