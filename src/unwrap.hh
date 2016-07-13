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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/unwrap.hh,v $
 * $Revision: 3.8 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Declaration of routines for unwrapping 
 ****************************************************************/



#ifndef UNWRAP_H                        // guard
#define UNWRAP_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "constants.hh"                 // typedefs
#include "readinput.hh"                 // input structs
#include "orbitbk.hh"                   // my orbit class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class



// ______ Use stanford software through unix calls _____
void unwraptreeframon(
        const input_gen     &generalinput,
        const input_unwrap  &unwrapinput,
        const productinfo   &interferogram);

// ______ Use snaphu software through unix calls _____
void snaphu_unwrap(
        const input_gen     &generalinput,
        const input_unwrap  &unwrapinput,
        const productinfo   &interferogram,
        const slcimage      &master,
        const slcimage      &slave,
              orbit         &masterorbit,
              orbit         &slaveorbit,
        const input_ell     &ellips);

#endif // UNWRAP_H



