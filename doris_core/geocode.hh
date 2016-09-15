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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/geocode.hh,v $
 * $Revision: 3.9 $
 * $Date: 2005/08/24 10:03:18 $
 * $Author: kampes $
 *
 * Declaration of routines for computation of endproducts (DEM, defo.map, )
 ****************************************************************/


#ifndef GEOCODE_H                       // guard
#define GEOCODE_H

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
#include "bk_baseline.hh"               // my 'baseline' class




// ______ Use schwabisch approx. method ______
void slant2hschwabisch(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage     &master,
        const slcimage     &slave,
        const productinfo     &interferogram,
        orbit               &masterorbit,
        orbit               &slaveorbit);


// ______ Use method ramon, derivative ______
void slant2hambiguity(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage     &master,
        const slcimage     &slave,
        const productinfo     &interferogram,
        orbit               &masterorbit,
        orbit               &slaveorbit, 
        const BASELINE      &baseline);


// ______ Use standard method? ______
void slant2hrodriguez(
        const input_gen     &generalinput,
        const input_slant2h &slant2hinput,
        const input_ell     &ellips,
        const slcimage     &master,
        const slcimage     &slave,
        const productinfo     &interferogram,
        const matrix<real8> &coeff_flatearth,
        orbit               &masterorbit,
        orbit               &slaveorbit,
        const BASELINE      &baseline);


// ______ Geocode after s2h ______
void geocode(
        const input_gen     &generalinput,
        const input_geocode &geocodeinput,
        const input_ell     &ellips,
        const slcimage     &master,
        const productinfo     &interferogram,
        orbit               &masterorbit);


#endif // GEOCODE_H



