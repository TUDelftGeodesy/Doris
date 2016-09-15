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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/products.hh,v $
 * $Revision: 3.9 $
 * $Date: 2006/05/18 10:03:18 $
 * $Author: kampes $
 *
 * Declaration of routines for computation of products (interferogram etc.)
 ****************************************************************/


#ifndef PRODUCTS_H                      // guard
#define PRODUCTS_H

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


// ______ Simulate amplitude ______ [MA]
void sim_amplitude(
        const input_gen        &generalinput,
        const input_ell        &ellips,
        const input_simamp     &simamp,
        const slcimage         &master,
        orbit                  &masterorbit);

// ______ Compute the (complex) interferogram ______
void compinterfero(
        const slcimage          &master,
        const slcimage          &slave,
        const input_gen         &input_general,
        const input_interfero   &input_i_interfero);

// ______ Subtract ref.pha from complex interferogram ______
// ______ evaluate polynomial from comprefpha ______
void subtrrefpha(
        const slcimage          &master,
        const productinfo       &interferogram,
        const input_gen         &input_general,
        const input_subtrrefpha &input_i_subtrrefpha,
        const matrix<real8>     &coeff_refpha,
        const matrix<real8>     &coeff_h2ph); // added by FvL
// ______ evaluate ref. pha foreach pixel ______
void subtrrefpha(
        const input_ell         &ellips,
        const slcimage          &master,
        const slcimage          &slave,
        const productinfo       &interferogram,
        const input_gen         &input_general,
        const input_subtrrefpha &input_i_subtrrefpha,
              orbit             &masterorbit,
              orbit             &slaveorbit);


// ______ Compute the (complex) coherence image ______
void compcoherence(
        const slcimage          &master,
        const slcimage          &slave,
        const input_gen         &input_general,
        const input_coherence   &input_i_coherence,
        const matrix<real8>     &coeff_flatearth);

// ______ Compute the (complex) coherence image [NEW method] ______ added by DON
void compcoherence(
        const slcimage          &master,
        const slcimage          &slave,
        const productinfo       &radarcoderefdem,
        const input_gen         &input_general,
        const input_coherence   &input_i_coherence,
        const matrix<real8>     &coeff_flatearth);

// ______ Subtract ref.dem from complex interferogram ______
void subtrrefdem(
        const productinfo       &interferogram,         // filename, size info
        const productinfo       &radarcoderefdem,       // filename, size info
        const input_gen         &input_general,         // overwrite, memory
        const input_subtrrefdem &input_i_subtrrefdem);


// ______ Scale defo-interferogram with ratio of Bperp ______
void dinsar(
        const input_gen         &input_general,         // overwrite, memory
        const input_dinsar     &dinsarinput,
        const input_ell         &ellips,
        const slcimage          &master,
              orbit             &masterorbit,
        const slcimage          &defoslave,
              orbit             &defoorbit,
        const productinfo       &defointerferogram);


#endif // PRODUCTS_H



