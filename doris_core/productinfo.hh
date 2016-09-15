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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/productinfo.hh,v $        *
 * $Revision: 3.10 $                                            *
 * $Date: 2005/08/24 10:03:18 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * The productinfo class contains the definition of the data    *
 * and functions for 'products' (i.e. not slc images, but the   *
 * interferogram, coherence image, DEM etc.)                    *
 * Data mainly public because this used to be a struct and I    *
 * did not want to change the other code.                       *
 * It also consists of functions reading files etc.             *
 #%// BK 25-Aug-2000
 ****************************************************************/


#ifndef PRODUCTINFO_H
#define PRODUCTINFO_H

using namespace std;                    // BK 29-Mar-2003, new compiler?

// Jia defined this for compilation under windows
// Bert Kampes, 24-Aug-2005
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "constants.hh"                 // typedefs
#include <cstring>                      // strcpy, req. on some systems



// ====== Define template functions (no member no friend) ======
// ______ (matrix class is declared way below) ______
template <class Type> class matrix;




// ====== Struct slcimage: information on master/slave ======
class productinfo                               // info on 'products'
  {
  public:
    char          file[EIGHTY];                   // current filename
    // ______ window / multilook factors ______
    window        win;                            // current window, line(1:N) etc
    uint          multilookL;                     // multilookfactor in line (azi) dir.
    uint          multilookP;                     // multilookfactor in pixel (ra) dir.
    // ______ file format ______
    int16         formatflag;                     // current read formatflag


    // ______ Public function in struct ______
    // ______ constructor ______
    productinfo()               
      {
      formatflag = -1;// undefined
      multilookL =  1;
      multilookP =  1;
      } // rest ==0

    // ______ fill it from info in resultfiles ______
    void fillproductinfo(const char *file, const char *iden);

    // ______ assignment operator ______
    productinfo& operator = (productinfo X)
      {
      if (this != &X)
        {
        strcpy(file,X.file);
        win        = X.win;
        multilookL = X.multilookL;
        multilookP = X.multilookP;
        formatflag = X.formatflag;
        }
      return *this;
      };

    // ______ show content ______
    inline void showdata() const                  // show content
      {DEBUG << "\ncurrent file: \t" << file
             << "\nformatflag:   \t" << formatflag
             << "\nmultilook:    \t" << multilookL << " " << multilookP
             << "\nwindow:       \t" << win.linelo << " " << win.linehi
                              << " " << win.pixlo  << " " << win.pixhi;
       DEBUG.print();
      }

    // ______ read data from file ______
    matrix<real4> readphase(window win) const;

    // ______ read data from file ______
    matrix<complr4> readdata(window win) const;
    matrix<real4> readdatar4(window win) const; // [MA]

  }; // END class productinfo


#endif // PRODUCTINFO_H


