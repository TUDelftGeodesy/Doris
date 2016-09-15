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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/products.cc,v $   *
 * $Revision: 3.21 $                                            *
 * $Date: 2006/05/18 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * implementation of computation of products.                   *
 * -complex interferogram                                       *
 * -interferogram (minus reference fase)                        *
 * -magnitude image -> coherence image                          *
 * -differential insar.                                         *
 ****************************************************************/


#include "constants.hh"                 // global constants
#include "matrixbk.hh"                  // my matrix class
#include "orbitbk.hh"                   // my orbit class
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class
#include "products.hh"                  // header file
#include "utilities.hh"                 // utils
#include "ioroutines.hh"                // 
#include "exceptions.hh"                // my exceptions class
#include "bk_baseline.hh"               // my baseline class
#include "referencephase.hh"            // getcorners, griddatalinear for simamp ; // getcorners --> should go to utilities.hh I guess [MA] TODO

#include <algorithm>                    // generic max

#ifdef WIN32
// Jia defined this here, i did this in constants.hh
// Bert Kampes, 24-Aug-2005
//#define max _MAX
//#define min _MIN
#endif

#ifdef WIN32 // [MA] following added due to ntohl requirments in master_timing simulation
  // Jia defined this.
  // Bert Kampes, 24-Aug-2005
  #include "winsock2.h"
#else
  #include <netinet/in.h>                 // ntohl byteorder x86-HP unix
#endif


/****************************************************************
 *    compinterfero                                             *
 *                                                              *
 * Compute products:                                            *
 *  - (compex) interferogram                                    *
 * note: master-slave                                           *
 * Assumed that slave.currentwin is in master coord. system     *
 * and is smaller than or equal to maste.currentwin.            *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, filenames                                *
 * Output:                                                      *
 *  - files on disk                                             *
 *    Bert Kampes, 07-Apr-1999                                  *
 *                                                              *
 * bugfix computations, subtract reference phase                *
 * for all points before multilooking.                          *
 *    Bert Kampes, 06-Oct-1999                                  *
 *                                                              *
 * Change log: - FvL, 23-Jul-2014, removed unwanted automatic   *
 *   removal of pre-calculated reference phase                  *
 ****************************************************************/
void compinterfero(
        const slcimage         &master,
        const slcimage         &slave,
        const input_gen        &input_general,
        const input_interfero  &input_i_interfero)
  {
  TRACE_FUNCTION("compinterfero (BK 06-Oct-1999)");
  INFO << "INTERFERO: master input file: " << master.file;
  INFO.print();
  INFO << "INTERFERO: slave input file:  " << slave.file;
  INFO.print();

#ifdef __DEBUG
  // ______This should be checked before, impossible______
  if (slave.currentwindow.linelo < master.currentwindow.linelo ||
      slave.currentwindow.linehi > master.currentwindow.linehi ||
      slave.currentwindow.pixlo  < master.currentwindow.pixlo  ||
      slave.currentwindow.pixhi  > master.currentwindow.pixhi    )
    {
    slave.currentwindow.disp();
    master.currentwindow.disp();
    slave.currentwindow.disp();
    PRINT_ERROR("Panic, interferogram.win smaller than master win.")
    throw(file_error);
    }
#endif

  INFO << "INTERFERO: interferogram (master coord.): "
       << slave.currentwindow.linelo << ":" << slave.currentwindow.linehi << "; "
       << slave.currentwindow.pixlo  << ":" << slave.currentwindow.pixhi;
  INFO.print();
  
  // ___ Check what to do ___
  //const uint  BUFFERMEMSIZE = input_general.memory;
  const real8  BUFFERMEMSIZE = input_general.memory;
  const int32 multiL        = input_i_interfero.multilookL;
  const int32 multiP        = input_i_interfero.multilookP;
  bool nocint               = true;                     // output complex phase image
  bool noint                = true;                     // no output real phase image
  if (specified(input_i_interfero.focint))
    nocint  = false;
  if (specified(input_i_interfero.foint))
    noint  = false;

  // ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;

  // ====== Open output files ======
  ofstream ofilecint;
  if (!nocint)
    {
    openfstream(ofilecint,input_i_interfero.focint,input_general.overwrit);
    bk_assert(ofilecint,input_i_interfero.focint,__FILE__,__LINE__);
    }
  ofstream ofileint;
  if (!noint)
    {
    openfstream(ofileint,input_i_interfero.foint,input_general.overwrit);
    bk_assert(ofileint,input_i_interfero.foint,__FILE__,__LINE__);
    }


  // ====== allocate matrices ======
  const int32 numpixels      = (slave.currentwindow.pixels());
  const int32 bytesperline   = numpixels * sizeof(complr4);
  const real4 numbigmatrices = 3.2;             // M, S
  //int32 numlines             = int32((BUFFERMEMSIZE/numbigmatrices)/bytesperline);// lines in buffer
  int32 numlines             = int32(ceil( (BUFFERMEMSIZE/numbigmatrices) / bytesperline )); // lines in buffer [MA] see reference.cc BUFFERMEMSIZE
  while (numlines%multiL)                               // correct numlines to multiple of multiL
    numlines -= 1;
  if (numlines%multiL)
    {
    PRINT_ERROR("panic : totally impossible on HP aCC compiler.")
    throw(input_error);
    }
  int32 nummllines  = numlines/multiL;                  // exact
  int32 nummlpixels = numpixels/multiP;                 // floor...
  if (nummllines < 1)
    {
    PRINT_ERROR("Please increase memory (MEMORY card) or decrease multiL (INT_MULTILOOK card).")
    throw(input_error);
    }
  DEBUG << "Number of lines per buffer: " << numlines; 
  DEBUG.print();

  // ______ Number of lines on file ______
  // Ifile contains lines of resampled slave ?
  //const uint Mfilelines = master.currentwindow.lines();
  const uint Ifilelines = slave.currentwindow.lines();

  // ______ Window in master system to load from file ______
  window winfile(slave.currentwindow.linelo,
                 slave.currentwindow.linelo + numlines - 1,
                 slave.currentwindow.pixlo,
                 slave.currentwindow.pixhi);


  // ====== Loop over all totally filled buffers ======
  // ______ Mis-use Master to store complex interferogram (also multilooked)
  //register int32 i,blocks;
  const int32 numfullbuffers = Ifilelines/numlines;     // fully filled buffers
  const int32 numrestlines   = Ifilelines%numlines;     // restlines total
  const int32 nummlrestlines = numrestlines/multiL;     // floor...
  const int32 EXTRABUFFER    = nummlrestlines ? 1 : 0;

  for (int32 blocks=1; blocks<=numfullbuffers+EXTRABUFFER; blocks++)
    {
    // ______ Progress info ______
    PROGRESS << "INTERFERO: "
         << int32(100*real4(blocks-1)/(real4(Ifilelines)/real4(numlines))+0.5)
         << "%";
    PROGRESS.print();
  
    // ______ Check last block ______
    if (blocks == (numfullbuffers+1))           // there is an extra (smaller) block
      {
      numlines       = multiL * nummlrestlines;
      winfile.linehi = numlines + winfile.linelo - 1;
      }

    // ______ Fill buffers master/slave from disk ______
    matrix<complr4> MASTER = master.readdata(winfile);
    matrix<complr4> SLAVE  = slave.readdata(winfile);

    // ______ Compute M = M* conj(S) ______
    MASTER *= conj(SLAVE);              // ?better SLAVE.conj(); for speed and memory

    // ====== Multilook if appropriate ======
    matrix<complr4> ML = multilook(MASTER,multiL,multiP);

    // ====== Write (partial) products to file ======  
    if (!nocint)                        // complex interferogram
      ofilecint << ML;
    if (!noint)                         // phase image
      ofileint << fast_angle(ML);
  
    // ______ Update window from file ______
    winfile.linelo = winfile.linehi + 1;
    winfile.linehi += numlines;
    } // loop blocks

  // ====== Write results to file ======
  ofstream scratchlogfile("scratchloginterfero", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"compinterfero: scratchloginterfero",__FILE__,__LINE__);

  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* COMPUTATION OF INTERFEROGRAM"
    << "\n*******************************************************************"
    << "\nInput file master (format): \t\t\t"
    <<  master.file << " " << master.formatflag
    << "\nInput file slave (format): \t\t\t"
    <<  slave.file << " " << slave.formatflag
    << "\ncomplex interferogram: \t\t"
    <<  input_i_interfero.focint
    << "\ninterferogram: \t\t"
    <<  input_i_interfero.foint
    << "\nNumber of lines (multilooked): \t"
    <<  Ifilelines / multiL
    << "\nNumber of pixels (multilooked): \t"
    <<  nummlpixels
    << "\nMultilookfactor in line direction: \t"
    <<  multiL
    << "\nMultilookfactor in pixel direction: \t"
    <<  multiP
    << "\nTip: dismph " << input_i_interfero.focint 
    << " " << nummlpixels << " 1 500"
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchresinterfero", ios::out | ios::trunc);
  bk_assert(scratchresfile,"compinterfero: scratchresinterfero",__FILE__,__LINE__);

  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_interfero]
                 << "\n*******************************************************************";
  if (!nocint)
    {
    scratchresfile 
                 << "\nData_output_file: \t\t\t"
                 <<  input_i_interfero.focint
                 << "\nData_output_format: \t\t\t"
                 << "complex_real4";
    }
  if (!noint)
    {
    scratchresfile 
                 << "\nData_output_file_real_interferogram: \t"
                 <<  input_i_interfero.foint
                 << "\nData_output_format_real_interferogram: \t\t"
                 << "real4";
    }
  scratchresfile 
                 << "\nFirst_line (w.r.t. original_master): \t"
                 <<  slave.currentwindow.linelo
                 << "\nLast_line (w.r.t. original_master): \t"
                 <<  slave.currentwindow.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixlo
                 << "\nLast_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixhi
                 << "\nMultilookfactor_azimuth_direction: \t"
                 <<  multiL
                 << "\nMultilookfactor_range_direction: \t"
                 <<  multiP
                 << "\nNumber of lines (multilooked): \t\t"
                 <<  Ifilelines / multiL
                 << "\nNumber of pixels (multilooked): \t"
                 <<  nummlpixels
                 << "\n*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_interfero] << "_NORMAL"
                 << "\n*******************************************************************\n";
  scratchresfile.close();

  // ______Tidy up______
  PROGRESS.print("finished compinterfero.");
  if (!nocint)
    ofilecint.close();
  if (!noint)
    ofileint.close();
  } // END compinterfero



/****************************************************************
 *    compcoherence                                             *
 *                                                              *
 * Compute products:                                            *
 *  - (compex) coherence.                                       *
 * note: master-slave-ref                                       *
 * Assumed that slave.currentwin is in master coord. system     *
 * and is smaller than or equal to maste.currentwin.            *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, filenames                                *
 * Output:                                                      *
 *  - files on disk                                             *
 *                                                              *
 *    Bert Kampes, 16-Apr-1999                                  *
 ****************************************************************/
void compcoherence(
        const slcimage        &master,
        const slcimage        &slave,
        const input_gen       &input_general,
        const input_coherence &input_i_coherence,
        const matrix<real8>   &coeff_flatearth)
  {
  TRACE_FUNCTION("compcoherence (BK 16-Apr-1999)")
#ifdef __DEBUG
  // ______This should be checked before, impossible______
  if (slave.currentwindow.linelo < master.currentwindow.linelo ||
      slave.currentwindow.linehi > master.currentwindow.linehi ||
      slave.currentwindow.pixlo  < master.currentwindow.pixlo  ||
      slave.currentwindow.pixhi  > master.currentwindow.pixhi    )
    {
    PRINT_ERROR("Panic, impossible 3333.")
    throw(input_error);
    }
#endif

  bool nocoh  = true;                                   // no output real coherence
  bool noccoh = true;                                   // no output complex coherence
  bool noflatearthcorrection = false;                   // do correction
  string METHOD = "REFPHASE_ONLY";                      // updated when necessary
  if (specified(input_i_coherence.focoh))
    nocoh   = false;
  if (specified(input_i_coherence.foccoh))
    noccoh  = false;
  if (coeff_flatearth.size() == 0)                      // step flatearth not done
    {                                                   //  (not in result file)
    WARNING.print("Computation of coherence without subtracting of reference image.");
    noflatearthcorrection = true;
    METHOD = "NO_REFPHASE_REMOVED";
    }

// ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;
  INFO << "compcoherence: polynomial normalized by factors: "
       << minL << " " << maxL << " " << minP << " " << maxP
       << " to [-2,2]";
  INFO.print();

// ______ Some other variables ______
  //const uint  BUFFERMEMSIZE = input_general.memory;
  const real8  BUFFERMEMSIZE = input_general.memory;
  const int32 multiL    = input_i_coherence.multilookL;
  const int32 multiP    = input_i_coherence.multilookP;
  //const int32 multiLP         = multiL*multiP;
  const int32 winsizeL  = input_i_coherence.cohsizeL;
  const int32 winsizeP  = input_i_coherence.cohsizeP;

// ====== Open output files ======
  ofstream ofileccoh;
  if (!noccoh)
    {
    openfstream(ofileccoh,input_i_coherence.foccoh,input_general.overwrit);
    bk_assert(ofileccoh,input_i_coherence.foccoh,__FILE__,__LINE__);
    }
  ofstream ofilecoh;
  if (!nocoh)
    {
    openfstream(ofilecoh,input_i_coherence.focoh,input_general.overwrit);
    bk_assert(ofilecoh,input_i_coherence.focoh,__FILE__,__LINE__);
    }


// ====== allocate matrices ======
  const int32 numpixels = (slave.currentwindow.pixhi-slave.currentwindow.pixlo+1);
  const int32 bytesperline = numpixels * sizeof(complr4);
  int32 numlines;                                       // lines in buffer

// ______ estimate numlines for correct memory usage ______
  if (noflatearthcorrection)
    //numlines = int32((BUFFERMEMSIZE / 3.0) / bytesperline);
    numlines = int32(ceil( (BUFFERMEMSIZE / 3.0) / bytesperline ));  // [MA] BUFFERMEMSIZE as real8
  else
    numlines = int32(ceil( (BUFFERMEMSIZE / 3.5) / bytesperline ));

  while (numlines%multiL)       // correct numlines to multiple of multiL
    numlines -= 1;
  if (numlines%multiL)
    {
    PRINT_ERROR("PANIC: multilookwindow not exactly in numlines.")
    throw(input_error);
    }

  int32 nummllines  = numlines/multiL;  // exact zero fraction
  int32 nummlpixels = numpixels/multiP; // i.e., floor()
  if (nummllines < 1)
    {
    PRINT_ERROR("increase memory (MEMORY card) or reduce multiL (COH_MULTILOOK card).")
    throw(input_error);
    }

// ______ adapt numlines to extra for window ______
  int32 numlines2 = numlines;           // number of lines output of coherence
  numlines += (winsizeL - 1);           // number of lines from file

// ______ Number of lines on file ______
  //const uint Mfilelines = master.currentwindow.lines();
  const uint Ifilelines = slave.currentwindow.lines();


// ______ Window in master system to load from file ______
//    uint firstline = slave.currentwindow.linelo;
//    uint lastline  = firstline + numlines - 1;
// do little different cause of firstblock
  //register int32 i,j;
  uint  zerolinesstart = (winsizeL-1)/2;// overlap previous block?
  uint  trailinglines  = (winsizeL)/2;//   overlap next block?
  int32 extrazerolines = 0;// larger last block

/* 
// ____ ? why is this here? it seems not required ___
// ____ ? e.g., coh win=64, ml win=10 fails if this is done
// ____ ? but I guess there was a situation, e.g., when 
// ____ ? the slave is resampled larger than the master?
// ____ ? or when the ML window is much larger than the coh win.
// ____ ? BK, sep.2004: commented out.

  uint  writeallzero   = zerolinesstart/multiL;// ??
  int32 zerorestlines  = zerolinesstart%multiL;// ??
  INFO << "coherence: starting with zero lines: " << writeallzero;
  INFO.print();
  if (writeallzero != 0)
    {
    if (!noccoh)
      {
      matrix<complr4> zeroline(1,nummlpixels);
      for (i=0; i<writeallzero; i++)
        ofileccoh << zeroline;
      }
    if (!nocoh)
      {
      matrix<real4> zeroline(1,nummlpixels);
      for (i=0; i<writeallzero; i++)
        ofilecoh << zeroline;
      }
    }
*/

  const int32 BUF= 1+winsizeL/multiL;           // lines
  window winfile(slave.currentwindow.linelo,
                 slave.currentwindow.linelo + BUF*multiL + trailinglines - 1,
                 slave.currentwindow.pixlo,
                 slave.currentwindow.pixhi);

// ______ Buffers ______
  matrix<complr4> SLAVE(winfile.linehi-winfile.linelo+1,numpixels);
  matrix<complr4> MASTER;

// ______ Axis for reference phase ______
  matrix<real4> p_axis;
  if (!noflatearthcorrection)
    {
    p_axis.resize(numpixels,1);
    for (int32 i=0; i<numpixels; i++)
      p_axis(i,0) = winfile.pixlo + i;
    // ______ Normalize data ______
    normalize(p_axis,minP,maxP);
    }


// ====== Loop over all totally filled buffers ======
// ______ Misuse Slave to store complex interferogram (also multilooked)
  bool firstblock = true;
  bool lastblock  = false;
  for (register int32 blocks=0; blocks<100000000; blocks++)     // for ever
    {
    INFO << "coherence block: " << blocks;
    INFO.print();
    // ====== Initialize for (smaller) last block, if approriate ======
    if (lastblock)                                      // lastblock already done 
      break;
    if (winfile.linehi > slave.currentwindow.linehi)    // more then there is on file
      {
      lastblock = true;
      // ______ find out if execution of this buffer is required ______
      int32 restlines = Ifilelines%multiL;              // extra lines on file
      if (winfile.linelo+zerolinesstart >= slave.currentwindow.linehi-restlines+1)
        break;

      // ______ Resize matrices for last loop (1) ______
      // ______ (not that well programmed...) ______
      int32 lastlinerequired = slave.currentwindow.linehi - restlines + trailinglines;
      if (lastlinerequired > int32(slave.currentwindow.linehi))
        {
        winfile.linehi = slave.currentwindow.linehi;
        extrazerolines = lastlinerequired - slave.currentwindow.linehi;
        }
      else
        {
        winfile.linehi = lastlinerequired;
        }
      INFO << "coherence: extra zero lines last buffer: " << extrazerolines;
      INFO.print();
      MASTER.resize(winfile.linehi-winfile.linelo+1,numpixels);
      SLAVE.resize (winfile.linehi-winfile.linelo+1,numpixels);
      }

    else // not lastblock: give approx PROGRESS
      {
      // ______ Progress info: not totally accurate ______
      PROGRESS << "COHERENCE: "
           << int32(100.0*real4(blocks)/(real4(Ifilelines)/real4(numlines2))+0.5)
           << "%";
      PROGRESS.print();
      }

// ______ Fill buffers master/slave from disk ______
      INFO << "coherence winfile: [" << winfile.linelo << ":" << winfile.linehi
           << ", " << winfile.pixlo << ":" << winfile.pixhi << "]";
      INFO.print();
      MASTER = master.readdata(winfile);
      SLAVE  = slave.readdata(winfile);

// ______ Add zero lines if firstblock ______
    if (firstblock == true)
      {
      matrix<complr4> TMP = SLAVE;
      SLAVE.resize(BUF*multiL+winsizeL-1,TMP.pixels());
      const window wintmp(SLAVE.lines()-TMP.lines(), SLAVE.lines()-1,
                    0, SLAVE.pixels()-1);
      const window windef(0, 0, 0, 0);
      SLAVE.setdata(wintmp,TMP,windef);
      TMP = MASTER;
      MASTER.resize(BUF*multiL+winsizeL-1,TMP.pixels());
      MASTER.setdata(wintmp,TMP,windef);
      INFO << "coherence: increase lines for first block: "
           << SLAVE.lines()-TMP.lines();
      INFO.print();
      }

// ______ Resize matrices for last loop (2) ______
// ______ (not that well programmed...) ______
    if (extrazerolines != 0)
      {
      const window wintmp(0, MASTER.lines()-1, 0, MASTER.pixels()-1);
      const window windef(0, 0, 0, 0);
      matrix<complr4> TMP = SLAVE;
      SLAVE.resize(TMP.lines()+extrazerolines,TMP.pixels());
      SLAVE.setdata(wintmp,TMP,windef);
      TMP = MASTER;
      MASTER.resize(TMP.lines()+extrazerolines,TMP.pixels());
      MASTER.setdata(wintmp,TMP,windef);
      INFO << "coherence: increase lines for last block: "
           << SLAVE.lines()-TMP.lines();
      INFO.print();
      }

    // ______ Compute reference phase ______
    uint sizeL = SLAVE.lines();
    uint sizeP = SLAVE.pixels();
    if (!noflatearthcorrection)
      {
      matrix<real4> l_axis(sizeL,1);
      for (uint i=0; i<sizeL; i++)
        l_axis(i,0) = winfile.linelo + i;
      // ______ Normalize data ______
      normalize(l_axis,minL,maxL);
      matrix<real4> REFPHASE = polyval<real4>(l_axis, p_axis, coeff_flatearth);
      // ______ Complex interferogram in master, norms in slave ______
      for (uint i=0; i<sizeL; i++)
        {
        for (uint j=0; j<sizeP; j++)
          {
          real4 tmp    = norm(MASTER(i,j));
          MASTER(i,j) *= (conj(SLAVE(i,j)) *            // store cmplx interf.
                         complr4(fast_cos(REFPHASE(i,j)),fast_min_sin(REFPHASE(i,j))));
          SLAVE(i,j)   = complr4(norm(SLAVE(i,j)),tmp);         // store norm
          }
        }
      }

    // ______ Complex interferogram in master, norms in slave ______
    else        // no reference phase available
      {
      for (uint i=0; i<sizeL; i++)
        {
        for (uint j=0; j<sizeP; j++)
          {
          real4 tmp    = norm(MASTER(i,j));
          MASTER(i,j) *= conj(SLAVE(i,j));
          SLAVE(i,j)   = complr4(norm(SLAVE(i,j)),tmp);         // store norms
          }
        }
      }

// ______ Compute (complex) coherence ______
// ______ multilook and write to file ______
    INFO << "block: " << blocks << "; num output lines: "
         << real8(SLAVE.lines()-winsizeL+1.0)/real8(multiL);// should be exact x.0
    INFO.print();
    INFO << "SLAVE number of lines this buffer block: " << SLAVE.lines();
    INFO.print();
    if (!noccoh)        // complex coherence requested
      {
      if (!nocoh)       // and coherence
        {
        matrix<complr4> CCOHERENCE = multilook(
                coherence(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
        ofileccoh << CCOHERENCE;
        ofilecoh  << magnitude(CCOHERENCE);
        }
      else              // thus only complex coherence
        {
        ofileccoh << multilook(
                coherence(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
        }
      }

    else if (!nocoh)    // thus only (real) coherence requested
      {
      //ofilecoh  << multilook(
      //coherence2(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
      // test: see size of blocks by using tmp matrix...
      matrix<real4> COHERENCE = 
                coherence2(MASTER,SLAVE,winsizeL,winsizeP);
      INFO << "block: " << blocks << "; size COHERENCE matrix: "
           << COHERENCE.lines();// multiple of multiL ??
      INFO.print();
      ofilecoh  << multilook(COHERENCE, multiL,multiP);
      }

    // ______ impossible request, is checked before. ______
    else
      {
      PRINT_ERROR("2212 panic impossible.")
      throw(unhandled_case_error);
      }


// ______ update windows/matrix size ______
    winfile.linelo = winfile.linehi - trailinglines + 1 - zerolinesstart;
    winfile.linehi = winfile.linelo + numlines - 1;

    if (firstblock)
      {
      firstblock = false;// next time
      MASTER.resize(winfile.linehi-winfile.linelo+1,numpixels);
      SLAVE.resize(winfile.linehi-winfile.linelo+1,numpixels);
      }
    } // loop over all blocks



// ====== Write results to file ======
  ofstream scratchlogfile("scratchlogcoherence", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"compcoherence: scratchlogcoherence",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* COMPUTATION OF COHERENCE"
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\t"  << "REFPHASE_ONLY"
    << "\nInput file master (format):     \t"
    <<  master.file << " " << master.formatflag
    << "\nInput file slave (format):      \t"
    <<  slave.file << " " << slave.formatflag
    << "\ncomplex coherence:              \t"
    <<  input_i_coherence.foccoh
    << "\ncoherence:                      \t"
    <<  input_i_coherence.focoh
    << "\nNumber of lines (multilooked):  \t"
    <<  Ifilelines / multiL
    << "\nNumber of pixels (multilooked): \t"
    <<  nummlpixels
    << "\nMultilookfactor in line direction: \t"
    <<  multiL
    << "\nMultilookfactor in pixel direction: \t"
    <<  multiP
    << "\nNumber of lines window for coherence estimation: " 
    <<  winsizeL
    << "\nNumber of pixels window for coherence estimation: " 
    <<  winsizeP
    << "\nNumber of ml lines per buffer during computation: "
    <<  BUF
    << "\nTip: disfloat " << input_i_coherence.focoh
    << " " << nummlpixels << " 1 500"
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchrescoherence", ios::out | ios::trunc);
  bk_assert(scratchresfile,"compcoherence: scratchrescoherence",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 //<< "\n*_Start_coherence"
                 << "\n*_Start_" << processcontrol[pr_i_coherence]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\t\t"  << METHOD;
  if (!nocoh)
    {
    scratchresfile 
                 << "\nData_output_file: \t\t\t"
                 <<  input_i_coherence.focoh
                 << "\nData_output_format: \t\t\t"
                 << "real4";
    }
  if (!noccoh)
    {
    scratchresfile 
                 << "\nData_output_file_complex_coherence: "
                 <<  input_i_coherence.foccoh
                 << "\nData_output_format_complex_coherence: \t\t"
                 << "complex_real4";
    }
  scratchresfile 
                 << "\nFirst_line (w.r.t. original_master): \t" // updateproductinfo greps these
                 <<  slave.currentwindow.linelo
                 << "\nLast_line (w.r.t. original_master): \t"
                 <<  slave.currentwindow.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixlo
                 << "\nLast_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixhi
                 << "\nMultilookfactor_azimuth_direction: \t"
                 <<  multiL
                 << "\nMultilookfactor_range_direction: \t"
                 <<  multiP
                 << "\nNumber of lines (multilooked): \t\t"
                 <<  Ifilelines / multiL
                 << "\nNumber of pixels (multilooked): \t"
                 <<  nummlpixels
                 << "\n*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_coherence] << "_NORMAL"
                 << "\n*******************************************************************\n";
  scratchresfile.close();

// ______Tidy up______
  PROGRESS.print("finished compcoherence.");
  if (!noccoh)
    ofileccoh.close();
  if (!nocoh)
    ofilecoh.close();

  } // END compcoherence


/****************************************************************
 *    compcoherence                                             *
 *                                                              *
 * Compute products:                                            *
 *  - (complex) coherence.                                      *
 * Note: master-slave-ref.phase-demphase                        *
 * Assumptions:                                                 *
 *  - slave.currentwin is in master coord. system and           *
 *    is smaller than or equal to maste.currentwin.             *
 *  - DEM ref.phase is not multilooked.                         *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, filenames                                *
 * Output:                                                      *
 *  - files on disk                                             *
 *                                                              *
 *    Bert Kampes, 16-Apr-1999                                  *
 *    Davide Oscar Nitti, 14-Nov-2008 (removal of topo slope)   *
 *    Mahmut Arikan,      09-Jan-2009 (nonML refdem update  )   *
 ****************************************************************/
void compcoherence(
        const slcimage        &master,
        const slcimage        &slave,
        const productinfo     &radarcodedrefdem,
        const input_gen       &input_general,
        const input_coherence &input_i_coherence,
        const matrix<real8>   &coeff_flatearth)
  {
  TRACE_FUNCTION("compcoherence (don 14-Nov-2008)")
#ifdef __DEBUG
  // ______This should be checked before, impossible______
  if (slave.currentwindow.linelo < master.currentwindow.linelo ||
      slave.currentwindow.linehi > master.currentwindow.linehi ||
      slave.currentwindow.pixlo  < master.currentwindow.pixlo  ||
      slave.currentwindow.pixhi  > master.currentwindow.pixhi    )
    {
    PRINT_ERROR("Panic, impossible 3333.")
    throw(input_error);
    }
#endif

  // _____ start added by MA _____
  string refdemfilenoML = string(radarcodedrefdem.file) + ".noML" ; // filename for no-mlooked refdem phase
  const uint &refdemmlL   = radarcodedrefdem.multilookL;
  const uint &refdemmlP   = radarcodedrefdem.multilookP;

  if ( !existed(refdemfilenoML.c_str()) &&  refdemmlL != 1)    // [MA] multilooking of the refdem must match with the master and slave. 
    {
    ERROR << "Coherence: missing non-multilook "
          << refdemfilenoML << " file."; 
    PRINT_ERROR(ERROR.get_str())
    ERROR.print("and MultilookfactorL ref. dem not equal to 1.");
    throw(input_error);
    }
  else if (!existed(refdemfilenoML.c_str()) && refdemmlP != 1)
    {
    ERROR << "Coherence: missing non-multilook "
          << refdemfilenoML << " file."; 
    PRINT_ERROR(ERROR.get_str())
    ERROR.print("and MultilookfactorP ref. dem not equal to 1.");
    throw(input_error);
    }
  // _____ end added by MA _____

  bool nocoh  = true;                                   // no output real coherence
  bool noccoh = true;                                   // no output complex coherence
  bool noflatearthcorrection = false;                   // do correction
  char checkrefdemIncludeFE[20];
  readres(checkrefdemIncludeFE,sizeof(checkrefdemIncludeFE),
          input_general.i_resfile,"Include_flatearth:");
  if (specified(input_i_coherence.focoh))
    nocoh   = false;
  if (specified(input_i_coherence.foccoh))
    noccoh  = false;
  if (!strcmp(checkrefdemIncludeFE,"Yes"))
    {
    noflatearthcorrection = true; // flatearth phase included in refdem_phase
    }
  else
    {
    if (coeff_flatearth.size() == 0)                    // step flatearth not done
      {                                                 //  (not in result file)
      noflatearthcorrection = true;
      PRINT_ERROR("step flatearth not done")
      throw(input_error);
      }
    }

// ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;
  INFO << "compcoherence: polynomial normalized by factors: "
       << minL << " " << maxL << " " << minP << " " << maxP
       << " to [-2,2]";
  INFO.print();

// ______ Some other variables ______
  //const uint  BUFFERMEMSIZE = input_general.memory;
  const real8  BUFFERMEMSIZE = input_general.memory;
  const int32 multiL    = input_i_coherence.multilookL;
  const int32 multiP    = input_i_coherence.multilookP;
  //const int32 multiLP         = multiL*multiP;
  const int32 winsizeL  = input_i_coherence.cohsizeL;
  const int32 winsizeP  = input_i_coherence.cohsizeP;

// ====== Open output files ======
  ofstream ofileccoh;
  if (!noccoh)
    {
    openfstream(ofileccoh,input_i_coherence.foccoh,input_general.overwrit);
    bk_assert(ofileccoh,input_i_coherence.foccoh,__FILE__,__LINE__);
    }
  ofstream ofilecoh;
  if (!nocoh)
    {
    openfstream(ofilecoh,input_i_coherence.focoh,input_general.overwrit);
    bk_assert(ofilecoh,input_i_coherence.focoh,__FILE__,__LINE__);
    }


// ====== allocate matrices ======
  const int32 numpixels = (slave.currentwindow.pixhi-slave.currentwindow.pixlo+1);
  const int32 bytesperline = numpixels * sizeof(complr4);
  int32 numlines;                                       // lines in buffer

// ______ estimate numlines for correct memory usage ______
  if (noflatearthcorrection)
    //numlines = int32((BUFFERMEMSIZE / 3.0) / bytesperline);
    numlines = int32(ceil( (BUFFERMEMSIZE / 3.0) / bytesperline ));
  else
    numlines = int32(ceil( (BUFFERMEMSIZE / 3.5) / bytesperline ));

  while (numlines%multiL)       // correct numlines to multiple of multiL
    numlines -= 1;
  if (numlines%multiL)
    {
    PRINT_ERROR("PANIC: multilookwindow not exactly in numlines.")
    throw(input_error);
    }

  int32 nummllines  = numlines/multiL;  // exact zero fraction
  int32 nummlpixels = numpixels/multiP; // i.e., floor()
  if (nummllines < 1)
    {
    PRINT_ERROR("increase memory (MEMORY card) or reduce multiL (COH_MULTILOOK card).")
    throw(input_error);
    }

// ______ adapt numlines to extra for window ______
  int32 numlines2 = numlines;           // number of lines output of coherence
  numlines += (winsizeL - 1);           // number of lines from file

// ______ Number of lines on file ______
  //const uint Mfilelines = master.currentwindow.lines();
  const uint Ifilelines = slave.currentwindow.lines();


// ______ Window in master system to load from file ______
//    uint firstline = slave.currentwindow.linelo;
//    uint lastline  = firstline + numlines - 1;
// do little different cause of firstblock
  //register int32 i,j;
  uint  zerolinesstart = (winsizeL-1)/2;// overlap previous block?
  uint  trailinglines  = (winsizeL)/2;//   overlap next block?
  int32 extrazerolines = 0;// larger last block

/* 
// ____ ? why is this here? it seems not required ___
// ____ ? e.g., coh win=64, ml win=10 fails if this is done
// ____ ? but I guess there was a situation, e.g., when 
// ____ ? the slave is resampled larger than the master?
// ____ ? or when the ML window is much larger than the coh win.
// ____ ? BK, sep.2004: commented out.

  uint  writeallzero   = zerolinesstart/multiL;// ??
  int32 zerorestlines  = zerolinesstart%multiL;// ??
  INFO << "coherence: starting with zero lines: " << writeallzero;
  INFO.print();
  if (writeallzero != 0)
    {
    if (!noccoh)
      {
      matrix<complr4> zeroline(1,nummlpixels);
      for (i=0; i<writeallzero; i++)
        ofileccoh << zeroline;
      }
    if (!nocoh)
      {
      matrix<real4> zeroline(1,nummlpixels);
      for (i=0; i<writeallzero; i++)
        ofilecoh << zeroline;
      }
    }
*/

  const int32 BUF= 1+winsizeL/multiL;           // lines
  window winfile(slave.currentwindow.linelo,
                 slave.currentwindow.linelo + BUF*multiL + trailinglines - 1,
                 slave.currentwindow.pixlo,
                 slave.currentwindow.pixhi);

// ______ Buffers ______
  matrix<complr4> SLAVE(winfile.linehi-winfile.linelo+1,numpixels);
  matrix<real4> REFDEMPHA(winfile.linehi-winfile.linelo+1,numpixels);
  matrix<complr4> MASTER;

// ______ Axis for reference phase ______
  matrix<real4> p_axis;
  if (!noflatearthcorrection)
    {
    p_axis.resize(numpixels,1);
    for (int32 i=0; i<numpixels; i++)
      p_axis(i,0) = winfile.pixlo + i;
    // ______ Normalize data ______
    normalize(p_axis,minP,maxP);
    }


// ====== Loop over all totally filled buffers ======
// ______ Misuse Slave to store complex interferogram (also multilooked)
  bool firstblock = true;
  bool lastblock  = false;
  for (register int32 blocks=0; blocks<100000000; blocks++)     // for ever
    {
    INFO << "coherence block: " << blocks;
    INFO.print();
    // ====== Initialize for (smaller) last block, if approriate ======
    if (lastblock)                                      // lastblock already done 
      break;
    if (winfile.linehi > slave.currentwindow.linehi)    // more then there is on file
      {
      lastblock = true;
      // ______ find out if execution of this buffer is required ______
      int32 restlines = Ifilelines%multiL;              // extra lines on file
      if (winfile.linelo+zerolinesstart >= slave.currentwindow.linehi-restlines+1)
        break;

      // ______ Resize matrices for last loop (1) ______
      // ______ (not that well programmed...) ______
      int32 lastlinerequired = slave.currentwindow.linehi - restlines + trailinglines;
      if (lastlinerequired > int32(slave.currentwindow.linehi))
        {
        winfile.linehi = slave.currentwindow.linehi;
        extrazerolines = lastlinerequired - slave.currentwindow.linehi;
        }
      else
        {
        winfile.linehi = lastlinerequired;
        }
      INFO << "coherence: extra zero lines last buffer: " << extrazerolines;
      INFO.print();
      MASTER.resize(winfile.linehi-winfile.linelo+1,numpixels);
      SLAVE.resize (winfile.linehi-winfile.linelo+1,numpixels);
      REFDEMPHA.resize (winfile.linehi-winfile.linelo+1,numpixels);
      }

    else // not lastblock: give approx PROGRESS
      {
      // ______ Progress info: not totally accurate ______
      PROGRESS << "COHERENCE: "
           << int32(100.0*real4(blocks)/(real4(Ifilelines)/real4(numlines2))+0.5)
           << "%";
      PROGRESS.print();
      }

// ______ Fill buffers master/slave from disk ______
      INFO << "coherence winfile: [" << winfile.linelo << ":" << winfile.linehi
           << ", " << winfile.pixlo << ":" << winfile.pixhi << "]";
      INFO.print();
      MASTER = master.readdata(winfile);
      SLAVE  = slave.readdata(winfile);

      // _____ start added by MA _____
      if ( refdemmlL == 1 && refdemmlP == 1  ) // [MA] give precedence to regular refdem phase if both files exists
                                               // use the regular file since it's not multilooked. 
        {
         INFO << "topo refphase file: [" << radarcodedrefdem.file << "is used." ;
         INFO.print();
         readfile(REFDEMPHA,radarcodedrefdem.file,
                  slave.currentwindow.lines(),winfile,slave.currentwindow);
        }
      else if ( existed(refdemfilenoML.c_str()) )                                 // use mas_sla.demphase.noML file
        {
         INFO << "topo refphase file: [" << refdemfilenoML << "is used." ;
         INFO.print();
         readfile(REFDEMPHA,refdemfilenoML.c_str(),
                  slave.currentwindow.lines(),winfile,slave.currentwindow);
        }
      // _____ end added by MA _____


// ______ Add zero lines if firstblock ______
    if (firstblock == true)
      {
      matrix<complr4> TMP = SLAVE;
      SLAVE.resize(BUF*multiL+winsizeL-1,TMP.pixels());
      const window wintmp(SLAVE.lines()-TMP.lines(), SLAVE.lines()-1,
                    0, SLAVE.pixels()-1);
      const window windef(0, 0, 0, 0);
      SLAVE.setdata(wintmp,TMP,windef);
      TMP = MASTER;
      MASTER.resize(BUF*multiL+winsizeL-1,TMP.pixels());
      MASTER.setdata(wintmp,TMP,windef);
      INFO << "coherence: increase lines for first block: "
           << SLAVE.lines()-TMP.lines();
      INFO.print();
      matrix<real4> TMPr4 = REFDEMPHA;
      REFDEMPHA.resize(BUF*multiL+winsizeL-1,TMPr4.pixels());
      REFDEMPHA.setdata(wintmp,TMPr4,windef);
      }

// ______ Resize matrices for last loop (2) ______
// ______ (not that well programmed...) ______
    if (extrazerolines != 0)
      {
      const window wintmp(0, MASTER.lines()-1, 0, MASTER.pixels()-1);
      const window windef(0, 0, 0, 0);
      matrix<complr4> TMP = SLAVE;
      SLAVE.resize(TMP.lines()+extrazerolines,TMP.pixels());
      SLAVE.setdata(wintmp,TMP,windef);
      TMP = MASTER;
      MASTER.resize(TMP.lines()+extrazerolines,TMP.pixels());
      MASTER.setdata(wintmp,TMP,windef);
      INFO << "coherence: increase lines for last block: "
           << SLAVE.lines()-TMP.lines();
      INFO.print();
      matrix<real4> TMPr4 = REFDEMPHA;
      REFDEMPHA.resize(TMPr4.lines()+extrazerolines,TMPr4.pixels());
      REFDEMPHA.setdata(wintmp,TMPr4,windef);
      }

    // ______ Compute reference phase ______
    uint sizeL = SLAVE.lines();
    uint sizeP = SLAVE.pixels();
    if (!noflatearthcorrection)
      {
      matrix<real4> l_axis(sizeL,1);
      for (uint i=0; i<sizeL; i++)
        l_axis(i,0) = winfile.linelo + i;
      // ______ Normalize data ______
      normalize(l_axis,minL,maxL);
      matrix<real4> REFPHASE = polyval<real4>(l_axis, p_axis, coeff_flatearth);
      // ______ Complex interferogram in master, norms in slave ______
      for (uint i=0; i<sizeL; i++)
        {
        for (uint j=0; j<sizeP; j++)
          {
          real4 tmp    = norm(MASTER(i,j));
          MASTER(i,j) *= (conj(SLAVE(i,j)) *            // store cmplx interf.
                         complr4(fast_cos(REFPHASE(i,j)),fast_min_sin(REFPHASE(i,j))) *
                         complr4(fast_cos(REFDEMPHA(i,j)),fast_min_sin(REFDEMPHA(i,j))));
          SLAVE(i,j)   = complr4(norm(SLAVE(i,j)),tmp);         // store norm
          }
        }
      }

    // ______ Complex interferogram in master, norms in slave ______
    else        // reference phase included in refdem phase
      {
      for (uint i=0; i<sizeL; i++)
        {
        for (uint j=0; j<sizeP; j++)
          {
          real4 tmp    = norm(MASTER(i,j));
          MASTER(i,j) *= (conj(SLAVE(i,j)) * 
                         complr4(fast_cos(REFDEMPHA(i,j)),fast_min_sin(REFDEMPHA(i,j))));
          SLAVE(i,j)   = complr4(norm(SLAVE(i,j)),tmp);         // store norms
          }
        }
      }

// ______ Compute (complex) coherence ______
// ______ multilook and write to file ______
    INFO << "block: " << blocks << "; num output lines: "
         << real8(SLAVE.lines()-winsizeL+1.0)/real8(multiL);// should be exact x.0
    INFO.print();
    INFO << "SLAVE number of lines this buffer block: " << SLAVE.lines();
    INFO.print();
    if (!noccoh)        // complex coherence requested
      {
      if (!nocoh)       // and coherence
        {
        matrix<complr4> CCOHERENCE = multilook(
                coherence(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
        ofileccoh << CCOHERENCE;
        ofilecoh  << magnitude(CCOHERENCE);
        }
      else              // thus only complex coherence
        {
        ofileccoh << multilook(
                coherence(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
        }
      }

    else if (!nocoh)    // thus only (real) coherence requested
      {
      //ofilecoh  << multilook(
      //coherence2(MASTER,SLAVE,winsizeL,winsizeP),multiL,multiP);
      // test: see size of blocks by using tmp matrix...
      matrix<real4> COHERENCE = 
                coherence2(MASTER,SLAVE,winsizeL,winsizeP);
      INFO << "block: " << blocks << "; size COHERENCE matrix: "
           << COHERENCE.lines();// multiple of multiL ??
      INFO.print();
      ofilecoh  << multilook(COHERENCE, multiL,multiP);
      }

    // ______ impossible request, is checked before. ______
    else
      {
      PRINT_ERROR("2212 panic impossible.")
      throw(unhandled_case_error);
      }


// ______ update windows/matrix size ______
    winfile.linelo = winfile.linehi - trailinglines + 1 - zerolinesstart;
    winfile.linehi = winfile.linelo + numlines - 1;

    if (firstblock)
      {
      firstblock = false;// next time
      MASTER.resize(winfile.linehi-winfile.linelo+1,numpixels);
      SLAVE.resize(winfile.linehi-winfile.linelo+1,numpixels);
      REFDEMPHA.resize(winfile.linehi-winfile.linelo+1,numpixels);
      }
    } // loop over all blocks


// ====== Write results to file ======
  ofstream scratchlogfile("scratchlogcoherence", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"compcoherence: scratchlogcoherence",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* COMPUTATION OF COHERENCE"
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\t"  << "INCLUDE_REFDEM"
    << "\nInput file master (format):     \t"
    <<  master.file << " " << master.formatflag
    << "\nInput file slave (format):      \t"
    <<  slave.file << " " << slave.formatflag
    << "\ncomplex coherence:              \t"
    <<  input_i_coherence.foccoh
    << "\ncoherence:                      \t"
    <<  input_i_coherence.focoh
    << "\nNumber of lines (multilooked):  \t"
    <<  Ifilelines / multiL
    << "\nNumber of pixels (multilooked): \t"
    <<  nummlpixels
    << "\nMultilookfactor in line direction: \t"
    <<  multiL
    << "\nMultilookfactor in pixel direction: \t"
    <<  multiP
    << "\nNumber of lines window for coherence estimation: " 
    <<  winsizeL
    << "\nNumber of pixels window for coherence estimation: " 
    <<  winsizeP
    << "\nNumber of ml lines per buffer during computation: "
    <<  BUF
    << "\nTip: disfloat " << input_i_coherence.focoh
    << " " << nummlpixels << " 1 500"
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchrescoherence", ios::out | ios::trunc);
  bk_assert(scratchresfile,"compcoherence: scratchrescoherence",__FILE__,__LINE__);
  scratchresfile << "\n\n*******************************************************************"
                 << "\n*_Start_" << processcontrol[pr_i_coherence]
                 << "\n*******************************************************************"
                 << "\nMethod: \t\t\t\t"  << "INCLUDE_REFDEM";
  if (!nocoh)
    {
    scratchresfile 
                 << "\nData_output_file: \t\t\t"
                 <<  input_i_coherence.focoh
                 << "\nData_output_format: \t\t\t"
                 << "real4";
    }
  if (!noccoh)
    {
    scratchresfile 
                 << "\nData_output_file_complex_coherence: "
                 <<  input_i_coherence.foccoh
                 << "\nData_output_format_complex_coherence: \t\t"
                 << "complex_real4";
    }
  scratchresfile 
                 << "\nFirst_line (w.r.t. original_master): \t" // updateproductinfo greps these
                 <<  slave.currentwindow.linelo
                 << "\nLast_line (w.r.t. original_master): \t"
                 <<  slave.currentwindow.linehi
                 << "\nFirst_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixlo
                 << "\nLast_pixel (w.r.t. original_master): \t"
                 <<  slave.currentwindow.pixhi
                 << "\nMultilookfactor_azimuth_direction: \t"
                 <<  multiL
                 << "\nMultilookfactor_range_direction: \t"
                 <<  multiP
                 << "\nNumber of lines (multilooked): \t\t"
                 <<  Ifilelines / multiL
                 << "\nNumber of pixels (multilooked): \t"
                 <<  nummlpixels
                 << "\n*******************************************************************"
                 << "\n* End_" << processcontrol[pr_i_coherence] << "_NORMAL"
                 << "\n*******************************************************************\n";
  scratchresfile.close();

// ______Tidy up______
  PROGRESS.print("finished compcoherence.");
  if (!noccoh)
    ofileccoh.close();
  if (!nocoh)
    ofilecoh.close();

  } // END compcoherence


/****************************************************************
 *    subtrrefpha (polynomial)                                  *
 *                                                              *
 * Compute complex_interferogram .* conj(REF.PHA)               *
 * note: master-slave-ref.phase                                 *
 * ref.phase is 2d-polynomial model of ellipsoid phase          *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, multilook, filenames                     *
 * Output:                                                      *
 *  - new files on disk                                         *
 *    Bert Kampes, 09-Feb-2000                                  *
 ****************************************************************/
void subtrrefpha(
        const slcimage          &master,        // normalization factor polynomial
        const productinfo       &interferogram,
        const input_gen         &input_general,
        const input_subtrrefpha &input_i_subtrrefpha,
        const matrix<real8>     &coeff_flatearth,
        const matrix<real8>     &coeff_h2ph) // added by FvL
  {
  TRACE_FUNCTION("subtrrefpha (polynomial) (BK 09-Feb-2000)")
  // ====== Input options ======
  const int32 multiL    = input_i_subtrrefpha.multilookL;
  const int32 multiP    = input_i_subtrrefpha.multilookP;

  // ====== Get number of buffers etc. ======
  const int32 mldiskL   = interferogram.multilookL;     // cint on disk
  const int32 mldiskP   = interferogram.multilookP;     // cint on disk

  const int32 numlinesdisk    = interferogram.win.lines()/mldiskL;
  const int32 numpixelsdisk   = interferogram.win.pixels()/mldiskP;
  // ______ In original master radar coordinate system ______
  const real4 veryfirstline   = real4(interferogram.win.linelo) +
                                (real4(mldiskL) - 1.) / 2.;
  const real4 firstpixel      = real4(interferogram.win.pixlo)  +
                                (real4(mldiskP) - 1.) / 2.;

  const int32 totalmlL        = mldiskL*multiL;
  const int32 totalmlP        = mldiskP*multiP;
  const int32 numlinesoutput  = numlinesdisk/multiL;    // floor
  const int32 numpixelsoutput = numpixelsdisk/multiP;   // floor
 // const int32 lastlineoutput  = interferogram.win.linelo        // [MA] These were used be dumped to the res file.
 //                               + totalmlL*numlinesoutput - 1;  // Now, They are replaced by interferogram.win.linehi
 // const int32 lastpixeloutput = interferogram.win.pixlo         // and interferogram.win.pixhi.
 //                               + totalmlP*numpixelsoutput - 1; // see: comprefdem
  const bool outputh2ph       =  specified(input_i_subtrrefpha.foh2ph); // if spec. then output, added by FvL

// ______ Normalize data for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;
  INFO << "subtrrefpha: polynomial ref.phase normalized by factors: "
       << minL << " " << maxL << " " << minP << " " << maxP
       << " to [-2,2]";
  INFO.print();

// ====== Open files ======
// ______ If requested, dump ref.pha and do nothing else ______
  ifstream ifcint;                              // input file complex interf.
  ofstream ofilecint;                           // output file complex interf.
  ofstream ofrefpha;                            // output file ref. phase
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    openfstream(ofrefpha,input_i_subtrrefpha.forefpha,input_general.overwrit);
    bk_assert(ofrefpha,input_i_subtrrefpha.forefpha,__FILE__,__LINE__);
    }
  else  // compute cint and dump this...  
    {
    openfstream(ofilecint,input_i_subtrrefpha.focint,input_general.overwrit);
    bk_assert(ofilecint,input_i_subtrrefpha.focint,__FILE__,__LINE__);
    
    if (interferogram.formatflag != FORMATCR4)
      {
      PRINT_ERROR("code: .. Complex interferogram on disk must be complex real4.")
      throw(unhandled_case_error);
      }

    openfstream(ifcint,interferogram.file);
    bk_assert(ifcint,interferogram.file,__FILE__,__LINE__);
    }
  
//______________________________________________________________________
//__________ added by FvL ______________________________________________
//______________________________________________________________________   
  ofstream ofh2ph;
  if (outputh2ph==true)
    {
    openfstream(ofh2ph,input_i_subtrrefpha.foh2ph,input_general.overwrit);
    bk_assert(ofh2ph,input_i_subtrrefpha.foh2ph,__FILE__,__LINE__);
    }   
//______________________________________________________________________
//__________ end added by FvL ______________________________________________
//______________________________________________________________________   

// ====== Loop parameters ======
  //const uint  BUFFERMEMSIZE = input_general.memory;
  const real8  BUFFERMEMSIZE = input_general.memory;
  const int32 bytesperline = 2*sizeof(real4)*numpixelsdisk;


  const real4 nummatrices = 3;                          // CINT and REFPHA
  //int32 bufferlines = int32((BUFFERMEMSIZE / nummatrices) / bytesperline);        // lines in buffer
  int32 bufferlines = int32(ceil((BUFFERMEMSIZE / nummatrices) / bytesperline ));        // lines in buffer
  while (bufferlines%multiL)       // correct bufferlines to multiple of multiL
    {
    bufferlines -= 1;
    }
  DEBUG << "bufferlines per buffer: " << bufferlines;
  DEBUG.print();
  if (bufferlines%multiL != 0)
    {
    PRINT_ERROR("panic: refphase bufferlines totally impossible on HP.")
    throw(some_error);
    }

  const int32 numfullbuffers = numlinesdisk/bufferlines;        // floor
  int32 restlines            = numlinesdisk%bufferlines;
  while (restlines%multiL)            // correct bufferlines to multiple of multiL
    {
    restlines -= 1;
    }
  DEBUG << "restlines last buffer: " << restlines
       << "; and: " << numfullbuffers
       << " buffers of #lines: " << bufferlines;
  DEBUG.print();

  if (restlines%multiL != 0)
    {
    PRINT_ERROR("panic: i thought restlines is always exactly fitted this way??.")
    throw(some_error);
    }
  const int32 EXTRABUFFER = (restlines) ? 1 : 0;

  // ______ Axis for polyval ______
  register int32 i;
  matrix<real4> p_axis(numpixelsdisk,1);        // must be standing
  for (i=0; i<numpixelsdisk; i++)
    {
    p_axis(i,0) = firstpixel + i*mldiskP;       // grid on disk
    }
  normalize(p_axis,minP,maxP);                  // normalize

// ====== Loop over buffers ======
  for (register int32 buffer=0; buffer<numfullbuffers+EXTRABUFFER; ++buffer)
    {
    const real4 firstline = veryfirstline + buffer*bufferlines*mldiskL;
    if (buffer == numfullbuffers)               // i.e., last smaller buffer
      bufferlines = restlines;                  // ==multiple of multiL // [MA] this isn't nec. fixed in v4.02

    // ______ Progress information ______
    PROGRESS << "SUBTRREFPHA: buffer: "   << buffer+1   << "; "
                      << "line: " << firstline << " : "
                      << firstline+bufferlines-1;
    PROGRESS.print();

    // ====== Evaluate reference phase for (large) grid ======
    // ______ suspect that real4 is not enough for normalization.
    matrix<real4> l_axis(bufferlines,1);
    for (i=0; i<bufferlines; i++)
      l_axis(i,0) = firstline + i*mldiskL;      // grid on disk

    normalize(l_axis,minL,maxL);                // normalize data

    matrix<real8> REFPHASE = polyval<real8>(l_axis, p_axis, coeff_flatearth); // computation done in real8

    // ______ Dump this and continue if requested ______
    if (outputh2ph==true)                // dump h2ph, added by FvL
      {
      matrix<real8> H2PH = polyval<real8>(l_axis, p_axis, coeff_h2ph); // added by FvL

      matrix<real4> output_layer(H2PH.lines(),H2PH.pixels());
      convert_type(H2PH, output_layer); // ex: real8 --> real4
      ofh2ph << multilook(output_layer, multiL, multiP); // H2PH
      H2PH(1,1); //deallocate
      }
    if (input_i_subtrrefpha.dumponlyrefpha)     // dump complex ref.phase for checking
      {
      //matrix<complr4> REFPHASECR4 = fast_angle2cmplx(REFPHASE);
      //ofrefpha << REFPHASECR4;
      
      matrix<complr8> REFPHASECR8 = fast_angle2cmplx(REFPHASE);
      
      //matrix<real4> output_layer(REFPHASE.lines(),REFPHASE.pixels());
      
      matrix<complr4> output_layer(REFPHASECR8.lines(),REFPHASECR8.pixels());
      convert_type(REFPHASECR8, output_layer); // [MA] real8 --> real4

      ofrefpha << multilook(output_layer, multiL, multiP); // unwrapped REFPHASE 
      continue;                                 // next buffer
      }

    // ====== Read in buffer of complex interferogram ======
    matrix<complr4> CINT(bufferlines,numpixelsdisk);
    ifcint >> CINT;             // readin: assumed format complr4

    // ====== Subtract phase by complex multiplication (Euler) ======
    // CINT *= conj(angle2cmplx(REFPHASE));     // pointwise multiplication with conj.
    // ______ changes CINT ______
    fast_dotmultconjphase(CINT,REFPHASE);               // pointwise multiplication with conj.

    // ====== Write multilooked output to new file ======
    ofilecint << multilook(CINT, multiL, multiP);
    } // buffer loop
  ifcint.close();
  ofilecint.close();

  if (input_i_subtrrefpha.dumponlyrefpha)
    ofrefpha.close();
  if (outputh2ph==true) // added by FvL
    ofh2ph.close();


// ====== Log info/results ======
  ofstream scratchlogfile("scratchlogsubtrrefpha", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"subtrrefpha: scratchlogsubtrrefpha",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* SUBTRREFPHA" 
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\t\tpolynomial";
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    scratchlogfile
      << "\nOnly dump of reference phase, no subtraction."
      << "\nData_output_file_ref.phase: \t\t\t"
      <<  input_i_subtrrefpha.forefpha;
    }
  else
    {
    scratchlogfile
      << "\nInput_file_complex_interferogram: \t\t"
      <<  interferogram.file
      << "\nData_output_file_complex_interferogram: \t"
      <<  input_i_subtrrefpha.focint
      << "\nData_output_format_complex_interferogram: \t"
      << "complex_real4";
    }
  scratchlogfile
    << "\nFirst_line (w.r.t. original_master): \t\t" // updateproductinfo greps these
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master): \t\t"
    <<  interferogram.win.linehi                            //lastlineoutput
    << "\nFirst_pixel (w.r.t. original_master): \t\t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master): \t\t"
    <<  interferogram.win.pixhi                             //lastpixeloutput
    << "\nMultilookfactor_azimuth_direction: \t\t"
    <<  totalmlL
    << "\nMultilookfactor_range_direction: \t\t"
    <<  totalmlP
    << "\nNumber of lines (multilooked): \t\t\t"
    <<  numlinesoutput
    << "\nNumber of pixels (multilooked): \t\t"
    <<  numpixelsoutput
    << "\n\nMultilookfactors input complex interferogram: \t"
    <<  mldiskL << " " << mldiskP
    << "\nMultilookfactors requested in this step: \t"
    <<  multiL << " " << multiP
    << "\nNumber of buffers used (size): \t\t\t"
    <<  numfullbuffers << " (" << bufferlines << "," << numpixelsdisk << ")"
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();


  ofstream scratchresfile("scratchressubtrrefpha", ios::out | ios::trunc);
  bk_assert(scratchresfile,"subtrrefpha: scratchressubtrrefpha",__FILE__,__LINE__);
  scratchresfile
    // ______ Updateproductinfo greps these ______
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_subtrrefpha]
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\tpolynomial"
    << "\nData_output_file: \t\t\t";
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    scratchresfile 
      << "NO_OUTPUT_ONLY_DUMPING_REF_PHA"
      << "\nFile_name of ref.phase: \t\t"
      << input_i_subtrrefpha.forefpha
      << "\nData_output_format: \t\t\t"
      << "complex_real4";
    }
  else
    {
    scratchresfile
      <<  input_i_subtrrefpha.focint
      << "\nData_output_format: \t\t\t"
      << "complex_real4";
    }
  scratchresfile
    << "\nFirst_line (w.r.t. original_master): \t"
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master): \t"
    <<  interferogram.win.linehi                            //lastlineoutput
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixhi                             //lastpixeloutput
    << "\nMultilookfactor_azimuth_direction: \t"
    <<  totalmlL
    << "\nMultilookfactor_range_direction: \t"
    <<  totalmlP
    << "\nNumber of lines (multilooked): \t\t"
    <<  numlinesoutput
    << "\nNumber of pixels (multilooked): \t"
    <<  numpixelsoutput
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_subtrrefpha] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();

// ====== Tidy up ======
  // ??
  } // END subtrrefpha polynomial



/****************************************************************
 *    subtrrefpha  (exact)                                      *
 *                                                              *
 * Compute complex_interferogram .* conj(REF.PHA)               *
 * note: master-slave-ref.phase                                 *
 * ref.phase is computed here by evaluating system of 3 eq.     *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, multilook, filenames                     *
 * Output:                                                      *
 *  - new files on disk                                         *
 *    Bert Kampes, 03-Jul-2000                                  *
 ****************************************************************/
void subtrrefpha(
        const input_ell    &ellips,
        const slcimage     &master,
        const slcimage     &slave,
        const productinfo  &interferogram,
        const input_gen    &input_general,
        const input_subtrrefpha &input_i_subtrrefpha,
        orbit              &masterorbit,
        orbit              &slaveorbit)
  {
  TRACE_FUNCTION("subtrrefpha (BK 03-Jul-2000)")
  const int32 MAXITER   = 10;                           // number iterations
  const real8 CRITERPOS = 1e-6;                         // meters
  const real8 CRITERTIM = 1e-10;                        // seconds
  INFO << "SUBTRREFPHA: MAXITER: "   << MAXITER   << "; "
                    << "CRITERPOS: " << CRITERPOS << " m; "
                    << "CRITERTIM: " << CRITERTIM << " s";
  INFO.print();

// ====== Input options ======
  const int32 multiL    = input_i_subtrrefpha.multilookL;
  const int32 multiP    = input_i_subtrrefpha.multilookP;

// ====== Get some variables ======
  const int32 mldiskL   = interferogram.multilookL;     // cint on disk
  const int32 mldiskP   = interferogram.multilookP;     // cint on disk

  const int32 numlinesdisk    = interferogram.win.lines()/mldiskL;
  const int32 numpixelsdisk   = interferogram.win.pixels()/mldiskP;
  // ______ In original master radar coordinate system ______
  const real4 veryfirstline   = real4(interferogram.win.linelo) +
                                (real4(mldiskL) - 1.) / 2.;
  const real4 firstpixel      = real4(interferogram.win.pixlo)  +
                                (real4(mldiskP) - 1.) / 2.;

  const int32 totalmlL        = mldiskL*multiL;
  const int32 totalmlP        = mldiskP*multiP;
  const int32 numlinesoutput  = numlinesdisk/multiL;    // floor
  const int32 numpixelsoutput = numpixelsdisk/multiP;   // floor
 // const int32 lastlineoutput  = interferogram.win.linelo        // [MA] These were used be dumped to the res file.
 //                               + totalmlL*numlinesoutput - 1;  // Now, They are replaced by interferogram.win.linehi
 // const int32 lastpixeloutput = interferogram.win.pixlo         // and interferogram.win.pixhi.
 //                               + totalmlP*numpixelsoutput - 1; // see: comprefdem
  const bool outputh2ph       =  specified(input_i_subtrrefpha.foh2ph); // if spec. then output, added by FvL

// ====== Open files ======
// ______ If requested, dump ref.pha and do nothing else ______
  ifstream ifcint;                              // input file complex interf.
  ofstream ofilecint;                           // output file complex interf.
  ofstream ofrefpha;                            // output file ref. phase
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    openfstream(ofrefpha,input_i_subtrrefpha.forefpha,input_general.overwrit);
    bk_assert(ofrefpha,input_i_subtrrefpha.forefpha,__FILE__,__LINE__);
    }
  else  // compute cint and dump this...  
    {
    openfstream(ofilecint,input_i_subtrrefpha.focint,input_general.overwrit);
    bk_assert(ofilecint,input_i_subtrrefpha.focint,__FILE__,__LINE__);

    if (interferogram.formatflag != FORMATCR4)
      {
      PRINT_ERROR("code: .. Complex interferogram on disk must be complex real4.")
      throw(unhandled_case_error);
      }
    openfstream(ifcint,interferogram.file);
    bk_assert(ifcint,interferogram.file,__FILE__,__LINE__);
    }

//______________________________________________________________________
//__________ added by FvL ______________________________________________
//______________________________________________________________________   
  ofstream ofh2ph;
  if (outputh2ph==true)
    {
    openfstream(ofh2ph,input_i_subtrrefpha.foh2ph,input_general.overwrit);
    bk_assert(ofh2ph,input_i_subtrrefpha.foh2ph,__FILE__,__LINE__);
    }   
//______________________________________________________________________
//__________ end added by FvL ______________________________________________
//______________________________________________________________________   

// ====== Loop parameters ======
  matrix <real8> REFPHASE(multiL,numpixelsdisk);          // computed, read multiple of multiL
  matrix <real8> H2PH(multiL,numpixelsdisk);            // computed, added by FvL
  int32 tenpercent = int32(floor(numlinesoutput/10.));  // better round this
  if (tenpercent==0) tenpercent = 1000;
  int32 percentage = 0;
  real8 line = veryfirstline;                   // master coord. system, buffer0
  const real8 m_minpi4cdivlam = (-4*PI*SOL)/master.wavelength;
  const real8 s_minpi4cdivlam = (-4*PI*SOL)/slave.wavelength;

// ====== Compute delta range for all points ======
// ______ Read in numlinesoutput buffers of size multiL ______
  for (register int32 buffer=0; buffer<numlinesoutput; ++buffer)  // [MA] TODO not exactly a buffer operation, 
    {                                                             // but we can decrease disk access time by using large buffers as usual
    // ______ Give progress info each ten percent ______
    if (buffer%tenpercent==0)
      {
      PROGRESS << "SUBTRREFPHA: " << setw(3) << percentage << "%";
      PROGRESS.print();
      percentage += 10;
      }

    // ______ Compute ref. phase for this buffer ______
    for (register int32 i=0; i<multiL; ++i)
      {
      for (register int32 j=0; j<numpixelsdisk; ++j)
        {
        real8 pixel = firstpixel + j*mldiskP;           // master coord. system

        // ______ Compute range time for this pixel ______
        //const real8 m_trange = pix2tr(pixel,master.t_range1,master.rsr2x);
        const real8 m_trange = master.pix2tr(pixel);
        const real8 m_tazi = master.line2ta(line); // added by FvL

        // ______ Compute xyz of this point P from position in image ______
        cn P;                                       // point, returned by lp2xyz
        lp2xyz(line,pixel,ellips,master,masterorbit,
               P,MAXITER,CRITERPOS);

        // ______ Compute xyz for slave satellite from P ______
        real8 s_tazi;                               // returned, not used
        real8 s_trange;                             // returned
        xyz2t(s_tazi,s_trange,slave,
              slaveorbit,
              P,MAXITER,CRITERTIM);


        // ______Compute delta range ~= phase______
        // ______ real8 dr = dist(m_possat,pospoint) - dist(s_possat,pospoint);
        // ______ real8 phase = -pi4*(dr/LAMBDA);
        // ______  dr    == M-S         want if no flatearth M-S - flatearth = M-S-(M-S)=0
        // ______  phase == -4pi*dr/lambda == 4pi*(S-M)/lambda
        // BK: 24-9: actually defined as: phi = +pi4/lambda * (r1-r2) ???
        // real8 phase = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
        //y(i,0) = pi4*((dist(s_possat,pospoint)-dist(m_possat,pospoint))/LAMBDA);
        //y(i,0) = pi4divlam*(s_possat.dist(pospoint)-m_possat.dist(pospoint));
        //REFPHASE(i,j) = minpi4cdivlam*(m_trange - s_trange);
        REFPHASE(i,j) = m_minpi4cdivlam*m_trange -
                        s_minpi4cdivlam*s_trange;
                      
        if (outputh2ph==true)
        {
        // ____________________________________________________________________________________
        // _____________ Vector with h2ph factors for random number of points by FvL __________
        //_____________________________________________________________________________________

        cn Psat_master = masterorbit.getxyz(m_tazi);
        cn Psat_slave = slaveorbit.getxyz(s_tazi);
        real8 B    = Psat_master.dist(Psat_slave);                      // abs. value
                // const real8 Bpar = P.dist(M) - P.dist(S);            // sign ok
        real8 Bpar = SOL*(m_trange-s_trange);           // sign ok

                // ______ if (MP>SP) then S is to the right of slant line, then B perp is positive.
        cn r1 = Psat_master.min(P);
        cn r2 = Psat_slave.min(P);
        // real8 theta = Psat_master.angle(r1); // return look angle
        real8 theta       = P.angle(r1);   // incidence angle
        real8 theta_slave = P.angle(r2);   // incidence angle slave
        real8 Bperp = (theta > theta_slave) ?                   // sign ok
           sqrt(sqr(B)-sqr(Bpar)) :
          -sqrt(sqr(B)-sqr(Bpar));

        H2PH(i,j) = Bperp/(m_trange*SOL*sin(theta));
        }
        // ____________________________________________________________________________________
        // _____________ End added part by FvL ________________________________________________
        //_____________________________________________________________________________________
        } // pixels

      line += mldiskL;                                  // update here
      } // multilines loop

    // ______ Either dump refphase or subtract it ______
    if (input_i_subtrrefpha.dumponlyrefpha)
      {
      //matrix<complr4> REFPHASECR4 = fast_angle2cmplx(REFPHASE);
      //ofrefpha << multilook(REFPHASECR4,multiL,multiP);
      
      matrix<complr8> REFPHASECR8 = fast_angle2cmplx(REFPHASE); 
      matrix<complr4> output_layer(REFPHASECR8.lines(),REFPHASECR8.pixels());

      convert_type(REFPHASECR8, output_layer);
      ofrefpha << multilook(output_layer, multiL, multiP);

      }
    else
      { // subtrrefphase
      matrix <complr4> CINT(multiL,numpixelsdisk);      // read from disk as real4
      ifcint >> CINT;                           // read next buffer, filepointer++
      // ______ Changes CINT ______
      fast_dotmultconjphase(CINT,REFPHASE);             // pointwise multiplication with conj.
      ofilecint << multilook(CINT,multiL,multiP);               // write next line
      REFPHASE(1,1);                                    //deallocate
      }
    if (outputh2ph==true) // dump h2ph, added by FvL
      {
      matrix<real4> output_layer(H2PH.lines(),H2PH.pixels());
      convert_type(H2PH, output_layer);                   // [MA] ex: real8 --> real4
      ofh2ph << multilook(output_layer, multiL, multiP);  // H2PH --> outfile
      H2PH(1,1);                                          //deallocate
      }
    } // azimuth buffers


// ______ Close files ______
  if (input_i_subtrrefpha.dumponlyrefpha)
    ofrefpha.close();
  else
    {
    ofilecint.close();
    ifcint.close();
    }
  if (outputh2ph==true) // added by FvL
    ofh2ph.close();


// ====== Log info/results ======         // TODO put date and time to log as well 
  ofstream scratchlogfile("scratchlogsubtrrefpha", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"subtrrefpha: scratchlogsubtrrefpha",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* SUBTRREFPHA" 
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\t\texact";                                     
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    scratchlogfile
      << "\nOnly dump of reference phase, no subtraction."
      << "\nData_output_file_ref.phase: \t\t\t"
      <<  input_i_subtrrefpha.forefpha;
    }
  else
    {
    scratchlogfile
      << "\nInput file complex interferogram: \t\t"
      <<  interferogram.file
      << "\nData_output_file_complex_interferogram: \t"
      <<  input_i_subtrrefpha.focint
      << "\nData_output_format_complex_interferogram: \t"
      << "complex_real4";
    }
  scratchlogfile
    << "\nFirst_line (w.r.t. original_master): \t\t" // updateproductinfo greps these
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master): \t\t"
    <<  interferogram.win.linehi                            //lastlineoutput
    << "\nFirst_pixel (w.r.t. original_master): \t\t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master): \t\t"
    <<  interferogram.win.pixhi                             //lastpixeloutput
    << "\nMultilookfactor_azimuth_direction: \t\t"
    <<  totalmlL
    << "\nMultilookfactor_range_direction: \t\t"
    <<  totalmlP
    << "\nNumber of lines (multilooked): \t\t\t"
    <<  numlinesoutput
    << "\nNumber of pixels (multilooked): \t\t"
    <<  numpixelsoutput
    << "\n\nMultilookfactors input complex interferogram: \t"
    <<  mldiskL << " " << mldiskP
    << "\nMultilookfactors requested in this step: \t"
    <<  multiL << " " << multiP
    << "\nNumber of buffers used (size): \t\t\t" 
    //<<  numfullbuffers << "(" << bufferlines << "," << numpixelsdisk << ")"
    <<  numlinesoutput/multiL  << "(" << multiL << "," << numpixelsoutput << ")"
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();


  ofstream scratchresfile("scratchressubtrrefpha", ios::out | ios::trunc);
  bk_assert(scratchresfile,"subtrrefpha: scratchressubtrrefpha",__FILE__,__LINE__);
  scratchresfile
    // ______ Updateproductinfo greps these ______
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_subtrrefpha]
    << "\n*******************************************************************"
    << "\nMethod: \t\t\t\texact"
    << "\nData_output_file: \t\t\t";
  if (input_i_subtrrefpha.dumponlyrefpha)
    {
    scratchresfile 
      << "NO_OUTPUT_ONLY_DUMPING_REF_PHA"
      << "\nFile_name of ref.phase: \t\t"
      << input_i_subtrrefpha.forefpha
      << "\nData_output_format: \t\t\t"
      << "complex_real4";
    }
  else
    {
    scratchresfile
      <<  input_i_subtrrefpha.focint
      << "\nData_output_format: \t\t\t"
      << "complex_real4";
    }
  scratchresfile
    << "\nFirst_line (w.r.t. original_master): \t"
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master): \t"
    <<  interferogram.win.linehi                            //lastlineoutput
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixhi                             //lastpixeloutput
    << "\nMultilookfactor_azimuth_direction: \t"
    <<  totalmlL
    << "\nMultilookfactor_range_direction: \t"
    <<  totalmlP
    << "\nNumber of lines (multilooked): \t\t"
    <<  numlinesoutput
    << "\nNumber of pixels (multilooked): \t"
    <<  numpixelsoutput
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_subtrrefpha] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();


// ====== Tidy up ======
// ??
  } // END subtrrefpha exact




/****************************************************************
 *    subtrrefdem                                               *
 *                                                              *
 * Compute complex_interferogram .* conj(REF.DEM)               *
 * note: master-slave-ref.phase-demphase                        *
 * DEM ref.phase is a file containing float phase values        *
 * Optionally multilook again (not expected)                    *
 * Subtract in small buffers, assume size is same as cint.      *
 *                                                              *
 * Input:                                                       *
 *  - input arguments, filenames                                *
 * Output:                                                      *
 *  - new complex interferogram on disk                         *
 *    Bert Kampes, 10-Feb-2000                                  *
 * added offset if specified by user, not automatically computed*
 * Card SRD_OFFSET l p                                          *
 #%// BK 24-Apr-2002                                            *
 ****************************************************************/
void subtrrefdem(
        const productinfo       &interferogram,
        const productinfo       &radarcodedrefdem,
        const input_gen         &input_general,
        const input_subtrrefdem &input_i_subtrrefdem)
  {
  TRACE_FUNCTION("subtrrefdem (BK 10-Feb-2000)")
  const int32 additional_offsetL = input_i_subtrrefdem.offsetL;
  const int32 additional_offsetP = input_i_subtrrefdem.offsetP;

// ====== Handle input ======
// do  not multilook in this step..., neither cutout
//  const int32 mlL = input_i_subtrrefdem.mlL;
//  const int32 mlP = input_i_subtrrefdem.mlP;

// ______ Add offset from correlation to these windows ______
  //window cint   = interferogram.wininterfero;
  window refdem = radarcodedrefdem.win;

// IF INPUT METHOD = CORRELATE
// FIRST COMPUTE BEST SHIFT DEM AT PIXEL LEVEL
// THEN ADD THIS SHIFT TO WIN.REFDEM ...
// Do this by correlation of the phase image at lot of patches?
// or for total image within +-4 pixels.


//WARNING.print("BERT TESTING ADDITIONAL OFFSET SPECIFIED BY USER IN INPUT");
 refdem.linelo += (additional_offsetL * radarcodedrefdem.multilookL);
 refdem.linehi += (additional_offsetL * radarcodedrefdem.multilookL);
 refdem.pixlo  += (additional_offsetP * radarcodedrefdem.multilookP);
 refdem.pixhi  += (additional_offsetP * radarcodedrefdem.multilookP);
// then DEM is shifted to the right wrt. cint


// ====== Compute overlap interferogram, ref.dem ======
// ______ Output cint always same size input cint ______
  if (interferogram.multilookL != radarcodedrefdem.multilookL)
    {
    PRINT_ERROR("MultilookfactorL complex interferogram, ref. dem not equal.")
    throw(file_error);
    }
  if (interferogram.multilookP != radarcodedrefdem.multilookP)
    {
    PRINT_ERROR("MultilookfactorP complex interferogram, ref. dem not equal.")
    throw(file_error);
    }
  if ((interferogram.win.linelo-refdem.linelo)%radarcodedrefdem.multilookL != 0)
    WARNING.print("Seems reference phase DEM does not lie at same grid as complex interferogram.");
  if ((interferogram.win.pixlo-refdem.pixlo)%radarcodedrefdem.multilookP != 0)
    WARNING.print("Seems reference phase DEM does not lie at same grid as complex interferogram.");

  const int32 cintfilelines  =
    interferogram.win.lines() /interferogram.multilookL;
  const int32 cintfilepixels =
    interferogram.win.pixels()/interferogram.multilookP;
  const int32 demfilelines   =
    radarcodedrefdem.win.lines() /radarcodedrefdem.multilookL;
  const int32 demfilepixels  =
    radarcodedrefdem.win.pixels()/radarcodedrefdem.multilookP;

// YOU FORGOT TO TAKE PRIOR MULTILOOKING IN ACCOUNT !!!!!!!
// CHECK THIS (BK 11-feb-2000)
  const int32 offsetL        =
    (int32(refdem.linelo)-int32(interferogram.win.linelo))/
     int32(interferogram.multilookL);
  const int32 skiplinesstart = max(0,offsetL);  // number of ml.lines no overlap
  const int32 skiplinesend   = max(0,cintfilelines-demfilelines-offsetL);
  const int32 numpixoverlap  =
    int32(min(int32(interferogram.win.pixhi),int32(refdem.pixhi)) -
     max(int32(interferogram.win.pixlo),int32(refdem.pixlo)) + 1 ) /
     int32(radarcodedrefdem.multilookP);
  const int32 startcintP     =          // startindex in CINT array
    max(0,(int32(refdem.pixlo-interferogram.win.pixlo))/
           int32(radarcodedrefdem.multilookP));
  const int32 startrefdemP   =          // startindex in REFDEM array
    max(0,(int32(interferogram.win.pixlo-refdem.pixlo))/
           int32(radarcodedrefdem.multilookP));

  DEBUG << " skiplinesstart: " << skiplinesstart << " skiplinesend: " << skiplinesend
       << " offsetL: " << offsetL << " numpixoverlap: " << numpixoverlap
       << " startcintP: " << startcintP << " startrefdemP: " << startrefdemP;
  DEBUG.print();


// ====== Open files ======
  if (interferogram.formatflag != FORMATCR4)
    {
    PRINT_ERROR("code ..: Complex interferogram on disk assumed to be complex real4.")
    throw(file_error);
    }
  if (radarcodedrefdem.formatflag != FORMATR4)
    {
    PRINT_ERROR("code ..: Reference phase DEM on disk assumed to be real4.")
    throw(file_error);
    }

  ifstream ifcint;
  openfstream(ifcint,interferogram.file);
  bk_assert(ifcint,interferogram.file,__FILE__,__LINE__);

  ifstream ifrefdem;
  openfstream(ifrefdem,radarcodedrefdem.file);
  bk_assert(ifrefdem,radarcodedrefdem.file,__FILE__,__LINE__);

  ofstream ofilecint;
  openfstream(ofilecint,input_i_subtrrefdem.focint,input_general.overwrit);
  bk_assert(ofilecint,input_i_subtrrefdem.focint,__FILE__,__LINE__);


// ====== Loop over complex interferogram per line ======
  matrix<complr4> CINT   (1,cintfilepixels);    // read one line at a time
  matrix<real4>   REFDEM (1,demfilepixels);     // read one line at a time
  register int32 i,line;

#define COMP_COHER
#ifdef COMP_COHER
  real8 coher = 0;// BK 24-Apr-2002
#endif

  for (line=0; line<skiplinesstart; ++line)     // no overlap cint,dem
    {
    DEBUG.print("Copying line of interferogram since no overlap at start.");
    ifcint    >> CINT;
    ofilecint << CINT;
    }
  for (line=0; line<-offsetL; ++line)           // set file pointer refdem
    ifrefdem >> REFDEM;                         // do nothing ...
  for (line=0; line<cintfilelines-skiplinesstart-skiplinesend; ++line)
    {

    // ______ Read full line complex interferogram/ref.dem phase ______
    ifcint   >> CINT;
    ifrefdem >> REFDEM;

    // ______ Multiply by complex conjugated of phase (subtraction) ______
    // might be faster in matrix notation (DEM2=getpartdem; cint*=dem2)
    for (i=0; i<numpixoverlap; ++i)
      CINT[0][i+startcintP] *= complr4(fast_cos(REFDEM[0][i+startrefdemP]),
                                       fast_min_sin(REFDEM[0][i+startrefdemP]));
#ifdef COMP_COHER
    // ______ Coherence is how much CINT and DEM are alike _____
    // ______ I give the measure here as the sum over the pixels ______
    // coh=abs[sum_k=1:K{abs(CINT)*exp(i CINT)exp(-i DEM)}/sqrt(abs(CINT))]
    complr4 coher_line = 0.;
    real4   sum_line   = 0.;
    for (i=0; i<numpixoverlap; ++i)
      {
      coher_line   += CINT[0][i+startcintP];
      sum_line     += abs(CINT[0][i+startcintP]);
      }
    coher += abs(coher_line)/sum_line;
#endif

    // ______ Write line complex interferogram ______
    ofilecint << CINT;
    }
  ifrefdem.close();
#ifdef COMP_COHER
  coher /= (cintfilelines-skiplinesstart-skiplinesend);
  INFO << "coherence interferogram_synthetic phase = " << coher;
  INFO.print();
#endif
  for (line=0; line<skiplinesend; ++line)       // no overlap cint,dem
    {
    DEBUG.print("Copying line of interferogram since no overlap at end.");
    ifcint    >> CINT;
    ofilecint << CINT;
    } // loop line
  ifcint.close();
  ofilecint.close();


// ====== Log info/results ======
  ofstream scratchlogfile("scratchlogsubtrrefdem", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"subtrrefdem: scratchlogsubtrrefdem",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* SUBTRREFDEM PHASE" 
    << "\n*******************************************************************"
    << "\nInput file complex interferogram: \t\t\t"
    <<  interferogram.file
    << "\nInput file reference dem: \t\t\t"
    <<  radarcodedrefdem.file
    << "\nAdditional_azimuth_shift:             \t"
    << additional_offsetL
    << "\nAdditional_range_shift:               \t"
    << additional_offsetP
    << "\nCoherence IFG DEMPHASE:               \t"
    << coher
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();
 
  ofstream scratchresfile("scratchressubtrrefdem", ios::out | ios::trunc);
  bk_assert(scratchresfile,"subtrrefdem: scratchressubtrrefdem",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_subtrrefdem]
    << "\n*******************************************************************"
    << "\nMethod:                               \t"
    << "NOT_USED"
    << "\nAdditional_azimuth_shift:             \t"
    << additional_offsetL
    << "\nAdditional_range_shift:               \t"
    << additional_offsetP
    << "\nData_output_file:                     \t"
    <<  input_i_subtrrefdem.focint
    << "\nData_output_format:                   \t"
    << "complex_real4"
    << "\nFirst_line (w.r.t. original_master):  \t"
    <<  interferogram.win.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  interferogram.win.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  interferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  interferogram.win.pixhi
    << "\nMultilookfactor_azimuth_direction:    \t"
    <<  interferogram.multilookL
    << "\nMultilookfactor_range_direction:      \t"
    <<  interferogram.multilookP
    << "\nNumber of lines (multilooked):        \t"
    <<  cintfilelines
    << "\nNumber of pixels (multilooked):       \t"
    <<  cintfilepixels
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_subtrrefdem] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  } // END subtrrefdem



/****************************************************************
 *    dinsar                                                    *
 * Differential insar with an unwrapped topo interferogram      *
 * (hgt or real4 format) and a wrapped(!) defo interf.          * 
 * if r4 then NaN==-999 is problem with unwrapping, else hgt    *
 * if ampl. =0 then problem flagged with unwrapping.            *
 * The topography is removed from the deformation interferogram *
 * by the formula (prime ' denotes defo pair):                  *
 * dr       = lambda\4pi * [phi' - phi(Bperp'/Bperp)]           *
 * phi_diff = phi(Bperp'/Bperp) - phi'                          *
 * where Bperp is the perpendicular baseline for points on the  *
 * ellipsoid (and not the true one)!                            *
 * I implemented this by computing the baseline for a number    *
 * of points for topo and defo, and then modeling the ratio     *
 * as a 2D polynomial of degree 1 for the image.                *
 * Then evaluating this to compute the new phase (defo only).   *
 * I assume the interferogram files are coregistered on each    *
 * other and have the same dimensions.                          *
 *                                                              *
 * If TOPOMASTER file is empty (" "), then use current master   *
 * res file for master orbit (== 3pass), else get orbit         *
 * (==4pass).                                                   *
 *                                                              *
 * Input:                                                       *
 *  -input parameters                                           *
 *  -orbits                                                     *
 *  -info on input files                                        *
 *  -                                                           *
 * Output:                                                      *
 *  -complex float file with differential phase.                *
 *   (set to (0,0) for not ok unwrapped parts)                  *
 *                                                              *
 * See also Zebker, 1994.                                       *
 #%// BK 22-Sep-2000                                            *
 ****************************************************************/
void dinsar(
        const input_gen         &input_general,
        const input_dinsar      &dinsarinput,
        const input_ell         &ellips,
        const slcimage          &master,
              orbit             &masterorbit,
        const slcimage          &defoslave,
              orbit             &defoorbit,
        const productinfo       &defointerferogram
        )
  {
  TRACE_FUNCTION("dinsar (BK 22-Sep-2000)")
  // ====== Get input, check file dimensions ======
  char        difffile[2*ONE27];                  // output file name
  strcpy(difffile,dinsarinput.fodinsar);

  // ______ Fill these from resultfiles topo-pair processing ______
  // ______ Processing topo pair has to be until unwrap ______ 
  slcimage    topomaster;               // only fill if 4 pass
  slcimage    toposlave;                // info on slave image
  productinfo topounwrappedinterf;      // interferogram
  orbit       topomasterorbit;          // only fill if 4 pass
  orbit       toposlaveorbit;           // always fill

  // ______ Check 4 pass if topomaster is specified (diff than m_res) ______
  bool FOURPASS = true;                 // assume 4pass
  if (!specified(dinsarinput.topomasterresfile))
    FOURPASS = false;
  if (!strcmp(dinsarinput.topomasterresfile,input_general.m_resfile))
    FOURPASS = false;

  if (FOURPASS==true)
    {
    INFO.print("Using 4 pass differential interferometry (card IN_TOPOMASTER).");
    INFO.print("Reading primary (topography pair) master parameters.");
    topomaster.fillslcimage(dinsarinput.topomasterresfile);     // fill wavelength etc.
    INFO.print("Modelling primary (topography pair) master orbit.");
    topomasterorbit.initialize(dinsarinput.topomasterresfile);  // fill interp. coeff.
    }

  INFO.print("Reading primary (topography pair) slave parameters.");
  toposlave.fillslcimage(dinsarinput.toposlaveresfile); // fill wavelength etc.
  INFO.print("Reading (topography pair) unwrapped section.");
  char   SECTIONID[ONE27];
  strcpy(SECTIONID,"*_Start_");
  strcat(SECTIONID,processcontrol[pr_i_unwrap]);
  topounwrappedinterf.fillproductinfo(dinsarinput.topointresfile,SECTIONID); // fill info
  INFO.print("Modelling primary (topography pair) slave orbit.");
  toposlaveorbit.initialize(dinsarinput.toposlaveresfile);      // fill interp. coeff.


  // ______ Check dimensions ______
  if (defointerferogram.multilookL != topounwrappedinterf.multilookL)
    {
    WARNING << "multilookfactor defo != factor topo: " << defointerferogram.multilookL
         << " != " << topounwrappedinterf.multilookL;
    WARNING.print();
    }
  if (defointerferogram.multilookP != topounwrappedinterf.multilookP)
    {
    WARNING << "multilookfactor defo != factor topo: " << defointerferogram.multilookP
         << " != " << topounwrappedinterf.multilookP;
    WARNING.print();
    }
  if (defointerferogram.win != topounwrappedinterf.win)
    WARNING.print("window defo pair != window topo pair");

  if (defointerferogram.formatflag != FORMATCR4)
    {
    PRINT_ERROR("format defo interferogram not complr4.")
    throw(file_error);
    }

  // not overlap, multilooked, center of pixel, trouble
  //const window overlap = getoverlap(topounwrappedinterf.win,
  //                                defointerferogram.win);


  // ______ Normalization factors for polynomial ______
  const real8 minL = master.originalwindow.linelo;
  const real8 maxL = master.originalwindow.linehi;
  const real8 minP = master.originalwindow.pixlo;
  const real8 maxP = master.originalwindow.pixhi;
  INFO << "dinsar: polynomial for ratio normalized by: "
          << minL << " " << maxL << " " << minP << " " << maxP
          << " to [-2,2]";
  INFO.print();

/*  // TODO
  BASELINE topo_baseline, defo_baseline;

  if (FOURPASS==true)
    {
  topomasterorbit
    topo_baseline.model_parameters(topomaster,toposlave,topomasterorbit,toposlaveorbit,ellips);
    defo_baseline.model_parameters(defomaster,defoslave,defomasterorbit,defoslaveorbit,ellips);
    }
  else
    {
    }
  const real8 Bperp = topo_baseline.get_bperp(line,pixel);
  const real8 Bperp = defo_baseline.get_bperp(line,pixel);
  ...
*/

  // ====== Model perpendicular baseline for master and slave ======
  // ______ compute B on grid every 500 lines, 100 pixels 
  // ______ in window for topo/defo ______ 
  const int32 numpointsL = 20;                                  // grid for modelling
  const int32 numpointsP = 10;                                  // grid for modelling
  real8 dlines  = (topounwrappedinterf.win.linehi-topounwrappedinterf.win.linelo) /
                  (numpointsL-1);
  real8 dpixels = (topounwrappedinterf.win.pixhi -topounwrappedinterf.win.pixlo)  /
                  (numpointsP-1);
  INFO << "Computing baseline on grid (" 
       << topounwrappedinterf.win.linelo << ":" << dlines << ":"
       << topounwrappedinterf.win.linehi << ","
       << topounwrappedinterf.win.pixlo << ":" << dpixels << ":"
       << topounwrappedinterf.win.pixhi
       << ") = (" << numpointsL << "x" << numpointsP << ")";
  INFO.print();

  matrix<real8> LINENUMBER(numpointsL*numpointsP,1);
  matrix<real8> PIXELNUMBER(numpointsL*numpointsP,1);
  int32 i,j;
  int32 k=0;
  for (i=0; i<numpointsL; ++i)
    {
    for (j=0; j<numpointsP; ++j)
      {
      LINENUMBER(k,0)  = topounwrappedinterf.win.linelo + i*dlines;     // line coordinate
      PIXELNUMBER(k,0) = topounwrappedinterf.win.pixlo + j*dpixels;     // pixel coordinate
      ++k;
      }
    }

  // ______ compute master, point on ellips, position for these lines ______
  cn masterpos[numpointsL*numpointsP];
  cn pointpos[numpointsL*numpointsP];
  cn defoslavepos[numpointsL*numpointsP];
  cn topomasterpos[numpointsL*numpointsP];      // 4 pass then diff from masterpos
  cn toposlavepos[numpointsL*numpointsP];

  real8 lastline = -1.0;
  real8 B,Bpar,Bperp;
  matrix<real8> Bperptopo(LINENUMBER.lines(),1);
  matrix<real8> Bperpdefo(LINENUMBER.lines(),1);
  for (int32 i=0; i<int32(LINENUMBER.lines()); ++i)
    {
    real8 line  = LINENUMBER(i,0);
    real8 pixel = PIXELNUMBER(i,0);
    if (line == lastline)
      {
      masterpos[i] = masterpos[i-1];                            // same position
      }
    else
      {
      lastline = line;
      const real8 m_tazi = master.line2ta(line);
      masterpos[i]       = masterorbit.getxyz(m_tazi);
      }
    // ______ Do it the *slow* way, get 3d positions slaves ______
    lp2xyz(line,pixel,ellips,master,masterorbit,pointpos[i]);           // fill pointpos
    xyz2orb(defoslavepos[i],defoslave,defoorbit,pointpos[i]);           // fill defopos
    xyz2orb(toposlavepos[i],toposlave,toposlaveorbit,pointpos[i]);      // fill toposlavepos
    if (FOURPASS==true) // fill topomasterpos
      xyz2orb(topomasterpos[i],topomaster,topomasterorbit,pointpos[i]);
    else // 3 pass, same topomaster as defomaster...
      topomasterpos[i] = masterpos[i];

    // ______ Do it the *slow* way, based on 3 (4) 3d positions ______
    BBparBperp(B,Bpar,Bperp,topomasterpos[i],pointpos[i],toposlavepos[i]);// fill Bperp
    Bperptopo(i,0) = Bperp;
    BBparBperp(B,Bpar,Bperp,masterpos[i],pointpos[i],defoslavepos[i]);  // fill Bperp
    Bperpdefo(i,0) = Bperp;
    }

  // ______ Now model ratio Bperpdefo/Bperptopo as linear ______
  // ______ r(l,p) = a00 + a10*l + a01*p ______
  // ______ give stats on max. error ______
  matrix <real8> Ratio = Bperpdefo/Bperptopo;

  #ifdef __DEBUG
  DEBUG.print("data dump to files: LINENUMBER, PIXELNUMBER, Bperptopo, Bperpdefo, Ratio");
  dumpasc("LINENUMBER",LINENUMBER);
  dumpasc("PIXELNUMBER",PIXELNUMBER);
  dumpasc("Bperptopo",Bperptopo);
  dumpasc("Bperpdefo",Bperpdefo);
  dumpasc("Ratio",Ratio);
  #endif

  // old, better but changes info i want later.
  //normalize(LINENUMBER,minL,maxL);
  //normalize(PIXELNUMBER,minP,maxP);
  //A.setcolumn(0,1.);
  //A.setcolumn(1,LINENUMBER);
  //A.setcolumn(2,PIXELNUMBER);

  // ______ Set designmatrix, compute normalmatrix, righthandside ______
  matrix<real8> A(Ratio.lines(),3);
  for (uint i=0; i<A.lines(); ++i)
    {
    A(i,0) = 1.0;
    A(i,1) = normalize(LINENUMBER(i,0),minL,maxL);
    A(i,2) = normalize(PIXELNUMBER(i,0),minP,maxP);
    }
  matrix<real8> N   = matTxmat(A,A);
  matrix<real8> rhs = matTxmat(A,Ratio);

  // ______ Compute solution ______
  matrix<real8> Qx_hat = N;
  choles(Qx_hat);               // Cholesky factorisation normalmatrix
  solvechol(Qx_hat,rhs);        // Estimate unknowns (a00,a10,a01) in rhs
  invertchol(Qx_hat);           // Covariance matrix of unknowns


  // ______ Test inverse (thus stability cholesky) ______
  for (uint i=0; i<Qx_hat.lines(); i++)
    for (uint j=0; j<i; j++)
      Qx_hat(j,i) = Qx_hat(i,j);// repiar only stored Lower tri
  const real8 maxdev = max(abs(N*Qx_hat-eye(real8(Qx_hat.lines()))));
  INFO << "dinsar: max(abs(N*inv(N)-I)) = " << maxdev;
  INFO.print();
  if (maxdev > 0.01)
    {
    ERROR << ". Too large, normalization factors <-> crop?";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  else if (maxdev > 0.001)
    {
    WARNING.print("Deviation quite large.  careful!");
    }


  // ______ Some other stuff for logfile ______
  //  matrix<real8> Qy_hat        = A * (matxmatT(Qx_hat,A));
  matrix<real8> y_hat   = A * rhs;
  matrix<real8> e_hat   = Ratio - y_hat;

  uint pos,dummy;
  const real8 maxerrorratio = max(abs(e_hat),pos,dummy);
  const real8 maxrelerror   = 100.0 * maxerrorratio / Ratio(pos,0);
  INFO << "maximum error for l,p: " << LINENUMBER(pos,0) << ","
       << PIXELNUMBER(pos,0) << "; Ratio=" << Ratio(pos,0)
       << " estimate=" << y_hat(pos,0) << "; rel. err=" << maxrelerror << "%. ";
  INFO.print();
  if (maxrelerror < 5.0 )
    {
    INFO.print("max err OK");
    }
  else
    {
    WARNING.print("max err quite large");
    WARNING.print("Error in deformation vector larger than 5% due to mismodeling baseline!");
    }



  // ====== Per 100 lines, read in topo and defo interf. ======
  ofstream ofcint,ofscaledtopo;
  openfstream(ofcint,dinsarinput.fodinsar,input_general.overwrit);
  bk_assert(ofcint,dinsarinput.fodinsar,__FILE__,__LINE__);

  bool writescaledtopo=false;
  if (specified(dinsarinput.foscaleduint))
    {
    INFO.print("writing scaled version of unwrapped topo interferogram in real4 format.");
    writescaledtopo=true;
    openfstream(ofscaledtopo,dinsarinput.foscaleduint,input_general.overwrit);
    bk_assert(ofscaledtopo,dinsarinput.foscaleduint,__FILE__,__LINE__);
    }

  const int32 numlines  = defointerferogram.win.lines()  / defointerferogram.multilookL;
  const int32 numpixels = defointerferogram.win.pixels() / defointerferogram.multilookP;
  const real4 firstline = defointerferogram.win.linelo +
                        (real8(defointerferogram.multilookL-1.)/2.);
  const real4 firstpixel= defointerferogram.win.pixlo +
                        (real8(defointerferogram.multilookP-1.)/2.);

  matrix<real4> ratioline(1,numpixels);
  for (int32 i=0; i<numpixels; ++i)
    ratioline(0,i) = firstpixel+i*defointerferogram.multilookP;
  normalize(ratioline,real4(minP),real4(maxP));
  ratioline *= real4(rhs(2,0));         //     a01*p
  ratioline += real4(rhs(0,0));         // a00+a01*p

  // ______ read in matrices line by line, correct phase ______
  ifstream ifdefocint;
  openfstream(ifdefocint,defointerferogram.file);
  bk_assert(ifdefocint,defointerferogram.file,__FILE__,__LINE__);
  matrix<complr4> DEFO(1,numpixels);    // buffer

  // test if reading phase was ok...
  //cerr << "test: writing hgt phase and cint.\n";
  //ofstream oftest1("TOPO.raw", ios::out | ios::binary | ios::trunc);
  //ofstream oftest2("CINT.raw", ios::out | ios::binary | ios::trunc);

  int32 tenpercent = int32(floor(numlines/10.));
  if (tenpercent==0) tenpercent = 1000;
  int32 percent = 0;
  for (int32 i=0; i<numlines; ++i)
    {
    if (i%tenpercent==0)
      {
      PROGRESS << "DINSAR: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    // ______ ratio=a00+a10*l+a01*p ______
    const real4 line    = firstline + i*defointerferogram.multilookL;
    matrix<real4> ratio = ratioline +
                          real4(rhs(1,0))*normalize(line,real4(minL),real4(maxL));

    // ______ read from file, correct, write to file ______
    ifdefocint >> DEFO;                         // read next full line
    const window filewin(i+1,i+1,1,numpixels);
    matrix<real4> TOPO = topounwrappedinterf.readphase(filewin);

    // seems faster, but how to check for unwrapping?
    //TOPO *= ratio;    // scaled topo to defo baseline
    //DEFO *= complr4(cos(TOPO),-sin(TOPO));
    // better matrix <int32> index = TOPO.find(NaN); later reset??
    // and topo=topo*ratio; and defo(index)=(0,0);
    // but how to implement this best in matrixclass?
    // BK 24-Oct-2000
    for (int32 j=0; j<numpixels; ++j)
      {
      (TOPO(0,j)==NaN) ?                // if unwrapping ok, then subtract scaled phase
        DEFO(0,j)  = complr4(0.0, 0.0) :
        DEFO(0,j) *= complr4(fast_cos(ratio(0,j)*TOPO(0,j)),fast_min_sin(ratio(0,j)*TOPO(0,j)));
      }
    ofcint << DEFO;             // now topo-corrected phase

    // ______ Slow, only debugging ______
    if (writescaledtopo)
      {
      if (!(i%1000))
        PROGRESS.print("Writing scaled topo interferogram to file.");
      TOPO *= ratio;    // scaled topo to defo baseline
      ofscaledtopo << TOPO;
      }
    } 
  ifdefocint.close();
  ofcint.close();
  if (writescaledtopo)
    ofscaledtopo.close();


  // ====== Write result file; tidy up ======
  ofstream scratchlogfile("scratchlogdinsar", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"dinsar: scratchlogdinsar",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* " << processcontrol[pr_i_dinsar]
    << "\n*******************************************************************"
    << "\nDegree 2d polynomial for baseline modeling: "
    << "1"
  // ______ Write estimated coefficients ______
    << "\nEstimated coefficients:\n"
    << "\nx_hat \tstd:\n";
  for (uint i=0; i<rhs.lines(); ++i)
    scratchlogfile <<  rhs(i,0) << " \t" << sqrt(Qx_hat(i,i)) << endl;
  // ______ Write covariance matrix ______
  scratchlogfile << "\nQx_hat:\n";
  for (uint i=0; i<Qx_hat.lines(); i++)
    {
    for (uint j=0; j<Qx_hat.pixels(); j++)
      {
      scratchlogfile <<  Qx_hat(i,j) << " ";
      }
    scratchlogfile << endl;
    }
  scratchlogfile << "\nMaximum deviation N*inv(N):"
                 << maxdev
    << "\nBaseline for pixel: "
    << LINENUMBER(0,0) << ", " << PIXELNUMBER(0,0) << ":"
    << "\n  Bperp_defo:  " << Bperpdefo(0,0)
    << "\n  Bperp_topo:  " << Bperptopo(0,0) << endl
    << "\nSome more info for each observation:"
    << "\nline \t pixel \t obs \t obs_hat \t err_hat\n";
  for (uint i=0; i<LINENUMBER.lines(); i++)
    scratchlogfile <<  LINENUMBER(i,0) << "\t" <<  PIXELNUMBER(i,0) << "\t"
                   << Ratio(i,0) << "\t" << y_hat(i,0) << "\t" << e_hat(i,0) << "\n";
  scratchlogfile
    << "\nMaximum absolute error: \t"
    <<  maxerrorratio << " = " << maxrelerror << "%"
    << "\ninput filename topo interferogram: "
    <<  defointerferogram.file
    << "\ninput filename defo interferogram: "
    <<  topounwrappedinterf.file
    << "\noutput filename topocorrected defo: "
    <<  dinsarinput.fodinsar;
  if (writescaledtopo)
    scratchlogfile
      << "\noutput filename scaled topo interferogram: "
      <<  dinsarinput.foscaleduint;
  scratchlogfile
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  ofstream scratchresfile("scratchresdinsar", ios::out | ios::trunc);
  bk_assert(scratchresfile,"dinsar: scratchresdinsar",__FILE__,__LINE__);
  scratchresfile.setf(ios::scientific, ios::floatfield);
  scratchresfile.setf(ios::right, ios::adjustfield);
  scratchresfile.precision(8);
  scratchresfile.width(18);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_dinsar]
    << "\n*******************************************************************"
    << "\nData_output_file:                     \t"
    <<  dinsarinput.fodinsar
    << "\nData_output_format:                   \t"
    << "complex_real4"
    << "\nFirst_line (w.r.t. original_master):  \t"
    <<  defointerferogram.win.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  defointerferogram.win.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  defointerferogram.win.pixlo
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  defointerferogram.win.pixhi
    << "\nMultilookfactor_azimuth_direction:    \t"
    <<  defointerferogram.multilookL
    << "\nMultilookfactor_range_direction:      \t"
    <<  defointerferogram.multilookP
    << "\nNumber of lines (multilooked):        \t"
    <<  numlines
    << "\nNumber of pixels (multilooked):       \t"
    <<  numpixels
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_dinsar] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  } // END dinsar



/****************************************************************
 *    Simulate amplitude                                        *
 *                                                              *
 * Compute synthetic amplitude of a SAR image based on an       *
 * external DEM.                                                *
 *                                                              *
 * Batuhan Osmanoglu <batu@rsmas.miami.edu>  (via phase slope)  *
 * Freek van Leijen  <f.j.vanleijen@tudelft.nl> (topo slope)    *
 * Mahmut Arikan     <m.arikan@tudelft.nl>                      *
 *                                                              *
 * See paper by Eineder 2003                                    *
 *                                                              *
 ****************************************************************/
void sim_amplitude(
        const input_gen        &generalinput,
        const input_ell        &ellips,
        const input_simamp     &simamp,
        const slcimage         &master,
        orbit                  &masterorbit)
  {
  TRACE_FUNCTION("sim_amplitude (MA, FvL, BO 31-OCT-2008)")
  
  const string STEP="SAM: ";
  const string STEPlog="SAM";
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;

  const real8 lat0file    = simamp.demlatleftupper;// first pix on disk w02090
  const real8 lon0file    = simamp.demlonleftupper;// first pix on disk
  const real8 DEMdeltalat = simamp.demdeltalat;    // in radians
  const real8 DEMdeltalon = simamp.demdeltalon;    // in radians
  const int32 numberoflonpixels = simamp.demcols;  // NCOLS on file
  const int32 numberoflatpixels = simamp.demrows;  // NROWS on file
  const real8 NODATA      =  simamp.demnodata;     // (BK 4 may 2001)
  //const bool outputdemi   =  specified(simamp.fodemi);// if spec. then output
  const bool outputrefdemhei =  specified(simamp.fodemlp);   // MA SPT
  const bool outputthetalp   =  specified(simamp.fothetalp);  // MA SPT

  //const int32 mlL = 1/master.ovs_az;   // multilookfactor azimuth
  //const int32 mlP = 1/master.ovs_rg;   // multilookfactor range this is ifgmetric convention of Doris
  const int32 mlL = 1;   // multilookfactor azimuth
  const int32 mlP = 1;   // multilookfactor range SLCs should be 1 even if oversampled.
  //INFO << "master multilook_az = " << mlL;
  //INFO.print();
  //INFO << "master multilook_rg = " << mlP;
  //INFO.print();
  //DEBUG << "master wavelength = " << master.wavelength;
  //DEBUG.print();
  //DEBUG << "slave  wavelength = " << slave.wavelength;
  //DEBUG.print();

  const real8 latNfile = lat0file-DEMdeltalat*(numberoflatpixels-1); // upper=max. lat value
  const real8 lonNfile = lon0file+DEMdeltalon*(numberoflonpixels-1); // left=min. lon value

  // ______ Extra info ______
  INFO << "DEM input: w/e/s/n:          \t"
       << rad2deg(lon0file) << "/" << rad2deg(lonNfile) << "/"
       << rad2deg(latNfile) << "/" << rad2deg(lat0file);
  INFO.print();

  // ______ Get corners of master (approx) to select DEM ______
  // ______ in radians (if height were zero)______
  real8 extralat  = (1.5*DEMdeltalat + deg2rad(4.0/25.0));
  real8 extralong = (1.5*DEMdeltalon + deg2rad(4.0/25.0));

  real8 phimin;
  real8 phimax;
  real8 lambdamin;
  real8 lambdamax;
  int32 indexphi0DEM;
  int32 indexphiNDEM;
  int32 indexlambda0DEM;
  int32 indexlambdaNDEM;
  getcorners(master.currentwindow.linelo,master.currentwindow.linehi,
             master.currentwindow.pixlo,master.currentwindow.pixhi,
             extralat,extralong,lat0file,lon0file,
             DEMdeltalat,DEMdeltalon,numberoflatpixels,numberoflonpixels,
             ellips,master,masterorbit,phimin,phimax,lambdamin,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

  // ______ Extra info ______
  INFO << "DEM input required: w/e/s/n: \t"
       << rad2deg(lambdamin) << "/" << rad2deg(lambdamax) << "/"
       << rad2deg(phimin)    << "/" << rad2deg(phimax);
  INFO.print();
  INFO << "DEM input window (l0,lN,p0,pN):    \t"
       << indexphi0DEM << " " << indexphiNDEM << " "
       << indexlambda0DEM  << " " << indexlambdaNDEM ;
  INFO.print();
  INFO << "For master window (l0,lN,p0,pN):    \t"
       << master.currentwindow.linelo << " "
       << master.currentwindow.linehi << " "
       << master.currentwindow.pixlo << " "
       << master.currentwindow.pixhi;
  INFO.print();


  // ______ Check corners of DEM ______
  // check if DEM is appropriate for master
  // DEM should at least partially cover master
  // note: phi is [90:-90]
  if (phimax <= latNfile)// DEM is more north than master
    {
    ERROR << "master outside DEM: most South latitude: " << rad2deg(latNfile)
         << " [deg]; master requires: " << rad2deg(phimax) 
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  // DEM is more south than master
  if (phimin >= lat0file)// largest latitude at first line of file
    {
    ERROR << "master outside DEM: most North latitude: " << rad2deg(lat0file)
         << " [deg]; master requires: " << rad2deg(phimax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  if (lambdamax <= lon0file)
    {
    ERROR << "master outside DEM: most West longitude: " << rad2deg(lon0file)
         << " [deg]; master window requires: " << rad2deg(lambdamax)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }
  if (lambdamin >= lonNfile)
    {
    ERROR << "master outside DEM: most East longitude: " << rad2deg(lonNfile)
         << " [deg]; master window requires: " << rad2deg(lambdamin)
         << " [deg]";
    PRINT_ERROR(ERROR.get_str())
    throw(some_error);
    }


  //===================================================================
  //============ First loop: radarcode DEM ============================
  //============ (DEM geometry)            ============================
  //===================================================================

  int32 numvalid        = 0;// number of good values, not NODATA in buffer
  int32 numNODATA       = 0;// number of NODATA values in buffer
  real8 meancroppedDEM  = 0.0;// to detect byte order problems, formats
  real8 min_input_dem   =  100000.0;// stats
  real8 max_input_dem   = -100000.0;// stats

  // ______ Compute buffer size radarcoding DEM______
  const real8 BUFFERMEMSIZE = generalinput.memory;// Bytes
  const int32 NcolsDEM = indexlambdaNDEM-indexlambda0DEM+1;
  const int32 NrowsDEM = indexphiNDEM-indexphi0DEM+1;
  const real8 Nrows_possible_DEM = BUFFERMEMSIZE / (5*8*NcolsDEM);
  int32 bufferlines = int32(ceil(Nrows_possible_DEM));
  if (bufferlines>NrowsDEM) bufferlines=NrowsDEM;
  int32 numfullbuffers = NrowsDEM / bufferlines;
  int32 restlines      = NrowsDEM % bufferlines;
  int32 extrabuffer = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "DEM output total pixels: " << NcolsDEM;
  INFO.print();
  INFO << "DEM output total lines : " << NrowsDEM;
  INFO.print();
  INFO << "Radar coding of DEM in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();
  
  const string fsuccess = "sam_m_theta.success";  // [MA] a file to check successfull run of first loop


  //if(!existed("sam_m_theta.temp")) // The condition to skip first loop
  if(!existed(fsuccess.c_str())) // The condition to skip first loop
    {

  // ______ Open (temporary) output files ______
  // DEM heights 
  ofstream demofile;
  openfstream(demofile,simamp.fodem,generalinput.overwrit);
  bk_assert(demofile,simamp.fodem,__FILE__,__LINE__);

  // theta (incidence angle) (in DEM coords)
  ofstream thetaoutfile("sam_m_theta.temp", ios::out | ios::trunc); 
  bk_assert(thetaoutfile,"sam_m_theta.temp",__FILE__,__LINE__);

  // master line coordinates of DEM
  ofstream masterdemlineoutfile("sam_m_demline.temp", ios::out | ios::trunc); 
  bk_assert(masterdemlineoutfile,"sam_m_demline.temp",__FILE__,__LINE__);
  
  // master pixel coordinates of DEM
  ofstream masterdempixeloutfile("sam_m_dempixel.temp", ios::out | ios::trunc); 
  bk_assert(masterdempixeloutfile,"sam_m_dempixel.temp",__FILE__,__LINE__);


  // ______ DEM loop per buffer ______
  register int32 j,i;// DEM index grid counter, register j first to ensure allocation
  for (register int32 buffer=0; buffer<numfullbuffers+extrabuffer; ++buffer)
    {

     // Determine indices for buffer
    const int32 indexphi0BUFFER = indexphi0DEM+buffer*bufferlines;
    const int32 blines          = (buffer == numfullbuffers) ? restlines : bufferlines;
    const int32 indexphiNBUFFER = indexphi0BUFFER+(blines-1);
    matrix<real4> DEM(blines,NcolsDEM);

    // ______ Extra info ______
    PROGRESS << STEP << "Buffer# [l0:lN, p0:pN]: " << buffer+1 << " ["
         << indexphi0BUFFER  << ": " << indexphiNBUFFER  << ", "
         << indexlambda0DEM << ": " << indexlambdaNDEM << "]";
    PROGRESS.print();

    // ______ lat/lon for first pixel in matrix read from file ______
    // ______ upper is max. latitude, left is min. longitude ______
    const real8 upperleftphi    = lat0file-indexphi0BUFFER*DEMdeltalat;
    const real8 upperleftlambda = lon0file+indexlambda0DEM*DEMdeltalon;

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0BUFFER,indexphiNBUFFER,
                        indexlambda0DEM,indexlambdaNDEM);

    // ______ Read in grdfile of DEM in matrix R4 (raw data, no header) _______
    // ______ added formats (BK 4-May-2001) ______
    PROGRESS << STEP << "Reading crop of DEM for buffer: " << buffer+1;
    PROGRESS.print();
    DEBUG.print("Reading input DEM into real4 matrix (buffer).");
    switch (simamp.iformatflag)                         // TODO [MA] potential external function
      {
      // ______ Read as short BE, then convert to host order ______
      case FORMATI2_BIGENDIAN:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,simamp.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = real4(ntohs(DEMi2(iii,jjj)));// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATI2:
        {
        matrix<int16> DEMi2(blines,NcolsDEM);
        readfile(DEMi2,simamp.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMi2(iii,jjj);// cast to real4
        DEMi2.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: SHORT SIGNED INT.");
        break;
        }

      case FORMATR4:
        readfile(DEM,simamp.firefdem,numberoflatpixels,winfromfile,zerooffset);
        INFO.print("Read crop of input DEM: format: REAL4.");
        break;
      case FORMATR8:
        {
        matrix<real8> DEMr8(blines,NcolsDEM);
        readfile(DEMr8,simamp.firefdem,numberoflatpixels,winfromfile,zerooffset);
        for (int32 iii=0; iii<DEM.lines(); ++iii)
          for (int32 jjj=0; jjj<DEM.pixels(); ++jjj)
            DEM(iii,jjj) = DEMr8(iii,jjj);// cast to real4
        DEMr8.resize(1,1);// dealloc...
        INFO.print("Read crop of input DEM: format: REAL8.");
        break;
        }
      default:
        PRINT_ERROR("totally impossible, checked input.")
        //throw(unhandled_case_error);
      }


    // ----- Loop over DEM for stats ------------------------
    real8 min_dem_buffer =  100000.0;
    real8 max_dem_buffer = -100000.0;
    for (i=0; i<DEM.lines(); ++i)
      {
      // ----- Loop over oversampled matrix in x ------
      for (j=0; j<DEM.pixels(); ++j)
        {
        if(DEM(i,j)!=NODATA)
          {
          numvalid++;
          meancroppedDEM += DEM(i,j);// divide by numvalid later
          if (DEM(i,j)<min_dem_buffer) min_dem_buffer=DEM(i,j);//buffer
          if (DEM(i,j)>max_dem_buffer) max_dem_buffer=DEM(i,j);// stats
          }
        else
          {
          numNODATA++;
          }
        }//loop dem for stats
      }//loop dem for stats
    min_input_dem = min(min_input_dem,min_dem_buffer);//global stats
    max_input_dem = max(max_input_dem,max_dem_buffer);//global stats


    // ====== Radarcoding DEM ==============================
    // ______ DEM contains values from leftupper with ______
    // ______ spacing (DEMdeltalat,DEMdeltalon) ______
    // ______ Transform DEM to l,p,refphase ______
    PROGRESS.print("Converting DEM to radar system for this buffer.");
    const int32 NpointsDEM = DEM.size();
    const int32 NpixelsDEM = DEM.pixels();
    // ______ Extra info ______
    INFO << "Number of points in DEM: "
         << NpointsDEM;
    INFO.print();

    matrix<real8> masterDEMline(DEM.lines(),DEM.pixels());
    matrix<real8> masterDEMpixel(DEM.lines(),DEM.pixels());
    matrix<real8> theta_array(DEM.lines(),DEM.pixels());

    // --- Loop DEM ---
    cn P;
    real8 phi,lambda,height,l,p,theta; //,sim_amp; already a matrix see below [MA]


    phi = upperleftphi;
    for (i=0; i<DEM.lines(); ++i)
      {
      if ((i%100)==0)
        {
        // ______ Extra info ______
        PROGRESS << STEP << "Radarcoding buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " DEM line: " << i << " ("
             << floor(.5+(100.*real8(i)/real8(DEM.lines())))
             << "%)";
        PROGRESS.print();
        }

      lambda = upperleftlambda;
      for (j=0; j<DEM.pixels(); ++j)
        {
          
          height = DEM(i,j);
          ell2lp(l,p,ellips,master,masterorbit,phi,lambda,height,MAXITER,CRITERTIM);
          masterDEMline(i,j) = l; 
          masterDEMpixel(i,j) = p;
          P      = ellips.ell2xyz(phi,lambda,height);     // returns P(x,y,z)
          
          real8 t_range_master;
          real8 t_azi_master;
          xyz2t(t_azi_master,t_range_master,
              master,masterorbit,
              P,MAXITER,CRITERTIM);
          cn Psat_master = masterorbit.getxyz(t_azi_master);
          cn r1 = Psat_master.min(P);
          //theta = Psat_master.angle(r1); // theta look angle
          theta = P.angle(r1);           // theta incidence angle
         
// DEBUG << "Psat_master: " << Psat_master.x << " " << Psat_master.y << " " << Psat_master.z  
//          << " r1 (x,y,z): " << r1.x << " " << r1.y << " " << r1.z << " theta: " << theta ;
// DEBUG.print();
 
          theta_array(i,j) = theta;

          lambda += DEMdeltalon;
        } // loop DEM pixels
      
      // ______ update latitude of next line ______
      phi -= DEMdeltalat;           // upper left is max. value
      } // loop DEM lines
    

    // Write results to output files 
    PROGRESS << STEP << "Writing radar coded DEM and the incidence angle (theta) to files, buffer: " << buffer+1 << " of " << numfullbuffers + extrabuffer ;
    PROGRESS.print();
    demofile << DEM;
    masterdemlineoutfile << masterDEMline;
    masterdempixeloutfile << masterDEMpixel;
    thetaoutfile << theta_array;

    masterDEMline.resize(1,1); //deallocate
    masterDEMpixel.resize(1,1); //deallocate
    DEM.resize(1,1); //deallocate
    theta_array(1,1); //deallocate
    } // buffer loop

  demofile.close();
  masterdemlineoutfile.close();
  masterdempixeloutfile.close();
  thetaoutfile.close();
    }
  else
    {
    INFO << "Found completed sam_m_theta.temp!";
    INFO.print();
    INFO << "Using the existing file sam*.temp from the Radar coding step... ";
    INFO.print();
    }

  //===================================================================
  //============ End first loop: radarcode DEM ========================
  //============ (DEM geometry)            ============================
  //===================================================================


  //===================================================================
  //============ Second loop: interpolation               =============
  //============ (radar geometry)                         =============
  //===================================================================
  
  PROGRESS << STEP << "Start interpolation...";
  PROGRESS.print();

  // ______ Line/pixel of first point in original master coordinates ______
  // ______ maybe this should be changed to be x+(ml/2) ?? but depends on
  // ______ definition of range_to_first_bin is to center or start..
  // Bert Kampes, 08-Apr-2005: chose center by adding ml/2
  const int32 Nlinesml    = master.currentwindow.lines()  / mlL;
  const int32 Npixelsml   = master.currentwindow.pixels() / mlP;
  const real8 offset = 0;
  
  const real8 veryfirstline = real8(master.currentwindow.linelo) + (real8(mlL)-1.0)/2.0;
  const real8 verylastline = veryfirstline + real8((Nlinesml-1)*mlL);
  const real8 firstpixel = real8(master.currentwindow.pixlo) + (real8(mlP)-1.0)/2.0;
  const real8 lastpixel = firstpixel + real8((Npixelsml-1)*mlP);


  //Determine range-azimuth spacing ratio, needed for proper triangulation
  cn P1, P2 , P3, P4;
  lp2xyz(veryfirstline,firstpixel,ellips,master,masterorbit,
           P1,MAXITER,CRITERPOS);
  lp2xyz(veryfirstline,lastpixel,ellips,master,masterorbit,
           P2,MAXITER,CRITERPOS);
  lp2xyz(verylastline,firstpixel,ellips,master,masterorbit,
           P3,MAXITER,CRITERPOS);
  lp2xyz(verylastline,lastpixel,ellips,master,masterorbit,
           P4,MAXITER,CRITERPOS);

  const real8 r_spacing  = ( (P1.min(P2)).norm() + (P3.min(P4)).norm() ) / 2 /(lastpixel - firstpixel) ;
  const real8 az_spacing = ( (P1.min(P3)).norm() + (P2.min(P4)).norm() ) /2 /(verylastline - veryfirstline); 
  const real8 r_az_ratio = r_spacing/az_spacing;

  INFO << "Master azimuth spacing: " << az_spacing;
  INFO.print();
  INFO << "Master range spacing: " << r_spacing;
  INFO.print();
  INFO << "Range-azimuth spacing ratio: " << r_az_ratio;
  INFO.print();

  // ______ Compute buffer size interpolation______
  const real8 Nlinesml_possible = BUFFERMEMSIZE / (6*8*Npixelsml);
  bufferlines = int32(ceil(Nlinesml_possible));
  if (bufferlines > Nlinesml) bufferlines=Nlinesml;
  numfullbuffers = Nlinesml / bufferlines;
  restlines      = Nlinesml % bufferlines;
  extrabuffer = (restlines == 0) ? 0 : 1;

  // ______ Extra info ______
  INFO << "Interpolation in: " << numfullbuffers << " buffers of " 
       << bufferlines << " lines and " << extrabuffer << " extra buffer of "
       << restlines << " lines.";
  INFO.print();


  // ______ Open output files ______
  ofstream simampofile;
  openfstream(simampofile,simamp.fosimamp,generalinput.overwrit);
  bk_assert(simampofile,simamp.fosimamp,__FILE__,__LINE__);

  // if request for height in radar coordinates l,p
  ofstream samheiofile;                    // MA
  if (outputrefdemhei==true)
    {
    openfstream(samheiofile,simamp.fodemlp,generalinput.overwrit);
    bk_assert(samheiofile,simamp.fodemlp,__FILE__,__LINE__);
    }

  // if request for theta in radar coordinates l,p
  ofstream samthetaofile;                    // MA
  if (outputrefdemhei==true)
    {
    openfstream(samthetaofile,simamp.fothetalp,generalinput.overwrit);
    bk_assert(samthetaofile,simamp.fothetalp,__FILE__,__LINE__);
    }

  // ______ interpolation loop per buffer ______
  
  for (register int32 buffer = 0; buffer < numfullbuffers + extrabuffer; ++buffer)
    {

    // Determine indices for buffer
    const int32 blines = (buffer == numfullbuffers) ? restlines : bufferlines;
    const real8 firstline_buffer = veryfirstline+buffer*bufferlines*mlL;
    const real8 lastline_buffer = firstline_buffer+(blines-1)*mlL;

    // ______ Extra info ______
    PROGRESS << STEP << "Interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer << " [l0:lN, p0:pN]: " << " ["
         << firstline_buffer << ": " << lastline_buffer << ", "
         << firstpixel << ": " << lastpixel << "]";
    PROGRESS.print();

    // Get corners of buffer
    real8 phimin_az;
    real8 phimax_az;
    real8 lambdamin_az;
    real8 lambdamax_az;
    getcorners(firstline_buffer+offset,lastline_buffer+offset,
             firstpixel+offset,lastpixel+offset,
             extralat,extralong,phimax,lambdamin,
             DEMdeltalat,DEMdeltalon,NrowsDEM,NcolsDEM,
             ellips,master,masterorbit,phimin_az,phimax_az,lambdamin_az,lambdamax,
             indexphi0DEM,indexphiNDEM,indexlambda0DEM,indexlambdaNDEM);

    window zerooffset  (0,0,0,0);
    window winfromfile (indexphi0DEM,indexphiNDEM,
                        indexlambda0DEM,indexlambdaNDEM);
    const int32 NrowsDEM_buffer = indexphiNDEM-indexphi0DEM+1;
    const int32 NcolsDEM_buffer = indexlambdaNDEM-indexlambda0DEM+1;

    PROGRESS << STEP << "Reading input for interpolation buffer: " << buffer+1 << "of" << numfullbuffers + extrabuffer;
    PROGRESS.print();

    // read x,y
    matrix<real8> DEMline_buffer(NrowsDEM_buffer,NcolsDEM_buffer);
    matrix<real8> DEMpixel_buffer(NrowsDEM_buffer,NcolsDEM_buffer);

    readfile(DEMline_buffer,"sam_m_demline.temp",NrowsDEM,winfromfile,zerooffset);
    readfile(DEMpixel_buffer,"sam_m_dempixel.temp",NrowsDEM,winfromfile,zerooffset);
    
    // read z (multiple, number can easily be increased, e.g. simulated intensity)
    int32 Nz = 2;                                       //number of z, matrix depth
    matrix<real8> input_buffer(NrowsDEM_buffer *Nz ,NcolsDEM_buffer);

    matrix<real8> temp_input_buffer(NrowsDEM_buffer,NcolsDEM_buffer);

    // Read in real4 dem heights and convert to real8 for computation [MA] TODO potential external function
    matrix<real4> dem_input(NrowsDEM_buffer,NcolsDEM_buffer);
    readfile(dem_input,simamp.fodem,NrowsDEM,winfromfile,zerooffset); // TODO modify this so it can handle conversion as well
    for (register int32 i =0 ; i < NrowsDEM_buffer ; i++)
      for(register int32 j = 0; j < NcolsDEM_buffer; j++)
        temp_input_buffer(i,j) = real8(dem_input(i,j));

    dem_input.resize(1,1); // deallocation

//    readfile(temp_input_buffer,simamp.fodem,NrowsDEM,winfromfile,zerooffset);
    input_buffer.setdata(0, 0, temp_input_buffer);

    readfile(temp_input_buffer,"sam_m_theta.temp",NrowsDEM,winfromfile,zerooffset);
    input_buffer.setdata(NrowsDEM_buffer * (Nz-1) , 0, temp_input_buffer);

    DEBUG << "NrowsDEM: " << NrowsDEM << " window " << indexphi0DEM << " " << indexphiNDEM << " " << indexlambda0DEM << " " << indexlambdaNDEM ;
    DEBUG.print();

    
    // initialize output array
    Nz = 2;  
    matrix<real8> output_buffer(blines * Nz, Npixelsml);

    // interpolation
    DEBUG << "fl_buf: " << firstline_buffer << " ll_buf: " << lastline_buffer 
             << " p: "  << firstpixel << " P: " <<  lastpixel
             << " mlL: "  <<  mlL << " mlP: " << mlP << " r_a_ratio: " << r_az_ratio 
             << " offset: " << offset << " nodata: "  << NODATA ;
    DEBUG.print();
    DEBUG << "theta: " << input_buffer(NrowsDEM_buffer * (Nz-1),0) << " =? " << temp_input_buffer(0,0) << " " << NrowsDEM << " " << simamp.fodem ;
    DEBUG.print();

    griddatalinear(DEMline_buffer,DEMpixel_buffer,input_buffer,
                   firstline_buffer,lastline_buffer,firstpixel,lastpixel,
                   mlL,mlP,r_az_ratio,offset,NODATA,output_buffer);

    DEBUG << "temp_input_buffer size: "  << temp_input_buffer.size()  << " l: " << temp_input_buffer.lines()  << " p: " << temp_input_buffer.pixels() << "\n\t:"
             << " input_buffer size    : "  << input_buffer.size()  << " l: " << input_buffer.lines()  << " p: " << input_buffer.pixels() << "\n\t:" 
             << " output_buffer size   : " << output_buffer.size() << " l: " << output_buffer.lines() << " p: " << output_buffer.pixels() ;
    DEBUG.print();

    DEBUG << "NrowsDEM_buffer,NcolsDEM_buffer: " <<  NrowsDEM_buffer  << "," << NcolsDEM_buffer ;
    DEBUG.print();
    DEBUG << "blines, Npixelsml: " <<  blines << "," << Npixelsml ;
    DEBUG.print();

    if (outputrefdemhei==true)
      {
        Nz = 1; // radarcoded topo
        matrix<real4> topo(blines,Npixelsml); //topo matrix in master coords. in local scope
        for (register int32 i = 0 ; i < blines ; i++)
          for(register int32 j = 0; j < Npixelsml; j++)
            topo(i,j)= real4(output_buffer((Nz-1)*blines+i,j)); // real8 --> real4
                                                   // could use window() but here needed real4 conversion
                                                   // samheiofile << output_buffer(window((Nz-1) * blines,Nz * blines - 1, 0, Npixelsml -1 ));
        samheiofile << topo;
      }

    if (outputthetalp==true)  // theta matrix dump
      {
        Nz = 2; // radarcoded theta
        matrix<real4> theta(blines,Npixelsml); //theta matrix in master coords. in local scope
        for (register int32 i = 0 ; i < blines ; i++)
          for(register int32 j = 0; j < Npixelsml; j++)
            theta(i,j)= real4(output_buffer((Nz-1)*blines+i,j)); // real8 --> real4
        samthetaofile << theta;
      }

//    real8 theta = output_buffer(blines+int(blines/2),int(Npixelsml/2)); // is it constant if so const

    const real8 slant_range = SOL/(2*master.rsr2x);
    real8 grad, alpha, local_inc_angle, theta, ground_range;

    matrix<real4> sim_amp(blines,Npixelsml); //sim_amp matrix in master coords.

    Nz=2;
    for (register int16 i = 0 ; i < blines ; i ++)
      for(register int16 j = 0; j < Npixelsml; j++)
        if (j < Npixelsml-1)
          {
            theta = output_buffer((Nz-1)*blines+i,j);           // incidence angle
            ground_range = slant_range/sin(theta);
            grad = output_buffer(i,j+1) - output_buffer(i,j);   // hei gradient
            alpha = atan(grad/ground_range);                    // slope
            local_inc_angle = theta-alpha;
            //sim_amp(i,j) = (real4) alpha;
            //sim_amp(i,j) = (real4) local_inc_angle;
            //sim_amp(i,j) = (real4) cos(local_inc_angle)/sin(local_inc_angle); // see eineder2003
            // sim_amp(i,j) = (real4) sin(-local_inc_angle);   // used
            //sim_amp(i,j) = (real4) sin(-local_inc_angle)+1;   // used, +1 shift to positive range of values.
            //sim_amp(i,j) = (real4) pow(sim_amp(i,j),2);      // intensity to be tested
            if ( output_buffer(i,j) == 0 )        // search for sea where hei ~ 0; [TODO] introduce external mask
              {
               sim_amp(i,j) = 0;
              }
            else
              {
               sim_amp(i,j) = (real4) sin(-local_inc_angle)+1;   // used, +1 shift to positive range of values.
               //sim_amp(i,j) = (real4) pow((sin(-local_inc_angle)+1),2);      // intensity 
              }
          }
        else
          {
            sim_amp(i,j) = sim_amp(i,j-1);  // cp previous coln as last coln
          }

    // Write results to output files
    PROGRESS << STEP << "Writing radar coded SIM_AMPLITUDE to file, buffer#: " << buffer+1 << " of " << numfullbuffers + extrabuffer ;
    PROGRESS.print();
    simampofile << sim_amp;

  
    DEMline_buffer.resize(1,1);      //deallocate
    DEMpixel_buffer.resize(1,1);
    input_buffer.resize(1,1);
    temp_input_buffer.resize(1,1);
    output_buffer.resize(1,1);
    sim_amp.resize(1,1);
  
    } //end loop 

  INFO << "Closing output files";
  INFO.print();

  simampofile.close();
  if (outputrefdemhei==true)   samheiofile.close();
  if (outputthetalp==true) samthetaofile.close();

  //if(!remove(fsuccess.c_str()))  // when success report removed.
  //  {
  //  DEBUG << STEP << "Removed existing " << fsuccess << "file.";
  //  DEBUG.print();
  //  } ;

  const string assertmsg = STEP + fsuccess; 
  ofstream scratchfile(fsuccess.c_str(), ios::out | ios::trunc); // [MA]
  //bk_assert(scratchfile, strcat("SAM: ", fsuccess.c_str()) ,__FILE__,__LINE__); // this does not work
  bk_assert(scratchfile, assertmsg.c_str() ,__FILE__,__LINE__);
  //bk_assert(scratchfile, fsuccess.c_str() ,__FILE__,__LINE__);
  scratchfile << "sam_m_theta.temp succeeded!" << endl;
  scratchfile.close();
 

  // ====== Write output information ======
  //*char croppeddemi[ONE27];
  //*strcpy(croppeddemi,"NO output requested");
  //*if (outputdemi) strcpy(croppeddemi,simamp.fodemi);
  INFO << "Min. value of input DEM covering sim. amplitude (master): " << min_input_dem; 
  INFO.print();
  INFO << "Max. value of input DEM covering sim. amplitude (master): " << max_input_dem; 
  INFO.print();

  ofstream scratchlogfile("scratchlogsimamp", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"simamp: scratchlogsimamp",__FILE__,__LINE__);
  scratchlogfile
    << "\n*******************************************************************"
    << "\n* " << processcontrol[pr_m_simamp]
    << "\n*******************************************************************"
    << "\n1) DEM source file:                   \t" <<  simamp.firefdem
    << "\nFormat:                               \t";
    switch (simamp.iformatflag)
      {
      case FORMATI2:
        {
        scratchlogfile << "SHORT SIGNED INTEGER (HOST ENDIANNESS)";
        break;
        }
      case FORMATI2_BIGENDIAN:
        {
        scratchlogfile << "SHORT SIGNED INTEGER, BIG ENDIAN";
        break;
        }
      case FORMATR4:
        {
        scratchlogfile << "REAL4 SIGNED FLOAT";
        break;
        }
      case FORMATR8:
        {
        scratchlogfile << "REAL8 SIGNED DOUBLE";
        break;
        }
      default:
        {
        scratchlogfile << "UNKNOWN? IMPOSSIBLE...";
        break;
        }
      }
  scratchlogfile
    << "\nByte order:                           \t" <<  "check it yourself..."
    << "\nNumber of lines:                      \t" <<  numberoflatpixels
    << "\nNumber of pixels:                     \t" <<  numberoflonpixels
    << "\nResolution latitude:                   \t" << rad2deg(DEMdeltalat) << " [deg]"
    << "\nResolution longitude:                  \t" << rad2deg(DEMdeltalon) << " [deg]"
    << "\nMost West point in input DEM:         \t" << rad2deg(lon0file)
    << "\nMost East point in input DEM:         \t" << rad2deg(lonNfile)
    << "\nMost South point in input DEM:        \t" << rad2deg(latNfile)
    << "\nMost North point in input DEM:        \t" << rad2deg(lat0file)
    << "\nMin. value of input DEM covering master: " << min_input_dem
    << "\nMax. value of input DEM covering master: " << max_input_dem
    << "\n2) Output file cropped DEM:           \t" <<  simamp.fodem 
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):         \t" << NrowsDEM 
    << "\nNumber of pixels (multilooked):        \t" << NcolsDEM
    << "\nDEM extend w/e/s/n            :      \t" << rad2deg(lambdamin) << "/" 
    << rad2deg(lambdamax) << "/" << rad2deg(phimin)    << "/" << rad2deg(phimax)
//    << "\nMean value:                           \t" <<  meancroppedDEM
//    << "\n3) Output file interpolated crop DEM: \t" <<  croppeddemi
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):        \t" << Nlinesml
    << "\nNumber of pixels (multilooked):       \t" << Npixelsml
    << "\n3) Output file synthetic amplitude:   \t" <<  simamp.fosimamp
    << "\nFormat:                               \t" << "REAL4"
    << "\nByte order:                           \t" << "(same as host)"
    << "\nNumber of lines (multilooked):        \t" <<  Nlinesml
    << "\nNumber of pixels (multilooked):       \t" <<  Npixelsml;
    if (outputrefdemhei==true)
      {
      scratchlogfile
        << "\n4) Output file radarcoded DEM:        \t" <<  simamp.fodemlp
        << "\nFormat:                               \t" << "REAL4"
        << "\nByte order:                           \t" << "(same as host)"
        << "\nNumber of lines (multilooked):        \t" <<  Nlinesml
        << "\nNumber of pixels (multilooked):       \t" <<  Npixelsml;
      }
    if (outputthetalp==true) 
      {
      scratchlogfile
        << "\n5) Output file local_inc angle theta: \t" <<  simamp.fothetalp
        << "\nFormat:                               \t" << "REAL4"
        << "\nByte order:                           \t" << "(same as host)"
        << "\nNumber of lines (multilooked):        \t" <<  Nlinesml
        << "\nNumber of pixels (multilooked):       \t" <<  Npixelsml;
      }
    scratchlogfile
//    << "\nDeltaline_slave00_dem:                    \t" << deltaline_slave00  // [MA no gridlinear so no corners]
//    << "\nDeltapixel_slave00_dem:                   \t" << deltapixel_slave00 // [? put getcorners instead]
//    << "\nDeltaline_slave0N_dem:                    \t" << deltaline_slave0N
//    << "\nDeltapixel_slave0N_dem:                   \t" << deltapixel_slave0N
//    << "\nDeltaline_slaveN0_dem:                    \t" << deltaline_slaveN0
//    << "\nDeltapixel_slaveN0_dem:                   \t" << deltapixel_slaveN0
//    << "\nDeltaline_slaveNN_dem:                    \t" << deltaline_slaveNN
//    << "\nDeltapixel_slaveNN_dem:                   \t" << deltapixel_slaveNN
    << "\n*******************************************************************\n\n";
  scratchlogfile.close();


  ofstream scratchresfile("scratchressimamp", ios::out | ios::trunc);
  bk_assert(scratchresfile,"simamp: scratchressimamp",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_m_simamp]
    << "\n*******************************************************************";
  scratchresfile
    << "\nDEM source file:                      \t" << simamp.firefdem
    << "\nMin. of input DEM:                    \t" << min_input_dem
    << "\nMax. of input DEM:                    \t" << max_input_dem
    << "\nData_output_file:                     \t" << simamp.fosimamp
    << "\nData_output_format:                   \t" << "real4";
    if (outputrefdemhei==true)
      scratchresfile << "\n" << STEPlog << "_outdem_file:                     \t" <<  simamp.fodemlp;
    if (outputthetalp==true) 
      scratchresfile << "\n" << STEPlog << "_theta_file:                   \t" <<  simamp.fothetalp;
  scratchresfile
    << "\nFirst_line (w.r.t. original_master):  \t"
    <<  master.currentwindow.linelo
    << "\nLast_line (w.r.t. original_master):   \t"
    <<  master.currentwindow.linehi
    << "\nFirst_pixel (w.r.t. original_master): \t"
    <<  master.currentwindow.pixlo 
    << "\nLast_pixel (w.r.t. original_master):  \t"
    <<  master.currentwindow.pixhi
    << "\nMultilookfactor_azimuth_direction:    \t" <<  mlL
    << "\nMultilookfactor_range_direction:      \t" <<  mlP
    << "\nNumber of lines (multilooked):        \t" <<  Nlinesml
    << "\nNumber of pixels (multilooked):       \t" <<  Npixelsml
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_m_simamp] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();


  } // END simulate_amp;
