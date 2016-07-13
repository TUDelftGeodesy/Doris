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
 * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/filtering.cc,v $  *
 * $Revision: 3.19 $                                            *
 * $Date: 2005/10/06 11:09:20 $                                 *
 * $Author: kampes $                                            *
 *                                                              *
 * Implementation of filtering routines:                        *
 * -rangefilter (adaptive, after resampling)                    *
 * -rangefilter (based on baseline, points on ellips)           *
 * -phasefilter: goldstein                                      *
 * -phasefilter: spatial convolution                            *
 * -phasefilter: spectral method                                *
 * -azimuthfilter: based on fDC polynomials                     *
 ****************************************************************/

// #define _LIB_VERSION _IEEE_

#include "matrixbk.hh"
#include "constants.hh"                 // global constants
#include "filtering.hh"                 // header file
#include "utilities.hh"                 // myhamming
#include "slcimage.hh"                  // my slc image class
#include "productinfo.hh"               // my 'products' class
#include "orbitbk.hh"                   // my orbit class
#include "exceptions.hh"                 // my exceptions class

#include <strstream>                    // for memory stream
#include <iomanip>                      // setp

#ifdef WIN32
  // Jia defined min/max here, I did it in constants.hh
  // Bert Kampes, 24-Aug-2005
  //#define max _MAX
  //#define min _MIN
#endif






/****************************************************************
 *    rangefilter                                               *
 * adaptive range filtering. To be performed after resampling.  *
 * peak estimation of local fringe frequency in spectral        *
 * domain of oversampled([2] or 4) complex interferogram.       *
 * then master and slave spectra are filtered, optionally       *
 * de- and re-weight with hamming.                              *
 *                                                              *
 * Input:                                                       *
 *  - files for master/slave, input options                     *
 * Output:                                                      *
 *  - filtered files for master/slave                           *
 *                                                              *
 *    Bert Kampes, 29-Mar-2000                                  *
 ****************************************************************/
void rangefilter(
        const input_gen       &generalinput,
        const slcimage        &master,
        const slcimage        &slave,
        const productinfo     &interferogram,
        const input_filtrange &filtrangeinput)
  {
  TRACE_FUNCTION("rangefilter (BK 29-Mar-2000)");
  INFO << "RANGEFILT: Master input file: " << master.file;
  INFO.print();
  INFO << "RANGEFILT: Slave input file:  " << slave.file;
  INFO.print();

  // ====== Handle input ======
  const real8 hamming      = filtrangeinput.hammingalpha;// alpha=1: const
  const int32 nlmean       = filtrangeinput.nlmean;      // odd number range lines to average
  const int32 fftlength    = filtrangeinput.fftlength;   // length of adaptive
  const int32 OVERLAP      = filtrangeinput.overlap;     // overlap input blocks
  const real8 SNRthreshold = filtrangeinput.SNRthreshold;// criterium to filter
  const bool doweightcorrel= filtrangeinput.doweightcorrel;
  const int32 oversamplefactor = filtrangeinput.oversample;

  // RSR,RBW in Hz
  const real8 RSR = 0.5*master.rsr2x;           // range sampling rate fr 18.96MHz
  const real8 RBW = master.rbw;                 // range band width fs 15.55MHz


// ______ Number of lines on file ______
// ______ Assume resampling has been done and interferowin contains master coord.
  //const int32 Mfilelines        = master.currentwindow.lines();
  //const int32 Ifilelines        = interferogram.win.lines();
  const int32 numpixels         = interferogram.win.pixels();

  // ______ Memory buffers ______
  //const uint  BUFFERMEMSIZE  = generalinput.memory;
  const real8  BUFFERMEMSIZE  = generalinput.memory;
  const int32 bytesperline   = numpixels * sizeof(complr4);
  const real4 numbigmatrices = 4.2;                                               // in/out Master/Slave +rest
  int32 bufferlines = int32(ceil( (BUFFERMEMSIZE/numbigmatrices)/bytesperline )); // numlines in buffer 


// ====== Open output files ======
  ofstream ofilemaster;
  openfstream(ofilemaster,filtrangeinput.fomaster,generalinput.overwrit);
  bk_assert(ofilemaster,filtrangeinput.fomaster,__FILE__,__LINE__);

  ofstream ofileslave;
  openfstream(ofileslave, filtrangeinput.foslave,generalinput.overwrit);
  bk_assert(ofileslave, filtrangeinput.foslave, __FILE__,__LINE__);


// ====== Compute indices, counters, etc. ======
  const int32 numrangeblocks  = int32(numpixels/(fftlength-2*OVERLAP)); // approx
  real8 overallSNR         = 0.;                        // statistics
  real8 overallnotfiltered = 0.;                        // statistics

  // ______ Buffer counters ______
  int32 azimuthbuffer      = 1;                         // loop counter
  bool azimuthdone         = false;
  int32 startlinethisblock = interferogram.win.linelo;
  window filteredpart      (0, bufferlines-uint((nlmean+1)/2),
                              0, numpixels-1);                  // output to disk

  // ______ Give info ______
  INFO << "overlap master slave: ["
       << interferogram.win.linelo << ":" << interferogram.win.linehi
       << "; " << interferogram.win.pixlo << ":" << interferogram.win.pixhi
       << "] (master system)";
  INFO.print();

  // ====== Loop over azimuth buffers over whole width ======
  // ______ loop forever, break when finished ______
  for (azimuthbuffer=1; azimuthbuffer<1e5; ++azimuthbuffer)
    {
    PROGRESS << "azimuthbuffer (" << bufferlines << " lines): " 
         << azimuthbuffer << ".";
    PROGRESS.print();

    // ====== Read data master slave for this block ======
    int32 endlinethisblock = startlinethisblock+bufferlines-1;
    // ______ check slightly larger last block ______
    if (endlinethisblock>=int32(interferogram.win.linehi)-nlmean-int32(bufferlines/10))
      {
      azimuthdone         = true;                               // finish this one
      endlinethisblock    = interferogram.win.linehi;           // last line on file
      bufferlines         = endlinethisblock-startlinethisblock+1;
      filteredpart.linehi = bufferlines - 1;    // write all at end
      }

    // ______ Allocate matrices and read window from file ______
    const window winfile(startlinethisblock, endlinethisblock,
                         interferogram.win.pixlo,
                         interferogram.win.pixhi);
    // ______ Input buffer ______
    const matrix<complr4> MASTER_ORIG = master.readdata(winfile);
    const matrix<complr4> SLAVE_ORIG  = slave.readdata(winfile);
    // ______ Output buffer ______
    matrix<complr4> MASTER_FILT(MASTER_ORIG.lines(),MASTER_ORIG.pixels());
    matrix<complr4> SLAVE_FILT (SLAVE_ORIG.lines(), SLAVE_ORIG.pixels());

    // ====== Loop over range buffers dividing width ======
    // ______ Indices are in inbuffer local coordinates [0:numpixs] ______
    int32 startpixel    = 0;                    // valid for first range block
    int32 filteredpixlo = 0;                    // valid for first range block
    int32 partpixlo     = 0;                    // local system in PART
    int32 partpixhi     = fftlength-1-OVERLAP;  // local system in PART
    bool  rangedone     = false;                // not finished yet
    for (;;)// forever
      {
      int32 endpixel      = startpixel+fftlength-1;     // BUFFER -> PART
      int32 filteredpixhi = endpixel - OVERLAP;         // PART   -> BUFFER
      if (endpixel>=int32(MASTER_ORIG.pixels())-1)      // last block back-shifted
        {
        rangedone     = true;
        endpixel      = MASTER_ORIG.pixels()-1;
        startpixel    = endpixel - fftlength + 1;
        filteredpixhi = endpixel;               // write all for last block
        partpixhi     = fftlength-1;            // write all for last block
        partpixlo     = partpixhi - (filteredpixhi-filteredpixlo);
        }

      // ______ Some info ______
      TRACE << "range block [" << startlinethisblock << ":" << endlinethisblock
           << "; " << startpixel+interferogram.win.pixlo << ":"
                   << endpixel+interferogram.win.pixlo << "]";
      TRACE.print();

      // ______ Cut out master/slave for adaptive range filtering ______
      const window win(0,MASTER_ORIG.lines()-1,startpixel,endpixel);
      matrix<complr4> PARTMASTER(win,MASTER_ORIG);              // construct as part
      matrix<complr4> PARTSLAVE (win,SLAVE_ORIG);               // construct as part

      // ______ Do the actual range filtering ______
      // ______ filtered lines: index[0:L-1]: [(nlmean-1)/2:L-1-(nlmean-1)/2] ______
      real8 SNRmean, percentnotfiltered;
      rfilterblock(PARTMASTER,PARTSLAVE,                        // returned
                   nlmean,SNRthreshold,RSR,RBW,                 // parameters
                   hamming,oversamplefactor,doweightcorrel,     // control algorithm
                   SNRmean, percentnotfiltered);                // statistics returned
      overallSNR         += SNRmean;
      overallnotfiltered += percentnotfiltered;

      // ______ Put filtered part in outputbuffer ______
      const window winpart(win.linelo,win.linehi, partpixlo,partpixhi);
      const window wintoset(win.linelo,win.linehi, filteredpixlo,filteredpixhi);
      MASTER_FILT.setdata(wintoset,PARTMASTER,winpart);
      SLAVE_FILT.setdata (wintoset,PARTSLAVE, winpart);

      // ______ Update counters and check if it is already time to go ______
      if (rangedone==true) break;
      partpixlo     = OVERLAP;                  // valid for all next blocks
      startpixel    = endpixel - 2*OVERLAP + 1; // get PART from inputbuffer
      filteredpixlo = startpixel + OVERLAP;     // valid for all blocks except first
      } // loop over range blocks

    // ====== Write filtered parts master,slave to output files ======
    switch (filtrangeinput.oformatflag)
      {
      case FORMATCR4:
        {
        writefile(ofilemaster,MASTER_FILT,filteredpart);
        writefile(ofileslave, SLAVE_FILT, filteredpart);
        break;
        }
      // ______ Convert first to ci2 before writing to file ______
      case FORMATCI2:
        {
        matrix <compli16> TMP(MASTER_FILT.lines(),MASTER_FILT.pixels());
        for (uint ii=filteredpart.linelo; ii<=filteredpart.linehi; ++ii)
          for (uint jj=filteredpart.pixlo; jj<=filteredpart.pixhi; ++jj)
            TMP(ii,jj) = cr4toci2(MASTER_FILT(ii,jj));
        writefile(ofilemaster,TMP,filteredpart);
        for (uint ii=filteredpart.linelo; ii<=filteredpart.linehi; ++ii)
          for (uint jj=filteredpart.pixlo; jj<=filteredpart.pixhi; ++jj)
            TMP(ii,jj) = cr4toci2(SLAVE_FILT(ii,jj));
        writefile(ofileslave, TMP, filteredpart);
        break;
        }
      default:
        {
        PRINT_ERROR("Totally impossible, checked input.")
        throw(unhandled_case_error);
        }
      }
    // ______ Update counters, check finished ______
    if (azimuthdone==true) break;
    filteredpart.linelo = uint((nlmean-1.)/2.); // same for all buffers except first
    startlinethisblock  = endlinethisblock - nlmean + 2; 
    } // loop over azimuth buffers

  // _______ Some stats, better dump SUN rasterfile all shifts ______
  overallSNR         /= (azimuthbuffer*numrangeblocks);
  overallnotfiltered /= (azimuthbuffer*numrangeblocks);
  INFO << "Mean SNR for all blocks (approx): " << overallSNR;
  INFO.print();
  INFO << "Filtered (approx): " << setprecision(3) 
       << 100.00-overallnotfiltered << "%";
  if (overallnotfiltered<60.) 
    {
    INFO.print();
    }
  else
    {
    WARNING.print(INFO.get_str());
    INFO.rewind();
    }


// ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltrange", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtrange: scratchlogfiltrange",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* RANGE FILTERING ADAPTIVE FOR MASTER AND SLAVE..."
    << "\n*******************************************************************"
    << "\nInput file master (format): \t\t\t"
    <<  master.file << " " << master.formatflag
    << "\nInput file slave (format): \t\t\t"
    <<  slave.file << " " << slave.formatflag
    << "\nOutput file filtered master (format): \t\t\t"
    <<  filtrangeinput.fomaster << " ";
    if (filtrangeinput.oformatflag==FORMATCR4)
      scratchlogfile << "complex r4";
    if (filtrangeinput.oformatflag==FORMATCI2)
      scratchlogfile << "complex short int";
  scratchlogfile
    << "\nOutput file filtered slave (format): \t\t\t"
    <<  filtrangeinput.foslave << " " << FORMATCR4
    << "\n..."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  for (int32 fileID=1; fileID<=2; ++fileID)
    {
    char oresfile[EIGHTY];
    char odatafile[EIGHTY];
    char odataformat[EIGHTY];
    char processcf[EIGHTY];
    if (filtrangeinput.oformatflag==FORMATCR4)
      strcpy(odataformat,"complex_real4");
    else if (filtrangeinput.oformatflag==FORMATCI2)
      strcpy(odataformat,"complex_short");
    else
      {
      PRINT_ERROR("case problem.")
      throw(unhandled_case_error);
      }
    switch (fileID)
      {
      case 1:                                                   // master
        strcpy(oresfile,"scratchresMfiltrange");
        strcpy(odatafile,filtrangeinput.fomaster);
        strcpy(processcf,processcontrol[pr_m_filtrange]);       // control flag
        break;
      case 2:                                                   // slave
        strcpy(oresfile,"scratchresSfiltrange");
        strcpy(odatafile,filtrangeinput.foslave);
        strcpy(processcf,processcontrol[pr_s_filtrange]);       // control flag
        break;
      default:
        {
        PRINT_ERROR("panic: ID!={1,2}")
        throw(unhandled_case_error);
        }
      }

    // ______ updateproductinfo greps info from this file ______
    ofstream scratchresfile(oresfile, ios::out | ios::trunc);
    bk_assert(scratchresfile,oresfile,__FILE__,__LINE__);
    scratchresfile
      << "\n\n*******************************************************************"
      << "\n*_Start_" << processcf
      << "\n*******************************************************************"
      << "\nMethod_rangefilt:                     \t"
      << "adaptive"
      << "\nData_output_file:                     \t"
      <<  odatafile
      << "\nData_output_format:                   \t"
      <<  odataformat
      << "\nFirst_line (w.r.t. original_master):  \t"
      <<  interferogram.win.linelo
      << "\nLast_line (w.r.t. original_master):   \t"
      <<  interferogram.win.linehi
      << "\nFirst_pixel (w.r.t. original_master): \t"
      <<  interferogram.win.pixlo
      << "\nLast_pixel (w.r.t. original_master):  \t"
      <<  interferogram.win.pixhi
      << "\n*******************************************************************"
      << "\n* End_" << processcf << "_NORMAL"
      << "\n*******************************************************************\n";
    scratchresfile.close();
    } // for master and slave
  // ______Tidy up______
  } // END rangefilter



/****************************************************************
 *    rfilterblock                                              *
 * Computes powerspectrum of complex interferogram.             *
 * (product oversampled master. conj oversampled slave in range)*
 * A peak in this spectrum corresponds to frequency shift.      *
 * The master and slave are LPF filtered for this shift.        *
 *                                                              *
 * Optionally the oversampling can be turned off, since no use  *
 * if only small baseline, and flat terrain.                    *
 * The powerspectrum can be weighted to give more influence to  *
 * higher frequencies (conv. of 2 blocks, should be hamming).   *
 * The peak is detected by taking a mean of nlmean lines (odd). *
 * Filtering is applied if the SNR (N*power peak / power rest)  *
 * is above a user supplied threshold.                          *
 * At LPF filtering of the master/slave a hamming window may    *
 * be applied first to deweight, then to re-weight the spectrum *
 *                                                              *
 * Should filter based on zero terrain slope if below SNR,      *
 * but this requried knowledge of orbits, pixel,line coordinate *
 #%// BK 13-Nov-2000                                            *
 *                                                              *
 * Input:                                                       *
 *  - MASTER: block of master, that will be filtered            *
 *  - SLAVE:  block of slave, that will be filtered             *
 * Output:                                                      *
 *  - MASTER (SLAVE): filtered from indeces[0:numl-1]           *
 *    (nlmean-1)/2 to numlines-(nlmean-1)/2-1                   *
 *                                                              *
 *    Bert Kampes, 29-Mar-2000                                  *
 ****************************************************************/
void rfilterblock(
        matrix<complr4> &master,                // updated
        matrix<complr4> &slave,                 // updated
        int32 nlmean,                           // number of lines to take mean over
        real8 SNRthreshold,                     // 
        real8 RSR,                              // in MHz
        real8 RBW,                              // in MHz
        real8 alphahamming,                     // parameter hamming filtering [0,1]
        int32  osfactor,                        // oversampling factor interf. gen.
        bool  doweightcorrel,                   // correct correl for #elem.
        real8 &meanSNR,                         // returned
        real8 &percentnotfiltered)              // returned
  {
  const int32 numlines    = master.lines();
  const int32 numpixs     = master.pixels();
  const int32 outputlines = numlines - nlmean + 1;
  const int32 firstline   = int32((nlmean-1)/2);        // indices in matrix system
  const int32 lastline    = firstline + outputlines - 1;
  const bool dohamming    = (alphahamming<0.9999) ? true : false;
  // use oversampling before int. gen.
  const bool dooversample = (osfactor!=1) ? true : false;
  int32 notfiltered=0;                                  // counter


#ifdef __DEBUG
  TRACE_FUNCTION("rfilterblock (BK 29-Mar-2000)")
//  if (!isodd(nlmean))                          // MA error check is moved to readinput.cc # rm lines for major release
//    {
//    PRINT_ERROR("nlmean has to be odd.")
//    throw(argument_error);
//    }
  if (!ispower2(numpixs))
    {
    PRINT_ERROR("numpixels (FFT) has to be power of 2.")
    throw(argument_error);
    }
  if (!ispower2(osfactor))
    {
    PRINT_ERROR("oversample factor (FFT) has to be power of 2.")
    throw(argument_error);
    }
  if (slave.lines()!=numlines)
    {
    PRINT_ERROR("slave not same size as master.")
    throw(argument_error);
    }
  if (slave.pixels()!=numpixs)
    {
    PRINT_ERROR("slave not same size as master.")
    throw(argument_error);
    }
#endif
  if (outputlines < 1)
    {
    WARNING.print("no outputlines, continuing."); 
    return;
    }


  // ______ Shift parameters ______
  register int32 i,j;
  const real4 deltaf =  RSR/real4(numpixs);
  const real4 fr     = -RSR/2.;
  matrix<real4> freqaxis(1,numpixs);
  for (i=0; i<numpixs; ++i)
    freqaxis(0,i) = fr+real4(i)*deltaf;
 
  matrix<real4> inversehamming;
  if (dohamming)
    {
    inversehamming = myhamming(freqaxis,RBW,RSR,alphahamming);
    for (i=0; i<numpixs; ++i)
      if (inversehamming(0,i)!=0.)
        inversehamming(0,i)= 1./inversehamming(0,i);
    }

// ______ Compute complex interferogram -> power ______
  matrix<complr4> cint;
  if (dooversample)
    cint = dotmult(oversample(master,1,osfactor),oversample(conj(slave),1,osfactor));
  else
    cint = dotmult(master,conj(slave));
  const int32 fftlength = cint.pixels();

  DEBUG.print("is real4 accurate enough?");// seems so
  fft(cint,2);                                  // cint=fft over rows
  matrix<real4> power   = intensity(cint);      // power=cint.*conj(cint);


  // ______ Use weighted correlation due to bias in normal definition ______
  // ______ Actually better deweight with autoconvoluted hamming.
  // ______ No use a triangle for #points used for correlation estimation
  // ______ not in combination with dooversample...
  if (doweightcorrel)
    {
    INFO.print("De-weighting power spectrum of convolution.");
    //matrix<real4> weighting(1,fftlength);
    //for (i=0; i<fftlength; ++i)
    //  weighting(0,i) = abs(i)...;
    //power *= weightingfunction...== inv.triangle

    // weigth = numpoints in spectral convolution for fft squared for power...
    const int32 indexnopeak = int32((1.-(RBW/RSR))*real4(numpixs));
    for (j=0; j<fftlength; ++j)
      {
      const int32 npnts = abs(numpixs-j);
      const real4 weight =                      // oversample: numpixs==fftlength/2
        (npnts<indexnopeak) ? sqr(numpixs) : sqr(npnts);        // ==zero
      for (i=0; i<numlines; ++i)
        {
        power(i,j) /= weight;
        }
      }
    }


  // ______ Average power to reduce noise ______
  fft(master,2);                                // master=fft over rows
  fft(slave,2);                                 // slave=fft over rows
  TRACE.print("Took FFT over rows of master, slave.");
  matrix<real4> nlmeanpower = sum(power(0,nlmean-1, 0,fftlength-1),1);

  uint shift      = 0;                          // returned by max
  uint dummy      = 0;                          // returned by max
  meanSNR         = 0.;
  real8 meanSHIFT = 0.;

  // ______ Start actual filtering ______
  // ______ oline is index in matrix system ______
  for (register int32 oline=firstline; oline<=lastline; ++oline)
    {
    matrix<real4> totalp   = sum(nlmeanpower,2);        // 1x1 matrix ...
    const real4 totalpower = totalp(0,0);
    const real4 maxvalue   = max(nlmeanpower,dummy,shift);      // shift returned
    uint lastshift         = shift;     // use this if current shift not ok.
    const real8 SNR        = fftlength*(maxvalue/(totalpower-maxvalue));
    meanSNR               += SNR;

    // ______ Check for negative shift ______
    bool negshift = false;
    if (shift > uint(fftlength/2))
      {
      // ERR: uint! shift    = uint(abs(real8(shift-fftlength)));
      shift     = uint(fftlength)-shift;
      lastshift = shift;// use this if current shift not OK.
      negshift  = true;
      }

    // ______ Do actual filtering ______
    if (SNR<SNRthreshold)
      {
      notfiltered++;                                    // update counter
      shift = lastshift;// but use this anyway.
      //WARNING.print("using last shift for filter");
      }
    meanSHIFT += shift;
    matrix<real4> filter;
    if (dohamming)
        {
        // ______ Newhamming is scaled and centered around new mean ______
        filter  = myhamming(freqaxis-real4(.5*shift*deltaf),
                            RBW-real8(shift*deltaf),
                            RSR,alphahamming);          // fftshifted
        filter *= inversehamming;
        }
      else                                              // no weighting of spectra
        {
        filter = myrect((freqaxis-real4(.5*shift*deltaf)) /
                          (real4(RBW-shift*deltaf)));   // fftshifted
        }
      // ______ Use freq. as returned by fft ______
      // ______ Note that filter_s = fliplr(filter_m) ______
      // ______ and that this is also valid after ifftshift ______
      ifftshift(filter);                                // fftsh works on data!

      // ====== Actual spectral filtering ======
      // ______ Decide which side to filter, may be dependent on ______
      // ______ definition of FFT, this is ok for VECLIB ______
      // ______ if you have trouble with this step, either check your FFT
      // ______ (or use uinternal one, or add a card for changing false to true below
      if (negshift==false)
        {
        dotmult(master[oline],filter,1); 
        filter.fliplr();                                // works on data!
        dotmult(slave[oline],filter,1);
        }
      else
        {
        dotmult(slave[oline],filter,1);
        filter.fliplr();                                // works on data!
        dotmult(master[oline],filter,1);
        }
// following is removed, we now always filter with last know spectral offset
// Bert Kampes, 13-Sep-2004
//      } // SNR>threshold
//    else
//      {
//      notfiltered++;                                  // update counter
//      }

    // ______ Update 'walking' mean ______
    if (oline!=lastline)                        // then breaks
      {
      matrix<real4> line1 = power.getrow(oline-firstline);
      matrix<real4> lineN = power.getrow(oline-firstline+nlmean);
      nlmeanpower += (lineN-line1);
      }
    } // loop over outputlines

  // ______ IFFT of spectrally filtered data, and return these ______
  ifft(master,2);                               // master=ifft over rows
  ifft(slave,2);                                // slave=ifft over rows

  // ______ Return these to main ______
  meanSHIFT          /= (outputlines-notfiltered);
  meanSNR            /= outputlines;
  percentnotfiltered  = 100. * (real4(notfiltered)/real4(outputlines));


  // ______ Some info for this block ______
  const real8 meanfrfreq = meanSHIFT*deltaf;    // Hz?
  DEBUG << "mean SHIFT for block"
       << ": "  << meanSHIFT
       << " = " << meanfrfreq/1e6 << " MHz (fringe freq.).";
  DEBUG.print();
  DEBUG << "mean SNR for block"
       << ": " << meanSNR;
  DEBUG.print();
  DEBUG << "filtered for block"
       << ": " << setprecision(3) << 100.00-percentnotfiltered << "%" << ends;
  DEBUG.print();
  if (percentnotfiltered>60.0) 
    {
    WARNING.print(DEBUG.get_str());
    DEBUG.reset();
    }
  } // END rfilterblock



/****************************************************************
 *    phasefilter                                               *
 * goldsteins method, see routine goldstein and smooth.         *
 * After Goldstein and Werner, Radar interferogram filtering    *
 * for geophysical applications. GRL 25-21 pp 4035-4038, 1998.  *
 * and: ESA Florence 1997, vol2, pp969-972, Goldstein & Werner  *
 * "Radar ice motion interferometry".                           *
 #%// BK 24-Oct-2000                                            *
 ****************************************************************/
void phasefilter(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput,
        const productinfo     &coherence)
  {
  TRACE_FUNCTION("phasefilter (BK 24-Oct-2000)")

  // ====== Handle input ======
  bool doexternalfile = false;
  if (specified(filtphaseinput.fifiltphase))
    doexternalfile = true;
  char infile[EIGHTY];                          // file 2b filtered
  char cohfile[EIGHTY];				// coherence file for modgoldstein
  strcpy(infile,interferogram.file);
  if (filtphaseinput.method == fp_modgoldstein )
    strcpy(cohfile,coherence.file);

  int32 numlinesinput  = int32(interferogram.win.lines()/interferogram.multilookL);
  int32 numpixelsinput = int32(interferogram.win.pixels()/interferogram.multilookP);
  if (doexternalfile)                           // overwrite default file
    {
    numlinesinput  = filtphaseinput.finumlines;
    strcpy(infile,filtphaseinput.fifiltphase);
    ifstream ifile;
    openfstream(ifile,infile);
    bk_assert(ifile,infile,__FILE__,__LINE__);
    ifile.seekg(0,ios::end);                    // pointer to end...
    const streamoff  sizefile     = ifile.tellg();   // opened ate, [MA] 64 bit pointer on a 32-bit system
    ifile.close();
    numpixelsinput = sizefile/sizeof(complr4)/numlinesinput;// floor
    if (streamoff(numlinesinput)*streamoff(numpixelsinput)*streamoff(sizeof(complr4)) != sizefile)
      WARNING.print("Format infile not CR4, or numlines not correct.");
    INFO.print("Using input file for phase filtering, not default.");
    }
  INFO << "phasefilter: inputfile: (#l,#p): " << infile << " ("
       << numlinesinput << "," << numpixelsinput << ")";
  INFO.print();

  // ______ Set variables ______
  const real8 ALPHA    = filtphaseinput.alpha;          // 0<a<1;
  const int32 SIZE     = filtphaseinput.blocksize;      // power of 2
  const int32 OVERLAP  = filtphaseinput.overlap;        // half of the overlap

  // ______ Open output file ______
  ofstream ofile;
  openfstream(ofile,filtphaseinput.fofiltphase,generalinput.overwrit);
  bk_assert(ofile,filtphaseinput.fofiltphase,__FILE__,__LINE__);


  // ______ initialize indices ______
  const int32 numout   = SIZE-(2*OVERLAP);      // number of output lines
  int32 cintlinelo     = 0;                     // index in CINT to get 1st buffer
  int32 cintlinehi     = SIZE-1;                // index in CINT to get 1st buffer
  int32 filteredlinelo = 0;                     // index in FILTERED to write (1st)
  int32 filteredlinehi = SIZE-OVERLAP-1;        // index in FILTERED to write (1st)
  bool  lastbuffer     = false;                 // not yet, just starting...

  if (!doexternalfile)
    if (interferogram.formatflag!=FORMATCR4)
      {
      PRINT_ERROR("Sorry, for phasefilter, format of interferogram must be CR4")
      throw(argument_error);
      }
  matrix<complr4> CINT(SIZE,numpixelsinput);    // type must be cr4!
//  if (filtphaseinput.method == fp_modgoldstein )
//    {
    matrix<real4> COH(SIZE,numpixelsinput);    // type must be r4!
//    }

  int32 lineswritten     = OVERLAP;             // first buffer extra OVERLAP lines
  const int32 numbuffers = numlinesinput/numout;        // approx?
  int32 tenpercent = int32((numbuffers/10)+.5); // round
  if (tenpercent==0) tenpercent = 1000;
  PROGRESS.print("FILTPHASE:  0%");
  int32 percent          = 10;

  // ====== Start filter loop buffers of BLOCKSIZE ======
  for (int32 buffer=1; buffer<1e6; ++buffer)
    {
    // ______ Check if this will be the last buffer ______
    if (cintlinehi >= numlinesinput-1)          // -1 since in matrix system
      {
      lastbuffer     = true;
      cintlinehi     = numlinesinput-1;
      cintlinelo     = cintlinehi-SIZE+1;       // make sure SIZE lines are read
      filteredlinehi = SIZE-1;                  // write upto lastline for last buffer
      const int32 lines2bwritten = numlinesinput-lineswritten;
      filteredlinelo = filteredlinehi - lines2bwritten + 1;
      if (lines2bwritten<1) 
        WARNING.print("PANIC: this will crash, lines2bwritten<1?");
      }

    // ______ Read in buffers of BLOCKSIZE lines complex interferogram ______
    const window windummy (0,0,0,0);            // no offset in readfile
    const window wincint  (cintlinelo,cintlinehi,0,numpixelsinput-1);
    readfile(CINT,infile,numlinesinput,wincint,windummy);
    matrix<complr4> FILTERED(SIZE,numpixelsinput);
    if (filtphaseinput.method == fp_goldstein )
      {
      // ______ Filter data in buffer of BLOCKSIZE lines ______
      // ______ output has same size, but write only part to disk ______
      FILTERED = goldstein(CINT,ALPHA,OVERLAP,filtphaseinput.kernel);
      }
    else // filtphaseinput.method == fp_modgoldstein    
      {
      readfile(COH,cohfile,numlinesinput,wincint,windummy);
      FILTERED = modgoldstein(CINT,COH,OVERLAP,filtphaseinput.kernel);
      }
    // ______ Write filtered data ______
    const window winfiltered(filteredlinelo,filteredlinehi,0,numpixelsinput-1);
    writefile(ofile,FILTERED,winfiltered);

    // ______ Check if all done ______
    if (lastbuffer)
      break;

    // ______ Update indexes in matrices, will be corrected for last block ______
    lineswritten  += numout;            // first block overlap added
    cintlinelo    += numout;            // next buffer
    cintlinehi    += numout;            // next buffer
    filteredlinelo = OVERLAP;           // index in buffer,
                                        // valid for all middle buffers, but not last
    // ______ Give progress message ______
    if (buffer%tenpercent==0)
      {
      PROGRESS << "FILTPHASE: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    } // end loop buffers
  ofile.close();

  // ______ exit if only external file was desired ______
  if (doexternalfile)
    {
    cerr  << "\n Finished external file phasefilter, Exiting\n";
    PROGRESS.print("\n Finished external file phasefilter, Exiting\n");
    exit(0);
    }

  // ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltphase", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtphase: scratchlogfiltphase",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* PHASE FILTER COMPLEX INTERFEROGRAM: \t"
    << "\n*******************************************************************"
    <<  infile
    << "\nOutput file filtered master (format): \t\t\t"
    <<  filtphaseinput.fofiltphase << " "
    << "\n..."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  // ______ updateproductinfo greps info from this file ______
  ofstream scratchresfile("scratchresfiltphase", ios::out | ios::trunc);
  bk_assert(scratchresfile,"scratchresfiltphase",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_filtphase]
    << "\n*******************************************************************";
  if (filtphaseinput.method == fp_goldstein ){   
    scratchresfile 
    << "\nMethod_phasefilt: goldstein: size, alpha, overlap: \t"
    <<  SIZE << " " << ALPHA << " " << OVERLAP;
    }
  else{
    scratchresfile 
    << "\nMethod_phasefilt: modgoldstein: size, overlap: \t"    
    <<  SIZE << " " << OVERLAP;
    }
  scratchresfile << "\n1D Smoothing kernel for |spectrum|:   \t";
  for (uint ii=0; ii<filtphaseinput.kernel.pixels(); ++ii)
    scratchresfile << " " << filtphaseinput.kernel(0,ii);
  scratchresfile
    << "\nInput_file:                           \t"
    <<  infile
    << "\nData_output_file:                     \t"
    <<  filtphaseinput.fofiltphase
    << "\nData_output_format:                   \t"
    <<  "complex_real4"
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
    <<  numlinesinput
    << "\nNumber of pixels (multilooked):       \t"
    <<  numpixelsinput
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_filtphase] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  // ====== Tidy up ======
  } // END phasefilter


/****************************************************************
 *    phasefilter goldstein                                     *
 * Input is matrix of SIZE (e.g. 32) lines, and N range pixels. *
 * Filtered OUTPUT is same size as input block.                 *
 * Because of overlap, only write to disk in calling routine    *
 * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
 *                                                              *
 * Smoothing of the amplitude of the spectrum is performed by   *
 * spatial convolution with a block kernel of size 2*SMOOTH+1.  *
 * (Which is done by FFT's). e.g. a spatial moving average with *
 * kernel (1d) k=[1 1 1 1 1]/5; kernel2d = transpose(k)*k.      *
 * Blocks in range direction.                                   * 
 *                                                              *
 * After Goldstein and Werner, Radar interferogram filtering    *
 * for geophysical applications. GRL 25-21 pp 4035-4038, 1998.  *
 * and: ESA Florence 1997, vol2, pp969-972, Goldstein & Werner  *
 * "Radar ice motion interferometry".                           *
 #%// BK 25-Oct-2000                                            *
 ****************************************************************/
matrix<complr4> goldstein(
        const matrix<complr4> &CINT,
        const real4            ALPHA,
        const int32            OVERLAP,
        const matrix<real4>   &smoothkernel)    // lying down
  {
  TRACE_FUNCTION("ROUTINE: goldstein (BK 25-Oct-2000)")
  //#define CHECKINDICESCALLING
  #ifdef CHECKINDICESCALLING
    return CINT;
  #else
  // ______ Allocate output matrix ______
  const int32 SIZE = CINT.lines();
  const int32 NPIX = CINT.pixels();
  matrix<complr4> FILTERED(SIZE,NPIX);          // output

  // ______ Get block from buffer ______
  const int32 numout  = SIZE-(2*OVERLAP);       // number of output pixels
  int32 cintpixlo     = 0;                      // index in CINT to get 1st block
  int32 cintpixhi     = SIZE-1;                 // index in CINT to get 1st block
  int32 outblockpixlo = 0;                      // index in BLOCK (only 1st block)
  int32 outblockpixhi = SIZE-1-OVERLAP;         // index in BLOCK (except last block)
  int32 outpixlo      = outblockpixlo;          // index in FILTERED (1st block)
  int32 outpixhi      = outblockpixhi;          // index in FILTERED
  bool  lastblockdone = false;                  // only just started...
  // note that int32() floors division 
  const int32 SMOOTH  = int32(smoothkernel.pixels())/2; // half block size, odd kernel
  const bool dosmooth = (SMOOTH==0) ? false : true;
  DEBUG << "SMOOTH: " << SMOOTH;// problem with uint<0 index in smoothkernel
  DEBUG.print();

  // ______ use FFT's for convolution with smoothkernel ______
  // ______ this could also be done static, or in the calling routine ______
  // ______ KERNEL2D is FFT2 of even kernel (no imag part after fft!) ______
  matrix<complr4> KERNEL2D;
  if (dosmooth==true)
    {
    matrix<complr4> kernel(1,SIZE);             // init to zeros
    for (register int32 ii=-SMOOTH; ii<=SMOOTH; ++ii)// 1d kernel function of block
      {
      //kernel(0,(ii+SIZE)%SIZE) = smoothkernel(0,ii-SMOOTH);
      // BK 07-Apr-2003
      // Cygwin fails on index very large number, i.e., uint problem somewhere
      // e.g.: [30,31,0,1,2] <--> [0,1,2,3,4]
      int32 tmp1 = (ii+SIZE)%SIZE;
      int32 tmp2 = ii+SMOOTH;// used to be ii-SMOOTH: wrong
      DEBUG << "tmp1: " << tmp1 << "; tmp2: " << tmp2;
      DEBUG.print();
      kernel(0,tmp1) = complr4(smoothkernel(0,tmp2),real4(0.0));
      }
    KERNEL2D = matTxmat(kernel,kernel);
    fft2d(KERNEL2D);                            // should be real sinc
    }
  DEBUG.print("kernel created for smoothing spectrum");

  // ====== Loop forever, stop after lastblockdone ======
  for (;;)      //forever
    {
    if (cintpixhi>=NPIX-1)                      // check if we are doing the last block
      {
      lastblockdone = true;
      cintpixhi     = NPIX-1;                   // prevent reading after file
      cintpixlo     = cintpixhi-SIZE+1;         // but make sure SIZE pixels are read
      outpixhi      = cintpixhi;                // index in FILTERED 2b written
      outblockpixhi = SIZE-1;                   // write all to the end
      outblockpixlo = outblockpixhi - (outpixhi-outpixlo+1) + 1;
      }
    const window wincint     (0,SIZE-1,cintpixlo,cintpixhi);
    const window winblock    (0,SIZE-1,outblockpixlo,outblockpixhi);
    const window winfiltered (0,SIZE-1,outpixlo,outpixhi);

    // ______ Construct BLOCK as part of CINT ______
    matrix<complr4> BLOCK(wincint,CINT);

    //#define CHECKINDEXONLY
    #ifndef CHECKINDEXONLY
    // ______ Get spectrum/amplitude/smooth/filter ______
    fft2d(BLOCK);
    matrix<real4> AMPLITUDE = magnitude(BLOCK);

    // ______ use FFT's for convolution with rect ______
    if (dosmooth==true)
      AMPLITUDE = smooth(AMPLITUDE,KERNEL2D);

    //dumpasc("A",AMPLITUDE);
    //AMPLITUDE = smooth(AMPLITUDE,SMOOTH);
    //dumpasc("As",AMPLITUDE); exit(1);
    const real4 maxamplitude = max(AMPLITUDE);
    if (maxamplitude>1e-20) //?
      {
      AMPLITUDE /= maxamplitude;
      AMPLITUDE.mypow(ALPHA);
      #ifdef WIN32
      BLOCK = timesCxR(BLOCK,AMPLITUDE);// weight spectrum
      #else
      BLOCK     *= AMPLITUDE;           // weight spectrum
      #endif
      }
    else
      {
      WARNING.print("no filtering, maxamplitude<1e-20, zeros in this block?");
      }
    ifft2d(BLOCK);
    #endif // check index blocks
    // ______ Set correct part that is filtered in output matrix ______
    FILTERED.setdata(winfiltered,BLOCK,winblock);

    // ______ Exit if finished ______
    if (lastblockdone)
      return FILTERED;                  // return

    // ______ Update indexes in matrices, will be corrected for last block ______
    cintpixlo    += numout;             // next block
    cintpixhi    += numout;             // next block
    outblockpixlo = OVERLAP;            // index in block, valid for all middle blocks
    outpixlo      = outpixhi+1;         // index in FILTERED, next range line
    outpixhi      = outpixlo+numout-1;  // index in FILTERED
    } // for all blocks in this buffer
  #endif // check index calling routine
  } // END goldstein phase filter

/****************************************************************
 *    phasefilter modgoldstein                                  *
 * Input is matrix of SIZE (e.g. 32) lines, and N range pixels. *
 * Filtered OUTPUT is same size as input block.                 *
 * Because of overlap, only write to disk in calling routine    *
 * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
 *                                                              *
 * Smoothing of the amplitude of the spectrum is performed by   *
 * spatial convolution with a block kernel of size 2*SMOOTH+1.  *
 * (Which is done by FFT's). e.g. a spatial moving average with *
 * kernel (1d) k=[1 1 1 1 1]/5; kernel2d = transpose(k)*k.      *
 * Blocks in range direction.                                   * 
 *                                                              *
 * After Baran et al. 2003, A Modification to the Goldstein     *
 * Radar Interferogram Filter, IEEE Trans. GRS, v41/9, Sep.2003 *
 #%// Batu 2011 11 03                                           *
 #%// BK 25-Oct-2000                                            *
 ****************************************************************/
matrix<complr4> modgoldstein(
        const matrix<complr4> &CINT,
        const matrix<real4>   &COH,
        const int32            OVERLAP,
        const matrix<real4>   &smoothkernel)    // lying down
  {
  TRACE_FUNCTION("ROUTINE: modgoldstein (BK 03-Nov-2011)")
  //#define CHECKINDICESCALLING
  #ifdef CHECKINDICESCALLING
    return CINT;
  #else
  // ______ Allocate output matrix ______
  const int32 SIZE = CINT.lines();
  const int32 NPIX = CINT.pixels();
  matrix<complr4> FILTERED(SIZE,NPIX);          // output

  // ______ Get block from buffer ______
  const int32 numout  = SIZE-(2*OVERLAP);       // number of output pixels
  int32 cintpixlo     = 0;                      // index in CINT to get 1st block
  int32 cintpixhi     = SIZE-1;                 // index in CINT to get 1st block
  int32 outblockpixlo = 0;                      // index in BLOCK (only 1st block)
  int32 outblockpixhi = SIZE-1-OVERLAP;         // index in BLOCK (except last block)
  int32 outpixlo      = outblockpixlo;          // index in FILTERED (1st block)
  int32 outpixhi      = outblockpixhi;          // index in FILTERED
  bool  lastblockdone = false;                  // only just started...
  // note that int32() floors division 
  const int32 SMOOTH  = int32(smoothkernel.pixels())/2; // half block size, odd kernel
  const bool dosmooth = (SMOOTH==0) ? false : true;
  DEBUG << "SMOOTH: " << SMOOTH;// problem with uint<0 index in smoothkernel
  DEBUG.print();

  // ______ use FFT's for convolution with smoothkernel ______
  // ______ this could also be done static, or in the calling routine ______
  // ______ KERNEL2D is FFT2 of even kernel (no imag part after fft!) ______
  matrix<complr4> KERNEL2D;
  if (dosmooth==true)
    {
    matrix<complr4> kernel(1,SIZE);             // init to zeros
    for (register int32 ii=-SMOOTH; ii<=SMOOTH; ++ii)// 1d kernel function of block
      {
      //kernel(0,(ii+SIZE)%SIZE) = smoothkernel(0,ii-SMOOTH);
      // BK 07-Apr-2003
      // Cygwin fails on index very large number, i.e., uint problem somewhere
      // e.g.: [30,31,0,1,2] <--> [0,1,2,3,4]
      int32 tmp1 = (ii+SIZE)%SIZE;
      int32 tmp2 = ii+SMOOTH;// used to be ii-SMOOTH: wrong
      DEBUG << "tmp1: " << tmp1 << "; tmp2: " << tmp2;
      DEBUG.print();
      kernel(0,tmp1) = complr4(smoothkernel(0,tmp2),real4(0.0));
      }
    KERNEL2D = matTxmat(kernel,kernel);
    fft2d(KERNEL2D);                            // should be real sinc
    }
  DEBUG.print("kernel created for smoothing spectrum");

  // ====== Loop forever, stop after lastblockdone ======
  for (;;)      //forever
    {
    if (cintpixhi>=NPIX-1)                      // check if we are doing the last block
      {
      lastblockdone = true;
      cintpixhi     = NPIX-1;                   // prevent reading after file
      cintpixlo     = cintpixhi-SIZE+1;         // but make sure SIZE pixels are read
      outpixhi      = cintpixhi;                // index in FILTERED 2b written
      outblockpixhi = SIZE-1;                   // write all to the end
      outblockpixlo = outblockpixhi - (outpixhi-outpixlo+1) + 1;
      }
    const window wincint     (0,SIZE-1,cintpixlo,cintpixhi);
    const window winblock    (0,SIZE-1,outblockpixlo,outblockpixhi);
    const window winfiltered (0,SIZE-1,outpixlo,outpixhi);

    // ______ Construct BLOCK as part of CINT ______
    matrix<complr4> BLOCK(wincint,CINT);
    matrix<real4> BLKCH(wincint,COH);

    //#define CHECKINDEXONLY
    #ifndef CHECKINDEXONLY
    // ______ Get spectrum/amplitude/smooth/filter ______
    fft2d(BLOCK);
    matrix<real4> AMPLITUDE = magnitude(BLOCK);

    // ______ use FFT's for convolution with rect ______
    if (dosmooth==true)
      AMPLITUDE = smooth(AMPLITUDE,KERNEL2D);

    //dumpasc("A",AMPLITUDE);
    //AMPLITUDE = smooth(AMPLITUDE,SMOOTH);
    //dumpasc("As",AMPLITUDE); exit(1);
    const real4 maxamplitude = max(AMPLITUDE);
    if (maxamplitude>1e-20) //?
      {
      AMPLITUDE /= maxamplitude;
      real4 ALPHA=1-mean(BLKCH); // calculate alpha based on 1-coherence
      AMPLITUDE.mypow(ALPHA);
      #ifdef WIN32
      BLOCK = timesCxR(BLOCK,AMPLITUDE);// weight spectrum
      #else
      BLOCK     *= AMPLITUDE;           // weight spectrum
      #endif
      }
    else
      {
      WARNING.print("no filtering, maxamplitude<1e-20, zeros in this block?");
      }
    ifft2d(BLOCK);
    #endif // check index blocks
    // ______ Set correct part that is filtered in output matrix ______
    FILTERED.setdata(winfiltered,BLOCK,winblock);

    // ______ Exit if finished ______
    if (lastblockdone)
      return FILTERED;                  // return

    // ______ Update indexes in matrices, will be corrected for last block ______
    cintpixlo    += numout;             // next block
    cintpixhi    += numout;             // next block
    outblockpixlo = OVERLAP;            // index in block, valid for all middle blocks
    outpixlo      = outpixhi+1;         // index in FILTERED, next range line
    outpixhi      = outpixlo+numout-1;  // index in FILTERED
    } // for all blocks in this buffer
  #endif // check index calling routine
  } // END goldstein phase filter


/****************************************************************
 * B = smooth(A,KERNEL)                                         *
 * (circular) spatial moving average with a (2N+1,2N+1) block.  *
 * See also matlab script smooth.m for some tests.              *
 * implementation as convolution with FFT's                     *
 * input: KERNEL is the FFT of the kernel (block)               *
 #%// BK 26-Oct-2000                                            *
 ****************************************************************/
matrix<real4> smooth(
        const matrix<real4> &A,
        const matrix<complr4> &KERNEL2D)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("smooth (BK 26-Oct-2000)");
  #endif
  matrix<complr4> DATA = mat2cr4(A);            // or define fft(R4)
  fft2d(DATA);                                  // or define fft(R4)
  // ______ create kernel in calling routine, e.g., like ______
  // ______ Kernel has to be even! ______
  //const int32 L = A.lines();
  //const int32 P = A.pixels();
  //matrix<complr4> kernel(1,L);                        // init to zeros
  //for (register int32 ii=-N; ii<=N; ++ii)     // 1d kernel function of block
  //  kernel(0,(ii+L)%L) = 1./(2*N+1);
  //matrix<complr4> KERNEL2D = matTxmat(kernel,kernel);
  //fft2d(KERNEL2D);                            // should be real sinc
  DATA *= KERNEL2D;                             // no need for conj. with real fft...
  ifft2d(DATA);                                 // convolution, but still complex...
  return real(DATA);                            // you know it is real only...
  } // END smooth



/****************************************************************
 * B = smooth(A,blocksize)                                      *
 * (circular) spatial moving average with a (2N+1,2N+1) block.  *
 * See also matlab script smooth.m for some tests.              *
 #%// BK 26-Oct-2000                                            *
 ****************************************************************/
matrix<real4> smooth(
        const matrix<real4> &A,
        int32 N)
  {
  #ifdef __DEBUGMAT2
  matDEBUG.print("ROUTINE: smooth (BK 26-Oct-2000)");
  #endif
  // ______ Check if smoothing desired ______
  if (N==0) return A;

  // define one of SPACEDOMAIN and SPECTRALDOMAIN to select algorithm
  // #define SPACEDOMAIN
  #define SPECTRALDOMAIN

  // ====== In space domain ======
  #ifdef SPACEDOMAIN
  const int32 L = A.lines();
  const int32 P = A.pixels();
  matrix<real4> SMOOTH(L,P);            // init to zero...
  register real8 sum = 0.;
  register int32 indexii;
  const real8 Nsmooth = (2*N+1)*(2*N+1);
  for (register int32 i=0; i<L; ++i)
    {
    for (register int32 j=0; j<P; ++j)
      {
      // ______ Smooth this pixel ______
      for (register int32 ii=-N; ii<=N; ++ii)
        {
        indexii = (i+ii+L)%L;
        for (register int32 jj=-N; jj<=N; ++jj)
          {
          sum += A(indexii,(j+jj+P)%P);
          }
        }
      SMOOTH(i,j) = real4(sum/Nsmooth);
      sum = 0.;
      }
    }
  return SMOOTH;
  #endif

  // ====== Do the same but faster ======
  // ______ FFT's, but overhead due to conversion r4<->cr4 ______
  #ifdef SPECTRALDOMAIN
  const int32 L = A.lines();
  //const int32 P = A.pixels();
  matrix<complr4> DATA = mat2cr4(A);            // or define fft(R4)
  fft2d(DATA);                                  // or define fft(R4)
  matrix<complr4> kernel(1,L);                  // init to zeros
  for (register int32 ii=-N; ii<=N; ++ii)       // 1d kernel function of block
    kernel(0,(ii+L)%L) = complr4(1.0/(2*N+1),0.0);// BK 07-Apr-2003
  matrix<complr4> KERNEL2D = matTxmat(kernel,kernel);
  fft2d(KERNEL2D);                              // should be real sinc
  DATA *= KERNEL2D;                             // no need for conj. with real fft...
  ifft2d(DATA);                                 // convolution, but still complex...
  return real(DATA);                            // you know it is real only...
  #endif
  } // END smooth



/****************************************************************
 *    spatialphasefilt                                          *
 * For the first block the part [0:OVERLAP-1] is set to 0.      *
 * For the last block the part [NPIX-1-OVERLAP:NPIX-1] is 0.    *
 #%// BK 30-Oct-2000                                            *
 ****************************************************************/
void spatialphasefilt(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput)
  {
  TRACE_FUNCTION("spatialphasefilt (BK 30-Oct-2000)")

  // ====== Handle input ======
  bool doexternalfile = false;
  char infile[EIGHTY];                          // file 2b filtered
  strcpy(infile,interferogram.file);
  int32 numlinesinput  = int32(interferogram.win.lines()/interferogram.multilookL);
  int32 numpixelsinput = int32(interferogram.win.pixels()/interferogram.multilookP);

  if (specified(filtphaseinput.fifiltphase))
    {
    doexternalfile = true;
    numlinesinput  = filtphaseinput.finumlines;
    strcpy(infile,filtphaseinput.fifiltphase);
    ifstream ifile;
    openfstream(ifile,infile);
    bk_assert(ifile,infile,__FILE__,__LINE__);
    ifile.seekg(0,ios::end);                    // pointer to end...
//  const int32 sizefile     = ifile.tellg();   // opened ate
    const streamoff sizefile     = ifile.tellg();   // opened ate, MA 64 bit pointer
    ifile.close();
    numpixelsinput = sizefile/sizeof(complr4)/numlinesinput;// floor
    // if (numlinesinput*numpixelsinput*sizeof(complr4) != uint(sizefile))
    if (streamoff(numlinesinput)*streamoff(numpixelsinput)*streamoff(sizeof(complr4)) != sizefile)
      WARNING.print("Format infile not CR4, or numlines not correct.");
    INFO.print("Using input file for phase filtering, not default.");
    }
  INFO << "phasefilter: inputfile: (#l,#p): " << infile << " ("
       << numlinesinput << "," << numpixelsinput << ")";
  INFO.print();

  // ______ Set variables ______
  // ______ edge of image not filtered ______
  //const int32 SIZE     = 512;                         // the larger the better
  int32 SIZE       = 1024;                              // blocksize power of 2
  if (generalinput.memory<100e6) SIZE/=2;               // check memory size 512
  if (generalinput.memory<20e6)  SIZE/=2;               // check memory size 256
  while (numlinesinput  < SIZE ) SIZE/=2;
  while (numpixelsinput < SIZE ) SIZE/=2;
  int32 KERNELSIZE_LINES;
  int32 KERNELSIZE_PIXELS;                              // ==L for 1D
  int32 OVERLAP_LINES;
  int32 OVERLAP_PIXELS;                                 // ==L for 1D

  // ______ use FFT's for convolution with rect ______
  // ______ this could also be done static, or in the calling routine ______
  // ______ KERNEL2D is FFT2 of even kernel (no imag part after fft!) ______
  // ______ one could also use ascii input file for 2d kernel ______
  matrix<complr4> KERNEL2D;
  if (!specified(filtphaseinput.fikernel2d))
    {
    KERNELSIZE_LINES  = filtphaseinput.kernel.pixels();
    KERNELSIZE_PIXELS = KERNELSIZE_LINES;
    OVERLAP_LINES     = int32(floor(KERNELSIZE_LINES/2.));
    OVERLAP_PIXELS    = OVERLAP_LINES;

    // ______ 1d kernel function ______
    matrix<complr4> kernel(1,SIZE);             // init to zeros
    for (register int32 ii=-OVERLAP_LINES; ii<=OVERLAP_LINES; ++ii)
      kernel(0,(ii+SIZE)%SIZE) = 
                          complr4(filtphaseinput.kernel(0,ii+OVERLAP_LINES),0.0);// BK 07-Apr-2003
      //kernel(0,(ii+SIZE)%SIZE) = filtphaseinput.kernel(0,ii+OVERLAP_LINES);
    KERNEL2D = matTxmat(kernel,kernel);
    }
  else // use filename and blocksize
    {
    KERNEL2D.resize(SIZE,SIZE);
    INFO.print("Reading 2d kernel from file (PF_IN_KERNEL2D card).");
    ifstream kernel2d(filtphaseinput.fikernel2d, ios::in);
    bk_assert(kernel2d,filtphaseinput.fikernel2d,__FILE__,__LINE__);
    real4 scalefactor;
    kernel2d >> KERNELSIZE_LINES >> KERNELSIZE_PIXELS >> scalefactor;
    if (abs(scalefactor)<=EPS) scalefactor=1.;
    INFO << "kernel2d: file: lines: rows: scale: "
         << filtphaseinput.fikernel2d << " " << KERNELSIZE_LINES << " "
         << KERNELSIZE_PIXELS << " " << scalefactor;
    INFO.print();
    if (!isodd(KERNELSIZE_LINES))
      {
      PRINT_ERROR("2D kernel must have odd number of rows!")
      throw(argument_error);
      }
    if (!isodd(KERNELSIZE_PIXELS))
      {
      PRINT_ERROR("2D kernel must have odd number of columns!")
      throw(argument_error);
      }
    if (KERNELSIZE_LINES >SIZE)
      {
      PRINT_ERROR("2D kernel has more rows than blocksize")
      throw(argument_error);
      }
    if (KERNELSIZE_PIXELS>SIZE)
      {
      PRINT_ERROR("2D kernel has more columns than blocksize")
      throw(argument_error);
      }
    // ______ DO NOT normalize to 1! ______
    OVERLAP_LINES  = int32((KERNELSIZE_LINES/2));// i.e., int32() floors
    OVERLAP_PIXELS = int32((KERNELSIZE_PIXELS/2));// i.e., int32() floors
    // ______ Shift kernel so centered in space domain around pixels ______
    for (int32 ii=-OVERLAP_LINES; ii<=OVERLAP_LINES; ++ii)
      {
      real4 tmpvalue;
      char dummyline[10*ONE27];                 // take care of very large kernel
      kernel2d.getline(dummyline,10*ONE27,'\n');
      const int32 indexii = (ii+SIZE)%SIZE;
      for (int32 jj=-OVERLAP_PIXELS; jj<=OVERLAP_PIXELS; ++jj)
        {
        const int32 indexjj = (jj+SIZE)%SIZE;
        kernel2d >> tmpvalue;
        KERNEL2D(indexii,indexjj) = complr4(scalefactor*tmpvalue);
        }
      }
    } // 2Dkernel

  // ______ Give some info ______
  INFO << "buffersize: (" << SIZE << ", " << SIZE << ")";
  INFO.print();
  INFO << "kernelsize: ("
       << KERNELSIZE_LINES << ", " << KERNELSIZE_PIXELS << ")";
  INFO.print();

  // ______ Prepare kernel for routine convolution ______
  fft2d(KERNEL2D);      // input to myconv function
  KERNEL2D.conj();      // input to myconv function, not required here, but correct.
                        // since fft of even function is real only...

  // ______ Open output file ______
  ofstream ofile;
  openfstream(ofile,filtphaseinput.fofiltphase,generalinput.overwrit);
  bk_assert(ofile,filtphaseinput.fofiltphase,__FILE__,__LINE__);


  // ______ initialize indices ______
  const int32 numout   = SIZE-(2*OVERLAP_LINES);// number of output lines per buffer
  int32 cintlinelo     = 0;                     // index in CINT to get 1st buffer
  int32 cintlinehi     = SIZE-1;                // index in CINT to get 1st buffer
  int32 filteredlinelo = 0;                     // index in FILTERED to write (1st)
  int32 filteredlinehi = SIZE-OVERLAP_LINES-1;  // index in FILTERED to write
  bool  lastbuffer     = false;                 // not yet, just starting...

  if (!doexternalfile)
    if (interferogram.formatflag!=FORMATCR4)
      {
      PRINT_ERROR("Sorry, for phasefilter, format of interferogram must be CR4")
      throw(argument_error);
      }
  matrix<complr4> CINT(SIZE,numpixelsinput);    // type must be cr4!

  int32 lineswritten     = OVERLAP_LINES;               // first buffer extra OVERLAP lines
  const int32 numbuffers = numlinesinput/numout;        // approx?
  int32 tenpercent       = int32((numbuffers/10)+.5); // round
  int32 percent          = 10;
  if (tenpercent==0) tenpercent = 1000;                 // avoid error x%0;
  PROGRESS.print("filtphase:  0%");

  // ====== Start filter loop buffers of BLOCKSIZE ======
  for (int32 buffer=1; buffer<1e6; ++buffer)
    {
    // ______ Check if this will be the last buffer ______
    if (cintlinehi >= numlinesinput-1)          // -1 since in matrix system
      {
      lastbuffer     = true;
      cintlinehi     = numlinesinput-1;
      cintlinelo     = cintlinehi-SIZE+1;       // make sure SIZE lines are read
      filteredlinehi = SIZE-1;                  // write 0's OVERLAP lines
      const int32 lines2bwritten = numlinesinput-lineswritten;
      filteredlinelo = filteredlinehi - lines2bwritten + 1;
      if (lines2bwritten<1) 
        WARNING.print("panic, this will crash, lines2bwritten<1?");
      }

    // ______ Read in buffers of BLOCKSIZE lines complex interferogram ______
    const window windummy (0,0,0,0);            // no offset in readfile
    const window wincint  (cintlinelo,cintlinehi,0,numpixelsinput-1);
    readfile(CINT,infile,numlinesinput,wincint,windummy);

    // ______ Filter data in buffer of BLOCKSIZE lines ______
    // ______ output has same size, but write only part to disk ______
    matrix<complr4> FILTERED = convbuffer(CINT,KERNEL2D,OVERLAP_PIXELS);

    // ______ Write filtered data (i.e. part of FILTERED) ______
    // ______ Write overlap zero lines for first buffer and lastbuffer (circular) ______
    if (filteredlinelo==0)                      // firstblock
      for (int ii=0; ii<OVERLAP_LINES; ++ii)
        FILTERED.setrow(ii,complr4(0.,0.));
    if (lastbuffer==true)
      for (int ii=SIZE-OVERLAP_LINES; ii<SIZE; ++ii)
        FILTERED.setrow(ii,complr4(0.,0.));
    const window winfiltered(filteredlinelo,filteredlinehi,0,numpixelsinput-1);
    writefile(ofile,FILTERED,winfiltered);

    // ______ Check if all done ______
    if (lastbuffer)
      break;

    // ______ Update indexes in matrices for all buffers except last ______
    lineswritten  += numout;            // first block overlap added
    cintlinelo    += numout;            // next buffer
    cintlinehi    += numout;            // next buffer
    filteredlinelo = OVERLAP_LINES;     // index in buffer for all restbuffers

    // ______ Give progress message ______
    if (buffer%tenpercent==0)
      {
      PROGRESS << "FILTPHASE: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    } // end loop buffers
  ofile.close();

  // ______ exit if only external file was desired ______
  if (doexternalfile)
    {
    cerr  << "\n Finished external file phasefilter, Exiting\n";
    PROGRESS.print("\n Finished external file phasefilter, Exiting\n");
    exit(0);
    }


  // ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltphase", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtphase: scratchlogfiltphase",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* PHASE FILTER COMPLEX INTERFEROGRAM: \t"
    << "\n*******************************************************************"
    <<  infile
    << "\nOutput file filtered master (format): \t\t\t"
    <<  filtphaseinput.fofiltphase << " "
    << "\n..."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  // ______ updateproductinfo greps info from this file ______
  ofstream scratchresfile("scratchresfiltphase", ios::out | ios::trunc);
  bk_assert(scratchresfile,"scratchresfiltphase",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_filtphase]
    << "\n*******************************************************************"
    << "\nMethod_phasefilt: spatial convolution: \t"
//    <<  kernel...
    << "\nInput_file:                           \t"
    <<  infile
    << "\nData_output_file:                     \t"
    <<  filtphaseinput.fofiltphase
    << "\nData_output_format:                   \t"
    <<  "complex_real4"
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
    <<  numlinesinput
    << "\nNumber of pixels (multilooked):       \t"
    <<  numpixelsinput
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_filtphase] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  // ====== Tidy up ======
  } // END spatialphasefilt



/****************************************************************
 *    phasefilter buffer by spatial conv. with kernel.          *
 * Input is matrix of SIZE (e.g. 256) lines, and N range pixels.*
 * Filtered OUTPUT is same size as input block.                 *
 * Because of overlap, only write to disk in calling routine    *
 * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
 * (in line direction)                                          *
 * spatial convolution with a kernel function, such as a block  *
 * function 111 (1D) (By FFT's).                                *
 * Processing is done in blocks in range direction.             * 
 * For the first block the part [0:OVERLAP-1] is set to 0.      *
 * For the last block the part [NPIX-1-OVERLAP:NPIX-1] is 0.    *
 *                                                              *
 * Input:                                                       *
 *  - matrix to be filtered of blocklines * numpixs             *
 *  - kernel2d: fft2 of 2d spatial kernel.                      *
 *  - overlap: half of the kernel size, e.g., 1 for 111.        *
 * Output:                                                      *
 *  - filtered matrix.                                          *
 *    ifft2d(BLOCK .* KERNEL2D) is returned, so if required for *
 *    non symmetrical kernel, offer the conj(KERNEL2D)!         *
 *                                                              *
 #%// BK 30-Oct-2000                                            *
 ****************************************************************/
matrix<complr4> convbuffer(
        const matrix<complr4> &CINT,
        const matrix<complr4> &KERNEL2D,
        const int32            OVERLAP)         // overlap in column direction
  {
  TRACE_FUNCTION("convbuffer (BK 30-Oct-2000)")

  //#define CHECKINDICESCALLING
  #ifdef CHECKINDICESCALLING
    return CINT;
  #else
  // ______ Allocate output matrix ______
  const int32 SIZE = CINT.lines();
  const int32 NPIX = CINT.pixels();
  matrix<complr4> FILTERED(SIZE,NPIX);          // allocate output (==0)

  // ______ Get block from buffer ______
  const int32 numout  = SIZE-(2*OVERLAP);       // number of output pixels per block
  int32 cintpixlo     = 0;                      // index in CINT to get 1st block
  int32 cintpixhi     = SIZE-1;                 // index in CINT to get 1st block
  //int32 outblockpixlo = 0;                    // index in BLOCK (only 1st block)
  int32 outblockpixlo = OVERLAP;                // index in block
  int32 outblockpixhi = SIZE-1-OVERLAP;         // index in BLOCK (except last block)
  int32 outpixlo      = outblockpixlo;          // index in FILTERED (1st block)
  int32 outpixhi      = outblockpixhi;          // index in FILTERED
  bool  lastblockdone = false;                  // only just started...


  // ====== Loop forever, stop after lastblockdone ======
  for (;;)      //forever
    {
    if (cintpixhi>=NPIX-1)                      // check if we are doing the last block
      {
      lastblockdone = true;
      cintpixhi     = NPIX-1;                   // prevent reading after file
      cintpixlo     = cintpixhi-SIZE+1;         // but make sure SIZE pixels are read
      // leave last few==0
      outpixhi      = NPIX-1-OVERLAP;           // index in FILTERED 2b written
      //outblockpixhi = SIZE-1;                 // write all to the end
      outblockpixlo = outblockpixhi - (outpixhi-outpixlo+1) + 1;
      }
    const window wincint     (0,SIZE-1,cintpixlo,cintpixhi);
    const window winblock    (0,SIZE-1,outblockpixlo,outblockpixhi);
    const window winfiltered (0,SIZE-1,outpixlo,outpixhi);

    // ______ Construct BLOCK as part of CINT ______
    matrix<complr4> BLOCK(wincint,CINT);

    //#define CHECKINDEXONLY
    #ifndef CHECKINDEXONLY
    // ______ use FFT's for convolution with kernel ______
    fft2d(BLOCK);
    BLOCK *= KERNEL2D;                  // no need for conj. cause kernel is even?
    ifft2d(BLOCK);                      // conv. in space domain
    #endif // check index blocks
    // ______ Set correct part that is filtered in output matrix ______
    FILTERED.setdata(winfiltered,BLOCK,winblock);

    // ______ Exit if finished ______
    if (lastblockdone)
      return FILTERED;                  // return

    // ______ Update indexes in matrices, will be corrected for last block ______
    cintpixlo    += numout;             // next block
    cintpixhi    += numout;             // next block
    outpixlo      = outpixhi+1;         // index in FILTERED, next range line
    outpixhi      = outpixlo+numout-1;  // index in FILTERED
    } // for all blocks in this buffer
  #endif // check index calling routine
  } // END convbuffer



/****************************************************************
 *    phasefilterspectral                                       *
 * loop over whole file and multiply spectrum of interferogram  *
 * with kernel specified in input file.                         *
 #%// BK 31-Oct-2000                                            *
 ****************************************************************/
void phasefilterspectral(
        const input_gen       &generalinput,
        const productinfo     &interferogram,
        const input_filtphase &filtphaseinput)
  {
  TRACE_FUNCTION("phasefilterspectral (BK 31-Oct-2000)")

  // ====== Handle input ======
  bool doexternalfile = false;
  if (specified(filtphaseinput.fifiltphase))
    doexternalfile = true;
  char infile[EIGHTY];                          // file 2b filtered
  strcpy(infile,interferogram.file);
  int32 numlinesinput  = int32(interferogram.win.lines()/interferogram.multilookL);
  int32 numpixelsinput = int32(interferogram.win.pixels()/interferogram.multilookP);
  if (doexternalfile)                           // overwrite default file
    {
    numlinesinput  = filtphaseinput.finumlines;
    strcpy(infile,filtphaseinput.fifiltphase);
    ifstream ifile;
    openfstream(ifile,infile);
    bk_assert(ifile,infile,__FILE__,__LINE__);
    ifile.seekg(0,ios::end);                    // pointer to end...
//  const int32 sizefile     = ifile.tellg();   // opened ate
    const streamoff sizefile     = ifile.tellg();   // opened ate, MA 64 bit pointer
    ifile.close();
    numpixelsinput = sizefile/sizeof(complr4)/numlinesinput;// floor
    //if (numlinesinput*numpixelsinput*sizeof(complr4) != uint(sizefile))
    if (streamoff(numlinesinput)*streamoff(numpixelsinput)*streamoff(sizeof(complr4)) != sizefile)
      WARNING.print("Format infile not CR4, or numlines not correct.");
    INFO.print("Using input file for phase filtering, not default.");
    }
  INFO << "phasefilter: inputfile: (#l,#p): " << infile << " ("
       << numlinesinput << "," << numpixelsinput << ")";
  INFO.print();

  // ______ Set variables ______
  const int32 SIZE     = filtphaseinput.blocksize;      // power of 2
  const int32 OVERLAP  = filtphaseinput.overlap;        // half of the overlap

  // ______ Open output file ______
  ofstream ofile;
  openfstream(ofile,filtphaseinput.fofiltphase,generalinput.overwrit);
  bk_assert(ofile,filtphaseinput.fofiltphase,__FILE__,__LINE__);


  if (!doexternalfile)
    if (interferogram.formatflag!=FORMATCR4)
      {
      PRINT_ERROR("Sorry, for phasefilter, format of interferogram must be CR4")
      throw(argument_error);
      }
  matrix<complr4> CINT(SIZE,numpixelsinput);    // type must be cr4!

  // ====== Obtain KERNEL2D to multiply spectrum by from input file ======
  INFO.print("Reading 2d kernel from file (PF_IN_KERNEL2D card).");
  ifstream kernel2d(filtphaseinput.fikernel2d, ios::in);
  bk_assert(kernel2d,filtphaseinput.fikernel2d,__FILE__,__LINE__);
  // ______ First read header of kernel2d file ______
  int32 KERNELSIZE_LINES  = 0;
  int32 KERNELSIZE_PIXELS = 0;
  real4 scalefactor       = 0.;
  kernel2d >> KERNELSIZE_LINES >> KERNELSIZE_PIXELS >> scalefactor;
  if (abs(scalefactor)<=EPS) scalefactor=1.;
  INFO << "kernel2d: file: lines: rows: scale: "
       << filtphaseinput.fikernel2d << " " << KERNELSIZE_LINES << " "
       << KERNELSIZE_PIXELS << " " << scalefactor;
  INFO.print();
  // why should it be odd for spectral filter...
  //if (!isodd(KERNELSIZE_LINES))  PRINT_ERROR("2D kernel must have odd number of rows!")
  //if (!isodd(KERNELSIZE_PIXELS)) PRINT_ERROR("2D kernel must have odd number of columns!");
  if (KERNELSIZE_LINES >SIZE)
    {
    PRINT_ERROR("2D kernel has more rows than blocksize");
    throw(argument_error);
    }
  if (KERNELSIZE_PIXELS>SIZE) 
    {
    PRINT_ERROR("2D kernel has more columns than blocksize");
    throw(argument_error);
    }

  // ______ Then read values of kernel2d ______
  // ______ Shift kernel so centered in spectral domain around zero freq. ______
  // ______ DO NOT normalize to 1 ______
  //const int32 HBS_L = int32(floor(KERNELSIZE_LINES/2));
  //const int32 HBS_P = int32(floor(KERNELSIZE_PIXELS/2));
  const int32 HBS_L = int32((KERNELSIZE_LINES/2));
  const int32 HBS_P = int32((KERNELSIZE_PIXELS/2));
  const int32 EXTRA_L = (iseven(KERNELSIZE_LINES))  ? 1:0;      // 1 less to fill
  const int32 EXTRA_P = (iseven(KERNELSIZE_PIXELS)) ? 1:0;      // 1 less to fill
  matrix<real4> KERNEL2D(SIZE,SIZE);                    // allocate THE matrix
  for (int32 ii=-HBS_L+EXTRA_L; ii<=HBS_L; ++ii)
    {
    char dummyline[10*ONE27];                           // prevent very large kernels
    kernel2d.getline(dummyline,10*ONE27,'\n');
    const int32 indexii = (ii+SIZE)%SIZE;
    for (int32 jj=-HBS_P+EXTRA_P; jj<=HBS_P; ++jj)
      {
      const int32 indexjj = (jj+SIZE)%SIZE;
      kernel2d >> KERNEL2D(indexii,indexjj);
      }
    }
  if (scalefactor!=1) KERNEL2D *= scalefactor;


  // ______ Initialize indices ______
  const int32 numout   = SIZE-(2*OVERLAP);      // number of output lines
  int32 cintlinelo     = 0;                     // index in CINT to get 1st buffer
  int32 cintlinehi     = SIZE-1;                // index in CINT to get 1st buffer
  int32 filteredlinelo = 0;                     // index in FILTERED to write (1st)
  int32 filteredlinehi = SIZE-OVERLAP-1;        // index in FILTERED to write (1st)
  bool  lastbuffer     = false;                 // not yet, just starting...
  int32 lineswritten   = OVERLAP;               // first buffer extra OVERLAP lines

  // ______ Progress messages ______
  int32 percent          = 10;
  //const int32 numbuffers = numlinesinput/numout;      // approx?
  const int32 numbuffers = numlinesinput/SIZE;  // approx for large mem?
  int32 tenpercent = int32((numbuffers/10)+.5); // round
  if (tenpercent==0) tenpercent = 1000;                 // avoid error: x%0
  PROGRESS.print("filtphase:  0%");

  // ====== Start filter loop buffers of BLOCKSIZE ======
  for (int32 buffer=1; buffer<1e6; ++buffer)
    {
    // ______ Check if this will be the last buffer ______
    if (cintlinehi >= numlinesinput-1)          // -1 since in matrix system
      {
      lastbuffer     = true;
      cintlinehi     = numlinesinput-1;
      cintlinelo     = cintlinehi-SIZE+1;       // make sure SIZE lines are read
      filteredlinehi = SIZE-1;                  // write upto lastline for last buffer
      const int32 lines2bwritten = numlinesinput-lineswritten;
      filteredlinelo = filteredlinehi - lines2bwritten + 1;
      if (lines2bwritten<1) 
        WARNING.print("panic, this will crash, lines2bwritten<1?");
      }

    // ______ Read in buffers of BLOCKSIZE lines complex interferogram ______
    const window windummy (0,0,0,0);            // no offset in readfile
    const window wincint  (cintlinelo,cintlinehi,0,numpixelsinput-1);
    readfile(CINT,infile,numlinesinput,wincint,windummy);

    // ______ Filter data in buffer of BLOCKSIZE lines ______
    // ______ output has same size, but write only part to disk ______
    const matrix<complr4> FILTERED = spectralfilt(CINT,KERNEL2D,OVERLAP);

    // ______ Write filtered data ______
    const window winfiltered (filteredlinelo,filteredlinehi,0,numpixelsinput-1);
    writefile(ofile,FILTERED,winfiltered);

    // ______ Check if all done ______
    if (lastbuffer)
      break;

    // ______ Update indexes in matrices, will be corrected for last block ______
    lineswritten  += numout;            // first block overlap added
    cintlinelo    += numout;            // next buffer
    cintlinehi    += numout;            // next buffer
    filteredlinelo = OVERLAP;           // index in buffer,
                                        // valid for all middle buffers, but not last
    // ______ Give progress message ______
    if (buffer%tenpercent==0)
      {
      PROGRESS << "FILTPHASE: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    } // end loop buffers
  ofile.close();

  // ______ exit if only external file was desired ______
  if (doexternalfile)
    {
    cerr << "\n  Finished external file phasefilter, Exiting\n";
    PROGRESS.print("\n Finished external file phasefilter, Exiting\n");
    exit(0);
    }

  // ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltphase", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtphase: scratchlogfiltphase",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* PHASE FILTER COMPLEX INTERFEROGRAM: \t"
    << "\n*******************************************************************"
    <<  infile
    << "\nOutput file filtered master (format): \t\t\t"
    <<  filtphaseinput.fofiltphase << " "
// TODO
    << "\n..."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  // ______ updateproductinfo greps info from this file ______
  ofstream scratchresfile("scratchresfiltphase", ios::out | ios::trunc);
  bk_assert(scratchresfile,"scratchresfiltphase",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_filtphase]
    << "\n*******************************************************************"
    << "\nMethod_phasefilt: spectral: size, overlap, kernelfile: \t"
// TODO
    <<  SIZE << " " << OVERLAP << " " << filtphaseinput.fikernel2d
    << "\nInput_file:                           \t"
    <<  infile
    << "\nData_output_file:                     \t"
    <<  filtphaseinput.fofiltphase
    << "\nData_output_format:                   \t"
    <<  "complex_real4"
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
    <<  numlinesinput
    << "\nNumber of pixels (multilooked):       \t"
    <<  numpixelsinput
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_filtphase] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  // ====== Tidy up ======
  } // END phasefilterspectral



/****************************************************************
 *    phasefilter spectral                                      *
 * Input is matrix of SIZE (e.g. 32) lines, and N range pixels. *
 * Filtered OUTPUT is same size as input block.                 *
 * Because of overlap, only write to disk in calling routine    *
 * part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]             *
 *                                                              *
 * Filtering is performed by pointwise multiplication of the    *
 * spectrum per block by the KERNEL2D (input).                  *
 * Blocks in range direction,                                   * 
 *                                                              *
 #%// BK 31-Oct-2000                                            *
 ****************************************************************/
matrix<complr4> spectralfilt(
        const matrix<complr4> &CINT,
        const matrix<real4>   &KERNEL2D,
        const int32            OVERLAP)
  {
  TRACE_FUNCTION("spectralfilt (BK 31-Oct-2000)")
  // ______ Allocate output matrix ______
  const int32 SIZE = CINT.lines();
  const int32 NPIX = CINT.pixels();
  matrix<complr4> FILTERED(SIZE,NPIX);          // output

  // ______ Get block from buffer ______
  const int32 numout  = SIZE-(2*OVERLAP);       // number of output pixels
  int32 cintpixlo     = 0;                      // index in CINT to get 1st block
  int32 cintpixhi     = SIZE-1;                 // index in CINT to get 1st block
  int32 outblockpixlo = 0;                      // index in BLOCK (only 1st block)
  int32 outblockpixhi = SIZE-1-OVERLAP;         // index in BLOCK (except last block)
  int32 outpixlo      = outblockpixlo;          // index in FILTERED (1st block)
  int32 outpixhi      = outblockpixhi;          // index in FILTERED
  bool  lastblockdone = false;                  // only just started...

  // ====== Loop forever, stop after lastblockdone ======
  for (;;)      //forever
    {
    if (cintpixhi>=NPIX-1)                      // check if we are doing the last block
      {
      lastblockdone = true;
      cintpixhi     = NPIX-1;                   // prevent reading after file
      cintpixlo     = cintpixhi-SIZE+1;         // but make sure SIZE pixels are read
      outpixhi      = cintpixhi;                // index in FILTERED 2b written
      outblockpixhi = SIZE-1;                   // write all to the end
      outblockpixlo = outblockpixhi - (outpixhi-outpixlo+1) + 1;
      }
    const window wincint     (0,SIZE-1,cintpixlo,cintpixhi);
    const window winblock    (0,SIZE-1,outblockpixlo,outblockpixhi);
    const window winfiltered (0,SIZE-1,outpixlo,outpixhi);

    // ______ Construct BLOCK as part of CINT ______
    matrix<complr4> BLOCK(wincint,CINT);

    // ______ Get spectrum/filter/ifft ______
    fft2d(BLOCK);
    #ifdef WIN32
    matrix<real4> temp = KERNEL2D;// Jia
    BLOCK = timesCxR(BLOCK,temp);// Jia
    #else
    BLOCK *= KERNEL2D;                  // the filter...
    #endif
    ifft2d(BLOCK);
    // ______ Set correct part that is filtered in output matrix ______
    FILTERED.setdata(winfiltered,BLOCK,winblock);

    // ______ Exit if finished ______
    if (lastblockdone)
      return FILTERED;                  // return

    // ______ Update indexes in matrices, will be corrected for last block ______
    cintpixlo    += numout;             // next block
    cintpixhi    += numout;             // next block
    outblockpixlo = OVERLAP;            // index in block, valid for all middle blocks
    outpixlo      = outpixhi+1;         // index in FILTERED, next range line
    outpixhi      = outpixlo+numout-1;  // index in FILTERED
    } // for all blocks in this buffer
  } // END spectralfilt



/****************************************************************
 *    azimuthfilter                                             *
 * Loop over whole master and slave image and filter out        *
 * part of the spectrum that is not common.                     *
 * Only do zero doppler freq. offset.                           *
 * do not use a polynomial from header for now.                 *
 * (next we will, but assume image are almost coreg. in range,  *
 *  so f_dc polynomial can be eval. same)                       * 
 * Per block in azimuth [1024] use a certain overlap with the   *
 * next block so that same data is partially used for spectrum  *
 * (not sure if this is requried).                              *
 * Filter is composed of: DE-hamming, RE-hamming (for correct   *
 * new size and center of the spectrum).                        *
 * Trick in processor.c: First call routine as:                 *
 *  (generalinput,filtaziinput,master,slave)                    *
 * in order to process the master, and then as:                 *
 *  (generalinput,filtaziinput,slave,master)                    *
 * to filter the slave slc image.                               *
 #%// BK 01-Nov-2000                                            *
 ****************************************************************/
void azimuthfilter(
        const input_gen       &generalinput,
        const input_filtazi   &filtaziinput,
        slcimage        &master,  // not const, fdc possibly reset here
        slcimage        &slave)   // not const, fdc possibly reset here
  {
  TRACE_FUNCTION("azimuthfilter (BK 01-Nov-2000)")
  // ====== Handle input ======
  char infile[EIGHTY];                          // file 2b filtered
  strcpy(infile,master.file);
  int32 numlinesinput  = master.currentwindow.lines();
  int32 numpixelsinput = master.currentwindow.pixels();
  INFO << "azimuthfilter: inputfile: (#l,#p): " << infile << " ("
       << numlinesinput << "," << numpixelsinput << ")";
  INFO.print();

  // ______ Set variables ______
  const int32 SIZE     = filtaziinput.fftlength;// power of 2
  const int32 OVERLAP  = filtaziinput.overlap;  // half of the overlap
  const real8 HAMMING  = filtaziinput.hammingalpha;
  matrix<real4> FILTER;                         // i.e., THE filter


  // ====== If master and slave both have almost constant Doppler, ======
  // ====== then reset polynomial to easy case ======
  // --- Check variation of Doppler for this crop; zeros in spectrum at correct place ---
  const real8 max_fdc_change  = 0.30*abs(master.prf - master.abw);// 100 Hz or so for ERS
  const real8 master_fdc_p0   = master.pix2fdc(master.currentwindow.pixlo);
  const real8 master_fdc_p05  = master.pix2fdc((master.currentwindow.pixhi-master.currentwindow.pixlo)/2);
  const real8 master_fdc_pN   = master.pix2fdc(master.currentwindow.pixhi);
  const real8 slave_fdc_p0    = slave.pix2fdc(slave.currentwindow.pixlo);
  const real8 slave_fdc_p05   = slave.pix2fdc((slave.currentwindow.pixhi-slave.currentwindow.pixlo)/2);
  const real8 slave_fdc_pN    = slave.pix2fdc(slave.currentwindow.pixhi);
  const real8 master_max_dfdc = max(abs(master_fdc_p0-master_fdc_p05),abs(master_fdc_p0-master_fdc_pN));
  const real8 slave_max_dfdc  = max(abs(slave_fdc_p0-slave_fdc_p05),abs(slave_fdc_p0-slave_fdc_pN));
  INFO << "Master: fDC Variation [Hz] = " << master_max_dfdc;
  INFO.print();
  INFO << "Slave: fDC Variation [Hz]  = " << slave_max_dfdc;
  INFO.print();
  if (master_max_dfdc<max_fdc_change &&
      slave_max_dfdc <max_fdc_change) 
    {
    INFO << "Variation of fDC for master and slave < " << max_fdc_change;
    INFO.print();
    INFO.print("Reseting Doppler polynomial to constant of center crop:");
    master.f_DC_a0 = master_fdc_p05;
    master.f_DC_a1 = 0.0;
    master.f_DC_a2 = 0.0;
    INFO << "Master: constant fDC [Hz] = " << master.f_DC_a0;
    INFO.print();
    slave.f_DC_a0  = slave_fdc_p05;
    slave.f_DC_a1  = 0.0;
    slave.f_DC_a2  = 0.0;
    INFO << "Slave: constant fDC [Hz] = " << slave.f_DC_a0;
    INFO.print();
    }

  // ====== Find out if same filter for all columns can be used ======
  bool samefd0 = false;                         // assume NOT same filter
  if (abs(master.f_DC_a1)<EPS &&
      abs(master.f_DC_a2)<EPS &&
      abs( slave.f_DC_a1)<EPS &&
      abs( slave.f_DC_a2)<EPS)          // true for eSARp'ed processed raw data.
    samefd0 = true;
  else
    DEBUG.print("Filtering data by computing fDC for each column."); 

  // ______ Compute filter here if same for all columns ______
  // ______ Else call routine filtazibuffer to compute filt per column ______
  // ______ based on polynomial for doppler centroid. Then assume ______
  // ______ images are almost alligned for 1st column in original system ______
  if (samefd0==true)
    {
    DEBUG.print("Filtering data by same fDC for each column."); 
    const bool dohamming = (HAMMING<0.9999) ? true : false;
    const real8 PRF      = master.prf;          // pulse repetition freq. [Hz]
    const real8 ABW      = master.abw;          // azimuth band width [Hz]
    const real8 fDC_m    = master.f_DC_a0;      // zero doppler freq. [Hz]
    const real8 fDC_s    = slave.f_DC_a0;       // zero doppler freq. [Hz]
    const real8 fDC_mean = 0.5*(fDC_m+fDC_s);   // mean doppler centroid freq.
    const real8 ABW_new  = max(1.0, 
      2.0*(0.5*ABW-abs(fDC_m-fDC_mean)));       // new bandwidth>1.0
    INFO << "New Azimuth Bandwidth: " << ABW_new << " [Hz]";
    INFO.print();
    INFO << "New central frequency: " << fDC_mean << " [Hz]";
    INFO.print();

    const real4 deltaf   =  PRF/real4(SIZE);
    const real4 fr       = -PRF/2.0;
    matrix<real4> freqaxis(1,SIZE);
    for (int32 i=0; i<SIZE; ++i)
      freqaxis(0,i) = fr+real4(i)*deltaf;       // [-fr:df:fr-df]
    if (dohamming)
      {
      // ______ NOT a good implementation for per col., cause wshift *AND* fftshift.
      // ______ DE-weight spectrum at centered at fDC_m ______
      // ______ spectrum should be periodic! (use wshift) ______
      matrix<real4> inversehamming = myhamming(freqaxis,ABW,PRF,HAMMING);
      for (int32 i=0; i<SIZE; ++i)
        if (inversehamming(0,i)!=0.0)// check no zero division
          inversehamming(0,i) = 1.0/inversehamming(0,i);
      // ______ Shift this circular by myshift pixels ______
      int32 myshift = int32(((SIZE*fDC_m/PRF) + 0.5));  // round
      wshift(inversehamming,-myshift);          // center at fDC_m

      // ______ Newhamming is scaled and centered around new mean ______
      myshift = int32(((SIZE*fDC_mean/PRF) + 0.5));             // round
      FILTER  = myhamming(freqaxis,
                          ABW_new,
                          PRF,HAMMING);         // fftshifted
      wshift(FILTER,-myshift);                  // center at fDC_mean
      //#define REALLYDEBUG
      #ifdef REALLYDEBUG
      WARNING.print("really debug switched on, dumping filters to ascii file.");
      dumpasc("INVHAM",inversehamming);
      dumpasc("NEWHAM",FILTER);
      #endif
      FILTER *= inversehamming; // filterslve=1/filter, or flipped in meanfdc...
      }
    else        // no weighting, but center at fDC_mean, size ABW_new
      {
      int32 myshift = int32(((SIZE*fDC_mean/PRF) + 0.5));       // round
      FILTER = myrect(freqaxis/real4(ABW_new)); // fftshifted
      wshift(FILTER,-myshift);                  // center at fDC_mean
      #ifdef REALLYDEBUG
      dumpasc("NEWHAM",FILTER);
      #endif
      }
    ifftshift(FILTER);                          // fftsh works on data!
    #ifdef REALLYDEBUG
    dumpasc("TOTshifted",FILTER);
    #endif
    }


  // ______ Open output file ______
  ofstream ofile;
  openfstream(ofile,filtaziinput.foname,generalinput.overwrit);
  bk_assert(ofile,filtaziinput.foname,__FILE__,__LINE__);

  // ______ Initialize indices ______
  const int32 numout   = SIZE-(2*OVERLAP);      // number of output lines
  int32 slclinelo      = 0;                     // index in SLCIMAGE for 1st buffer
  int32 slclinehi      = SIZE-1;                // index in SLCIMAGE for 1st buffer
  int32 filteredlinelo = 0;                     // index in FILTERED to write (1st)
  int32 filteredlinehi = SIZE-OVERLAP-1;        // index in FILTERED to write (1st)
  bool  lastbuffer     = false;                 // not yet, just starting...
  int32 lineswritten   = OVERLAP;               // first buffer extra OVERLAP lines

  // ______ Progress messages (bug report Mohanty Kamini Kanta 12/06/03) ______
  int32 percent          = 10;
  const int32 numbuffers = (numlinesinput+OVERLAP)/numout;      // approx?
  int32 tenpercent       = int32(rint(numbuffers/10.0));// round
  if (tenpercent==0) tenpercent = 1000;                 // avoid error: x%0
  PROGRESS.print("FILTAZI:  0%");

  // ====== Start filter loop buffers of BLOCKSIZE ======
  matrix<complr4> SLCIMAGE(SIZE,numpixelsinput);
  for (int32 buffer=1; buffer<1e6; ++buffer)
    {
    // ______ Check if this will be the last buffer ______
    if (slclinehi >= numlinesinput-1)           // -1 since in matrix system
      {
      lastbuffer     = true;
      slclinehi      = numlinesinput-1;
      slclinelo      = slclinehi-SIZE+1;        // make sure SIZE lines are read
      filteredlinehi = SIZE-1;                  // write upto lastline for last buffer
      // if there is 1 buffer but smaller than size=fftlength, it cannot work
      //if (slclinelo < 0) slclinelo = 0;     //KKM (not OK, BK: have to take fft)
      if (slclinelo < 0) 
        {
        WARNING.print("I think fftlength is too big (>size image).");       
        WARNING.print("...but i will continue and probably crash.");       
        WARNING.print("PLease run again with smaller fftlength.");       
        }

      const int32 lines2bwritten = numlinesinput-lineswritten;
      filteredlinelo = filteredlinehi - lines2bwritten + 1;
      if (lines2bwritten<1) 
        WARNING.print("panic, this will crash, lines2bwritten<1?");
      }

    // ______ Read buffer of BLOCKSIZE lines in SLCIMAGE ______
    const window windummy (0,0,0,0);            // no offset in readfile
    const window winslc   (slclinelo,slclinehi,0,numpixelsinput-1);
    switch (master.formatflag)
      {
      case FORMATCI2:
        {
        fileci2tomatcr4(SLCIMAGE,infile,numlinesinput,winslc,windummy);
        break;
        }
      case FORMATCR4:
        {
        readfile(SLCIMAGE,infile,numlinesinput,winslc,windummy);
        break;
        }
      default:
        {
        PRINT_ERROR("readdata::not correct format on file.")
        throw(unhandled_case_error);
        }
      }

    // ______ Filter data in buffer of BLOCKSIZE lines ______
    // ______ output has same size, but write only part to disk ______
    matrix<complr4> FILTERED;
    if (samefd0==true)
      {
      // ______ Same filter for all columns ______
      fft(SLCIMAGE,1);                          // fft foreach column
      FILTERED = diagxmat(FILTER,SLCIMAGE);     // filter each column
      ifft(FILTERED,1);                         // ifft foreach column
      }
    else
      {
      // ______ Filter depends on columns number ______
      FILTERED = blockazifilt(SLCIMAGE,master,slave,HAMMING);
      }

    // ______ Write filtered data ______
    const window winfiltered(filteredlinelo,filteredlinehi,0,numpixelsinput-1);
    switch (filtaziinput.oformatflag)
      {
      case FORMATCR4:
        {
        writefile(ofile,FILTERED,winfiltered);
        break;
        }
      // ______ Convert first to ci2 before writing to file ______
      case FORMATCI2:
        {
        matrix <compli16> TMP(FILTERED.lines(),FILTERED.pixels());
        for (uint ii=winfiltered.linelo; ii<=winfiltered.linehi; ++ii)
          for (uint jj=winfiltered.pixlo; jj<=winfiltered.pixhi; ++jj)
            TMP(ii,jj) = cr4toci2(FILTERED(ii,jj));
        writefile(ofile,TMP,winfiltered);
        break;
        }
      default:
        {
        PRINT_ERROR("Totally impossible, checked input.")
        throw(unhandled_case_error);
        }
      }

    // ______ Check if all done ______
    if (lastbuffer)
      break;

    // ______ Update indexes in matrices, will be corrected for last block ______
    lineswritten  += numout;            // first block overlap added
    slclinelo     += numout;            // next buffer
    slclinehi     += numout;            // next buffer
    filteredlinelo = OVERLAP;           // index in buffer,
                                        // valid for all middle buffers, but not last
    // ______ Give progress message ______
    if (buffer%tenpercent==0)
      {
      PROGRESS << "FILTAZI: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    } // end loop buffers
  ofile.close();


  // ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltazi", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtazi: scratchlogfiltazi",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* PHASE FILTER COMPLEX INTERFEROGRAM: \t"
    << "\n*******************************************************************"
    <<  infile
    << "\nOutput file filtered : \t\t\t"
    << "\n..."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  // ______ updateproductinfo greps info from this file ______
  ofstream scratchresfile("scratchresfiltazi", ios::out | ios::trunc);
  bk_assert(scratchresfile,"scratchresfiltazi",__FILE__,__LINE__);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_m_filtazi]
    << "\n*******************************************************************"
    << "\nInput_file:                           \t"
    <<  infile
    << "\nData_output_file:                     \t"
    <<  filtaziinput.foname
    << "\nData_output_format:                   \t";
  // ___ report of Mohanty Kamini Kanta i forgot this jun12/03___
  if (filtaziinput.oformatflag==FORMATCR4) scratchresfile <<  "complex_real4";
  if (filtaziinput.oformatflag==FORMATCI2) scratchresfile <<  "complex_short";
  scratchresfile
    << "\nFirst_line (w.r.t. original_image):   \t"
    <<  master.currentwindow.linelo
    << "\nLast_line (w.r.t. original_image):    \t"
    <<  master.currentwindow.linehi
    << "\nFirst_pixel (w.r.t. original_image):  \t"
    <<  master.currentwindow.pixlo
    << "\nLast_pixel (w.r.t. original_image):   \t"
    <<  master.currentwindow.pixhi
    << "\n*******************************************************************"
    << "\n* End_" << processcontrol[pr_m_filtazi] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();
  // ====== Tidy up ======
  } // END azimuthfilter



/****************************************************************
 *    azimuth filter per block                                  *
 * Input is matrix of SIZE (e.g. 1024) lines, and N range pixs. *
 * Input is SLC of master. slave_info gives fDC polynomial      *
 * for slave + coarse offset. HAMMING is alpha for myhamming f. *
 * Filtered OUTPUT is same size as input block.                 *
 * Because of overlap (azimuth), only write to disk in calling  *
 * routine part (in matrix coord.) [OVERLAP:SIZE-OVERLAP-1]     *
 * = SIZE-(2*OVERLAP);  // number of output pixels              *
 *                                                              *
 * Filtering is performed in the spectral domain                *
 * (1DFFT over azimuth for all columns at once)                 *
 * Filter is different for each column due to shift in fd_c     *
 * doppler centroid frequency.                                  *
 *                                                              *
 * ! It should still be studied if real4 matrices are accurate  *
 * enough, but I guess it is (BK).                              *
 *                                                              *
 #%// BK 06-Nov-2000                                            *
 ****************************************************************/
matrix<complr4> blockazifilt(
        const matrix<complr4> &SLCIMAGE,
        const slcimage        &master,          // PRF, BW, fd0
        const slcimage        &slave,           // PRF, BW, fd0
        const real8            HAMMING)
  {
  TRACE_FUNCTION("blockazifilt (BK 06-Nov-2000)")
  const uint SIZE  = SLCIMAGE.lines();  // fftlength
  const uint NCOLS = SLCIMAGE.pixels(); // width
  if (NCOLS != master.currentwindow.pixels())
    WARNING.print("this will crash, size input matrix not ok...");

  // ______ Compute fDC_master, fDC_slave for all columns ______
  // ______ Create axis to evaluate fDC polynomial for master/slave ______
  // ______ fDC(column) = fdc_a0 + fDC_a1*(col/RSR) + fDC_a2*(col/RSR)^2 ______
  // ______ fDC = y = Ax ______
  // ______ Capitals indicate matrices (FDC_M <-> fDC_m) ______
  DEBUG.print("Filtering data by evaluated polynomial fDC for each column."); 
  matrix<real8> xaxis(1,master.currentwindow.pixels());         // lying
  for (int32 i=master.currentwindow.pixlo; i<=int32(master.currentwindow.pixhi); ++i)
    xaxis(0,i-master.currentwindow.pixlo) = real8(i-1.0);
  xaxis /= (master.rsr2x/2.0);
  matrix<real8> FDC_M = (master.f_DC_a1 * xaxis);
  FDC_M += master.f_DC_a0;
  FDC_M += (master.f_DC_a2 * sqr(xaxis));// TODO: better use master.pix2fdc()
  // Bert Kampes, 20-Apr-2005


  // ______ fDC_slave for same(!) columns (coarse offset). ______
  // ______ offset defined as: cols=colm+offsetP ______
  for (int32 i=master.currentwindow.pixlo; i<=int32(master.currentwindow.pixhi); ++i)
    xaxis(0,i-master.currentwindow.pixlo) = real8(i-1.0) + real8(slave.coarseoffsetP);
  xaxis /= (slave.rsr2x/2.0);
  matrix<real8> FDC_S = (slave.f_DC_a1 * xaxis);
  FDC_S += slave.f_DC_a0;
  FDC_S += (slave.f_DC_a2 * sqr(xaxis));

  #ifdef __DEBUG
  DEBUG.print("Dumping matrices fDC_m, fDC_s (__DEBUG defined)");
  dumpasc("fDC_m",FDC_M);
  dumpasc("fDC_s",FDC_S);
  #endif
  
  // ______ Axis for filter in frequencies ______
// TODO check, rather shift, test matlab... or wshift,1D over dim1
// use fft properties to shift...

  const bool dohamming = (HAMMING<0.9999) ? true : false;
  const real8 PRF      = master.prf;            // pulse repetition freq. [Hz]
  const real8 ABW      = master.abw;            // azimuth band width [Hz]

  const real4 deltaf   =  PRF/real4(SIZE);
  const real4 fr       = -PRF/2.0;
  matrix<real4> freqaxis(1,SIZE);
  for (uint i=0; i<SIZE; ++i)
    freqaxis(0,i) = fr+real4(i)*deltaf; // [-fr:df:fr-df]

  matrix<real4> FILTER;                         // i.e., the filter per column
  matrix<real4> FILTERMAT(SIZE,NCOLS);          // i.e., THE filter
  for (uint i=0; i<NCOLS; ++i)
    {
    const real8 fDC_m    = FDC_M(0,i);          // zero doppler freq. [Hz]
    const real8 fDC_s    = FDC_S(0,i);          // zero doppler freq. [Hz]
    const real8 fDC_mean = 0.5*(fDC_m+fDC_s);   // mean doppler centroid freq.
    const real8 ABW_new  = max(1.0,
      2.0*(0.5*ABW-abs(fDC_m-fDC_mean)));       // new bandwidth > 1.0
    if (dohamming)
      {
      // ______ NOT a good implementation for per col., cause wshift AND fftshift.
      // ______ DE-weight spectrum at centered at fDC_m ______
      // ______ spectrum should be periodic! (use wshift) ______
      matrix<real4> inversehamming = myhamming(freqaxis,ABW,PRF,HAMMING);
      for (uint ii=0; ii<SIZE; ++ii)
        if (inversehamming(0,ii)!=0.0)// check zero division
              inversehamming(0,ii) = 1.0/inversehamming(0,ii);
      // ______ Shift this circular by myshift pixels ______
      int32 myshift = int32(rint((SIZE*fDC_m/PRF)));// round
      wshift(inversehamming,-myshift);          // center at fDC_m
  
      // ______ Newhamming is scaled and centered around new mean ______
      myshift = int32(rint((SIZE*fDC_mean/PRF)));// round
      FILTER  = myhamming(freqaxis,
                          ABW_new,
                          PRF,HAMMING);         // fftshifted
      wshift(FILTER,-myshift);                  // center at fDC_mean
      FILTER *= inversehamming;
      }
    else        // no weighting, but center at fDC_mean, size ABW_new
      {
      int32 myshift = int32(rint((SIZE*fDC_mean/PRF)));// round
      FILTER = myrect(freqaxis/real4(ABW_new)); // fftshifted
      wshift(FILTER,-myshift);                  // center at fDC_mean
      }
    ifftshift(FILTER);                          // fftsh works on data!
    FILTERMAT.setcolumn(i,FILTER);
    } // foreach column

  // ______ Filter slcdata ______
  matrix<complr4> FILTERED = SLCIMAGE;
  fft(FILTERED,1);                              // fft foreach column
  #ifdef WIN32
  FILTERED = timesCxR(FILTERED,FILTERMAT);// Jia
  #else
  FILTERED *= FILTERMAT;                        // filter each column
  #endif
  ifft(FILTERED,1);                             // ifft foreach column
  return FILTERED;
  } // END blockazifilt




/****************************************************************
 *    rangefiltporbits                                          *
 * -getoverlap approx between master/slave.                     *
 * -get next line master/slave.                                 *
 * -compute freq. shift due to different incidence angle.       *
 * -filter master and slave.                                    *
 * -write output.                                               *
 #%// BK 13-Nov-2000                                            *
 * BUG detected by Rens Swart (26-Apr-2001):                    *
 *  if fftlength < width/2 then output consists of replicas     *
 *  number of replicas = int(width/fftlength) - 1               *
 ****************************************************************/
void rangefiltporbits(
         const input_gen       &generalinput,
         const input_filtrange &filtrangeinput,
         const input_ell       &ellips,
         const slcimage        &master,
         const slcimage        &slave,
               orbit           &masterorbit,
               orbit           &slaveorbit)
  {
  TRACE_FUNCTION("rangefiltporbits (BK 13-Nov-2000)")

  // ______ RSR,RBW in Hz, SI ______
  const real8 HAMMINGA  = filtrangeinput.hammingalpha;
  const bool dohamming  = (HAMMINGA<0.9999) ? true : false;
  const real8 RSR       = 0.5*master.rsr2x;     // range sampling rate fr 18.96MHz
  const real8 RBW       = master.rbw;           // range band width fs 15.55MHz

  // ______ Use approximate overlap master slave ______
  // ______ s=m+offset --> m=s-offset ______
  // ______ prevent a-b uint (< 0) == very large ______ 
  window wintmp;
  int32 tmp     = slave.currentwindow.linelo - slave.coarseoffsetL;
  wintmp.linelo = (tmp<1) ? 1 : tmp;
  tmp           = slave.currentwindow.linehi - slave.coarseoffsetL;
  wintmp.linehi = (tmp<1) ? 1 : tmp;
  tmp           = slave.currentwindow.pixlo  - slave.coarseoffsetP;
  wintmp.pixlo  = (tmp<1) ? 1 : tmp;
  tmp           = slave.currentwindow.pixhi  - slave.coarseoffsetP;
  wintmp.pixhi  = (tmp<1) ? 1 : tmp;
  // ______ master.currentwin and wintmp are now in same system ______
  window overlap = getoverlap(master.currentwindow,wintmp);

  TRACE << "overlap window: " << overlap.linelo << " " << overlap.linehi
       << " " << overlap.pixlo << " " << overlap.pixhi;
  TRACE.print();

  // ______ Initialize loop parameters etc. ______
  int32 FFTLENGTH = filtrangeinput.fftlength;
  const int32 NCOLS     = overlap.pixels();             // width 2b processed
  while (FFTLENGTH>NCOLS)
    {
    WARNING << "FFTLENGTH>overlap: New length: "
         << FFTLENGTH << " --> " << FFTLENGTH/2;
    WARNING.print();
    FFTLENGTH /= 2;
    }
  const int32 NBLOCKS   = int32((NCOLS/FFTLENGTH));
  const int32 EXTRA     = ((NCOLS%FFTLENGTH)==0) ? 0:1; // last overlaps previous

  cn M;                 // master position
  cn S;                 // slave position
  cn P;                 // P on ellipsoid, perpendicular to orbits;
  const int32 MAXITER   = 10;
  const real8 CRITERPOS = 1e-6;
  const real8 CRITERTIM = 1e-10;
  INFO << "rangefiltporbits: MAXITER: "   << MAXITER   << "; "
       << "CRITERPOS: " << CRITERPOS << " m; "
       << "CRITERTIM: " << CRITERTIM << " s";
  INFO.print();

  // ====== Open output files ======
  ofstream of_m;
  openfstream(of_m,filtrangeinput.fomaster,generalinput.overwrit);
  bk_assert(of_m,filtrangeinput.fomaster,__FILE__,__LINE__);

  ofstream of_s;
  openfstream(of_s,filtrangeinput.foslave,generalinput.overwrit);
  bk_assert(of_s,filtrangeinput.foslave, __FILE__,__LINE__);

  // ______ Shift parameters ______
  //register int32 i;
  const real4 df =  RSR/real4(FFTLENGTH);
  const real4 fr = -RSR/2.;
  matrix<real4> freqaxis(1,FFTLENGTH);
  for (uint i=0; i<freqaxis.pixels(); ++i)
    freqaxis(0,i) = fr+real4(i)*df;
  matrix<real4> inversehamming;
  if (dohamming)
    {
    inversehamming = myhamming(freqaxis,RBW,RSR,HAMMINGA);
    for (uint i=0; i<inversehamming.pixels(); ++i)
      if (inversehamming(0,i)!=0.0)
        inversehamming(0,i)= 1.0/inversehamming(0,i);
    }

  // ______ Some extra info for debuggers ______
  DEBUG << "parameters used: NCOLS " << NCOLS   << " "
       << "FFTLENGTH: " << FFTLENGTH << " "
       << "NBLOCKS: " << NBLOCKS << " "
       << "EXTRA: " << EXTRA;
  DEBUG.print();
  DEBUG << "overlap window: " << overlap.linelo << " " << overlap.linehi
       << " " << overlap.pixlo << " " << overlap.pixhi;
  DEBUG.print();

  // ______ Progress messages ______
  int32 percent    = 10;
  int32 tenpercent = int32((overlap.lines()/10)+.5);    // round
  if (tenpercent==0) tenpercent = 1000;                 // avoid error: x%0
  PROGRESS.print("filtrange:  0%");

  real8 totalblocks = 0;        // to give some info at end
  real8 totaldeltaf = 0;        // to give some info at end

  // ====== Foreach line, read block, filter, write ======
  for (int32 line_m=overlap.linelo; line_m<=int32(overlap.linehi); ++line_m)    // in master coord.
    {
    //const real8 m_tazi = line2ta(line_m,master.t_azi1,master.prf);
    const real8 m_tazi = master.line2ta(line_m);
    M              = masterorbit.getxyz(m_tazi);
    int32 pixlo_m  = overlap.pixlo;
    int32 line_s   = line_m + slave.coarseoffsetL;              // in slave coord. system
    bool lastblock = false;
    // ______ For this line, filter all blocks ______
    for (int32 rangeblock=1; rangeblock<=NBLOCKS+EXTRA; ++rangeblock)
      {
      totalblocks++;
      int32 pixhi_m = pixlo_m + FFTLENGTH - 1;                  // i.e., fftlength pixels
      if (rangeblock==NBLOCKS+EXTRA)                            // only last one
        {
        pixhi_m   = overlap.pixhi;
        pixlo_m   = pixhi_m - FFTLENGTH + 1;
        lastblock = true;
        }

//cerr << "BERT: block, pixlo, hi: "
//<< rangeblock << " " << pixlo_m << " " << pixhi_m << endl;

      // ______ Compute frequency shift for this block ______
      const real8 middlepix = real8(pixlo_m+pixhi_m)/2.; // theta,Bperp here
      // better compute range with m_trange...
      lp2xyz(line_m,middlepix,ellips,master,masterorbit,P,MAXITER,CRITERPOS);
      // ______ Compute xyz slave satellite from P ______
      real8 s_trange,s_tazi;
      xyz2t(s_tazi,s_trange,slave,slaveorbit,P,MAXITER,CRITERTIM);
      // ______ Slave position ______
      S = slaveorbit.getxyz(s_tazi);
      cn MP = M.min(P);
      cn SP = S.min(P);
      const real8 incidence1 = MP.angle(P);             // theta1
      const real8 incidence2 = SP.angle(P);             // theta2
      const real8 deltatheta = incidence1-incidence2;   // counterclockwise
      // ______ (approx. in Hz, see also manual) ______
      const real4 deltaf = real4((-SOL*deltatheta) /
            (master.wavelength*tan(incidence1-filtrangeinput.terrainslope)));
      totaldeltaf += deltaf;

//cerr << "deltaf: " << deltaf << " Hz\n";
      // ______ Read block from disk ______
      const int32 pixlo_s    = pixlo_m + slave.coarseoffsetP;
      const int32 pixhi_s    = pixhi_m + slave.coarseoffsetP;
      const window win_m     (line_m,line_m,pixlo_m,pixhi_m);
      const window win_s     (line_s,line_s,pixlo_s,pixhi_s);

//  cout << "win_m:\n";
//  win_m.disp();
//  cout << "win_s:\n";
//  win_s.disp();

      matrix<complr4> MASTER = master.readdata(win_m);
      matrix<complr4> SLAVE  = slave.readdata(win_s);
      fft(MASTER,2);
      fft(SLAVE,2);

      // ______ Compose filter for master and slave, and filter both ______
      matrix<real4> FILTER;                             // for master, slave=fliplr(m)
      if (dohamming)
        {
        // ______ Newhamming is scaled and centered around new mean ______
        // Plus or min?, master or slave?
        FILTER  = myhamming(freqaxis-real4(.5*deltaf),// new center
                            RBW-deltaf,                 // new width
                            RSR,HAMMINGA);              // fftshifted
        FILTER *= inversehamming;                       // rect f. included
        }
      else                                              // no weighting of spectra
        {
        FILTER  = myrect((freqaxis-real4(.5*deltaf)) /
                          real4(RBW-deltaf));           // fftshifted
        }
      ifftshift(FILTER);        // fftsh works on data!

      // ______ Note that filter_s = fliplr(filter_m) ______
      // ______ and that this is also valid after ifftshift ______
      // ______ I tested what side had to be filtered: master first ______
      // ______ or slave first. Now it should always work. ______
      // ______ But comment out define if you think wrong side is filtered ______
      #define THISside
      #ifdef THISside
      #ifdef WIN32
        MASTER = timesCxR(MASTER,FILTER);// Jia
        FILTER.fliplr();
        SLAVE = timesCxR(SLAVE,FILTER);// Jia
      #else
        MASTER *= FILTER;       // filtered master range spectral domain
        FILTER.fliplr();
        SLAVE *= FILTER;                // filtered slave range spectral domain
        #endif
      #else
      WARNING.print("This is not as theory, sign Bperp defined otherwise?");
      SLAVE *= FILTER;          // filtered slave range spectral domain
      FILTER.fliplr();
      MASTER *= FILTER;         // filtered master range spectral domain
      #endif
      ifft(MASTER,2);           // back to space domain
      ifft(SLAVE,2);            // back to space domain

      // ______ Write to new output file for master and slave ______
      if (lastblock==false)
        {
        switch (filtrangeinput.oformatflag)
          {
          case FORMATCR4:
            {
            of_m << MASTER;
            of_s << SLAVE;
            break;
            }
          // ______ Convert first to ci2 before writing to file ______
          case FORMATCI2:
            {
            matrix <compli16> TMP(1,MASTER.pixels());
            for (uint ii=0; ii<TMP.pixels(); ++ii)
              TMP(0,ii) = cr4toci2(MASTER(0,ii));
            of_m << TMP;
            for (uint ii=0; ii<TMP.pixels(); ++ii)
              TMP(0,ii) = cr4toci2(SLAVE(0,ii));
            of_s << TMP;
            break;
            }
          default:
            {
            PRINT_ERROR("Totally impossible, checked input.")
            throw(unhandled_case_error);
            }
          } // switch
        }
      // ______ Filtered part F(NUMBLOCKS*FFTLENGTH:pixels()-1) ______
      else // only write extra part, which less than fftlength
        {
        const int32 lastwritten = NBLOCKS*FFTLENGTH;
        const int32 numtowrite  = overlap.pixels() - lastwritten;
        const window wintmp
          (0,0,MASTER.pixels()-numtowrite,MASTER.pixels()-1);
        matrix<complr4> TMP_M(wintmp,MASTER);   // construct as part
        matrix<complr4> TMP_S(wintmp,SLAVE);    // construct as part
        switch (filtrangeinput.oformatflag)
          {
          case FORMATCR4:
            {
            of_m << TMP_M;
            of_s << TMP_S;
            break;
            }
          // ______ Convert first to ci2 before writing to file ______
          case FORMATCI2:
            {
            matrix <compli16> TMP(1,TMP_M.pixels());
            for (uint ii=0; ii<TMP.pixels(); ++ii)
              TMP(0,ii) = cr4toci2(TMP_M(0,ii));
            of_m << TMP;
            for (uint ii=0; ii<TMP.pixels(); ++ii)
              TMP(0,ii) = cr4toci2(TMP_S(0,ii));
            of_s << TMP;
            break;
            }
          default:
            {
            PRINT_ERROR("Totally impossible, checked input.")
            throw(unhandled_case_error);
            }
          } // switch format
        } // write output

        // ______ Update block index ______
        // pixlo_m   = pixhi_m - FFTLENGTH + 1;
        pixlo_m   = pixhi_m + 1; // BUGFIX
      } // all blocks

    // ______ Give progress message ______
    if ((line_m-overlap.linelo+1%tenpercent)==0)
      {
      PROGRESS << "FILTRANGE: " << setw(3) << percent << "%";
      PROGRESS.print();
      percent += 10;
      }
    } // all lines of overlap

  // ====== Close files and finish up ======
  of_m.close();
  of_s.close();
  INFO << "Mean frequency shift over " << totalblocks << "blocks = "
       <<  (totaldeltaf/1e6)/totalblocks << " MHz.";
  INFO.print();


  // ====== Write results to file ======
  ofstream scratchlogfile("scratchlogfiltrange", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"filtrange: scratchlogfiltrange",__FILE__,__LINE__);
  scratchlogfile
    << "\n\n*******************************************************************"
    << "\n* RANGE FILTERING BASED ON ORBITS FOR MASTER AND SLAVE:"
    << "\n*******************************************************************"
    << "\nInput file master (format): \t\t\t"
    <<  master.file << " " << master.formatflag
    << "\nInput file slave (format): \t\t\t"
    <<  slave.file << " " << slave.formatflag
    << "\nOutput file filtered master (format): \t\t\t"
    <<  filtrangeinput.fomaster << " ";
    if (filtrangeinput.oformatflag==FORMATCR4)
      scratchlogfile << "complex r4";
    if (filtrangeinput.oformatflag==FORMATCI2)
      scratchlogfile << "complex short int";
  scratchlogfile
    << "\nOutput file filtered slave (format): \t\t\t"
    <<  filtrangeinput.foslave << " ";
    if (filtrangeinput.oformatflag==FORMATCR4)
      scratchlogfile << "complex r4";
    if (filtrangeinput.oformatflag==FORMATCI2)
      scratchlogfile << "complex short int";
  scratchlogfile
    << "\nmean shift: \t\t\t"
    <<  (totaldeltaf/1e6)/totalblocks << " MHz."
    << "\n*******************************************************************\n";
  scratchlogfile.close();

  for (int32 fileID=1; fileID<=2; ++fileID)
    {
    char oresfile[EIGHTY];
    char odatafile[EIGHTY];
    char odataformat[EIGHTY];
    char processcf[EIGHTY];
    int offsetL = 0;
    int offsetP = 0;
    if (filtrangeinput.oformatflag==FORMATCR4)
      strcpy(odataformat,"complex_real4");
    else if (filtrangeinput.oformatflag==FORMATCI2)
      strcpy(odataformat,"complex_short");
    else
      {
      PRINT_ERROR("case error")
      throw(unhandled_case_error);
      }
    switch (fileID)
      {
      case 1:                                                   // master
        strcpy(oresfile,"scratchresMfiltrange");
        strcpy(odatafile,filtrangeinput.fomaster);
        strcpy(processcf,processcontrol[pr_m_filtrange]);       // control flag
        break;
      case 2:                                                   // slave
        strcpy(oresfile,"scratchresSfiltrange");
        strcpy(odatafile,filtrangeinput.foslave);
        strcpy(processcf,processcontrol[pr_s_filtrange]);       // control flag
        offsetL = slave.coarseoffsetL;
        offsetP = slave.coarseoffsetP;
        break;
      default:
        {
        PRINT_ERROR("panic: ID!={1,2}")
        throw(unhandled_case_error);
        }
      }

    // ______ updateproductinfo greps info from this file ______
    ofstream scratchresfile(oresfile, ios::out | ios::trunc);
    bk_assert(scratchresfile,oresfile,__FILE__,__LINE__);
    scratchresfile
      << "\n\n*******************************************************************"
      << "\n*_Start_" << processcf
      << "\n*******************************************************************"
      << "\nMethod_rangefilt:                     \t"
      << "porbits"
      << "\nData_output_file:                     \t"
      <<  odatafile
      << "\nData_output_format:                   \t"
      <<  odataformat
      // ______ s=m+offset --> m=s-offset ______
      // ______ still in own system! ______
      << "\nFirst_line (w.r.t. original_image):   \t"
      <<  overlap.linelo + offsetL
      << "\nLast_line (w.r.t. original_image):    \t"
      <<  overlap.linehi + offsetL
      << "\nFirst_pixel (w.r.t. original_image):  \t"
      <<  overlap.pixlo  + offsetP
      << "\nLast_pixel (w.r.t. original_image):   \t"
      <<  overlap.pixhi  + offsetP
      << "\n*******************************************************************"
      << "\n* End_" << processcf << "_NORMAL"
      << "\n*******************************************************************\n";
    scratchresfile.close();
    } // for master and slave
  } // END rangefiltporbits




