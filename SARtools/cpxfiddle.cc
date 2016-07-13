/****************************************************************
 * cpxfiddle -- fiddle around with complex files.               *
 *      make -n cpxfiddle                                       *
 *      g++ -O -c -ocpxfiddle.o cpxfiddle.c                     *
 *      g++ -O cpxfiddle.o -o cpxfiddle                         *
 *                                                              *
 * TUDelft 10-May-2000                                          *
 * $Revision: 4.06 $  $Date: 2011/10/25 17:34:19 $              *
 *                                                              *
 * TODO:                                                        *
 *  - output binary -ofloat add -oshort to crop slc's           *
 #%// BK 01-Feb-2001 done.                                      *
 *  - -B (swap bytes)                                           *
 #%// BK 12-Mar-2001                                            *
 *  - band interleaved input formats                            *
 *  - -q(uiet)                                                  *
 *  - -o sunraster -q mag: rescale etc.                         *
 *  - handle input defaults better, no lot of options.          *
 #%// BK 20-Apr-2001                                            *
 *  - -S scalebar for sunraster.                                *
 *  - -version info.                                            *
 *  BK 02 APRIL 2003 ********************************************
 *  - removed nocreate ios bit since g++ v3.2 does not know it  *
 *    (not standard C++?)                                       *
 *  - removed cout.setf(ios::binary) for same reason, though    *
 *    this means cpxfiddle does not work properly.  How can we  *
 *    write binary to stdout??                                  *
 *   cout.write((char*)&x, sizeof(x));// seems to work...       *
 *   it actually seems fine, i already did it like this...      *
 *  -bugfix when byteswapping and multilooking                  *
 *  -changed multilooking to normalize by number of looks       *
 *  -for HGT input, or R4 disguished as cpx, added some things. *
 #%// BK 05-Jan-2001                                            *
 * added skip header bytes                                      *
 #%// BK 17-Oct-2003                                            *
 #%// PM        2007                                            *
 * Various fixes                                                *
 #%// MA 19-Jun-2008                                            *
 * added > 4GB filesize support                                 *
 #%// MA 14-Oct-2008                                            *
 * added -r min/max for scalling magnitude,normal,mixed outputs *
 #%// MA 24-Oct-2008                                            *
 * changed handling of float data for scaling computation       *
 * added --ignorenan only for scaling, see help                 *
 * added --fullstat for scaling, see help                       *
 * [TODO] proper scalling during high multilooking              *
 #%// MA 29-Aug-2009                                            *
 * fixed last column value when multilooked to float data       *
 #%// MA 23-MAR-2011                                            *
 * fixed long path name issue affecting other paramters         *
 #%// MA 25-NOV-2011                                            *
 * added missing input comments for -f                          *
 * ignore isinf values in --ignorenan option                    *
 ****************************************************************/

using namespace std;
#include <iostream>     // cerr etc
#include <iomanip>      // binary to stdout
#include <fstream>      // files
#include <cmath>        // atan2
#include <cstring>      // strcmp
#include <unistd.h>     // getopt
#include <getopt.h>     // getopt long options [MA]
#include <cstdlib>      // exit
#include <netinet/in.h> // ntohl htonl htons
#include <limits.h>  // with gcc 4.3.1 CHAR_MAX or <climits>
#include <math.h>       // NAN

// #include <typeinfo.h>   // short type run time id, does not exist??
// test...
// #include <sys/types.h>
//#include <machine/endian.h>


// ====== version id -h ======
//#define SWVERSION "v2.1 02-April-2003"
//#define SWVERSION "v2.2 16-Jun-2003"
//#define SWVERSION "v2.3 20-Jan-2004"
//#define SWVERSION "v2.4 19-Jun-2008"
//#define SWVERSION "v2.5 14-Oct-2008"
#define SWVERSION "v2.6 29-Aug-2009"


// ====== consts ======
typedef unsigned char uchar;

// ====== Typedefs for portability ====== // TODO update accordingly
typedef short int           int16;
typedef unsigned short int  uint16;
typedef int                 int32;  // sizeof (int) = sizeof (long) on 32bit data model; sizeof (int) < sizeof (long) on 64 bit.
typedef unsigned int        uint32;
typedef unsigned int        uint;
typedef long long           int64;  // new definition
typedef unsigned long long  uint64; // new
typedef float               real4;
typedef double              real8;

// ______ distinct flags ______
const int FORMATR8  = 108;      // i/oformat -o float numbytes=%10 for new input
const int FORMATR4  = 104;      // i/oformat -o float numbytes=%10 for new input
const int FORMATI2  = 102;      // i/oformat -o float numbytes=%10 for new input
const int FORMATUC1 = 101;      // i/oformat -o uchar numbytes=%10 for new input
const int SUNRASTER = 201;      // i/oformat -o uchar + header + colmap +extrawidth(?)
const int ASCII     = 999;      // oformat -o ascii
// ______ formats are distinct, but format%10 must equal bytes per elem ______
const int FORMATCUC1= 1;        // ==numbytes per real part of element %10
const int FORMATCC1 = 11;       // ==numbytes per real part of element %10
const int FORMATCI2 = 12;       // ==numbytes per real part of element %10
const int FORMATCI4 = 14;       // ==numbytes per real part of element %10
const int FORMATCR4 = 24;       // ==numbytes per real part of element %10
const int FORMATCR8 = 28;       // ==numbytes per real part of element %10
// ______ distinct flags ______
const int MAGNITUDE = 51;       // output magnitude
const int PHASE     = 52;       // output phase
const int REALPART  = 53;       // output real part
const int IMAGPART  = 54;       // output imag part
const int NORMAL    = 55;       // output as is
const int MIXED     = 56;       // output sunraster/uc1 mag/pha overlay
// ______ should be -1||1 used in loops ______
const int NOMIRROR  =  1;       // do no mirroring has to be  1 (used later)
const int DOMIRROR  = -1;       // do mirroring has to be -1 (used later)

// ====== Input variables ======
struct commandlineinput
  {
  //char  ifile[128];     // last argv  
  char  ifile[1024];     // last argv [Prabu reported the long path names problem : MA fixed]
  int   linelength;     // -w mandatory
  float exponent;       // -e [1.0]
  float scale;          // -s [1.0]
  int   firstline;      // -l [1]
  int   lastline;       // -L [all]
  int   firstpixel;     // -p [1]
  int   lastpixel;      // -P [all]
  int   sublines;       // -S [1/1]
  int   subpixels;      // -S [1/1]
  int   multilookL;     // -M [1/1]
  int   multilookP;     // -M [1/1]
  int   iformat;        // -f ci2 [cr4]
  int   oformat;        // -o float uchar sunraster short [ascii]
  int   dooutput;       // -q [mag] phase real imag
  int   verbose;        // -V
  int   mirrorX;        // -m XY YX X Y
  int   mirrorY;        // -m XY YX X Y
  int   headerlength;   // -H 32 number of bytes to skip before matrix starts
  //char  cmap[128];      // -c filename or identifier
  char  cmap[512];      // -c filename or identifier
  bool  dontohx;        // -B s swapbytes, big endian to host
  bool  dohtonx;        // -B b swapbytes, host to big endian.
  bool  scalebar;       // -b add scalebar for sunraster
  //bool  quiet;        // -q TODO, -V
  //
  // ______ Derived ______
  int   numlines;       // computed...
  int   bytesperelement;// sizeof real part
  int   bytesperpixel;  // sizeof complex pixel
  int   bytesperline;   // for seekg
  bool  realinput;      // for simplicity
  double   rangemin;    // -r [0.1/1.0]  minimum value of data, the same as clims [ 0.1 1.0 ] in Matlab [MA]
  double   rangemax;    // -r [0.1/1.0]  maximum value of data, used for scalling of data.
  bool     ignorenan ;  // ignore nan values during scalling
  bool     fullstat ;   // during scalling use fullstatistic generation
  };


// ====== Prototypes of functions in this file ======
void shortexpl();
void usage(char *programname);
void synopsis(char *programname);
bool handleinput(int argc, char* argv[], commandlineinput &input);
template <class Type>
void functie(Type, Type, const commandlineinput &input);
void rescale(uchar *OUT, float *IN, const int numin,
  const float MIN, const float MAX, const float NEWMIN, const float NEWMAX);
void rasterheader(unsigned int header[8], const int width, const int height);
void makecmap(unsigned char CMAP[3][256], const commandlineinput &input);
void makecmapmixed(unsigned char CMAP[3][256], const commandlineinput &input);




/****************************************************************
 * main                                                         *
 ****************************************************************/
int main(int argc,char* argv[])
  {
  // ====== Handle input ======
  // "what `which prog`"
  char ident[] = "@(#)Doris software, $Revision: 4.06 $, $Author: TUDelft $";
  //if (input.verbose)
  cerr << ident << endl;
  commandlineinput input;
  if (handleinput(argc,argv,input)==false) synopsis(argv[0]);

  // ______ Call template function to deal with different io ______
  char dummy = '1';
  switch (input.iformat)
    {
    case FORMATCC1:     // complex char
      functie(char(dummy),char(dummy),input);
      break;
    case FORMATUC1:     // byte, fall through
    case FORMATCUC1:    // unsigned complex char
      functie(uchar(dummy),uchar(dummy),input);
      break;
    case FORMATI2:      // short signed int, fall through
    case FORMATCI2:     // complex short integer
      functie(short(dummy),short(dummy),input);
      break;
    case FORMATCI4:     // complex integer
      functie(int(dummy),int(dummy),input);
      break;
    case FORMATR4:      // float, fall through
    case FORMATCR4:     // complex float
      functie(float(dummy),float(dummy),input);
      break;
    case FORMATR8:      // double, fall through
    case FORMATCR8:     // complex double
      functie(double(dummy),double(dummy),input);
      break;
    default: 
      cerr << argv[0] << ": PANIC: unknown format...\n";
    }

  // ====== Tidy up ======
  cerr << endl;
  return 0;
  } // END main



/****************************************************************
 * template read/comp/write function whole matrix               *
 * may be speeded up by defining a function readline, but       *
 * compiler can optimize this as well (?).                      *
 #%// BK 11-May-2000                                            *
 ****************************************************************/
template <class Type>
void functie(Type realpart, Type imagpart, const commandlineinput &input)
  {
    if (sizeof(Type) != input.bytesperelement) 
      {cerr << "PANIC: impossible...\n"; exit(1);}
    ifstream inf(input.ifile, ios::in | ios::binary);
    if (!inf)
      {
        cerr << "cpxfiddle: ERROR: opening input file: " 
             << input.ifile << "\n";
        exit(1);
      }
    
    // register int i,j,k,realindex,start;
    // MA larger file size
    // MA this should be signed especially for j see case Magnitude
    // it is better to rearrange types int or long for these variables
    // j is a potential int whereas realindex seem unsigned long long, check again! 200806
    register long long i,j,k,realindex; 
    register streamoff start;
    // ___ it seems this is not allowed in g++v3.2.  How to do it? ___
    // use new way if binary to stdout: cout.write((char*)&x, sizeof(x))
    //if (input.oformat!=ASCII) cout.write((char*)&x, sizeof(x));// binary to stdout!
    //if (input.oformat!=ASCII) cout.setf(ios::binary); // !! binary to stdout
    //if (input.oformat!=ASCII) cout.setf(ios::binary); // !! binary to stdout
    const float PI = 4.0*atan(1.0);
    
    // ______ read in the whole line (even if subsampling) 
    // ______ because simplicity for mirroring in X. (later multilooking?)
    // ______ fill output 0:numoutput as 2b written to stdout ______
    // ______ if MIXED: OUTPUT[1sthalf]==MAG, 2ndhalf==PHASE ______
    Type  *LINE;
    float *OUTPUT;                                      // 2* for NORMAL output
    LINE   = new Type[2*(input.lastpixel-input.firstpixel+1)];  // allocate memory
    OUTPUT = new float[2*(input.lastpixel-input.firstpixel+1)]; // allocate memory
    
    // ______ Allocate and set stuff for multilooking ______
    // ______ sublines = multilook set in 'handle input' ______
    Type  *LINEML;      // multilooked line
    if (input.multilookL != 1)  // ie. >1
      LINEML = new Type[2*(input.lastpixel-input.firstpixel+1)];    // allocate memory
    
    if (input.verbose)
      {
        int alloc = (int)(sizeof(LINE[0]) * 2 * (input.lastpixel - input.firstpixel+1) +
          sizeof(OUTPUT[0]) * 2 * (input.lastpixel-input.firstpixel+1));
        if (input.multilookL != 1)
          alloc += (int)(sizeof(LINEML[0]) * 2 * (input.lastpixel-input.firstpixel+1));
        cerr << "cpxfiddle: INFO: Allocating " << int(alloc/1000.) << "KB\n";
      }
    
    // ______ Indices to loop forward/backward over file ______
    const int starti = (input.mirrorY==NOMIRROR) ? input.firstline  :  input.lastline;
          int stopi  = (input.mirrorY==NOMIRROR) ? input.lastline   : -input.firstline;
    const int plusi  = (input.mirrorY==NOMIRROR) ? input.sublines   : -input.sublines;
    const int startj = (input.mirrorX==NOMIRROR) ? input.firstpixel :  input.lastpixel;
          int stopj  = (input.mirrorX==NOMIRROR) ? input.lastpixel  : -input.firstpixel;
    const int plusj  = (input.mirrorX==NOMIRROR) ? input.subpixels  : -input.subpixels;

    // ______ Correct index for multilooking, avoid reading to far ______
    if (input.multilookL!=1)
      stopi = int((stopi/plusi) * plusi);
    //stopi = int(floor(stopi/plusi) * plusi);
    if (input.multilookP!=1)
      stopj = int((stopj/plusj) * plusj);
    //stopj = int(floor(stopj/plusj) * plusj);
    
    // ______ number of elements to output per line ______
    int numoutput = int(ceil(float(input.lastpixel-input.firstpixel+1.)/
                             float(input.subpixels)));
    if (input.multilookP!=1)
      numoutput = int(floor(float(input.lastpixel-input.firstpixel+1.)/
                            float(input.subpixels)));
    if (input.dooutput==NORMAL)
      {
        if (input.realinput == true)
          {
            if (input.verbose==true)
              {
                cerr << "cpxfiddle: INFO: real input format detected.\n";
                cerr << "cpxfiddle: INFO: Number of (real) output pixels: "
                     << numoutput << "\n";
              }
          }
        else
          {
            if (input.verbose==true)
              {
                cerr << "cpxfiddle: INFO: Number of (complex) output pixels: "
                     << numoutput << "\n";
              }
            numoutput*=2;       // output real&imag part, size of array.
          }
      }
    else
      cerr << "cpxfiddle: INFO: Number of output pixels: " << numoutput << "\n";
    
    // ______ number of elements output height ______
    int height = int(ceil(float(input.lastline-input.firstline+1.)/
                          float(input.sublines)));
    if (input.multilookL!=1)
      height = int(floor(float(input.lastline-input.firstline+1.)/
                         float(input.sublines)));
    if (input.verbose==true)
      cerr << "cpxfiddle: INFO: Number of output lines:  " << height << "\n";
    
    
    // ______ Write sunraster header (and cmap) if requested ______
    double meanmag   = -999999999.999;// to rescale sunraster magnitude image
    double min_input = +999999999.999;
    double max_input = -999999999.999;
    if (input.oformat==SUNRASTER)
      {
        // ______ Header (cmap length==256*3) ______
        //      unsigned long int HEADER[8];    // 8 long ints
        unsigned int HEADER[8];                 // 8 long ints
        rasterheader(HEADER,numoutput,height);  // create header
        if (input.verbose==true)
          cerr << "Writing 32 byte header of SUNraster file.\n";
        /*
          << "\nmagic#:    " << ntohl(HEADER[0])
          << "\nwidth:     " << ntohl(HEADER[1])
          << "\nheight:    " << ntohl(HEADER[2])
          << "\ndepth:     " << ntohl(HEADER[3])
          << "\nlength:    " << ntohl(HEADER[4])
          << "\ntype:      " << ntohl(HEADER[5])
          << "\nmaptype:   " << ntohl(HEADER[6])
          << "\nmaplength: " << ntohl(HEADER[7])
          << endl;
        */
        cout.write((char*)&HEADER[0],32);               // binary to stdout
        // ______ Cmap ______
        // unsigned maplength = ntohl(HEADER[7]);       // 3 * length per color
        unsigned int maplength = 3*256;         // 3 * length per color
        unsigned char CMAP[3][256];                     // colormap see rasterheader
        if (input.dooutput==MIXED)
          makecmapmixed(CMAP,input);            // only [3][256]
        else
          makecmap(CMAP,input);                 // only [3][256]
        cout.write((char*)&CMAP[0],maplength);  // binary to stdout
        
        // ______ For magnitude output, get mean value approximately ______
        // ______ BK 4-apr-2001 ______
        string stattag;                                 // [MA] used in verbose
        (input.realinput == false) ? stattag = "magnitude" : stattag = "value"  ;
        if ( (input.dooutput==MAGNITUDE || input.dooutput==MIXED ||
            input.dooutput==NORMAL)  && 
            ( isnan(input.rangemin) || isnan(input.rangemax) ) ) // assume doris tricked for hgt.
          {
            meanmag         = 0.0;
            int numsamples  = 0;
            int samplesizeY = input.sublines;  // subsample per x lines...
            int samplesizeX = input.subpixels;
            int nX = (input.lastline-input.firstline+1)/samplesizeY;   // calculated for while loop
            int nY = (input.lastpixel-input.firstpixel+1)/samplesizeX;
//cerr << "fullstat tag: " << input.fullstat << endl;
            while (nX*nY>250000 && !input.fullstat ) // 200x200 window is ok in general // insert as !input.fullstat than use window
              {
                samplesizeX *= 2;
                samplesizeY *= 3; // reading file iritating...
                nX = (input.lastline-input.firstline+1)/samplesizeY;
                nY = (input.lastpixel-input.firstpixel+1)/samplesizeX;
              }
            for (i=input.firstline; i<=input.lastline; i+=samplesizeY)
              {
                start  = (streamoff) input.bytesperline*(i-1) +                  
                                    (input.firstpixel-1)*input.bytesperpixel;
                start += input.headerlength;// account for header.
                inf.seekg(start,ios::beg);
                inf.read((char*)&LINE[0],
                         (input.lastpixel-input.firstpixel+1)*input.bytesperpixel); 
                for (j=input.firstpixel; j<=input.lastpixel; j+=samplesizeX)
                  {
                    //realindex = 2*(j-input.firstpixel);               // index in LINE 
                    realindex = (input.realinput == false) ? 2*(j-input.firstpixel) : // false --> cpx [default]
                                                               (j-input.firstpixel);


                    if (input.ignorenan)  // [MA] ignore nan values during scaling
                      {
                      if  ( (input.realinput == false) && (isnan( LINE[realindex] ) || isnan( LINE[realindex+1] ) || 
                                                           isinf( LINE[realindex] ) || isinf( LINE[realindex+1] ) ) ) 
                        { //cerr << "cpxskip nan...\n"; 
                          continue; } // cpx 
                      else if  ( isnan( LINE[realindex] ) || isinf( LINE[realindex] )) 
                        {// cerr << "skip nan... line: "<< i << " pixel: " << j << "\n";
                         //cerr << "skip nan... \n";
                          continue; }   //  real 
                      }

                    numsamples++;
                    realpart  = LINE[realindex];
                    imagpart  = LINE[realindex+1];
                    // ______ Swap bytes if -B option specified ______
                    if (input.dontohx==true)
                      {
                        if (sizeof(Type)==sizeof(short))
                          {
                            realpart = ntohs(short(realpart));
                            imagpart = ntohs(short(imagpart));
                          }
                        else if (sizeof(Type)==sizeof(int))
                          {
                            realpart = ntohl(short(realpart));
                            imagpart = ntohl(short(imagpart));
                          }
                      }
                    if (input.dohtonx==true)
                      {
                        if (sizeof(Type)==sizeof(short))
                          {
                            realpart = htons(short(realpart));
                            imagpart = htons(short(imagpart));
                          }
                        else if (sizeof(Type)==sizeof(int))
                          {
                            realpart = htonl(short(realpart));
                            imagpart = htonl(short(imagpart));
                          }
                      }
                    if (input.exponent==1.)
                      {
                      meanmag  += ( input.realinput == false ) ? sqrt(float(realpart)*float(realpart)+
                                                                 float(imagpart)*float(imagpart)) 
                                                               : float(realpart);  // non-complex due usual mean
                     }
                    else
                      meanmag  += ( input.realinput == false ) ? pow(sqrt(float(realpart)*float(realpart)+
                                                                 float(imagpart)*float(imagpart)),
                                                                 input.exponent)
                                                               : pow(float(realpart),input.exponent);
                    // --- assume NORMAL input, unwrapped phase ---
                    if (realpart<min_input) min_input = realpart; 
                    if (input.realinput == false && imagpart<min_input) min_input = imagpart;
                    if (realpart>max_input) max_input = realpart;
                    if (input.realinput == false && imagpart>max_input) max_input = imagpart;
                  }
              }
            //meanmag = input.scale*meanmag/numsamples;// finally...
            meanmag = meanmag/numsamples;// finally...
            if (input.verbose==true || input.ignorenan)
              {
                cerr << "mean " << stattag << " based on " << numsamples
                     << " samples: " << meanmag << endl;
                cerr << "min  " << stattag << " based on " << numsamples
                     << " samples: " << min_input << endl;
                cerr << "max  " << stattag << " based on " << numsamples
                     << " samples: " << max_input << endl;
              }
          }
        else if ( input.dooutput==MAGNITUDE || input.dooutput==MIXED 
                                            || input.dooutput==NORMAL )      // [MA] then data rangemin,rangemax defined thus use them later to scale data range.
          {
            min_input = input.rangemin;
            max_input = input.rangemax;
            meanmag   = (max_input)/2.0;
            if (input.verbose==true || !isnan(input.rangemin) )
              {
                cerr << "mean " << stattag << " based on input" 
                     << " (rmin+rmax)/2: \t" << meanmag << endl;
                cerr << "min  " << stattag << "  based on input"
                     << " rangemin: \t" << min_input << endl;
                cerr << "max  " << stattag << "  based on input" 
                     << " rangemax: \t" << max_input << endl;
              }
               
          } 
      }
    
    // ______ Obtain information for placing scalebar for sunraster ______
    int outputlinecnt  = 0;
    int scalebarheight = int(height/25.);       // of total height in output lines
    int scalebarwidth  = int(numoutput/5.);     // of total width in output pixels
    int line0scalebar  = height - 2*scalebarheight;     // position lower
    int lineNscalebar  = line0scalebar+scalebarheight;  // position
    int pixel0scalebar = numoutput -scalebarheight-scalebarwidth;// position right
    
    // ______ Swap bytes to host order ______
    if (input.dontohx==true)
      if (input.verbose==true)
        cerr << "Swapping bytes (big endian to host type, function ntohs/ntohl).\n";
    if (input.dohtonx==true)
      if (input.verbose==true)
        cerr << "Swapping bytes (host type to big endian, function htons/htonl).\n";
    
    // ______ 10% should be multiple of plusi ______
    bool giveprogress = (numoutput > 1000) ? true : false;
    if (input.verbose==false) giveprogress = false;
    int onepercent = int(floor(abs(starti-stopi)/100.));
    //if (giveprogress) while ((tenpercent%plusi)) tenpercent--;
    if (giveprogress) while ((onepercent%plusi)) onepercent--;
    int percent = 0;
    
    // _____ Do normlize by number of looks, but not TYPE _______
    double numlooks = double(input.multilookL)*double(input.multilookP);
    if (input.verbose==true)
      cerr << "cpxfiddle: INFO: Number of looks:         " << numlooks << "\n";
    
    
    // ====== Read file, compute output, write to stdout ======
    // for (i=input.firstline; i<=input.lastline; i+=input.sublines)
    // ______ mirrorY in loop initialization ______
    for (i=starti; i*input.mirrorY<=stopi; i+=plusi)    // forwards or backwards
      {
        if (giveprogress==true)
          {
            if (!((i-starti)%onepercent))
              {
                if (input.verbose==true)
                cerr << "\rprocessed: " << percent << "%";
              percent += 1;
              //percent += 10;
            }
        }
      
      // ______ start is just before line i to be read ______
      start  = input.bytesperline*(i-1) + (input.firstpixel-1)*input.bytesperpixel;
      start += input.headerlength;// account for header.
      inf.seekg(start,ios::beg);
      inf.read((char*)&LINE[0],(input.lastpixel-input.firstpixel+1)*input.bytesperpixel); 
      
      // ______ Swap bytes if -B option specified ______
      if (input.dontohx==true)
        {
          if (sizeof(Type)==sizeof(short))
            {
              for (k=0; k<2*(input.lastpixel-input.firstpixel+1); ++k)
                {
                  LINE[k] = ntohs(short(LINE[k]));
                }
            }
          else if (sizeof(Type)==sizeof(int))
            {
              for (k=0; k<2*(input.lastpixel-input.firstpixel+1); ++k)
                {
                  LINE[k] = ntohl(short(LINE[k]));
                }
            }
          else
            cerr << "seems wrong, -B and not 2B/4B type, ignoring.\n";
        }
      if (input.dohtonx==true)
        {
          if (sizeof(Type)==sizeof(short))
            {
              for (k=0; k<2*(input.lastpixel-input.firstpixel+1); ++k)
                {
                  LINE[k] = htons(short(LINE[k]));
                }
            }
          else if (sizeof(Type)==sizeof(int))
            {
              for (k=0; k<2*(input.lastpixel-input.firstpixel+1); ++k)
                {
                  LINE[k] = htonl(short(LINE[k]));
                }
            }
          else
            cerr << "seems wrong, -B and not type 2B/4B, ignoring.\n";
        }
      
      
      // ______ Do multilooking in line direction here ______
      for (int mll=1; mll<input.multilookL; ++mll)
        {
          // ______ Read next/previous line and add it to LINE ______
          start  = input.bytesperline*(i-1+(mll*input.mirrorY)) +
            (input.firstpixel-1)*input.bytesperpixel;
          start += input.headerlength;// account for header.
          inf.seekg(start,ios::beg);
          inf.read((char*)&LINEML[0],
                   (input.lastpixel-input.firstpixel+1)*input.bytesperpixel); 
          for (int pix=0; pix<2*(input.lastpixel-input.firstpixel+1); ++pix)
            {
              // ______ Swap bytes if -B option specified ______
              if (input.dontohx==true)
                {
                  if (sizeof(Type)==sizeof(short))
                    LINEML[pix] = ntohs(short(LINEML[pix]));
                  else if (sizeof(Type)==sizeof(int))
                    LINEML[pix] = ntohl(short(LINEML[pix]));
                }
              if (input.dohtonx==true)
                {
                  if (sizeof(Type)==sizeof(short))
                    LINEML[pix] = htons(short(LINEML[pix]));
                  else if (sizeof(Type)==sizeof(int))
                    LINEML[pix] = htonl(short(LINEML[pix]));
                }
              // ______ Finally multilook ______
              LINE[pix] += LINEML[pix]; // add real and imag part
            }
        }
      
      // ______ Compute output, store in array OUTPUT ______
      // ______ multilook here in X (range) direction ______
      // ______ store in LINE[j,x,x,j,x,x,j,..]; j are multilooked values ______
      register int indexoutput = 0;                     // j is input index
      double dbl_real;// for normlaization with numlooks
      double dbl_imag;
      switch (input.dooutput)
        {
        case MAGNITUDE:
          // ______ j is index in pixel notation e[first,last] ______
          for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
            {
              realindex = 2*(j-input.firstpixel);               // index in LINE 
              realpart  = LINE[realindex];
              imagpart  = LINE[realindex+1];
              // ______ Multilook in X direction ______
              for (int mlp=1; mlp<input.multilookP; ++mlp)
                {
                  realpart += LINE[realindex+  2*(mlp*input.mirrorX)];
                  imagpart += LINE[realindex+1+2*(mlp*input.mirrorX)];
                }
              dbl_real = double(realpart) / numlooks;
              dbl_imag = double(imagpart) / numlooks;
              OUTPUT[indexoutput] = (input.exponent==1.) ?
                input.scale * sqrt(float(dbl_real*dbl_real+dbl_imag*dbl_imag)) :
                input.scale * pow(sqrt(float(dbl_real*dbl_real+dbl_imag*dbl_imag)),
                                  input.exponent);

              if (input.ignorenan)            // [MA] ignore nan values, assign 0 
                if  ( isnan(OUTPUT[indexoutput]) || isinf(OUTPUT[indexoutput])) 
                                                      OUTPUT[indexoutput] = 0;
                
              indexoutput++;
            }
          break; // case
          
      case PHASE:
        for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
          {
            realindex = 2*(j-input.firstpixel);
            realpart  = LINE[realindex];
            imagpart  = LINE[realindex+1];
            // ______ Multilook in X direction ______
            for (int mlp=1; mlp<input.multilookP; ++mlp)        // number of ml pixels
              {
                realpart += LINE[realindex+  2*(mlp*input.mirrorX)];
                imagpart += LINE[realindex+1+2*(mlp*input.mirrorX)];
              }
            dbl_real = double(realpart) / numlooks;
            dbl_imag = double(imagpart) / numlooks;
            OUTPUT[indexoutput] = (input.exponent==1.) ?
              input.scale * atan2(float(dbl_imag),float(dbl_real)) :
              input.scale * pow(atan2(float(dbl_imag),float(dbl_real)),
                                input.exponent);
            indexoutput++;
          }
        break; // case
        
        case MIXED:
          // ______ OUTPUT[0:numoutput-1] for MAG ______
          // ______ OUTPUT[numoutput:2*numoutput-1] for PHASE ______
          // ______ j is index in pixel notation e[first,last] ______
          for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
            {
              realindex = 2*(j-input.firstpixel);               // index in LINE 
              realpart  = LINE[realindex];
              imagpart  = LINE[realindex+1];
              // ______ Multilook in X direction ______
              for (int mlp=1; mlp<input.multilookP; ++mlp)
                {
                  realpart += LINE[realindex+  2*(mlp*input.mirrorX)];
                  imagpart += LINE[realindex+1+2*(mlp*input.mirrorX)];
                }
              dbl_real = double(realpart) / numlooks;
              dbl_imag = double(imagpart) / numlooks;
              OUTPUT[indexoutput] = (input.exponent==1.) ?
                input.scale * sqrt(float(dbl_real*dbl_real+dbl_imag*dbl_imag)) :
                input.scale * pow(sqrt(float(dbl_real*dbl_real+dbl_imag*dbl_imag)),
                                  input.exponent);
              // ______ PHASE in second half of output array (NO scale/exp) ______
              OUTPUT[indexoutput+numoutput] =
                atan2(float(dbl_imag),float(dbl_real));
              indexoutput++;
            }
          break; // case mixed
          
        case REALPART:
          for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
            {
              realindex = 2*(j-input.firstpixel);
              realpart  = LINE[realindex];
              // ______ Multilook in X direction ______
              for (int mlp=1; mlp<input.multilookP; ++mlp)      // number of ml pixels
                realpart += LINE[realindex + 2*(mlp*input.mirrorX)];
              
              dbl_real = double(realpart) / numlooks;
             
              OUTPUT[indexoutput] = (input.exponent==1.) ?
                input.scale * dbl_real :
                input.scale * pow(float(dbl_real),input.exponent);
              
              indexoutput++;
            }
          break; // case

      case IMAGPART:
        for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
          {
          realindex = 2*(j-input.firstpixel);
          imagpart  = LINE[realindex+1];
          // ______ Multilook in X direction ______
          for (int mlp=1; mlp<input.multilookP; ++mlp)  // number of ml pixels
            imagpart += LINE[realindex+1+2*(mlp*input.mirrorX)];
          dbl_imag = double(imagpart) / numlooks;
          OUTPUT[indexoutput] = (input.exponent==1.) ?
            input.scale * dbl_imag :
            input.scale * pow(float(dbl_imag),input.exponent);
          indexoutput++;
          }
        break; // case
      case NORMAL:
        {
        const uint16 f = (input.realinput == false) ? 2 :  // index factor: complex vs real
                                                      1 ;  // [MA] fixes last coln bug
        for (j=startj; j*input.mirrorX<=stopj; j+=plusj)  // forwards or backwards
          {
          // report and bugfix by Dietrich Wanke, for r4 mirroring in X #%// BK 30-Jul-2003
          //realindex = 2*(j-input.firstpixel);
          //realindex = (input.iformat != FORMATR4) ? 2*(j-input.firstpixel) : 
          //realindex = (input.realinput == false) ? 2*(j-input.firstpixel) : 
          //                                           (j-input.firstpixel);  // not nec since f factor is introduced
          realindex = f*(j-input.firstpixel);
          realpart  = LINE[realindex];
          imagpart  = LINE[realindex+1];
          // ______ Multilook in X direction ______
          for (int mlp=1; mlp<input.multilookP; ++mlp)  // number of ml pixels
            {

            //realpart += LINE[realindex+  2*(mlp*input.mirrorX)];
            //imagpart += LINE[realindex+1+2*(mlp*input.mirrorX)];
            realpart += LINE[realindex+  f*(mlp*input.mirrorX)];   // complex f==2, real f==1
            imagpart += LINE[realindex+1+f*(mlp*input.mirrorX)];
            }
          dbl_real = double(realpart) / numlooks;
          dbl_imag = double(imagpart) / numlooks;
          OUTPUT[indexoutput] = (input.exponent==1.) ?
            input.scale * dbl_real :
            input.scale * pow(float(dbl_real),input.exponent);
          indexoutput++;
          // made conditional followingbugfix by Dietrich Wanke #%// BK 30-Jul-2003
          //if (input.iformat != FORMATR4)// made conditional for r4 output
          if (input.realinput == false)// made conditional for r4 output
            {
            OUTPUT[indexoutput] = (input.exponent==1.) ?
              input.scale * dbl_imag :
              input.scale * pow(float(dbl_imag),input.exponent);
            indexoutput++;
            }
          }
        }
        break; // case
      default:
        cerr << "PANIC: impossible output request...\n";
      } // switch output type

    // ______ Ouptut to stdout, ______
    // ______ compiler moves switch outside loop for speed? ______
    outputlinecnt++;
    switch (input.oformat)
      {
      case ASCII:
        {
        for (k=0; k<numoutput; ++k)
          cout << OUTPUT[k] << " ";             // let compiler deal with format
        break;
        }
      case FORMATR4:
        {
        // ______ Write binary to stdout (redirect to file or pipe to next command) ______
        const int sizeoutput = sizeof(float);
        cout.write((char*)&OUTPUT[0],numoutput*sizeoutput);
        break;
        }
      case FORMATI2:
        {
        // ______ Write binary to stdout ______
        // ______ This is a cast, and not very smart, but works ______
        // BK 01-Feb-2001
        const int sizeoutput = sizeof(short);
        short *OUTPUTI2;
        OUTPUTI2 = new short[numoutput];        // allocate memory
          for (k=0; k<numoutput; ++k)
            OUTPUTI2[k] = short(OUTPUT[k]);
        cout.write((char*)&OUTPUTI2[0],numoutput*sizeoutput);
        break;
        }
      case SUNRASTER:
        {
        // ______ Header and colormap already written ______
        // ______ Fall through to UC1 ______
        }
      case FORMATUC1:
        {
        // ______ Rescale output [min,max] --> [0,255] ______
        // ______ only for -qphase for now ______
        const int sizeoutput = sizeof(char);            // output type is char
        uchar *OUT_UC;
        OUT_UC = new uchar[numoutput];
        switch (input.dooutput)
          {
          case NORMAL:
            {
            // ______ Assume someone is trying to plot a R4 or HGT file ______
            // ______ Fall through to PHASE for plotting______
            // ______But first adapt OUTPUT to range -pi,pi______
            // ______Because min/max may be unwrapped phase______
            for (int ii=0; ii<numoutput; ++ii)
              {
              if      (OUTPUT[ii]<min_input) OUTPUT[ii]=-PI;// avoid clipping
              else if (OUTPUT[ii]>max_input) OUTPUT[ii]= PI;// avoid clipping
              else OUTPUT[ii] = (OUTPUT[ii]-min_input)*
                                (2.0*PI/(max_input-min_input))-PI;
              }
            // ______ Fall through to PHASE for plotting______
            }
          case PHASE:
            rescale(OUT_UC,OUTPUT,numoutput,-PI,PI,0,255);          // fill uc
            // ______ Add scalebar if requested ______
            // ______ by changing the data in lower right corner ______
            if (input.scalebar)
              {
              if (outputlinecnt>=line0scalebar && outputlinecnt<=lineNscalebar)
                {
                // ______ change the data to [0:255] ______
                for (int ii=pixel0scalebar; ii<=pixel0scalebar+scalebarwidth; ++ii)
                  {
                  OUT_UC[ii]= uchar(float(ii-pixel0scalebar)*
                                    (255./float(scalebarwidth)));
                  }
                }
              }
            break;
          case MAGNITUDE:
            {
            // ______ First convert output to float 16:255, ______
            // ______ then make integers. scale influences this. ______
            // ______ so a larger scale factor thresholds more ______
            // ______ data to max (white) ______
            for (int ii=0; ii<numoutput; ++ii)
              {
              OUTPUT[ii] *= (150./meanmag);// scaled mag now has mean==150
              if      (OUTPUT[ii]<16)  OUTPUT[ii]=16;
              else if (OUTPUT[ii]>255) OUTPUT[ii]=255;
              }
            rescale(OUT_UC,OUTPUT,numoutput,0,255,0,255);    // fill uc
            break;
            }
          case MIXED:
            {
            // ______ We have 2 arrays, magnitude: OUPUT, phase: OUT_PHA
            // ______ create OUTPUT = numcolorsphase*MAG+PHASE ______
            float *OUT_PHA;                     // temp. store phase here as float
            OUT_PHA = new float[numoutput];
            for (int ii=0; ii<numoutput; ++ii)
              {
                OUT_PHA[ii] = OUTPUT[ii+numoutput];// phase is stored here
                OUTPUT[ii] *= (150./meanmag);// scaled mag now has mean==150
                if      (OUTPUT[ii]<16)  OUTPUT[ii]=16;
                else if (OUTPUT[ii]>255) OUTPUT[ii]=255;
              }
            uchar *OUT_UC_MAG;                          // temp store for mag.
            OUT_UC_MAG = new uchar[numoutput];
            rescale(OUT_UC_MAG,OUTPUT,numoutput,0,255,0,15);   // fill ucmag with mag
            rescale(OUT_UC,OUT_PHA,numoutput,-PI,PI,0,15);     // fill uc with phase
            // ______ transform data=16*mag+pha ______
            for (int ii=0; ii<numoutput; ++ii)
              OUT_UC[ii] += 16*OUT_UC_MAG[ii];//        e[0:255]
            // ______ Add scalebar if requested ______
            // ______ by changing the data in lower right corner ______
            if (input.scalebar)
              {
                if (outputlinecnt>=line0scalebar && outputlinecnt<=lineNscalebar)
                  {
                    // ______ change the data to [0:255] ______
                    for (int ii=pixel0scalebar; ii<=pixel0scalebar+scalebarwidth; ++ii)
                      {
                        OUT_UC[ii]= uchar(240+(float(ii-pixel0scalebar)*
                                               (16./float(scalebarwidth))));
                      }
                  }
              }
            break; // switch
            }
          default:
            cerr << "unknown option with sunraster output, uc1, assuming phase\n";
            rescale(OUT_UC,OUTPUT,numoutput,-PI,PI,0,255);          // fill uc
          }
        cout.write((char*)&OUT_UC[0],numoutput*sizeoutput); // binary to stdout
        break;
        }
      default:
        cerr << "PANIC: Impossible no oformat...\n";
      } // switch output format, write per output line
    
    // ______ Write some extra req. at end of line ______
    if (input.oformat==ASCII) cout << endl;
    
    // ______ Write extra zero after each line if odd width ______
    // ______ It seems format requires even width? (at least xv does) ______
    if (input.oformat==SUNRASTER)
      {
        if ((numoutput%2)!=0) // originally odd
          {
            // ______ Header and colormap already written ______
            // ______ But write extra 0 at end??? (why... bk) 
            // ______ This seems to be req. for xv?
            uchar sameaslast;
            rescale(&sameaslast,&OUTPUT[numoutput-1],1,-PI,PI,0,255);
            cout.write((char*)&sameaslast,1);
          }
      }
    
      } // loop over lines
    
    if (giveprogress==true) cerr << endl;
  } // END function



/****************************************************************
 * rescale                                                      *
 * rescale data in IN(min,max) --> out[MIN,MAX]                 *
 * o(x) = (i(x)-min) * (newrange/oldrange) + newmin             *
 #%// BK 17-Nov-2000                                            *
****************************************************************/
void rescale(
             uchar *OUT,
             float *IN,
             const int numin,
             const float MIN,
             const float MAX,
             const float NEWMIN,
             const float NEWMAX)
{
  if (MIN>MAX)
    {
      cerr << "wrong input in recale: min>max\n";
      exit (1);
    }
  if (NEWMIN>NEWMAX)
    {
      cerr << "wrong input in recale: newmin>newmax\n";
      exit (1);
    }
  if (NEWMIN<0 || NEWMAX>255)
    {
      cerr << "wrong input in recale: newmin<0 or newmax>255\n";
      exit (1);
    }
  //const float NEWMIN   = 0.;
  //const float NEWMAX   = 255.;
  const float NEWRANGE = NEWMAX-NEWMIN;
  const float RANGE    = MAX-MIN;
  const float SCALE    = NEWRANGE/RANGE;
  for (int ii=0; ii<numin; ++ii)
    {
    //if      (IN[ii]<MIN) OUT[ii]=uchar(NEWMIN);// avoid clipping
    //else if (IN[ii]>MAX) OUT[ii]=uchar(NEWMAX);// avoid clipping
    //else OUT[ii] = uchar(NEWMIN + (((IN[ii])-MIN) * SCALE));
    OUT[ii] = uchar(NEWMIN + (((IN[ii])-MIN) * SCALE));
    }
  } // END rescale



/****************************************************************
 * rasterheader                                                 *
 * or save as gif/jpg/bmp/tiff/etc?                             * 
 *                                                              *
 * ***** SUNraster header defined as 8 integers:                *
 * ras_magic:     [0x59a66a95] magic number                     *
 * ras_width:     [input]      width (pixels) of image          *
 * ras_height:    [input]      height (pixels) of image         *
 * ras_depth:     [8]          depth (1, 8, or 24 bits)         *
 * ras_length:    [computed]   length (bytes) of image          *
 * ras_type:      [1]          file; see RT_* below             *
 * ras_maptype:   [1]          colormap; see RMT_* below        *
 * ras_maplength: [3*256]      length (b) of following map      *
 *                                                              *
 * ***** Sun supported ras_type's RT_ ***                       *
 * 0       Raw pixrect image in 68000 byte order                *
 * 1       Raw pixrect image in 68000 byte order                *
 * 2       Run-length compression of bytes                      *
 * 3       XRGB or RGB instead of XBGR or BGR                   *
 * 4       tiff <-> standard rasterfile                         *
 * 5       iff (TAAC format) <-> standard rasterfile            *
 * 0xffff  Reserved for testing                                 *
 *                                                              *
 * ***** Sun registered ras_maptype's RMT_ ***                  *
 * 0       ras_maplength is expected to be 0                    *
 * 1       red[ras_maplength/3],green[],blue[]                  *
 * 2       RMT_RAW                                              *
 *                                                              *
 * It seems width must be even? (fix this)                      *
 #%// BK 17-Nov-2000                                            *
 ****************************************************************/
void rasterheader(unsigned int HEADER[8], const int width, const int height)
  {
    // ______ if not even width, make it even, and write last uchar twice ______
    const int W = ((width%2)==0) ? width : width+1;
    HEADER[0] = htonl(0x59a66a95);    // magic number, always same
    HEADER[1] = htonl(W);             // #columns, force even width
    HEADER[2] = htonl(height);        // #rows
    HEADER[3] = htonl(8);             // #planes in bits
    HEADER[4] = htonl(W*height);      // size in bytes
    HEADER[5] = htonl(1);             // rect type
    HEADER[6] = htonl(1);             // RMT_RGB follows, size = H[7]
    HEADER[7] = htonl(3*256);         // rgb, foreach possible pixelvalue
  } // END rasterheader


/****************************************************************
 * makecmap                                                     *
 * input option -c: only for sunraster output:                  *
 * generate a cmap from a input file, or from a name "gray", ?  *
 * input file consist of 256 ascii entries r g b                *
 #%// BK 22-Nov-2000                                            *
 ****************************************************************/
void makecmap(
              unsigned char CMAP[3][256],
              const commandlineinput &input)
{
  const int maplength = 256;                    // always
  if (input.oformat!=SUNRASTER)
    {
      cerr << "makecmap: oformat != sunraster.\n";
      exit(1);
    }
  
  // ______ Start internal maps, length 256 ______
  if (!strcmp(input.cmap,"gray"))
    {
      for (int color=0; color<3; ++color)
        {
          for (int ii=0; ii<maplength; ++ii)
            {
              CMAP[0][ii] = uchar(ii*256/maplength);    // red
              CMAP[1][ii] = uchar(ii*256/maplength);    // green
              CMAP[2][ii] = uchar(ii*256/maplength);    // blue
            }
        }
    }
  
  else if (!strcmp(input.cmap,"bert"))  // i.e. random for testing
    {
    for (int color=0; color<3; ++color)
      for (int ii=0; ii<maplength; ++ii)
        CMAP[color][ii] = uchar((color+1)*ii);  // rgb
    }
  
  else if (!strcmp(input.cmap,"jet"))
    {
      const double m = maplength-1;                             // 255
      const double n = floor(((m+1)/4)+0.5);
      double k       = 0.;
      for (int ii=0; ii<int(0.5*n); ++ii)
        {
          k++;
          double y    = (k+n/2)/n;
          CMAP[0][ii] = uchar(0);
          CMAP[1][ii] = uchar(0);
          CMAP[2][ii] = uchar(m*y);
        }
      k = 0.;
      for (int ii=int(0.5*n); ii<int(1.5*n); ++ii)
        {
          k++;
          double x    = k/n;
          CMAP[0][ii] = uchar(0);
          CMAP[1][ii] = uchar(m*x);
          CMAP[2][ii] = uchar(m);
        }
    k = 0.;
    double k2 = n;
    for (int ii=int(1.5*n); ii<int(2.5*n); ++ii)
      {
      k++;
      double x    = k/n;
      double x2   = k2/n;
      CMAP[0][ii] = uchar(m*x);
      CMAP[1][ii] = uchar(m);
      CMAP[2][ii] = uchar(m*x2);
      k2--;
      }
    k2 = n;
    for (int ii=int(2.5*n); ii<int(3.5*n); ++ii)
      {
      double x2   = k2/n;
      CMAP[0][ii] = uchar(m);
      CMAP[1][ii] = uchar(m*x2);
      CMAP[2][ii] = uchar(0);
      k2--;
      }
    k2 = n;
    for (int ii=int(3.5*n); ii<maplength; ++ii)
      {
      double y2   = (k2+n/2)/n;
      CMAP[0][ii] = uchar(m*y2);
      CMAP[1][ii] = uchar(0);
      CMAP[2][ii] = uchar(0);
      k2--;
      }
    } // jet

  else if (!strcmp(input.cmap,"cool"))
    {
    const double m = maplength-1;               // 255
    for (int ii=0; ii<maplength; ++ii)
      {
      CMAP[0][ii] = uchar(ii);                  // r
      CMAP[1][ii] = uchar(m-ii);                // g
      CMAP[2][ii] = uchar(m);                   // b
      }
    } // hot

  else if (!strcmp(input.cmap,"hot"))
    {
    const double m = maplength-1;                               // 255
    const double n = int(m*(3./8.));
    for (int ii=0; ii<=int(m); ++ii)
      {
      CMAP[0][ii] = (ii<n) ? uchar(m*ii/n):uchar(m);            // r
      if      (ii<n)   CMAP[1][ii] = uchar(0);                  // g
      else if (ii<2*n) CMAP[1][ii] = uchar(m*((ii-n)/n));       // g
      else             CMAP[1][ii] = uchar(m);                  // g
      CMAP[2][ii] = (ii<2*n) ? uchar(0):uchar(m*((ii-2*n)/(m-2*n)));  // b
      }
    } // hot

  else if (!strcmp(input.cmap,"hsv"))
    // ______ Created with matlab ph(256)*256, code with awk BK 04-Apr-2001 ______
    // >> map=round(255*ph(256));
    // >> fid=fopen('map.txt','w'); 
    // >> for ii=1:256
    // fprintf(fid,'%i %i %i\n',map(ii,:));
    // end
    // >> fclose(fid);
    // awk '{i++;print "CMAP[0]["i-1"]="$1";", "CMAP[1]["i-1"]="$2";",
    // "CMAP[2]["i-1"]="$3";"}' map.txt > map.code
    // vi: ":r map.code"
    // (spaces in vi: mark start with mb, end with me, then substitute
    // start of lnie with 4 spaces: ":'b,'e s/^/    /"
    {
    CMAP[0][0]=102; CMAP[1][0]=255; CMAP[2][0]=255;
    CMAP[0][1]=104; CMAP[1][1]=251; CMAP[2][1]=255;
    CMAP[0][2]=105; CMAP[1][2]=248; CMAP[2][2]=255;
    CMAP[0][3]=107; CMAP[1][3]=244; CMAP[2][3]=255;
    CMAP[0][4]=108; CMAP[1][4]=241; CMAP[2][4]=255;
    CMAP[0][5]=110; CMAP[1][5]=238; CMAP[2][5]=255;
    CMAP[0][6]=111; CMAP[1][6]=235; CMAP[2][6]=255;
    CMAP[0][7]=113; CMAP[1][7]=231; CMAP[2][7]=255;
    CMAP[0][8]=114; CMAP[1][8]=228; CMAP[2][8]=255;
    CMAP[0][9]=116; CMAP[1][9]=225; CMAP[2][9]=255;
    CMAP[0][10]=117; CMAP[1][10]=222; CMAP[2][10]=255;
    CMAP[0][11]=119; CMAP[1][11]=220; CMAP[2][11]=255;
    CMAP[0][12]=120; CMAP[1][12]=217; CMAP[2][12]=255;
    CMAP[0][13]=122; CMAP[1][13]=214; CMAP[2][13]=255;
    CMAP[0][14]=123; CMAP[1][14]=211; CMAP[2][14]=255;
    CMAP[0][15]=125; CMAP[1][15]=209; CMAP[2][15]=255;
    CMAP[0][16]=126; CMAP[1][16]=206; CMAP[2][16]=255;
    CMAP[0][17]=128; CMAP[1][17]=204; CMAP[2][17]=255;
    CMAP[0][18]=129; CMAP[1][18]=201; CMAP[2][18]=255;
    CMAP[0][19]=131; CMAP[1][19]=199; CMAP[2][19]=255;
    CMAP[0][20]=132; CMAP[1][20]=197; CMAP[2][20]=255;
    CMAP[0][21]=134; CMAP[1][21]=195; CMAP[2][21]=255;
    CMAP[0][22]=135; CMAP[1][22]=193; CMAP[2][22]=255;
    CMAP[0][23]=137; CMAP[1][23]=191; CMAP[2][23]=255;
    CMAP[0][24]=138; CMAP[1][24]=189; CMAP[2][24]=255;
    CMAP[0][25]=140; CMAP[1][25]=187; CMAP[2][25]=255;
    CMAP[0][26]=141; CMAP[1][26]=185; CMAP[2][26]=255;
    CMAP[0][27]=143; CMAP[1][27]=183; CMAP[2][27]=255;
    CMAP[0][28]=144; CMAP[1][28]=182; CMAP[2][28]=255;
    CMAP[0][29]=146; CMAP[1][29]=180; CMAP[2][29]=255;
    CMAP[0][30]=147; CMAP[1][30]=178; CMAP[2][30]=255;
    CMAP[0][31]=149; CMAP[1][31]=177; CMAP[2][31]=255;
    CMAP[0][32]=150; CMAP[1][32]=176; CMAP[2][32]=255;
    CMAP[0][33]=152; CMAP[1][33]=174; CMAP[2][33]=255;
    CMAP[0][34]=153; CMAP[1][34]=173; CMAP[2][34]=255;
    CMAP[0][35]=155; CMAP[1][35]=172; CMAP[2][35]=255;
    CMAP[0][36]=156; CMAP[1][36]=171; CMAP[2][36]=255;
    CMAP[0][37]=158; CMAP[1][37]=170; CMAP[2][37]=255;
    CMAP[0][38]=159; CMAP[1][38]=169; CMAP[2][38]=255;
    CMAP[0][39]=161; CMAP[1][39]=168; CMAP[2][39]=255;
    CMAP[0][40]=162; CMAP[1][40]=167; CMAP[2][40]=255;
    CMAP[0][41]=164; CMAP[1][41]=166; CMAP[2][41]=255;
    CMAP[0][42]=165; CMAP[1][42]=166; CMAP[2][42]=255;
    CMAP[0][43]=166; CMAP[1][43]=165; CMAP[2][43]=255;
    CMAP[0][44]=167; CMAP[1][44]=164; CMAP[2][44]=255;
    CMAP[0][45]=168; CMAP[1][45]=162; CMAP[2][45]=255;
    CMAP[0][46]=169; CMAP[1][46]=161; CMAP[2][46]=255;
    CMAP[0][47]=170; CMAP[1][47]=159; CMAP[2][47]=255;
    CMAP[0][48]=171; CMAP[1][48]=158; CMAP[2][48]=255;
    CMAP[0][49]=172; CMAP[1][49]=156; CMAP[2][49]=255;
    CMAP[0][50]=173; CMAP[1][50]=155; CMAP[2][50]=255;
    CMAP[0][51]=174; CMAP[1][51]=153; CMAP[2][51]=255;
    CMAP[0][52]=175; CMAP[1][52]=152; CMAP[2][52]=255;
    CMAP[0][53]=176; CMAP[1][53]=150; CMAP[2][53]=255;
    CMAP[0][54]=178; CMAP[1][54]=149; CMAP[2][54]=255;
    CMAP[0][55]=179; CMAP[1][55]=147; CMAP[2][55]=255;
    CMAP[0][56]=181; CMAP[1][56]=146; CMAP[2][56]=255;
    CMAP[0][57]=182; CMAP[1][57]=144; CMAP[2][57]=255;
    CMAP[0][58]=184; CMAP[1][58]=143; CMAP[2][58]=255;
    CMAP[0][59]=186; CMAP[1][59]=141; CMAP[2][59]=255;
    CMAP[0][60]=188; CMAP[1][60]=140; CMAP[2][60]=255;
    CMAP[0][61]=190; CMAP[1][61]=138; CMAP[2][61]=255;
    CMAP[0][62]=192; CMAP[1][62]=137; CMAP[2][62]=255;
    CMAP[0][63]=194; CMAP[1][63]=135; CMAP[2][63]=255;
    CMAP[0][64]=196; CMAP[1][64]=134; CMAP[2][64]=255;
    CMAP[0][65]=198; CMAP[1][65]=132; CMAP[2][65]=255;
    CMAP[0][66]=200; CMAP[1][66]=131; CMAP[2][66]=255;
    CMAP[0][67]=202; CMAP[1][67]=129; CMAP[2][67]=255;
    CMAP[0][68]=205; CMAP[1][68]=128; CMAP[2][68]=255;
    CMAP[0][69]=207; CMAP[1][69]=126; CMAP[2][69]=255;
    CMAP[0][70]=210; CMAP[1][70]=124; CMAP[2][70]=255;
    CMAP[0][71]=212; CMAP[1][71]=123; CMAP[2][71]=255;
    CMAP[0][72]=215; CMAP[1][72]=122; CMAP[2][72]=255;
    CMAP[0][73]=218; CMAP[1][73]=120; CMAP[2][73]=255;
    CMAP[0][74]=221; CMAP[1][74]=119; CMAP[2][74]=255;
    CMAP[0][75]=223; CMAP[1][75]=117; CMAP[2][75]=255;
    CMAP[0][76]=226; CMAP[1][76]=116; CMAP[2][76]=255;
    CMAP[0][77]=229; CMAP[1][77]=114; CMAP[2][77]=255;
    CMAP[0][78]=233; CMAP[1][78]=113; CMAP[2][78]=255;
    CMAP[0][79]=236; CMAP[1][79]=111; CMAP[2][79]=255;
    CMAP[0][80]=239; CMAP[1][80]=110; CMAP[2][80]=255;
    CMAP[0][81]=242; CMAP[1][81]=108; CMAP[2][81]=255;
    CMAP[0][82]=246; CMAP[1][82]=106; CMAP[2][82]=255;
    CMAP[0][83]=249; CMAP[1][83]=105; CMAP[2][83]=255;
    CMAP[0][84]=253; CMAP[1][84]=104; CMAP[2][84]=255;
    CMAP[0][85]=255; CMAP[1][85]=102; CMAP[2][85]=254;
    CMAP[0][86]=255; CMAP[1][86]=104; CMAP[2][86]=250;
    CMAP[0][87]=255; CMAP[1][87]=105; CMAP[2][87]=247;
    CMAP[0][88]=255; CMAP[1][88]=106; CMAP[2][88]=243;
    CMAP[0][89]=255; CMAP[1][89]=108; CMAP[2][89]=240;
    CMAP[0][90]=255; CMAP[1][90]=110; CMAP[2][90]=237;
    CMAP[0][91]=255; CMAP[1][91]=111; CMAP[2][91]=233;
    CMAP[0][92]=255; CMAP[1][92]=113; CMAP[2][92]=230;
    CMAP[0][93]=255; CMAP[1][93]=114; CMAP[2][93]=227;
    CMAP[0][94]=255; CMAP[1][94]=116; CMAP[2][94]=224;
    CMAP[0][95]=255; CMAP[1][95]=117; CMAP[2][95]=221;
    CMAP[0][96]=255; CMAP[1][96]=119; CMAP[2][96]=218;
    CMAP[0][97]=255; CMAP[1][97]=120; CMAP[2][97]=216;
    CMAP[0][98]=255; CMAP[1][98]=122; CMAP[2][98]=213;
    CMAP[0][99]=255; CMAP[1][99]=123; CMAP[2][99]=210;
    CMAP[0][100]=255; CMAP[1][100]=125; CMAP[2][100]=208;
    CMAP[0][101]=255; CMAP[1][101]=126; CMAP[2][101]=205;
    CMAP[0][102]=255; CMAP[1][102]=128; CMAP[2][102]=203;
    CMAP[0][103]=255; CMAP[1][103]=129; CMAP[2][103]=200;
    CMAP[0][104]=255; CMAP[1][104]=131; CMAP[2][104]=198;
    CMAP[0][105]=255; CMAP[1][105]=132; CMAP[2][105]=196;
    CMAP[0][106]=255; CMAP[1][106]=134; CMAP[2][106]=194;
    CMAP[0][107]=255; CMAP[1][107]=135; CMAP[2][107]=192;
    CMAP[0][108]=255; CMAP[1][108]=137; CMAP[2][108]=190;
    CMAP[0][109]=255; CMAP[1][109]=138; CMAP[2][109]=188;
    CMAP[0][110]=255; CMAP[1][110]=140; CMAP[2][110]=186;
    CMAP[0][111]=255; CMAP[1][111]=141; CMAP[2][111]=184;
    CMAP[0][112]=255; CMAP[1][112]=143; CMAP[2][112]=182;
    CMAP[0][113]=255; CMAP[1][113]=144; CMAP[2][113]=181;
    CMAP[0][114]=255; CMAP[1][114]=146; CMAP[2][114]=179;
    CMAP[0][115]=255; CMAP[1][115]=147; CMAP[2][115]=178;
    CMAP[0][116]=255; CMAP[1][116]=149; CMAP[2][116]=176;
    CMAP[0][117]=255; CMAP[1][117]=150; CMAP[2][117]=175;
    CMAP[0][118]=255; CMAP[1][118]=152; CMAP[2][118]=174;
    CMAP[0][119]=255; CMAP[1][119]=153; CMAP[2][119]=172;
    CMAP[0][120]=255; CMAP[1][120]=155; CMAP[2][120]=171;
    CMAP[0][121]=255; CMAP[1][121]=156; CMAP[2][121]=170;
    CMAP[0][122]=255; CMAP[1][122]=158; CMAP[2][122]=169;
    CMAP[0][123]=255; CMAP[1][123]=159; CMAP[2][123]=168;
    CMAP[0][124]=255; CMAP[1][124]=161; CMAP[2][124]=167;
    CMAP[0][125]=255; CMAP[1][125]=162; CMAP[2][125]=166;
    CMAP[0][126]=255; CMAP[1][126]=164; CMAP[2][126]=166;
    CMAP[0][127]=255; CMAP[1][127]=165; CMAP[2][127]=165;
    CMAP[0][128]=255; CMAP[1][128]=165; CMAP[2][128]=165;
    CMAP[0][129]=255; CMAP[1][129]=166; CMAP[2][129]=164;
    CMAP[0][130]=255; CMAP[1][130]=166; CMAP[2][130]=162;
    CMAP[0][131]=255; CMAP[1][131]=167; CMAP[2][131]=161;
    CMAP[0][132]=255; CMAP[1][132]=168; CMAP[2][132]=159;
    CMAP[0][133]=255; CMAP[1][133]=169; CMAP[2][133]=157;
    CMAP[0][134]=255; CMAP[1][134]=170; CMAP[2][134]=156;
    CMAP[0][135]=255; CMAP[1][135]=171; CMAP[2][135]=155;
    CMAP[0][136]=255; CMAP[1][136]=172; CMAP[2][136]=153;
    CMAP[0][137]=255; CMAP[1][137]=174; CMAP[2][137]=152;
    CMAP[0][138]=255; CMAP[1][138]=175; CMAP[2][138]=150;
    CMAP[0][139]=255; CMAP[1][139]=176; CMAP[2][139]=149;
    CMAP[0][140]=255; CMAP[1][140]=178; CMAP[2][140]=147;
    CMAP[0][141]=255; CMAP[1][141]=179; CMAP[2][141]=145;
    CMAP[0][142]=255; CMAP[1][142]=181; CMAP[2][142]=144;
    CMAP[0][143]=255; CMAP[1][143]=182; CMAP[2][143]=142;
    CMAP[0][144]=255; CMAP[1][144]=184; CMAP[2][144]=141;
    CMAP[0][145]=255; CMAP[1][145]=186; CMAP[2][145]=139;
    CMAP[0][146]=255; CMAP[1][146]=188; CMAP[2][146]=138;
    CMAP[0][147]=255; CMAP[1][147]=190; CMAP[2][147]=137;
    CMAP[0][148]=255; CMAP[1][148]=192; CMAP[2][148]=135;
    CMAP[0][149]=255; CMAP[1][149]=194; CMAP[2][149]=134;
    CMAP[0][150]=255; CMAP[1][150]=196; CMAP[2][150]=132;
    CMAP[0][151]=255; CMAP[1][151]=198; CMAP[2][151]=131;
    CMAP[0][152]=255; CMAP[1][152]=200; CMAP[2][152]=129;
    CMAP[0][153]=255; CMAP[1][153]=203; CMAP[2][153]=128;
    CMAP[0][154]=255; CMAP[1][154]=205; CMAP[2][154]=126;
    CMAP[0][155]=255; CMAP[1][155]=208; CMAP[2][155]=125;
    CMAP[0][156]=255; CMAP[1][156]=210; CMAP[2][156]=123;
    CMAP[0][157]=255; CMAP[1][157]=213; CMAP[2][157]=122;
    CMAP[0][158]=255; CMAP[1][158]=216; CMAP[2][158]=120;
    CMAP[0][159]=255; CMAP[1][159]=218; CMAP[2][159]=119;
    CMAP[0][160]=255; CMAP[1][160]=221; CMAP[2][160]=117;
    CMAP[0][161]=255; CMAP[1][161]=224; CMAP[2][161]=116;
    CMAP[0][162]=255; CMAP[1][162]=227; CMAP[2][162]=114;
    CMAP[0][163]=255; CMAP[1][163]=230; CMAP[2][163]=113;
    CMAP[0][164]=255; CMAP[1][164]=233; CMAP[2][164]=111;
    CMAP[0][165]=255; CMAP[1][165]=237; CMAP[2][165]=109;
    CMAP[0][166]=255; CMAP[1][166]=240; CMAP[2][166]=108;
    CMAP[0][167]=255; CMAP[1][167]=243; CMAP[2][167]=106;
    CMAP[0][168]=255; CMAP[1][168]=247; CMAP[2][168]=105;
    CMAP[0][169]=255; CMAP[1][169]=250; CMAP[2][169]=104;
    CMAP[0][170]=255; CMAP[1][170]=254; CMAP[2][170]=102;
    CMAP[0][171]=253; CMAP[1][171]=255; CMAP[2][171]=104;
    CMAP[0][172]=249; CMAP[1][172]=255; CMAP[2][172]=105;
    CMAP[0][173]=246; CMAP[1][173]=255; CMAP[2][173]=106;
    CMAP[0][174]=242; CMAP[1][174]=255; CMAP[2][174]=108;
    CMAP[0][175]=239; CMAP[1][175]=255; CMAP[2][175]=109;
    CMAP[0][176]=236; CMAP[1][176]=255; CMAP[2][176]=111;
    CMAP[0][177]=233; CMAP[1][177]=255; CMAP[2][177]=112;
    CMAP[0][178]=229; CMAP[1][178]=255; CMAP[2][178]=114;
    CMAP[0][179]=226; CMAP[1][179]=255; CMAP[2][179]=115;
    CMAP[0][180]=223; CMAP[1][180]=255; CMAP[2][180]=117;
    CMAP[0][181]=221; CMAP[1][181]=255; CMAP[2][181]=119;
    CMAP[0][182]=218; CMAP[1][182]=255; CMAP[2][182]=120;
    CMAP[0][183]=215; CMAP[1][183]=255; CMAP[2][183]=122;
    CMAP[0][184]=212; CMAP[1][184]=255; CMAP[2][184]=123;
    CMAP[0][185]=210; CMAP[1][185]=255; CMAP[2][185]=124;
    CMAP[0][186]=207; CMAP[1][186]=255; CMAP[2][186]=126;
    CMAP[0][187]=205; CMAP[1][187]=255; CMAP[2][187]=127;
    CMAP[0][188]=202; CMAP[1][188]=255; CMAP[2][188]=129;
    CMAP[0][189]=200; CMAP[1][189]=255; CMAP[2][189]=131;
    CMAP[0][190]=198; CMAP[1][190]=255; CMAP[2][190]=132;
    CMAP[0][191]=196; CMAP[1][191]=255; CMAP[2][191]=134;
    CMAP[0][192]=194; CMAP[1][192]=255; CMAP[2][192]=135;
    CMAP[0][193]=192; CMAP[1][193]=255; CMAP[2][193]=137;
    CMAP[0][194]=190; CMAP[1][194]=255; CMAP[2][194]=138;
    CMAP[0][195]=188; CMAP[1][195]=255; CMAP[2][195]=139;
    CMAP[0][196]=186; CMAP[1][196]=255; CMAP[2][196]=141;
    CMAP[0][197]=184; CMAP[1][197]=255; CMAP[2][197]=143;
    CMAP[0][198]=182; CMAP[1][198]=255; CMAP[2][198]=144;
    CMAP[0][199]=181; CMAP[1][199]=255; CMAP[2][199]=146;
    CMAP[0][200]=179; CMAP[1][200]=255; CMAP[2][200]=147;
    CMAP[0][201]=178; CMAP[1][201]=255; CMAP[2][201]=149;
    CMAP[0][202]=176; CMAP[1][202]=255; CMAP[2][202]=150;
    CMAP[0][203]=175; CMAP[1][203]=255; CMAP[2][203]=152;
    CMAP[0][204]=174; CMAP[1][204]=255; CMAP[2][204]=153;
    CMAP[0][205]=173; CMAP[1][205]=255; CMAP[2][205]=155;
    CMAP[0][206]=172; CMAP[1][206]=255; CMAP[2][206]=156;
    CMAP[0][207]=171; CMAP[1][207]=255; CMAP[2][207]=158;
    CMAP[0][208]=170; CMAP[1][208]=255; CMAP[2][208]=159;
    CMAP[0][209]=169; CMAP[1][209]=255; CMAP[2][209]=161;
    CMAP[0][210]=168; CMAP[1][210]=255; CMAP[2][210]=162;
    CMAP[0][211]=167; CMAP[1][211]=255; CMAP[2][211]=164;
    CMAP[0][212]=166; CMAP[1][212]=255; CMAP[2][212]=165;
    CMAP[0][213]=165; CMAP[1][213]=255; CMAP[2][213]=166;
    CMAP[0][214]=164; CMAP[1][214]=255; CMAP[2][214]=166;
    CMAP[0][215]=162; CMAP[1][215]=255; CMAP[2][215]=167;
    CMAP[0][216]=161; CMAP[1][216]=255; CMAP[2][216]=168;
    CMAP[0][217]=159; CMAP[1][217]=255; CMAP[2][217]=169;
    CMAP[0][218]=158; CMAP[1][218]=255; CMAP[2][218]=170;
    CMAP[0][219]=156; CMAP[1][219]=255; CMAP[2][219]=171;
    CMAP[0][220]=155; CMAP[1][220]=255; CMAP[2][220]=172;
    CMAP[0][221]=153; CMAP[1][221]=255; CMAP[2][221]=173;
    CMAP[0][222]=152; CMAP[1][222]=255; CMAP[2][222]=174;
    CMAP[0][223]=150; CMAP[1][223]=255; CMAP[2][223]=176;
    CMAP[0][224]=149; CMAP[1][224]=255; CMAP[2][224]=177;
    CMAP[0][225]=147; CMAP[1][225]=255; CMAP[2][225]=178;
    CMAP[0][226]=146; CMAP[1][226]=255; CMAP[2][226]=180;
    CMAP[0][227]=144; CMAP[1][227]=255; CMAP[2][227]=182;
    CMAP[0][228]=142; CMAP[1][228]=255; CMAP[2][228]=183;
    CMAP[0][229]=141; CMAP[1][229]=255; CMAP[2][229]=185;
    CMAP[0][230]=139; CMAP[1][230]=255; CMAP[2][230]=187;
    CMAP[0][231]=138; CMAP[1][231]=255; CMAP[2][231]=189;
    CMAP[0][232]=137; CMAP[1][232]=255; CMAP[2][232]=191;
    CMAP[0][233]=135; CMAP[1][233]=255; CMAP[2][233]=193;
    CMAP[0][234]=134; CMAP[1][234]=255; CMAP[2][234]=195;
    CMAP[0][235]=132; CMAP[1][235]=255; CMAP[2][235]=197;
    CMAP[0][236]=131; CMAP[1][236]=255; CMAP[2][236]=199;
    CMAP[0][237]=129; CMAP[1][237]=255; CMAP[2][237]=201;
    CMAP[0][238]=128; CMAP[1][238]=255; CMAP[2][238]=204;
    CMAP[0][239]=126; CMAP[1][239]=255; CMAP[2][239]=206;
    CMAP[0][240]=125; CMAP[1][240]=255; CMAP[2][240]=209;
    CMAP[0][241]=123; CMAP[1][241]=255; CMAP[2][241]=211;
    CMAP[0][242]=122; CMAP[1][242]=255; CMAP[2][242]=214;
    CMAP[0][243]=120; CMAP[1][243]=255; CMAP[2][243]=217;
    CMAP[0][244]=119; CMAP[1][244]=255; CMAP[2][244]=220;
    CMAP[0][245]=117; CMAP[1][245]=255; CMAP[2][245]=222;
    CMAP[0][246]=116; CMAP[1][246]=255; CMAP[2][246]=225;
    CMAP[0][247]=114; CMAP[1][247]=255; CMAP[2][247]=228;
    CMAP[0][248]=113; CMAP[1][248]=255; CMAP[2][248]=231;
    CMAP[0][249]=111; CMAP[1][249]=255; CMAP[2][249]=235;
    CMAP[0][250]=110; CMAP[1][250]=255; CMAP[2][250]=238;
    CMAP[0][251]=108; CMAP[1][251]=255; CMAP[2][251]=241;
    CMAP[0][252]=106; CMAP[1][252]=255; CMAP[2][252]=244;
    CMAP[0][253]=105; CMAP[1][253]=255; CMAP[2][253]=248;
    CMAP[0][254]=104; CMAP[1][254]=255; CMAP[2][254]=251;
    CMAP[0][255]=102; CMAP[1][255]=255; CMAP[2][255]=255;
    }

  else // load from ascii file, generated by, e.g., makecpt
       //  line1: r g b  (e.g. 255 0 255)
       //  line2: r g b  (e.g. 255 0 240)
       //  ...    ...
    {
    //ifstream cmap(input.cmap, ios::in | ios::nocreate);
    ifstream cmap(input.cmap, ios::in);
    if (!cmap) 
      {
      cerr << "cpxfiddle: ERROR: Could not open colormap file: " << input.cmap << "\n";
      exit(1);
      }
    char dummyline[128*2];
    int value;
    for (int ii=0; ii<maplength; ++ii)
      {
      for (int color=0; color<3; ++color)
        {
        cmap >> value;
        CMAP[color][ii] = uchar(value);
        }
      cmap.getline(dummyline,128,'\n');
      }
    cmap.close();
    }

  // ______ Display colormap for testing ______
  if (input.verbose==true)
    {
    #define VERBOSE
    #ifdef VERBOSE
    cerr << "Color table for SUNraster (rgb):\n";
    for (int ii=0; ii<maplength; ++ii)
      cerr << int(CMAP[0][ii]) << " "
           << int(CMAP[1][ii]) << " "
           << int(CMAP[2][ii]) << "\n";
    #endif
    }
  } // END makecmap


/****************************************************************
 * makecmapmixed                                                *
 * input option -c: only for sunraster output:                  *
 * generate a cmap from a input file, or from a name "gray", ?  *
 * input file consist of 256 ascii entries r g b                *
 * This function creates a mixed colormap, so that data that is *
 * d = numcolorsphase*mag + phase                               *
 * is displayed as overlay phase over mag.                      *
 * This function uses makecmap to generate a map of length 256  *
 * that is composed of 16 repeated colormap for phase with      *
 * different brightness levels.                                 *
 #%// BK 04-Apr-2001                                            *
 ****************************************************************/
void makecmapmixed(
        unsigned char CMAP[3][256],
        const commandlineinput &input)
  {
  // do not change!
  const int maplength      = 256;                       // always
  const int numcolorsphase = 16;
  const int numcolorsmag   = 16;
  // (numcolorsphase*numcolorsmag==maplength)

  // ______ First create a cmap for phase with 16 colors (subsample) ______
  makecmap(CMAP,input);
  unsigned char CMAP_PHA16[3][16];
  int i,rgb;
  for (i=0; i<numcolorsphase; ++i)
    for (rgb=0; rgb<3; ++rgb)
      CMAP_PHA16[rgb][i] = CMAP[rgb][int(numcolorsphase/2)+i*numcolorsphase];

  if (input.verbose==true)
    {
    cerr << "\nColromap for phase (16):\n";
    for (int ii=0; ii<16; ++ii)
      cerr << int(CMAP_PHA16[0][ii]) << " "
           << int(CMAP_PHA16[1][ii]) << " "
           << int(CMAP_PHA16[2][ii]) << "\n";
    cerr << "\n";
    }

  for (i=0; i<256; ++i)
    {
    int magnitudelevel = int(i/numcolorsphase);         // e[0:15]
    float level        = float(magnitudelevel)/float(numcolorsmag-1);// e[0,1]
    // ______ And add the gray colormap for the magnitude ______
    for (rgb=0; rgb<3; ++rgb)
      CMAP[rgb][i] = uchar(level*CMAP_PHA16[rgb][i%numcolorsphase]);
    }
  // ______ Display colormap for testing ______
  if (input.verbose==true)
    {
    #define VERBOSE
    #ifdef VERBOSE
    cerr << "\nWriting color table for SUNraster (mixed mag/phase):\n";
    for (int ii=0; ii<maplength; ++ii)
      cerr << int(CMAP[0][ii]) << " "
           << int(CMAP[1][ii]) << " "
           << int(CMAP[2][ii]) << "\n";
    #endif
    }
  } // END makecmapmixed


/****************************************************************
 * handleinput                                                  *
 * return input options/arguments                               *
 #%// BK 11-May-2000                                            *
 #%// MA 19-Oct-2008                                            *
 ****************************************************************/
bool handleinput(
        int argc,
        char* argv[],
        commandlineinput &input)
  {
  // at least      "prog -h"
  // at least      "prog -w 123 file" or "prog -w123 file"
  // leaves argc = " 1    2  3   4  "    " 1     2     3 "
  // if (argc < 3) {cerr << "argc: " << argc << endl; return false;} // speed
  //if (argc==1) {cerr << "argc: " << argc << endl; return false;} // speed
  if (argc==1) return false; // speed

  bool allmandatory = true;
  extern char *optarg;
  //extern int optind,opterr,optopt;
  extern int optind;
  int c;

  // ______ Set defaults ______
  input.linelength  = 0;
  input.iformat     = FORMATCR4;        // default
  input.oformat     = ASCII;            // default readable output to screen 
  input.dooutput    = MAGNITUDE;        // default magnitude
  input.exponent    = 1.;               // default
  input.scale       = 1.;               // default
  input.firstline   = 1;                // default
  input.lastline    = -1;               // compute last line
  input.firstpixel  = 1;                // default
  input.lastpixel   = -1;               // compute last pixel
  input.sublines    = 1;
  input.subpixels   = 1;
  input.multilookL  = 1;
  input.multilookP  = 1;
  input.headerlength = 0;               // default
  input.mirrorX     = NOMIRROR;         // default
  input.mirrorY     = NOMIRROR;         // default
  strcpy(input.cmap,"default");         // default (mag. gray/ otherwise hsv)
  //strcpy(input.cmap,"gray");          // default (mag.)
  input.numlines    = 0;
  input.dontohx     = false;            // default
  input.dohtonx     = false;            // default
  input.verbose     = false;            // default
  input.scalebar    = false;            // default
  input.realinput   = false;            // default cpx files
  input.ignorenan   = false;            // don't ignore NaN values
  input.fullstat    = false;            // use fullstat to get min max for big images, for scaling
  input.rangemin    = NAN;                // color range min 
  input.rangemax    = NAN;                // color range max

/* For long options that have no equivalent short option, use a
   non-character as a pseudo short option, starting with CHAR_MAX + 1.  */
enum
{
  IGNORENAN_OPTION = CHAR_MAX + 1,
  FULLSTAT_OPTION
  // TEMP_OPTION
};

static struct option const long_options[] =
{
  // {"all", no_argument, NULL, 'a'},
  // {"width", required_argument, NULL, 'w'},
  {"ignorenan", 0, NULL, IGNORENAN_OPTION},
  {"fullstat", 0, NULL, FULLSTAT_OPTION},
//  {GETOPT_HELP_OPTION_DECL},
//  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

  // ______ Get commandline ______
  char OPTSTRING[40] = "w:e:f:l:L:m:o:p:P:q:s:S:M:c:B:H:r:bhV";
  // while ((c = getopt(argc, argv, OPTSTRING)) != -1)
    string strarg = argv[1];
    for(int idx=2;idx<argc;idx++)
    {
     strarg = strarg + " " +  string(argv[idx]);
    } 
    cerr << "# args: " << strarg << "\n"; // for debuging arguments [PD+MA]

  while (1)
    {
    int option_index= -1 ;
    c = getopt_long(argc, argv, OPTSTRING, long_options, &option_index);

    if (c == -1)
        break;

    switch (c)
      {
      case 'w':
        input.linelength  = atoi(optarg);
        break;
      case 'e':
        input.exponent    = atof(optarg);
        break;
      case 's':
        input.scale       = atof(optarg);
        break;
      case 'l':
        input.firstline   = atoi(optarg);
        break;
      case 'L':
        input.lastline    = atoi(optarg);
        break;
      case 'p':
        input.firstpixel  = atoi(optarg);
        break;
      case 'P':
        input.lastpixel   = atoi(optarg);
        break;
      case 'S':         // subsample -Mp/l == -Mx/y == -M1/5
        {
        int i=0;
        char subP[8];
        while (optarg[i] != '/' && optarg[i] != '\0')
          {
          subP[i] = optarg[i];
          ++i;
          }
        subP[i]='\0';
        char subL[8];
        int j=0;
        for (j=0; j<sizeof(optarg)-i; ++j)
          subL[j] = optarg[j+i+1];
        subL[j]='\0';
        input.sublines  = atoi(subL);
        input.subpixels = atoi(subP);
        if (input.sublines  < 1) input.sublines = 1;
        if (input.subpixels < 1) input.subpixels = 1;
        break;
        }
      case 'q':
        if (!strcmp(optarg,"mag"))
          input.dooutput = MAGNITUDE;
        else if (!strcmp(optarg,"phase"))
          input.dooutput = PHASE;
        else if (!strcmp(optarg,"real"))
          input.dooutput = REALPART;
        else if (!strcmp(optarg,"imag"))
          input.dooutput = IMAGPART;
        else if (!strcmp(optarg,"normal"))
          input.dooutput = NORMAL;
        else if (!strcmp(optarg,"mixed"))
          input.dooutput = MIXED;
        else
          {
          cerr << argv[0] << ": ERROR: -q mag|phase|real|imag|normal|mixed: ~= "
               << optarg << endl;
          synopsis(argv[0]);
          }
        break;
      case 'f':
        if (!strcmp(optarg,"cc1"))
          input.iformat = FORMATCC1;
        else if (!strcmp(optarg,"cuc1"))
          input.iformat = FORMATCUC1;
        else if (!strcmp(optarg,"ci2"))
          input.iformat = FORMATCI2;
        else if (!strcmp(optarg,"ci4"))
          input.iformat = FORMATCI4;
        else if (!strcmp(optarg,"cr4"))
          input.iformat = FORMATCR4;
        else if (!strcmp(optarg,"cr8"))
          input.iformat = FORMATCR8;
        // add char input, though it is not complex...
        else if (!strcmp(optarg,"c1"))
          input.iformat = FORMATUC1;
        // add i2 input, though it is not complex...
        else if (!strcmp(optarg,"i2"))
          input.iformat = FORMATI2;
        // add real4 input, though it is not complex...
        else if (!strcmp(optarg,"r4"))
          input.iformat = FORMATR4;
        // add real8 input, though it is not complex...
        else if (!strcmp(optarg,"r8"))
          input.iformat = FORMATR8;
        else
          {
          cerr << argv[0] << ": ERROR: -f cc1|cuc1|ci2|ci4|cr4|cr8|c1|i2|r4|r8: ~= "
               << optarg << endl;
          synopsis(argv[0]);
          }
        break;
      case 'o':
        if (!strcmp(optarg,"ascii"))
          input.oformat = ASCII;
        else if (!strcmp(optarg,"float"))
          input.oformat = FORMATR4;
        else if (!strcmp(optarg,"short"))
          input.oformat = FORMATI2;
        else if (!strcmp(optarg,"uchar"))
          input.oformat = FORMATUC1;
        else if (!strcmp(optarg,"sunraster") || !strcmp(optarg,"ras") )  // [MA]  ras extention
          input.oformat = SUNRASTER;
        else
          {
          cerr << argv[0] << ": ERROR: -o ascii|float|short|uchar|sunraster ~= "
          << optarg << endl;
          synopsis(argv[0]);
          }
        break;
      case 'm':
        if (!strcmp(optarg,"XY") || !strcmp(optarg,"YX"))
          {
          input.mirrorX = DOMIRROR;
          input.mirrorY = DOMIRROR;
          }
        else if (!strcmp(optarg,"X"))
          input.mirrorX = DOMIRROR;
        else if (!strcmp(optarg,"Y"))
          input.mirrorY = DOMIRROR;
        else
          {
          cerr << argv[0] << ": ERROR: -m XY|YX|X|Y: " << optarg << endl;
          synopsis(argv[0]);
          }
        break;
      case 'M':         // multilook -Mp/l == -Mx/y == -M1/5
        {
        int i=0;
        char mlP[8];
        while (optarg[i] != '/' && optarg[i] != '\0')
          {
          mlP[i] = optarg[i];
          ++i;
          }
        mlP[i]='\0';
        int j=0;
        char mlL[8];
        for (j=0; j<strlen(optarg)-i; ++j)
          mlL[j] = optarg[j+i+1];
        mlL[j]='\0';
        input.multilookL = atoi(mlL);
        input.multilookP = atoi(mlP);
        if (input.multilookL < 1) input.multilookL = 1;
        if (input.multilookP < 1) input.multilookP = 1;
        break;
        }
      case 'c':
        strcpy(input.cmap,optarg);
        break;
      case 'B':
        if (!strcmp(optarg,"s"))
          input.dontohx = true;
        else if (!strcmp(optarg,"b"))
          input.dohtonx = true;
        else
          cerr << argv[0] << ": option to -B not recognized, continuing w/o swapping.\n";
        break;
      case 'H':
        input.headerlength = atoi(optarg);
        break;
      case 'r': // [MA] if not used rmin/rmax is NaN
       {        // TODO warning if phase is used: "no use during phase calculation"
        sscanf(optarg, "%lf%*c%lf", &input.rangemin,&input.rangemax); // %*c ignore / char in the middle, see also strtod
        cerr.precision(10);
        cerr << endl << argv[0] << ": INFO   : Using range min: " << input.rangemin << " & max: " << input.rangemax << " instead.\n"; // no samplin to calc min,max,mean
        if (input.rangemin > input.rangemax) 
           {
            cerr << argv[0] << ": ERROR  : range min: " << input.rangemin << " is bigger than max: " << input.rangemax << "\n";
            exit(1);
           } 
        else if ( input.rangemin == input.rangemax || isnan(input.rangemin) || isnan(input.rangemax) )
           {
            cerr << argv[0] << ": ERROR  : Check range min: " << input.rangemin << " & max: " << input.rangemax << " !\n";
            exit(1);
           }
        break;
        }
      case 'b':
        input.scalebar = true;
        break;
      case 'V':
        input.verbose = true;
        break;
      case IGNORENAN_OPTION:
        input.ignorenan = true;
        cerr << "debug: ignorenan" << IGNORENAN_OPTION << "\n";
        // clean up extern may be neccessary
        break;
      case FULLSTAT_OPTION:
        input.fullstat = true;
        cerr << "debug: fullstat" << IGNORENAN_OPTION << "\n";
        // clean up extern may be neccessary
        break;
      case 'h':
        usage(argv[0]);
        break;
      case '?':
        synopsis(argv[0]);
      } // switch
    } // commandlineoptions

  // ______ Check input ______
  if (input.linelength == 0)
    {
    cerr << argv[0] << ": ERROR: No width specified.\n";
    return false;
    }

  // filename: last argument
  // cerr << "OPTARG: " << argv[argc-1] << endl;
  // cerr << "OPTARG: " << argv[optind] << endl;
  if (argv[optind]=='\0')
    {
    cerr << argv[0] << ": ERROR: No input file specified.\n";
    return false;
    }
  strcpy(input.ifile,argv[optind]);

  // ______ Check for: "cpxfiddle --help" ______
  if (!strcmp(input.ifile,"help"))
    usage(argv[0]);

  // ____ new option for real4 input, not intented usage _____
  // ____ outformat must be normal for this _____
  // BK 13-Apr-2003
  if (input.iformat == FORMATR8  || 
      input.iformat == FORMATR4  || 
      input.iformat == FORMATUC1 || 
      input.iformat == FORMATI2)
    input.realinput = true;
  if (input.realinput)
    {
    if (input.dooutput != NORMAL)
      {
      cerr << argv[0] << ": ERROR: for real input formats, use -q normal\n";
      usage(argv[0]);
      }
    }

  // ______ Compute numlines ______
  // rather either numlines or width and compute the other.
  //ifstream infile(input.ifile, ios::in | ios::nocreate);
  ifstream infile(input.ifile, ios::in);
  if (!infile)
    {
    cerr << argv[0] << ": ERROR: inputfile: " << input.ifile << "not found.\n";
    synopsis(argv[0]);
    }
  infile.seekg(0,ios::end);                    // filepointer at end
  //const int fsize = int(infile.tellg()) - input.headerlength;// account for e.g., 32 byte header
  //MA file offset > 4 GB , 200803
  const streamoff fsize = infile.tellg() - streamoff(input.headerlength); // account for e.g., 32 byte header
  // here 
  infile.close();
  cerr << "# fsize: " << fsize << "\n";
  input.bytesperelement = input.iformat%10;     // real/imag is read seperately
  input.bytesperpixel   = 2*input.bytesperelement;// per complex pixel
  // ____ adapt bytes per pixel for real input.  -q normal must be ____
  // BK 13-Apr-2003
  if (input.realinput)
    {
    cerr << argv[0] << ": WARNING: using real input, width in non-complex pixels\n";
    input.bytesperpixel   = input.bytesperelement;// hopefully this is we have to do
    }
  input.numlines        = fsize / input.bytesperpixel / input.linelength;
  input.bytesperline    = input.bytesperpixel * input.linelength;
  // ______ Set other defaults ______
  if (input.lastline==-1)  input.lastline  = input.numlines;
  if (input.lastpixel==-1) input.lastpixel = input.linelength;

  // ______ Check wrong... ______
  if (input.lastline>input.numlines)
    {
    cerr << "#" << argv[0]
      << ": WARNING: -L: lastline>numlines; continuing with numlines of file.\n"; 
    input.lastline = input.numlines;
    }
  if (input.lastpixel>input.linelength)
    {
    cerr << "#" << argv[0] 
      << ": WARNING: -P: lastpixel>numpixels; continuing with width of file.\n";
    input.lastpixel = input.linelength;
    }
  // if (input.numlines*input.linelength*input.bytesperpixel != fsize) // MA error check fix
  if (streamoff(input.numlines)*streamoff(input.linelength)*streamoff(input.bytesperpixel) != fsize)
    {
    cerr << argv[0] << ": ERROR: -w does not seem ok (checked with numlines*bytes)\n";
    cerr << argv[0] << ": ERROR: but continuing, see it yourself...\n\a\a";
    //synopsis(argv[0]);
    }
  if (input.lastline<input.firstline)
    {
     cerr << argv[0] << ": ERROR: -l -L: lastline [" << input.lastline << "] < firstline [" << input.firstline << "]\n";
     synopsis(argv[0]);
    }
  if (input.lastpixel<input.firstpixel)
    {
    cerr << argv[0] << ": ERROR: -p -P: lastpixel<firstpixel\n";
    synopsis(argv[0]);
    }

  // ______ Set default cmap if not specified ______
  if (input.oformat==SUNRASTER)
    if (!strcmp(input.cmap,"default"))
      (input.dooutput==MAGNITUDE) ? 
        strcpy(input.cmap,"gray") : strcpy(input.cmap,"hsv");

  if (input.dooutput==PHASE)
    {
    if (input.exponent != 1.) cerr << "#WARNING: -e: exponent for phase output?\n";
    if (input.scale    != 1.) cerr << "#WARNING: -s: scale for phase output?\n";
    if (input.scale   == -1.) cerr << "#-s -1: ok, flipping phase.\n";
    if (input.oformat==SUNRASTER)
      if (!strcmp(input.cmap,"gray"))
        cerr << "Tip -c changes colormap, e.g. -c hot\n";
    }
  if (input.dooutput==MIXED)
    {
    if (input.oformat!=SUNRASTER && input.oformat!=FORMATUC1)
      {
      cerr << "ERROR: -q MIXED only for -o sunraster/uc1\n";
      synopsis(argv[0]);
      }
    else
      {
      if (!strcmp(input.cmap,"gray"))
        {
        cerr << "#WARNING: -q mixed: changing phase colormap from gray to hot (-c)\n";
        strcpy(input.cmap,"hot");
        }
      }
    }
  if (input.scalebar==true && input.oformat!=SUNRASTER)
    {
    cerr << "ERROR: -b only for sunraster phase/mixed\n";
    synopsis(argv[0]);
    }
  if (input.multilookL!=1 && input.sublines!=1)
    {
    cerr << argv[0] << ": ERROR: -M -Y: no subsampling and multilooking allowed.\n";
    synopsis(argv[0]);
    }
  if (input.multilookP!=1 && input.subpixels!=1)
    {
    cerr << argv[0] 
         << ": ERROR: -M -X: no subsampling and multilooking allowed.\n";
    synopsis(argv[0]);
    }


  // ______ give tip ______
  if (input.oformat==FORMATUC1 && input.dooutput==PHASE)
    cerr << argv[0] 
    << ": TIP: display with (ImageMagick, FILE is redirected file):\ndisplay -size "
    << input.linelength << "x" << input.numlines << " gray:FILE\n";
    

  // ____ be verbose if required ______
  if (input.verbose==true)
    {
    cerr                   << argv[0]               << " variables:\n" <<
    " linelength:        " << input.linelength      << "\n" <<
    " ifile:             " << input.ifile           << "\n" <<
    " iformat:           " << input.iformat         << "\n" <<
    " exponent:          " << input.exponent        << "\n" <<
    " scale:             " << input.scale           << "\n" <<
    " firstline:         " << input.firstline       << "\n" <<
    " lastline:          " << input.lastline        << "\n" <<
    " firstpixel:        " << input.firstpixel      << "\n" <<
    " lastpixel:         " << input.lastpixel       << "\n" <<
    " sublines:          " << input.sublines        << "\n" <<
    " subpixels:         " << input.subpixels       << "\n" <<
    " dooutput:          " << input.dooutput        << "\n" <<
    " oformat:           " << input.oformat         << "\n" <<
    " dohtonx:           " << input.dohtonx         << "\n" <<
    " dontohx:           " << input.dontohx         << "\n" <<
    " scalebar:          " << input.scalebar        << "\n" <<
    " mirrorX:           " << input.mirrorX         << "\n" <<
    " mirrorY:           " << input.mirrorY         << "\n" <<
    " multilookX:        " << input.multilookP      << "\n" <<
    " multilookY:        " << input.multilookL      << "\n" <<
    " cmap:              " << input.cmap            << "\n" <<
    " skip header bytes: " << input.headerlength    << "\n" <<
    " numlines:          " << input.numlines        << "\n" <<
    " bytesperelement:   " << input.bytesperelement << "\n" <<
    " bytesperpixel:     " << input.bytesperpixel   << "\n" <<
    " bytesperline:      " << input.bytesperline    << "\n" <<
    " range:             " << "[ " << input.rangemin << " " << input.rangemax << " ]" << endl;
    }

  // ______ Set variables used in 'functie' to multilook ______
  if (input.multilookL != 1)
    input.sublines = input.multilookL;                  // use this to set filepointer
  if (input.multilookP != 1)
    input.subpixels = input.multilookP;                 // use this to set filepointer

  return allmandatory;
  } // END handleinput


/****************************************************************
 * shortexpl                                                    *
 * for synopsis and usage                                       *
 #%// BK 11-May-2000                                            *
 ****************************************************************/
void shortexpl()
  {
  cerr << "\n   Dump content of complex binary file to stdout,"
       << "\n    either: as is, magnitude, phase, real or imaginary part."
       << "\n   Input files can be almost any complex file"
       << "\n    though not (yet) band interleaved."
       << "\n   Output can be manipulated by:"
       << "\n   multilooking, subsampling, mirroring, scaling, etc."
       << "\n   This program is useful for cropping and displaying, in combination"
       << "\n   with e.g., GMT, ImageMagick, or xv."
       << "\n   Output format to stdout can be binary."
       << "\n   Careful! only pipe or redirect this."
       << "\n   (use: \"cpxfiddle -h |& more\"     in csh "
       << "\n      or \"cpxfiddle -h 2>&1 | less\" in bash for more help.)"
       << endl;
  } // END shortexpl

/****************************************************************
 * synopsis                                                     *
 #%// BK 11-May-2000                                            *
 ****************************************************************/
void synopsis(char *programname)
  {
  cerr << "\n  PROGRAM: " << programname 
       << " (version " << SWVERSION << ")\n";
  cerr << "\n  SYNOPSIS:\n\t" << programname
       << " -w width [-f informat] [-q output] [-o outformat]"
       << "\n\t[-e exp] [-s scale] [-l line]"
       << " [-L line] [-p pixel] [-P pixel]"
       << "\n\t[-S x/y] [-M x/y] [-m mirror] [-c file]"
       << " [-r rmin/rmax] [-B swap] "
       << "\n\t[-H bytes] [-V] [-b] [-h[elp]] "
       << "[--] [--ignorenan] [--fullstat] inputfile\n";
  shortexpl();
  cerr << endl;
  exit(-1);
  } // END synopsis



/****************************************************************
 * usage                                                        *
 #%// BK 11-May-2000                                            *
 ****************************************************************/
void usage(char *programname)
  {
  cerr << "\n  PROGRAM: " << programname
       << " (version " << SWVERSION << ")\n";

  shortexpl();
  cerr << "\n\n"
       << "\n\t   1   range (X) -->"
       << "\n\t 1  ---------------------"
       << "\n\t   |                     | "
       << "\n\t   |  <Complex file>     |"
       << "\n\t   | e.g. float (2x4B)   |  azimuth (Y)"
       << "\n\t   |  format             |    |"
       << "\n\t   |                     |    |"
       << "\n\t   | (major row order)   |    |"
       << "\n\t   | (pixel interleav.)  |    \\/"
       << "\n\t   |                     | "
       << "\n\t    ---------------------\n"

       << "\n   First a cutout is made,"
       << "\n   then the data is mirrored,"
       << "\n   and afterwards multilooked or subsampled,"
       << "\n   and finally the data is scaled and exp-ed."
       << "\n   This explains why a mirrored image is not equal to"
       << "\n   the original if the number of lines cannot be exactly"
       << "\n   divided by the multilook or subsampling factor.\n"

       << "\n\n  DESCRIPTION [default]:"
       << "\n   -w mandatory      Line length (width, rangepixels, X direction)"
       << "\n   -e [1.0]          Exponent: out=scale*output^exp"
       << "\n   -s [1.0]          Scale     out=scale*output^exp"
       << "\n   -l [1]            First (azimuth Y) line"
       << "\n   -L [all]          Last (azimuth Y) line"
       << "\n   -p [1]            First (range X) pixel"
       << "\n   -P [all]          Last (range X) pixel"
       << "\n   -M [1/1]          Multilook factor in X/Y direction (range/azimuth)"
       << "\n                      No subsampling can be used combined with this"
       << "\n                      option. (Multilooking takes complex sum over"
       << "\n                       window(X,Y), divides by number of looks)"
       << "\n                      Output size: [floor(P-p+1)/mlX; floor(L-l+1)/mlY]."
       << "\n   -S [1/1]          Subsample factor in X/Y (range/azimuth)"
       << "\n                      Output dimensionY = ceil((L-l+1)/Y)."
       << "\n                      The last line does not have to equal -L."
       << "\n                      Output dimensionX = ceil((P-p+1)/X)."
       << "\n                      The last pixel does not have to equal -P."
       << "\n   -q [mag]          What to output:"
       << "\n                     NORMAL | MAG | PHASE | MIXED | REAL | IMAG"
       << "\n                      normal    = (real, imag),"
       << "\n                      magnitude = sqrt(real^2 + imag^2),"
       << "\n                      phase     = atan2(imag,real),"
       << "\n                      mixed     = phase overlay, only with -o sunraster"
       << "\n                      real      = line[2*j],"
       << "\n                      imag      = line[2*j+1]."
       << "\n                     Normal option can be (mis)used to fiddle with, e.g.,"
       << "\n                      float files (though even linelength required?)."
       << "\n   -f [cr4]          Input format identifier:"
       << "\n                     CC1 | CUC1 | CI2 | CI4 | CR4 | CR8"
       << "\n                      for complex 2x1B char, unsigned char, "
       << "\n                      short integer, integer, float, double"
       << "\n                      (major row order pixel interleaved complex data.)"
       << "\n   -o [ascii]        Output format identifier (to stdout):"
       << "\n                     UCHAR | SUNRASTER | FLOAT | SHORT | ASCII"
       << "\n                      uchar, short, sunraster, and float options"
       << "\n                      write binary to stdout!"
       << "\n   -c [gray]         Colormap for \"-o sunraster\" option."
       << "\n                      FILENAME | GRAY | JET | HOT | COOL | BERT | HSV"
       << "\n                     If cmap is a filename, then a ascii file with"
       << "\n                     256 lines, containing r g b [0:255] is assumed."
       << "\n                     If cmap is a identifier, this colormap is used."
       << "\n   -m code           Flag to mirror file content in X or Y direction:"
       << "\n                     X | Y | XY | YX"
       << "\n                      If mirroring Y then the first output line is"
       << "\n                      line number P (default last in file). This means"
       << "\n                      that if -Y is specified, the last line not"
       << "\n                      necessarely is line -l (firstline default 1)"
       << "\n                      The same is true for mirroring in X direction."
       << "\n   -r [rmin/rmax]    uses given minimum and maximum range on data as scalling parameters. "
       << "\n                      Basically, it skips data sampling step for statistics (min,max,mean) computation"
       << "\n                      and defines (min,max,mean) values from given parameters."
       << "\n                      It only effects magnitude, normal and mixed mode outputs."
       << "\n                      For example, use it to scale coherence maps to  0-1 range."
       << "\n                      Tip: no need to use --ignorenan when you are using -r option."
       << "\n   --ignorenan       ignores nan values during calculation of min, max and mean for scaling of sunraster outputs. "
       << "\n                      It only effects magnitude, normal and mixed mode output types."
       << "\n                      Additionally, when binary output is done for magnitude images 0 is set instead of nan values. " 
       << "\n   --fullstat        compute min, max and mean evaluating every pixel for scaling of sunraster outputs. "
       << "\n   -B b|s            Swap bytes."
       << "\n                      if -Bs then swap to host order."
       << "\n                       (call ntohs or ntohl, thus swapping from big endian"
       << "\n                       to host type, which is small endian for X86 PC's.)"
       << "\n                      if -Bb then swap to network order."
       << "\n                       (call htons or htonl, thus swapping to big endian"
       << "\n                       from host type, which is small endian for X86 PC's.)"
       << "\n                      On a big endian platform both these functions are"
       << "\n                      defined as NULL macros, and do nothing."
       << "\n   -b                Add a scalebar in lower right corner for -o sunraster"
       << "\n                      Only for -q phase and mixed."
       << "\n   -H bytes          Skip Header bytes, e.g., for GENESIS SUNraster floats."
       << "\n   -V                Verbose."
       << "\n   -h[elp]           This help."

       << endl
       << "\n\n  EXAMPLES:"
       << "\n   Check first few values in complex short file (width 100, height 200):"
       << "\n      " << programname << " -w100 -fci2 -qnormal -oascii -L5 -P3 -V -- cpxfile"
       << endl
       << "\n   To generate a grd file for use with GMT, without creating a large"
       << "\n   tmpfile, use something like:"
       << "\n      " << programname << " -w100 -fci2 -qphase -ofloat -- cpxfile | \\"
       << "\n      " << " xyz2grd -Gfile.grd -R1/100/1/200 -Zf"
       << endl
       << "\n   To crop a binary complex float file of width 781 and height 501:"
       << "\n   (redirect, or your screen will be flooded with binary output!)"
       << "\n      " << programname << " -w781 -fcr4 -qnormal -ofloat -l101 -L200 \\"
       << "\n      " << " -p51 -P750 cpxfile > cpxfile.cropped"
       << endl
       << "\n   To convert the phase of the same file to SUNraster format:"
       << "\n      " << programname << " -w781 -fcr4 -qphase -osunraster \\"
       << "\n      " << " -ccmap.rainbow -- cpxfile > phase.ras"
       << "\n      display phase.ras"
       << "\n   Where the colortable file is generated by (using GMT):"
       << "\n      makecpt -Crainbow -T0/256/1 | \\"
       << "\n        awk '{if (NR>3&&NR<260){print $2, $3, $4}}' > cmap.rainbow"
       << endl
       << "\n   To create a sort of 8 bit pixmap, visualize with ImageMagick's"
       << "\n   display as gray format:"
       << "\n      " << programname << " -qphase -fcr4 -ouchar -w1 -- cpxfile > phase.uc"
       << "\n      display -size 781x501 gray:phase.uc"
       << "\n   The file phase.uc could also be converted to another format, e.g.:"
       << "\n      rawtopgm 781 501 phase.uc > phase.pgm"
       << "\n      xv phase.pgm"
       << endl
       << "\n   But easier would have been to use " << programname 
       << " to generate the"
       << "\n   SUNraster file. This can be done by:"
       << "\n      " << programname << " -w100 -qmag -osunraster -e0.3 -s1.1 file.cr4 > file.ras"
       << endl
       << "\n   Where file.cr4 is complex real4, -e is used to equalize the "
       << "\n   histogram, and -s to threshold the data (larger scale means"
       << "\n   more white). To overlay the phase as layer over the magnitude:"
       << "\n      " << programname << " -w100 -qmixed -osunraster -e0.3 -- file.cr4 > file.ras"
       << "\n      " << "xv file.ras"
       << endl
       << "\n\n   To read the header of a (non complex!) sunraster file "
       << "\n      (8 int32 values):"
       << "\n      " << programname << " -w1 -L4 -qnormal -fci4 -- file.ras"
       << "\n   This program can be tricked to crop float files (non complex)"
       << "\n   This program can be used for conversions of files, trick with -w1"
       << "\n   To crop a bigendian file on a linux PC, swap and crop like:"
       << "\n      " << programname << " -w100 -fci2 -qnormal -oshort -Bs -- cpxfile > newfile"

       << endl
       << "\n\n   To read the content of a (non complex!) float vector with 32 byte header"
       << "\n      " << programname << " -w1 -qnormal -fr4 -H32  -oascii -- file.ras"

       << endl
       << "\n\n   To scale the content of a float file: coherence map to [ 0 1 ] range."
       << "\n      " << programname << " -w 1000 -q normal -o sunraster -b -c gray -f r4 cohmap.r4 -r 0,1 > cohmap.ras "

       << endl
       << "\n\n  BUGS:"
       << "\n   -M with -m:   First mirrored then multilooked."
       << "\n   -B:           Not tested properly, but simply calls functions."
       << "\n   all:          If cpxfile is a symbolic link, filesize cannot be"
       << "\n                 determined, and therefor the number of lines."
       << "\n                 Solution: Use the full path or use a hard link."
       << "\n   (lot of options, binary to screen, but these are features.)"

       << endl
       << "\n\n  SEE ALSO:"
       << "\n   cpx2ps, display, rawtopgm, makecpt (GMT)"

       << "\n\n  Please send your comments to: TUDelft, doris_users@tudelft.nl"
       << "\n   "
       << "\n   (BTW: view this help by piping stderr: " << programname << " -h |& more      in csh"
       << "\n                                       or " << programname << " -h 2>&1 | less  in bash.)"
       << endl << endl;
  exit(-1);
  } // END usage


//#EOF
