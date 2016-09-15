// convert.c
// prog converts binary file formats in to format out.
// formats: short, int, float, double (complex)
// BK 28-Jan-00
// MA 19-Dec-08 : support > 4GB

using namespace std;
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>


// ====== Typedefs for portability ======
typedef long long           int64;
typedef unsigned long long  uint64;

const int ONE27=127*2; // [MA] now 255
struct sizeofel
  {
  int ifile;
  int ofile;
  };

struct input
  {
  char ifilename[ONE27];
  char ofilename[ONE27];
  int iformat;
  int oformat;
  sizeofel bytesperelement;
  };


/* ***** */
void usage(char *programname)
  {
  cerr << "\n  Program: " << programname
       << "\n\tconvert binary data to other format.\n\n"
       << "  SYNOPSIS:\n\t" << programname
       << " inputfile outputfile informat oformat\n\n"
       << "\tinputfile: \tname of binary input file [INFILE]\n"
       << "\toutputfile: \tname of binary output file [OUTFILE]\n"
       << "\tifmt:      \tinput format [2]\n"
       << "\tofmt:      \toutput format [4]\n"
       << "\n  FORMATS:\n"
       << "\t1   \tsigned   (complex) char\n"
       << "\t101 \tunsigned (complex) char\n"
       << "\t2   \tsigned   (complex) short\n"
       << "\t102 \tunsigned (complex) short\n"
       << "\t3   \tsigned   (complex) int\n"
       << "\t103 \tunsigned (complex) int\n"
       << "\t4   \tsigned   (complex) float\n"
//       << "\t104 \tunsigned (complex) float\n"
       << "\t5   \tsigned   (complex) double\n"
//       << "\t105 \tunsigned (complex) double\n"
       << "\nNote: no rounding is performed when converting from e.g. float to int.\n"
       << endl;
  exit(-1);
  }

/* ***** */
void handleinput(int argc, char* argv[], input &options)
  {
  // default short to float
  options.iformat = 2;
  options.oformat = 4;
  strcpy(options.ifilename,"INFILE");
  strcpy(options.ofilename,"OUTFILE");

  // Handle command line arguments
  switch (argc)
    {
    case 5:                                     // ofmt
      options.oformat = atoi(argv[4]);
      //--- fall through --- //
    case 4:                                     // ifmt
      options.iformat = atoi(argv[3]);
      //--- fall through --- //
    case 3:                                     // of
      strcpy(options.ofilename,argv[2]);
      //--- fall through --- //
    case 2:                                     // if
      strcpy(options.ifilename,argv[1]);
      break;

    default:
      usage(argv[0]);
    } // switch handle arguments

cerr << "TEST: "
<< options.ifilename << "; "
<< options.ofilename << "; "
<< options.iformat << "; "
<< options.oformat << "; "
<< endl;

  if      (options.iformat==1 || options.iformat==101)
    options.bytesperelement.ifile = sizeof(char);
  else if (options.iformat==2 || options.iformat==102)
    options.bytesperelement.ifile = sizeof(short);
  else if (options.iformat==3 || options.iformat==103)
    options.bytesperelement.ifile = sizeof(int);
  else if (options.iformat==4 || options.iformat==104)
    options.bytesperelement.ifile = sizeof(float);
  else if (options.iformat==5 || options.iformat==105)
    options.bytesperelement.ifile = sizeof(double);
  else
    usage(argv[0]);
  if      (options.oformat==1 || options.oformat==101)
    options.bytesperelement.ofile = sizeof(char);
  else if (options.oformat==2 || options.oformat==102)
    options.bytesperelement.ofile = sizeof(short);
  else if (options.oformat==3 || options.oformat==103)
    options.bytesperelement.ofile = sizeof(int);
  else if (options.oformat==4 || options.oformat==104)
    options.bytesperelement.ofile = sizeof(float);
  else if (options.oformat==5 || options.oformat==105)
    options.bytesperelement.ofile = sizeof(double);
  else
    usage(argv[0]);

  if (options.iformat == options.oformat)
    usage(argv[0]);
  if (options.bytesperelement.ofile < options.bytesperelement.ifile)
    cerr << "WARNING: writing to less bytes, no rounding of numbers.\n";
  } // handleinput

//  /* *****
//   * rcw: read, convert, write
//   * template to read numelements of TypeA from opened ifile 
//   * and to write them to opened ofile as TypeB
//   * A and B are dummy function arguments, used to identificate correct template function.
//  */
//  template <class TypeA, class TypeB>
//  void rcw(TypeA in, TypeB out, ifstream &ifile, ofstream &ofile, int numelements)
//    {
//    const int sizeofin  = sizeof(TypeA);
//    const int sizeofout = sizeof(TypeB);
//    for (int i=0; i<numelements; i++)
//      {
//      ifile.read((char*)&in,sizeofin);	// read
//      out = TypeB(in);			// convert
//      ofile.write((char*)&out,sizeofout);	// write
//      }
//    }

/* *****
 * writedata
 * template to write numelements of TypeA from array
 * to opened ofile as TypeB
 * A and B are dummy function arguments, used to identificate correct template function.
*/
template <class TypeA>
void writedata(TypeA *in, ofstream &ofile, int numelements, int oformat, int sizeofout)
  {
  int i;
  switch (oformat)
    {
    case 1:
      {
      char out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 101:
      {
      unsigned char out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 2:
      {
      short out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 102:
      {
      unsigned short out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 3:
      {
      int out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 103:
      {
      unsigned int out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }
    case 4:
      {
      float out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }

//      case 104:
//        {
//        unsigned float out;
//        for (i=0; i<numelements; ++i)
//          {
//          out = char(in[i]);			// convert
//          ofile.write((char*)&out,sizeofout);	// write
//          }
//        break;
//        }

    case 5:
      {
      double out;
      for (i=0; i<numelements; ++i)
        {
        out = char(in[i]);			// convert
        ofile.write((char*)&out,sizeofout);	// write
        }
      break;
      }

//      case 105:
//        {
//        unsigned double out;
//        for (i=0; i<numelements; ++i)
//          {
//          out = char(in[i]);			// convert
//          ofile.write((char*)&out,sizeofout);	// write
//          }
//        break;
//        }

    default:
	cerr << "impossible error with checked input.\n";
	exit(6);
    }
  }

/* ***** */
int main(int argc,char* argv[])
  {
  //char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  char ident[] = "@(#)bkconvert: Doris software, $Revision: 4.03 $, $Author: TUDelft $";
  cerr << ident << endl;
  input options;
  handleinput(argc,argv,options);

  //ifstream ifile(options.ifilename, ios::in | ios::binary | ios::nocreate);
  ifstream ifile(options.ifilename, ios::in | ios::binary);
  if (!ifile)
    {
    cerr << "input file: \"" << options.ifilename << "\" cannot be opened!\n";
    return 1;
    }
  ifile.seekg(0,ios::end);		// filepointer at end
  // const int numbytes = ifile.tellg();
  const streamoff numbytes = ifile.tellg(); // [MA]
  if (numbytes%options.bytesperelement.ifile != 0)
    {
    cerr << "file size or format not ok!, exiting\n";
    return 4;
    }
  const streamoff numelements = (streamoff)(numbytes/options.bytesperelement.ifile);
  ifile.seekg(0,ios::beg);		// rewind

  //ofstream ofile(options.ofilename, ios::out | ios::binary | ios::noreplace);
  ofstream ofile(options.ofilename, ios::out | ios::binary);
  if (!ofile)
    {
    cerr << "output file: " << options.ofilename << " cannot be opened!\n"
         << "Should I try again (overwriting existing)? [y/n]\n";
    char dummychar;
    cin >> dummychar;
    if (dummychar=='y' || dummychar=='Y')
      ofstream ofile(options.ofilename, ios::out | ios::binary | ios::trunc);
    else
      return 2;
    }

  int buffersize = 512;         // number of elements in data buffer
  const uint64 NUMBUFFERS   = numelements/buffersize;
  const uint64 restelements = numelements%buffersize;
  const int EXTRABUF     = (restelements!=0) ? 1 : 0;

//    // input buffer arrays
//    // should be done dynamically to save memory.
//    char            data_schar[buffersize];
//    unsigned char   data_uchar[buffersize];
//    short           data_sshort[buffersize];
//    unsigned short  data_ushort[buffersize];
//    int             data_sint[buffersize];
//    unsigned int    data_uint[buffersize];
//    float           data_sfloat[buffersize];
//    unsigned float  data_ufloat[buffersize];
//    double          data_sdouble[buffersize];
//    unsigned double data_udouble[buffersize];

  // input buffer arrays
  // should be done dynamically to save memory.
  char            data_schar[512];
  unsigned char   data_uchar[512];
  short           data_sshort[512];
  unsigned short  data_ushort[512];
  int             data_sint[512];
  unsigned int    data_uint[512];
  float           data_sfloat[512];
  double          data_sdouble[512];

  for (int i=1; i<=NUMBUFFERS+EXTRABUF; i++)
    {
    if (i==NUMBUFFERS+1)	// last smaller buffer
      buffersize=restelements;

    // read data in array
    switch (options.iformat)
      {
      case 1:
        ifile.read((char*)&data_schar[0],buffersize*options.bytesperelement.ifile);
        writedata(data_schar,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 101:
        ifile.read((char*)&data_uchar[0],buffersize*options.bytesperelement.ifile);
        writedata(data_uchar,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 2:
        ifile.read((char*)&data_sshort[0],buffersize*options.bytesperelement.ifile);
        writedata(data_sshort,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 102:
        ifile.read((char*)&data_ushort[0],buffersize*options.bytesperelement.ifile);
        writedata(data_ushort,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 3:
        ifile.read((char*)&data_sint[0],buffersize*options.bytesperelement.ifile);
        writedata(data_sint,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 103:
        ifile.read((char*)&data_uint[0],buffersize*options.bytesperelement.ifile);
        writedata(data_uint,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;
      case 4:
        ifile.read((char*)&data_sfloat[0],buffersize*options.bytesperelement.ifile);
        writedata(data_sfloat,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;

//        case 104:
//          ifile.read((char*)&data_ufloat[0],buffersize*options.bytesperelement.ifile);
//          writedata(data_ufloat,
//          ofile, buffersize, options.oformat, options.bytesperelement.ofile);
//          break;

      case 5:
        ifile.read((char*)&data_sdouble[0],buffersize*options.bytesperelement.ifile);
        writedata(data_sdouble,
        ofile, buffersize, options.oformat, options.bytesperelement.ofile);
        break;

//        case 105:
//          ifile.read((char*)&data_udouble[0],buffersize*options.bytesperelement.ifile);
//          writedata(data_udouble,
//          ofile, buffersize, options.oformat, options.bytesperelement.ofile);
//          break;

      default:
	cerr << "impossible error with checked input.\n";
	return 5;
      }
    }
  ifile.close();
  ofile.close();

  //rcw(gettypea(options.iformat), gettypeb(options.oformat),
  //    ifile, ofile, numelements);

  cout << "Normal termination.\n";
  cerr << "Output file: \"" << options.ofilename << "\"\n";
  return 0;
  }



