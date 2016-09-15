// wrap.c
// geos1::aCC wrap.c -o wrap
// usage: wrap file [a b [ofile]]
// wrap file to interval [a b)
// file: binary float pixel interleaved
// default interval [-pi pi)
// default output file: ifile.wrap
//
// exit levels:
//  0: successful exit;
// -1: wrong input;
//  1: input file cannot be opened;
//  2: output file cannot be opened;

// idea to read in diff. formats (float int etc.)
// and output the values in format char (1B per pixel)
//
//  RCS: // // // // // // // // // // // // // // // // // // // // //
// $Header: /users/kampes/DEVELOP/DORIS/SARtools/RCS/wrap.cc,v 3.6 2003/06/24 08:37:12 kampes Exp $ //
// // // // // // // // // // // // // // // // // // // // // // // //


using namespace std;
#include <iostream>				// cout
#include <fstream>				// file
#include <cstdlib>				// exit
#include <cmath>				// atan etc.
#include <cstring>      // strcat etc.


void usage(char *programname)
  {
  cerr << "\nProgram: " << programname 
       << " wraps float binary file to interval [a,b)\n"
       << "\n  USAGE:\n\t" << programname
       << " infile [a b [ofile]]\n\n"
       << "  EXAMPLE:\n\t" << programname
       << " interferogram.raw -4pi 4pi interf4.raw\n"
       << "\ndefault output file    == infile.wrap"
       << "\ndefault interval [a b) == [-pi pi)\n\n";
  exit(-1);
  }

int main(int argc, char* argv[])
  {
  //char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  char ident[] = "@(#)wrap: Doris software, $Revision: 3.6 $, $Author: kampes $";
  cerr << ident << endl;
  const int ONE27   = 127;
  const float PI    = 4.*atan(1.);		// pi
  char ifile[ONE27];				// input file name
  char ofile[ONE27] = " ";			// output filename == "ifile.ml"
  char dummy[ONE27];				// dummy string
  //char dummy2[ONE27];				// dummy string
  float a = -999.;				// start interval
  float b = -999.;				// stop interval
  register int i;				// general counter
  const int sizeofelement = sizeof(float);	// data in float file

// ====== Handle input ======
  switch (argc)	// wrap infile a b outfile
    {
    case 3:
      usage(argv[0]);				// 3 arg. not allowed.
      break;					// not required

    case 5:
      strcpy(ofile,argv[4]);			// output filename arg4
      //--- fall through ---//

    case 4:
      strcpy(dummy,argv[3]);			// interval: b
      i = (int)strlen(dummy);
      if (dummy[i-1] == 'i')			// likely to be pi
	{
	if (!strcmp(dummy,"pi"))
	  b = PI;
	else if (!strcmp(dummy,"-pi"))
	  b = -PI;
	else
	  {
	  dummy[i-2]='\0';
	  b = PI*atof(dummy);
	  }
	}
      else
	{
        b = atof(dummy);
	}

      strcpy(dummy,argv[2]);			// interval: a
      i = (int)strlen(dummy);
      if (dummy[i-1] == 'i')			// likely to be pi
	{
	if (!strcmp(dummy,"pi"))
	  a = PI;
	else if (!strcmp(dummy,"-pi"))
	  a = -PI;
	else
	  {
	  dummy[i-2]='\0';
	  a = PI*atof(dummy);
	  }
	}
      else
	{
        a = atof(dummy);
	}
      //--- fall through ---//

    case 2:
      strcpy(ifile,argv[1]);			// input filename arg1
      break; // ---      ---//

    default:
      usage(argv[0]);
    } // switch input

// ______ Set defaults if required _____
  if (abs(a+999.) < 1e-10)
    a = -PI;
  if (abs(b+999.) < 1e-10)
    b =  PI;
  if (!strcmp(ofile," "))			// nothing specified
    {
    strcpy(ofile,ifile);
    strcat(ofile,".wrap");
    }

// ______ Check input ______
  if (a>=b)
    {
    cerr << "interval (a,b) = (" << a << "," << b 
	 << "): a should be smaller than b.\n";
    usage(argv[0]);
    }

  if (abs(a+b)>1e-10 && abs(a)>1e-10)
    {
    cerr << "interval (a,b): should be either: a=-b or a=0\n";
    usage(argv[0]);
    }
  if (!strcmp(ofile,ifile))
    {
    cerr << "input file name same as output file name: "
	 << ifile << " = " << ofile << endl;
    usage(argv[0]);
    }

  cerr << "Program parameters:\n\t" << argv[0] << " " 
       << ifile << " " << a << " " << b << " " << ofile << endl;


// ====== Start wrapping ======
//#if __GNUC__ >= 3
//  ifstream image(ifile, ios::in | ios::nocreate | ios::binary);
//#else
  ifstream image(ifile, ios::in | ios::binary);
//#endif

  image.seekg(0,ios::end);                    // filepointer at end

  if (!image) cerr << "Problem opening file: " << ifile << endl, exit(1);
  const int totalbytes = image.tellg();
  const int numberofpixels = totalbytes/sizeofelement;
  image.seekg(0,ios::beg);			// start of file
  ofstream wrapped(ofile, ios::out | ios::binary | ios::trunc);
  if (!wrapped) cerr << "Problem opening file: " << ofile << endl, exit(2);

  register float phase;
  const float bereik = b-a;
  const float normal = (.5*bereik)/PI;
  bool ais0 = false;
  if (abs(a)<1e-10)
    ais0 = true;

  for (i=0; i<numberofpixels; ++i)
    {
    image.read((char*)&phase,sizeofelement);

// ______ c = cos(p)+ i*sin(p); p=atan2(im,re); ______
    phase /= normal;
    phase  = atan2(sin(phase),cos(phase));	// [-pi,pi]
    phase *= normal;				// [a,-a]
    if (ais0)
      if (phase < 0)
	phase += bereik;			// [0,b]

    wrapped.write((char*)&phase,sizeofelement);
    } // loop over all pixels

// ______ Tidy up ______
  image.close();
  wrapped.close();

  cout << "\n\nAll done!\n";
  return 0;
  } // END
