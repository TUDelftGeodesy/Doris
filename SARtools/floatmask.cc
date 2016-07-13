// floatmult.c
// baarda:: g++ -O floatmult.c -o floatmult
// usage: floatmult infile float.number
//
// exit levels:
//  0: successful exit;
// -1: wrong input;
//  1: input file cannot be opened;
//  2: output file cannot be opened;
//
// BK 18-Oct-2000
// $Revision: 3.5 $  $Date: 2003/04/14 06:25:19 $
//
// PM 30-May-2006
// $Revision: 3.6 $  $Date: 2006/05/29 14:00:00 $
//   - now sstream compatible, see bzr for details
//
// MA 19-Dec-2008
// $Revision: 3.7 $
//  - support file size > 4GB
//
// BO 21-Feb-2013
//  - Modified to generate masks

using namespace std;
#include <iostream>				// cout
#include <sstream>				// file
#include <fstream>				// file
//#include <strstream>				// file
#include <cstdlib>				// exit
#include <cmath>				// atan etc.
#include <complex>				// conj etc.
#include <cstring>				// strcat

// ====== Typedefs for portability ======
typedef long long           int64;
typedef unsigned long long  uint64;


void usage(char *programname)
{
  cerr << "\nProgram: " << programname 
       << " pixelwise mask generation from a float file.\n"
       << " To be used to scale float files.\n"
       << " see also: floatmult, cpxmult, flapjack, cpxfiddle... \n"
       << "\n  USAGE:\n\t" << programname
       << " infile1 gt|lt|eq [threshold] \n\n"
       << "  EXAMPLE:\n\t" << programname << " interfero.coh lt 0.4 \n"
       << "\n default threshold is 0.5\n\n"
       << "\noutput file == infile1.gt|lt|eq.threshold\n\n"
       << "exit levels: 0:ok; -1: wrong input; 1:ifile; 2:ofile.\n"
       << "please send comments to: batu@gi.alaska.edu\n\n\n";
  exit(-1);
  }

bool (*check)(float, float) = NULL;
bool gt(float a, float b){ return a>b; }  
bool lt(float a, float b){ return a<b; }  
bool eq(float a, float b){ return a==b; }  

int main(int argc, char* argv[])
{
  char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  const int ONE27 = 127;
  char ifile[ONE27];				// input file name
  const int sizeofelement = sizeof(float);	// data in file
  float factor = 0.5;				// default
  // ====== Handle input ======
  switch (argc)
    {
    case 4:
      factor = atof(argv[3]);
    case 3:
      if (strcmp(argv[2], "gt") == 0)
        check=&gt;
      else if (strcmp(argv[2], "lt") == 0)
        check=&lt;
      else if (strcmp(argv[2], "eq") == 0)
        check=&eq;
      else{
        usage(argv[0]);
        exit(-1);
      }  
      //--- fall through ---//
    case 2:
      strcpy(ifile,argv[1]);			// input filename arg1
      break; // ---      ---//
    default:
      usage(argv[0]);
    } // switch input
  
  // ______ Check / echo input ______
  cerr << "Program parameters:\n\t" << argv[0] << " " 
       << ifile << " "<< argv[2] << " " << factor << endl;
  //if (factor == 1) usage(argv[0]);
  
  // ______ Set defaults if required _____
  char ofile[ONE27];		// output filename == "ifile.flapjack"
  
  // [PM:29.05.2006]: start
  //  ostrstream omem(ofile,ONE27);
  ostringstream omem;
  //  omem << ifile << "." << argv[0] << factor << ends;
  omem << ifile << "." << argv[2] << factor << ends;
  strcpy(ofile,omem.str().c_str());

  //quick debugging
  //cout << "Test [omem_stream]: " << omem.str() << "\n";
  //cout << "Test [ofile]: " << ofile << "\n";
  
  // [PM:29.05.2006]: stop
  
  // ====== Start complex exp. ======
  //ifstream infile1(ifile, ios::in | ios::nocreate | ios::binary);
  ifstream infile1(ifile, ios::in | ios::binary);
  if (!infile1) 
    {
      cerr << "Problem opening file: " << ifile << ends;
      exit(1);
    }
  infile1.seekg(0,ios::end);              	// filepointer at end
  const streamoff totalbytes1 = infile1.tellg();
  infile1.seekg(0,ios::beg);			// start of file1
  
  const uint64 numberofpixels = totalbytes1/sizeofelement;
  cerr << "Total number of float pixels in \""
       << ifile << "\": " << numberofpixels << endl; 
  ofstream outfile(ofile, ios::out | ios::binary | ios::trunc);
  if (!outfile)
    {
      cerr << "Problem opening output file: " << ofile << endl;
    exit(2);
    }
  
  int tenpercent = int(floor(numberofpixels/10.));
  int percent = 0;
  register float value;
  for (register int i=0; i<numberofpixels; ++i)
    {
      infile1.read((char*)&value,sizeofelement);
      if ( (*check)(value, factor) )
        value=0.0;
      else
        value=1.0;
      outfile.write((char*)&value,sizeofelement);
      if (!(i%tenpercent))
	{
	  cerr << "\rprocessed: " << percent << "%";
	  percent += 10;
	}
    } // loop over all pixels
  cerr << endl;
  
  
  // ====== Tidy up ======
  infile1.close();
  outfile.close();
  cerr << "\nThank you for using " << argv[0] << "!\n";
  cerr << "Output in: " << ofile << "\n";
  return 0;
} // END


