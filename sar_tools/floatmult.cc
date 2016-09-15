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
       << " pixelwise float multiplication of a (complex) float complex file.\n"
       << " To be used to scale float files, or magnitude of complex files.\n"
       << " see also: cpxmult, flapjack, cpxfiddle... \n"
       << "\n  USAGE:\n\t" << programname
       << " infile1 [factor==2.]\n\n"
       << "  EXAMPLE:\n\t" << programname << " cint.raw 1.47"
       << "\noutput file == infile1." << programname << ".factor\n\n"
       << "exit levels: 0:ok; -1: wrong input; 1:ifile; 2:ofile.\n"
       << "please sent comments to: doris_users@tudelft.nl\n\n\n";
  exit(-1);
  }


int main(int argc, char* argv[])
{
  char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  const int ONE27 = 127;
  char ifile[ONE27];				// input file name
  const int sizeofelement = sizeof(float);	// data in file
  float factor = 2.;				// default
  
  // ====== Handle input ======
  switch (argc)
    {
    case 3:
      factor = atof(argv[2]);			// input filename arg1
      //--- fall through ---//
    case 2:
      strcpy(ifile,argv[1]);			// input filename arg1
      break; // ---      ---//
    default:
      usage(argv[0]);
    } // switch input
  
  // ______ Check / echo input ______
  cerr << "Program parameters:\n\t" << argv[0] << " " 
       << ifile << " " << factor << endl;
  if (factor == 1) usage(argv[0]);
  
  // ______ Set defaults if required _____
  char ofile[ONE27];		// output filename == "ifile.flapjack"
  
  // [PM:29.05.2006]: start
  //  ostrstream omem(ofile,ONE27);
  ostringstream omem;
  //  omem << ifile << "." << argv[0] << factor << ends;
  omem << ifile << "." << argv[0] << factor << endl;
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
      value *= factor;
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


