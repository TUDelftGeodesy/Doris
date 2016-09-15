// flapjack.c
// baarda:: g++ -O flapjack.c -o flapjack
// usage: flapjack infile int.number
//
// exit levels:
//  0: successful exit;
// -1: wrong input;
//  1: input file cannot be opened;
//  2: output file cannot be opened;
//
// BK 18-Oct-2000
// $Revision: 3.5 $  $Date: 2003/04/14 06:25:52 $
//
// PM 29-May-2006
// $Revision: 3.6 $  $Date: 2006/05/29 17:34:00 $
//   - now sstream compatible, see bzr for details
//
// MA 19-Dec-2008
// $Revision: 3.7 $
//  - support file size > 4GB


#include <iostream>				// cout
//#include <strstream>				// file
#include <sstream>				// file
#include <fstream>				// file
#include <cstdlib>				// exit
#include <cmath>				// atan etc.
#include <complex>				// conj etc.
#include <cstring>				// strcat
#include <string>				// strcat

// ====== Typedefs for portability ======
typedef long long           int64;
typedef unsigned long long  uint64;

using namespace std;


void usage(char *programname)
  {
  cerr << "\nProgram: " << programname 
       << " pixelwise complex integer multiplication of a float complex\n"
       << " file. To be used to make linear combinations. \n"
       << " See also: cpxmult... \n"
       << "\n  USAGE:\n\t" << programname
       << " infile1 [factor==2]\n\n"
       << "  EXAMPLE:\n\t" << programname << " cint.raw 3"
       << "\noutput file == infile1." << programname << ".factor\n\n"
       << "exit levels: 0:ok; -1: wrong input; 1:ifile; 2:ofile.\n"
       << "please sent comments to: doris_users@tudelft.nl\n\n\n";
  exit(-1);
  }


int main(int argc, char* argv[])
  {
  char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  cerr << ident << endl;
  const int ONE27   = 127;
  char ifile[ONE27];				// input file name
  const int sizeofelement = sizeof(complex<float>);	// data in file
  int factor = 2;				// default

// ====== Handle input ======
  switch (argc)
    {
    case 3:
      factor = atoi(argv[2]);			// input filename arg1
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
  if (factor < 2) usage(argv[0]);


  /*
    

  */

  // ______ Set defaults if required _____
  char ofile[ONE27];			// output filename == "ifile.flapjack"

  // [PM:29.05.2006]: start
  //  ostrstream omem(ofile,ONE27);
  ostringstream omem;
  //  omem << ifile << "." << argv[0] << factor << ends;
  omem << ifile << "." << argv[0] << factor << ends;
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
      cerr << "Problem opening file: " << ifile << endl;
      exit(1);
    }
  infile1.seekg(0,ios::end);              	// filepointer at end
  const streamoff totalbytes1 = infile1.tellg();
  infile1.seekg(0,ios::beg);			// start of file1
  
  const uint64 numberofpixels = totalbytes1/sizeofelement;
  cerr << "Total number of complex pixels in \""
       << ifile << "\": " << numberofpixels << endl; 
  ofstream outfile(ofile, ios::out | ios::binary | ios::trunc);
  if (!outfile)
    {
      cerr << "Problem opening file: " << ofile << endl;
      exit(2);
    }
  
  int tenpercent = int(floor(numberofpixels/10.));
  int percent = 0;
  register complex<float> value,result;
  register float mag;
  if (factor==2)
    {
      for (register int i=0; i<numberofpixels; ++i)
	{
	  infile1.read((char*)&value,sizeofelement);
	  mag    = sqrt(value.real()*value.real()+
			value.imag()*value.imag());
	  result = (value*value) / mag;
	  outfile.write((char*)&result,sizeofelement);
	  if (!(i%tenpercent))
	    {
	      cerr << "\rprocessed: " << percent << "%";
	percent += 10;
	    }
	} // loop over all pixels
    }
  else if (factor==3)
    {
      for (register int i=0; i<numberofpixels; ++i)
	{
	  infile1.read((char*)&value,sizeofelement);
	  mag    =      value.real()*value.real()+
	            value.imag()*value.imag();
	  result = (value*value*value) / mag;
	  outfile.write((char*)&result,sizeofelement);
	  if (!(i%tenpercent))
	    {
	      cerr << "\rprocessed: " << percent << "%";
	      percent += 10;
	    }
	} // loop over all pixels
    }
  else // factor > 3
    {
      for (register int i=0; i<numberofpixels; ++i)
	{
	  infile1.read((char*)&value,sizeofelement);
	  mag = sqrt(value.real()*value.real()+
		     value.imag()*value.imag());
	  result = value;
	  value /= mag;		// magnitude now 1.
	  for (register int dummy=1; dummy<factor; ++dummy)
	    result *= value;
	  //result /= pow(mag,factor-1);
	  outfile.write((char*)&result,sizeofelement);
	  if (!(i%tenpercent))
	{
	  cerr << "\rprocessed: " << percent << "%";
	  percent += 10;
	}
	} // loop over all pixels
    }
  cerr << endl;
  
  
  // ====== Tidy up ======
  infile1.close();
  outfile.close();
  cerr << "\nThank you for using " << argv[0] << "!\n";
  cerr << "Output in: " << ofile << "\n";
  return 0;
  } // END


