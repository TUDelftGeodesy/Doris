// cpxmult.c
// baarda::g++ -O cpxmult.c -o cpxmult
// usage: cpxmult infile infile2 [outfile [add]]
// cpxmult: pixelwise complex subtraction of 2 files (if add=1 then adds).
// files: major row order binary (complex) float pixel interleaved
//
// exit levels:
//  0: successful exit;
// -1: wrong input;
//  1: input file cannot be opened;
//  2: output file cannot be opened;
//
// $Revision: 3.5 $  $Date: 2003/04/14 06:25:52 $
// BK 19-Apr-2000
// MA 64-bit update 2009-10-13

using namespace std;
#include <iostream>				// cout
#include <fstream>				// file
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
       << " subtracts or adds phase of two complex float files\n"
       << "\n  USAGE:\n\t" << programname
       << " infile1 infile2 [outfile [add==0]]\n\n"
       << "  infile[12] contain complex values a+ib\n"
       << "  outfile contains (a1+ib1)*conj(a2+ib2)\n\n" 
       << "  If add=1 then phase is added (not conj).\n\n"
       << "  EXAMPLE:\n\t" << programname
       << " cint.raw cint2.raw subtract.raw"
       << "\ndefault output file == infile1.min.infile2\n\n"
       << "exit levels: 0:ok; -1: wrong input; 1:ifile; 2:ofile.\n\n\n";
  exit(-1);
}


int main(int argc, char* argv[])
{
  char ident[] = "@(#)Doris software, doris_users@tudelft.nl";
  const int ONE27   = 127;
  char ifile1[ONE27];				// input file name
  char ifile2[ONE27];				// input file name
  char ofile[ONE27] = " ";			// output filename == "ifile.ml"
  bool add = false;				// take conj by default
  const int sizeofelement = sizeof(complex<float>);	// data in file
  
  // ====== Handle input ======
  switch (argc)
    {
    case 5:
      add=true;					// should be checked...
      //--- fall through ---//
    case 4:
      strcpy(ofile,argv[3]);			// output filename arg3
      //--- fall through ---//
    case 3:
      strcpy(ifile2,argv[2]);			// input filename arg1
      //--- fall through ---//
    case 2:
      strcpy(ifile1,argv[1]);			// input filename arg1
      break; // ---      ---//
    default:
      usage(argv[0]);
    } // switch input
  
  // ______ Set defaults if required _____
  if (!strcmp(ofile," "))			// nothing specified
    {
      strcpy(ofile,ifile1);
      strcat(ofile,".min.");
      strcpy(ofile,ifile2);
    }
  
  // ______ Check / echo input ______
  cerr << "Program parameters:\n\t" << argv[0] << " " 
       << ifile1 << " " << ifile2 << " " << ofile << " " << add << endl;
  if (!strcmp(ofile,ifile1) || !strcmp(ofile,ifile2))
    {
      cerr << "input file name same as other one: "
	   << ifile1 << " =? " << ofile << " =? " << ifile2 << endl;
      usage(argv[0]);
    }
  
  
  // ====== Start complex multiply ======
  //ifstream infile1(ifile1, ios::in | ios::nocreate | ios::binary);

  ifstream infile1(ifile1, ios::in | ios::binary);

  if (!infile1) cerr << "Problem opening file: " << ifile1 << endl, exit(1);
  infile1.seekg(0,ios::end);              	// filepointer at end
  
  //const int totalbytes1 = infile1.tellg();
  const streamoff totalbytes1 = infile1.tellg(); // MA
  infile1.seekg(0,ios::beg);			// start of file1
  //ifstream infile2(ifile2, ios::in | ios::nocreate | ios::binary);
  ifstream infile2(ifile2, ios::in | ios::binary);
  
  if (!infile2) cerr << "Problem opening file: " << ifile2 << endl, exit(1);
  infile2.seekg(0,ios::end);              	// filepointer at end
  
  //const int totalbytes2 = infile2.tellg();
  const streamoff totalbytes2 = infile2.tellg(); // MA
  infile2.seekg(0,ios::beg);			// start of file2
  
  if (totalbytes1!=totalbytes2)
    {
      cerr <<  "size of file1 != size of file2: "
	   << totalbytes1 << "!=" << totalbytes2 << endl;
      exit(2);
    }
  
  const uint64 numberofpixels = totalbytes1/sizeofelement;
  ofstream outfile(ofile, ios::out | ios::binary | ios::trunc);
  
  if (!outfile) cerr << "Problem opening file: " << ofile << endl, exit(2);
  int tenpercent = int(floor(numberofpixels/10.));
  int percent = 0;

  register complex<float> value1;
  register complex<float> value2;
  // ______ Good compiler would get rid of if in for loop, but 2b sure.
  if (add)
    {
      for (register int i=0; i<numberofpixels; ++i)
	{
	  infile1.read((char*)&value1,sizeofelement);
	  infile2.read((char*)&value2,sizeofelement);
	  value1 *= value2;
	  outfile.write((char*)&value1,sizeofelement);
	  if (!(i%tenpercent))
	    {
	      cerr << "\rprocessed: " << percent << "%";
	      percent += 10;
	    }
	} // loop over all pixels
    }
  else
    {
      for (register int i=0; i<numberofpixels; ++i)
	{
	  infile1.read((char*)&value1,sizeofelement);
	  infile2.read((char*)&value2,sizeofelement);
	  value1 *= conj(value2);
	  outfile.write((char*)&value1,sizeofelement);
	  if (!(i%tenpercent))
	    {
	      cerr << "\rprocessed: " << percent << "%";
	      percent += 10;
	    }
	} // loop over all pixels
    }
  
  
// ====== Tidy up ======
  infile1.close();
  infile2.close();
  outfile.close();
  cout << "\n\nThank you for using " << argv[0] << "!\n";
  return 0;
} // END
