//---------------------------------------------------//
// g++ dateconv.cc -o dateconv
// utility to convert a date string from one format to another
// see man ctime, strptime, strftime
//  (c) Bert Kampes, 06-Aug-2004
//---------------------------------------------------//

using namespace std;
#include <iostream>         // cout
#include <cstdio>         // sprintf?
#include <cstdlib>         // exit() 
#include <ctime>           // strptime, strftime
#include <cstring>        // strcpy


/* ***** */
int main(int argc, char* argv[])
  {
  bool verbose = false; //true;
  if (verbose) 
    cerr << "argc: " << argc << endl;
  if (argc != 4) 
    {
    cerr << argv[0] << "  -- Format/convert dates\n\n"
         << "Usage:  " << argv[0] << " datestr  informat  outformat\n"
         << "  where    datestr:    string, e.g., 15-JUL-2000\n"
         << "           informat:   specifier of datestr, see man strptime\n"
         << "           outformat:  specifier for output format, see man strptime\n\n"
	 << "  examples:\n"
	 << "    to obtain the day number in the year:\n"
	 << argv[0] << " \"15-JUL-2000\" \"%d-%b-%Y\" \"%j\"\n\n"
	 << "    to convert week number, day to date:\n"
	 << argv[0] << " \"week 15, day 2\" \"week %W, day %w\" \"%a %d-%B-%Y\"\n\n\n"
	 << "  also see man ctime, man strptime, man strftime\n";
    exit(1);
    }


  // --- Read input ---------------------------------------------
  char indatestring[128];
  char informat[128];
  char outformat[128];
  strcpy(indatestring,argv[1]);
  strcpy(informat,argv[2]);
  strcpy(outformat,argv[3]);
  if (verbose) 
    cerr << "input: " << indatestring << " " << informat << " " << outformat << endl;

  // --- fill tm struct using daynumber with strptime, re-format it using strftime ---
  struct tm tm_tmp;
  strptime(indatestring,informat,&tm_tmp);

  // --- And use strftime to format the output string ---
  char outdatestring[128];// output
  int q = (int)(strftime(outdatestring,127,outformat,&tm_tmp));

  // --- man ctime --------------------------------------------
  //int  tm_sec;        /* seconds after the minute - [0, 61] */
  //int  tm_min;        /* minutes after the hour - [0, 59] */
  //int  tm_hour;       /* hour since midnight - [0, 23] */
  //int  tm_mday;       /* day of the month - [1, 31] */
  //int  tm_mon;        /* months since January - [0, 11] */
  //int  tm_year;       /* years since 1900 */
  //int  tm_wday;       /* days since Sunday - [0, 6] */
  //int  tm_yday;       /* days since January 1 - [0, 365] */
  //int  tm_isdst;      /* flag for alternate daylight savings time */
  if (verbose) 
    cerr << "tm_tmp.tm_sec:  "   << tm_tmp.tm_sec   << endl
         << "tm_tmp.tm_min:  "   << tm_tmp.tm_min   << endl
         << "tm_tmp.tm_hour:  "  << tm_tmp.tm_hour  << endl
         << "tm_tmp.tm_mday:  "  << tm_tmp.tm_mday  << endl
         << "tm_tmp.tm_mon:  "   << tm_tmp.tm_mon   << endl
         << "tm_tmp.tm_year:  "  << tm_tmp.tm_year  << endl
         << "tm_tmp.tm_wday:  "  << tm_tmp.tm_wday  << endl
         << "tm_tmp.tm_yday:  "  << tm_tmp.tm_yday  << endl
         << "tm_tmp.tm_isdst:  " << tm_tmp.tm_isdst << endl;
  // ---------------------------------------------------------
  cout << outdatestring << endl;
  exit(0);
  } //EOF

