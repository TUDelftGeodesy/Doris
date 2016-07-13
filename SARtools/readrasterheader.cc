// prog to read a header...
// g++ -O readrasterheader.c -o readrasterheader
// BK 21-Nov-2000

using namespace std;
#include <iostream>
#include <fstream>
#include <cstdlib>                      // exit()
//#include <cstdio>
#include <netinet/in.h>                 // ntohl byteorder x86-HP unix
typedef unsigned long ulong;
typedef unsigned char uchar;

// changed reading get,read,put, functions for g++v3.2



/* ************************************************************ */
#define RAS_MAGIC       0x59a66a95
struct header
  {
  ulong   magic;              /* magic number */
  ulong   width;              /* width (pixels) of image */
  ulong   height;             /* height (pixels) of image */
  ulong   depth;              /* depth (1, 8, or 24 bits) of pixel */
  ulong   length;             /* length (bytes) of image */
  ulong   type;               /* type of file; see RT_* below */
  ulong   maptype;            /* type of colormap; see RMT_* below */
  ulong   maplength;          /* length (bytes) of following map */
  /* color map follows for ras_maplength bytes, followed by image */
  };
        /* Sun supported ras_type's */
#define RT_OLD          0       /* Raw pixrect image in 68000 byte order */
#define RT_STANDARD     1       /* Raw pixrect image in 68000 byte order */
#define RT_BYTE_ENCODED 2       /* Run-length compression of bytes */
#define RT_FORMAT_RGB   3       /* XRGB or RGB instead of XBGR or BGR */
#define RT_FORMAT_TIFF  4       /* tiff <-> standard rasterfile */
#define RT_FORMAT_IFF   5       /* iff (TAAC format) <-> standard rasterfile */

/* Sun registered ras_maptype's */
#define RMT_RAW         2
/* Sun supported ras_maptype's */
#define RMT_NONE        0       /* ras_maplength is expected to be 0 */
#define RMT_EQUAL_RGB   1       /* red[ras_maplength/3],green[],blue[] */

typedef struct {
        int              type;
        int              length;
        unsigned char   *map[3];
} colormap;


/* ************************************************************ */
int main (int argc,char *argv[])
  {
  if (argc<2)
    {
    cerr << argv[0] << " -- display header of SUN raster file.\n";
    cerr << "usage:   " << argv[0] << " rasfile\n";
    cerr << "example: " << argv[0] << " file.ras\n";
    cerr << "ex2:     " << argv[0] << " file.ras |& head -n 9\n";
    exit(1);
    }

  // ______ Open file ______
  ifstream ifp(argv[1]);
  if (!ifp) {cerr << "infile " << argv[1] << " does not exist\n"; exit(1);}
  
  // ______ Read HEADER 32 bytes and convert to big endian ______
  header HEADER;
  //cerr << "1: " << &HEADER << endl;
  //cerr << "1: " << (uchar *) &HEADER << endl;
  ifp.read(((char *) &HEADER), sizeof(HEADER));
//  for(int i=0; i<int(sizeof(HEADER)); ++i)
//    {
//    char tmp;
//    ifp.get(tmp);
//    HEADER[i] = uchar(tmp);
//    // orig: ifp.get(*((unsigned char*)(&HEADER) + i));
//    }
  //
  HEADER.magic     = ntohl((unsigned int)HEADER.magic);
  HEADER.width     = ntohl((unsigned int)HEADER.width);
  HEADER.height    = ntohl((unsigned int)HEADER.height);
  HEADER.depth     = ntohl((unsigned int)HEADER.depth);
  HEADER.length    = ntohl((unsigned int)HEADER.length);
  HEADER.type      = ntohl((unsigned int)HEADER.type);
  HEADER.maptype   = ntohl((unsigned int)HEADER.maptype);
  HEADER.maplength = ntohl((unsigned int)HEADER.maplength);

  cerr << "Header for file: " << argv[1]
       << "\nmagic#:    " << HEADER.magic
       << "\nwidth:     " << HEADER.width
       << "\nheight:    " << HEADER.height
       << "\ndepth:     " << HEADER.depth
       << "\nlength:    " << HEADER.length
       << "\ntype:      " << HEADER.type
       << "\nmaptype:   " << HEADER.maptype
       << "\nmaplength: " << HEADER.maplength
       << endl;

  // ______ Colormap, which type/how to output ??? ______
  // original i had: unsigned char CMAP[3][HEADER.maplength/3];
  // but if.get not for unsigned? so cast or read differently?
  unsigned char CMAP[3][HEADER.maplength/3];
  for(int color=0; color<3; ++color)
    {
    for(int i=0; i<int(HEADER.maplength/3); ++i)
      {
      ifp.get(*((char*)(&CMAP[color][i])));
      //char tmp;
      //ifp.get(tmp);
      //CMAP[color][i] = uchar(tmp);// is this the same?
      }
      //ifp.get(*((unsigned char*)(&CMAP[color][i])));
    }

  // ______ Write cmap to screen ______
  for(int i=0; i<int(HEADER.maplength/3); ++i)
    cerr << "rgb: " << int(CMAP[0][i]) << " " 
	            << int(CMAP[1][i]) << " " 
                    << int(CMAP[2][i]) << "\n";

  // ______ Tidy tidy ______
  ifp.close();
  return 0;
  }
