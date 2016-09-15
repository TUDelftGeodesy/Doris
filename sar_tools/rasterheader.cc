/* gcc -O rasterheader.c  -o rasterheader
 prog to create a header...
*/

using namespace std;
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <netinet/in.h>                 // ntohl byteorder x86-HP unix

int main (int argc,char *argv[])
{
  FILE *ofp;
  int c;
  //int bflg, aflg, errflg;
  extern char *optarg;
  //extern int optind, optopt;

  int header[8];
  int width,height,depth,datatype,length,maplength,maptype;

  int wrong=0;
  int computelength=1;

  /* defaults */
  width       = 0;
  height      = 0;
  depth       = 8;
  length      = 0;
  datatype    = 0;
  maplength   = 0;
  maptype     = 0;

  while ((c = getopt(argc, argv, "w:h:d:t:l:L:m:")) != -1)
    {
    switch (c)
      {
      case 'w':
                       width  = atoi(optarg);
                       break;
      case 'h':
                       height = atoi(optarg);
                       break;
      case 'd':
                       depth  = atoi(optarg);
                       break;
      case 't':
                       datatype  = atoi(optarg);
                       break;
      case 'l':
                       maplength = atoi(optarg);
                       break;
      case 'L':
                       length = atoi(optarg);
                       computelength = 0;
                       break;
      case 'm':
                       maptype = atoi(optarg);
                       break;
      case '?':
		       wrong = 1;
                       break;
    }
  }
  if ( width  == 0 ) wrong=1;
  if ( height == 0 ) wrong=1;
  /* check input */
  if ( wrong == 1 ) 
    {
    fprintf(stderr, "PROG: write header for sunraster file\n");
    fprintf(stderr, "  can be used together with cat.\n");
    fprintf(stderr, "  e.g. 1) make header\n");
    fprintf(stderr, "       2) cat header datafile > sunrasterfile.ras\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SYNOPSIS: prog -w width -h height [-ddepth]\n");
    fprintf(stderr, "               [-Ldatalength] [-tdatatype] [-mmaptype]\n");
    fprintf(stderr, "               [-lmaplength]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "USAGE (a SUNraster header is defined as 8 integers):\n");
    fprintf(stderr, "    ras_magic:     [0x59a66a95] magic number\n");
    fprintf(stderr, " -w ras_width:     [0]          width (pixels) of image\n");
    fprintf(stderr, " -h ras_height:    [0]          height (pixels) of image\n");
    fprintf(stderr, " -d ras_depth:     [8]          depth (1, 8, or 24 bits) of pixel\n");
    fprintf(stderr, " -L ras_length:    [computed]   length (bytes) of image\n");
    fprintf(stderr, " -t ras_type:      [0]          type of file; see RT_* below\n");
    fprintf(stderr, " -m ras_maptype:   [0]          type of colormap; see RMT_* below\n");
    fprintf(stderr, " -l ras_maplength: [0]          length (bytes) of following map\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "   /* Sun supported ras_type's RT_ */\n");
    fprintf(stderr, "0       Raw pixrect image in 68000 byte order \n");
    fprintf(stderr, "1       Raw pixrect image in 68000 byte order \n");
    fprintf(stderr, "2       Run-length compression of bytes \n");
    fprintf(stderr, "3       XRGB or RGB instead of XBGR or BGR \n");
    fprintf(stderr, "4       tiff <-> standard rasterfile \n");
    fprintf(stderr, "5       iff (TAAC format) <-> standard rasterfile \n");
    fprintf(stderr, "0xffff  Reserved for testing \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   /* Sun registered ras_maptype's RMT_ */ \n");
    fprintf(stderr, "0       ras_maplength is expected to be 0 \n");
    fprintf(stderr, "1       red[ras_maplength/3],green[],blue[] \n");
    fprintf(stderr, "2       RMT_RAW \n");
    fprintf(stderr, "\n");

    exit(1);
    }

  if (computelength==1) length = width * height * (int)(depth/8);
  fprintf(stderr, "\n");
  fprintf(stderr, "arguments used:");
  fprintf(stderr, "\n  width:     %i", width);
  fprintf(stderr, "\n  height:    %i", height);
  fprintf(stderr, "\n  depth:     %i", depth);
  fprintf(stderr, "\n  length:    %i", length);
  fprintf(stderr, "\n  datatype:  %i", datatype);
  fprintf(stderr, "\n  maptype:   %i", maptype);
  fprintf(stderr, "\n  maplength: %i", maplength);
  fprintf(stderr, "\n");

  header[0] = ntohl(0x59a66a95);	/* magic number, always same(?) */
  header[1] = ntohl(width);		/* #columns */
  header[2] = ntohl(height);		/* #rows */
  header[3] = ntohl(depth);		/* in bits  email s suchandt */
  header[4] = ntohl(length);		/* in bytes */
  header[5] = ntohl(datatype);		/* email s suchandt */
  header[6] = ntohl(maptype);		/* ... */
  header[7] = ntohl(maplength);		/* ... */

  /* Write SUN raster file */
  ofp=fopen("sunrasterheader","w");
  fwrite(header, sizeof(int), 8, ofp);
  fclose(ofp);
  free(header);
  return 0;
}



