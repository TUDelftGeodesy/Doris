/* BK: simple program to dump binary data for ENVISAT 
   for alternating pol, VV channel.
   #%// BK 24-Nov-2003
   $Id: envisat_dump_VV.c,v 1.2 2004/05/13 18:13:39 kampes Exp $

THIS program is not really good, but a fast hack for Z.Perksi..
*/
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "epr_api-2.2/src/epr_api.h"
#if defined(WIN32) && defined(_DEBUG)
#include <crtdbg.h>
#endif /* if defined(WIN32) && defined(_DEBUG) */


int main(int argc, char** argv) 
  {
  const char*     product_file_path;
  const char*     outfile;
  const char*     record_name;
  FILE*           outstream;
  EPR_SProductId* product_id;
  EPR_SDatasetId* MAIN_PROC_PM_ID;
  EPR_SDatasetId* MDS1;
  EPR_SRecord*    rec1;
  EPR_SRecord*    rec5;
  //EPR_SField*     numlines_field;
  //EPR_SField*     numpixels_field;
  //EPR_SField*     line_field;
  EPR_SField     numlines_field;
  EPR_SField     numpixels_field;
  EPR_SField     line_field;
  //EPR_SField*     line_num_field; // not used
  int             status;
  ulong           line_num;
  ulong           numlines;
  ulong           numberoflines;
  ulong           numberoffields;
  ulong           numpixels;
  //ulong           l0;
  //ulong           lN;
  //ulong           p0;
  //ulong           pN;
  unsigned int           l0;  // MA due to difference between 32-bit/64-bit data structure: define typedef like Doris see constant.h
  unsigned int           lN;
  unsigned int           p0;
  unsigned int           pN;
  //ulong           cnt;
  //ulong           x,y;  /* loop counter go upto 25.000 or so */
  int             cnt;
  int             x,y;  /* loop counter go upto 25.000 or so */

  short           realpart_short;/* written to output file */
  short           imagpart_short;/* written to output file */

  /* handle input */
  printf("argc: %f\n", (float)argc);
  if (argc != 3 && argc !=7 && argc !=8)
    {
    printf("Usage: envisat_dump_VV envisat-product outputfile [l0 lN p0 pN]\n");
    printf("  where envisat-product is the input filename of the alternating pol product\n");
    printf("        outputfile      is the output filename\n");
    printf("        l0              is the first azimuth line (starting at 1)\n");
    printf("        lN              is the last azimuth line\n");
    printf("        p0              is the first range pixel (starting at 1)\n");
    printf("        pN              is the last range pixel\n");
    printf("Example:\n");
    printf("  envisat_dump_VV ASA_IMS_1PNDPA20021025_175208_000000162010_00356_03416_0005.N1 crop.out 0 10 0 100\n\n");
    exit(1);
    }
  product_file_path = argv[1];
  outfile           = argv[2];
  printf("infile:  %s\n",product_file_path);
  printf("outfile: %s\n",outfile);
  if (argc>=7)
    {
    status=sscanf(argv[3],"%u",&l0);
    status=sscanf(argv[4],"%u",&lN);
    status=sscanf(argv[5],"%u",&p0);
    status=sscanf(argv[6],"%u",&pN);
    printf("sscanf l0: %f\n", (float)l0);
    printf("sscanf lN: %f\n", (float)lN);
    printf("sscanf p0: %f\n", (float)p0);
    printf("sscanf pN: %f\n", (float)pN);
    }
  /*
  record_name = (argc==8) ? argv[8] : "MDS1";
  printf("record name: %s\n",record_name);
  */

  /* Initialize the API. Set log-level to DEBUG and use default log-output (stdout) */
  epr_init_api(e_log_debug, epr_log_message, NULL);

  /* Open the product; an argument is a path to product data file */
  /* PRODUCT dataset record field element */
  product_id      = epr_open_product(product_file_path);
  /* product DATASET record field element */
  MAIN_PROC_PM_ID = epr_get_dataset_id(product_id, "MAIN_PROCESSING_PARAMS_ADS");
  /* see product specs, record is called MDS2! */
  MDS1            = epr_get_dataset_id(product_id, "MDS2");
  /* product dataset RECORD field element */
  rec1 = epr_read_record(MAIN_PROC_PM_ID,         0, NULL);
  /* product dataset record FIELD element */
  //numlines_field  = epr_get_field(rec1, "num_output_lines");
  //numpixels_field = epr_get_field(rec1, "num_samples_per_line");
  numlines_field  = *(epr_get_field(rec1, "num_output_lines"));
  numpixels_field = *(epr_get_field(rec1, "num_samples_per_line"));
  /*
  epr_free_record(rec1);
  */
  //epr_print_field(numlines_field, stdout);
  //epr_print_field(numpixels_field, stdout);
  epr_print_field(&numlines_field, stdout);
  epr_print_field(&numpixels_field, stdout);
  //numlines        = epr_get_field_elem_as_ulong(numlines_field, 0);
  //numpixels       = epr_get_field_elem_as_ulong(numpixels_field, 0);
  numlines        = epr_get_field_elem_as_uint(&numlines_field, 0);
  numpixels       = epr_get_field_elem_as_uint(&numpixels_field, 0);
  if (argc == 3)
    {
    l0 = 1;
    lN = numlines;
    p0 = 1;
    pN = numpixels;
    }
  /* loop over data record to get data and dump it to file */
  numberoflines = epr_get_num_records(MDS1);
  printf("numberoflines: %f\n", (float)numberoflines);
  printf("numlines: %f\n", (float)numlines);
  printf("numpixels: %f\n", (float)numpixels);
  printf("l0: %f\n", (float)l0);
  printf("lN: %f\n", (float)lN);
  printf("p0: %f\n", (float)p0);
  printf("pN: %f\n", (float)pN);
  /* check if number of records indeed equals the number of lines */
  if (numberoflines != numlines)
    {
    printf("numlines not equal in check, ASAR format error?.");
    exit(1);
    }


  /* --- Check if input as acceptable ---------------------------- */
  if (l0 < 1)  
    {
    printf("l0<1 not allowed.  first line is 1 not 0.\n");
    exit(1);
    }
  if (p0 < 1)  
    {
    printf("p0<1 not allowed.  first line is 1 not 0.\n");
    exit(1);
    }
  if (lN > numlines)  
    {
    printf("lN>numlines not allowed.\n");
    exit(1);
    }
  if (pN > numpixels)  
    {
    printf("pN>numpixels not allowed.\n");
    exit(1);
    }


  /* read in whole line of cpx data */
  outstream = fopen(outfile,"wb");
  for (y=l0;y<=lN;y++)
    {
    rec5       = epr_read_record(MDS1, y-1, NULL);
    /*
    numberoffields = epr_get_num_fields(rec5);
    line_num_field = epr_get_field(rec5, "line_num");
    line_num       = epr_get_field_elem_as_ulong(line_num_field,0);
    printf("numberoffields: %i\n", numberoffields);
    printf("line number: %i\n", line_num);
    */

    //line_field = epr_get_field(rec5, "proc_data");
    line_field = *(epr_get_field(rec5, "proc_data"));
    cnt = 2*(p0-1);/* p0 starts at 1 for first element */
    /* write out selected pixels */
    for (x=p0;x<=pN;x++)
      {
     // realpart_short = epr_get_field_elem_as_short(line_field,cnt);
     // cnt = cnt + 1;
     // imagpart_short = epr_get_field_elem_as_short(line_field,cnt);
     // cnt = cnt + 1;
      realpart_short = epr_get_field_elem_as_short(&line_field,cnt);
      cnt++;
      imagpart_short = epr_get_field_elem_as_short(&line_field,cnt);
      cnt++;
/*
printf("%i,%i: realpart: %f  %f\n", x,y,(float)realpart_short,(float)imagpart_short);
*/
      status = fwrite(&realpart_short,2,1,outstream);
      if (status != 1) fprintf(stderr,"fwrite could not write to disk?");
      status = fwrite(&imagpart_short,2,1,outstream);
      if (status != 1) fprintf(stderr,"fwrite could not write to disk?");
      }
    /* this program seemed to fill all memory for some reason?  try to free it */
    epr_free_record(rec5);
    }
  fclose(outstream);
  epr_close_product(product_id);
  /* Closes product reader API, release all allocated resources */
  epr_close_api();
  return 0;
  }


