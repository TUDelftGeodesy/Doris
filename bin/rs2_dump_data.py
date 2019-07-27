#!/usr/bin/env python

#-----------------------------------------------------------------#
# A python code for extracting image matrix from RS2 geotiff file 
# using gdal_translate
#
# Author: TUDelft - 2010
# Maintainer: Mahmut Arikan
#
# Developed based on tsx_dump_data.py code.
#
#-----------------------------------------------------------------#
import os,sys,time
from os import path

import xml.etree.ElementTree as etree  # parameters required for cropping and flipping


codeRevision=1.1   # this code revision number
def usage():
    print('INFO    : @(#)Doris InSAR software, $Revision: %s $, $Author: TUDelft $' % codeRevision)
    print()
    print('Usage   : rs2_dump_data.py <inputfile> <outputfile> [l0 lN p0 pN] [-res RESFILE]')
    print()
    print('          inputfile         is the input Radarsat-2 geotiff filename : master.tif')
    print('          outputfile        is the output filename                   : master.slc')
    print('          l0                is the first azimuth line (starting at 1)')
    print('          lN                is the last azimuth line')
    print('          p0                is the first range pixel (starting at 1)')
    print('          pN                is the last range pixel')
    print('          RESFILE           DORIS result file that is to be updated for crop metadata (optional)')
    print()
    print('          This software is part of Doris InSAR software package.\n')
    print('(c) 1999-2010 Delft University of Technology, the Netherlands.\n')

try:
    inputFileName  = sys.argv[1] 
    outputFileName = sys.argv[2]
except:
    print('\nError   : Unrecognized input or missing arguments\n\n')
    usage()
    sys.exit(1)

inputFileName    = sys.argv[1]
outputFileName   = sys.argv[2]

resFile = None
for element in range(len(sys.argv)):
    option = sys.argv[element];
    if option == '-res':
        resFile = str(sys.argv[element+1])
#        print(resFile)
        del sys.argv[element+1]
        del sys.argv[element]
        break


if len(sys.argv) == 3:
    outputWinFirstLine = None
    outputWinLastLine  = None
    outputWinFirstPix = None
    outputWinLastPix  = None
elif len(sys.argv) > 3 and len(sys.argv) < 9:
    outputWinFirstLine = int(sys.argv[3])    # Firstline
    outputWinLastLine  = int(sys.argv[4])    # Lastline
    outputWinFirstPix  = int(sys.argv[5])    # Firstpix
    outputWinLastPix   = int(sys.argv[6])    # Lastpix
elif len(sys.argv) > 3 and len(sys.argv) < 7:
    print('\nError   : Unrecognized input or missing arguments\n\n')
    usage()
    sys.exit(1)
else:
    outputWinFirstLine = None
    outputWinLastLine  = None
    outputWinFirstPix = None
    outputWinLastPix  = None



# Extract some important parameters for coordinate manipulation
productfile  = path.join(path.dirname(inputFileName), 'product.xml')
tree         = etree.parse(productfile) # 
NS           = 'http://www.rsi.ca/rs2/prod/xml/schemas'
PASS         = tree.find('.//{%s}passDirection'     % NS).text     # Ascending or Descending
line_order   = tree.find('.//{%s}lineTimeOrdering'  % NS).text
pixl_order   = tree.find('.//{%s}pixelTimeOrdering' % NS).text # Increasing or Decreasing
Nlines       = int(tree.find('.//{%s}numberOfLines'          % NS).text)
Npixls       = int(tree.find('.//{%s}numberOfSamplesPerLine' % NS).text)
in_FirstLine = outputWinFirstLine
in_LastLine  = outputWinLastLine
in_FirstPix  = outputWinFirstPix
in_LastPix   = outputWinLastPix

if line_order == 'Increasing':
    print('INFO     : Detected a imagery %s line time order'  % line_order) 
else:
    print('INFO     : Detected a imagery %s line time order'  % line_order) 
    print('INFO     : Adjusting first line and last line ') 
    outputWinFirstLine = Nlines-in_LastLine+1
    outputWinLastLine  = Nlines-in_FirstLine+1
    print(Nlines, Npixls)
    print(outputWinFirstLine,outputWinLastLine,outputWinFirstPix,outputWinLastPix)

if pixl_order == 'Increasing':
    print('INFO     : Detected a imagery %s pixel time order' % pixl_order) 
else:
    print('INFO     : Detected a imagery %s pixel time order' % pixl_order) 
    print('INFO     : Adjusting first pixel and last pixel ') 
    outputWinFirstPix = Npixls-in_LastPix+1
    outputWinLastPix  = Npixls-in_FirstPix+1
    print(Nlines, Npixls)
    print(outputWinFirstLine,outputWinLastLine,outputWinFirstPix,outputWinLastPix)


# GDAL Extract image matrix using gdal_translate
# Ex: gdal_translate -of ENVI -ot Int16 -co INTERLEAVE=BIP imagery_HH.tif rs2_slc.raw
gdalCall = 'gdal_translate'

# GDAL parameters
outputDataFormat = 'ENVI'
#outputDataType   = 'CInt16'
outputDataType   = 'Int16'
outputDataType2   = 'ci2'  # for cpxfiddle
outputInterleave = 'INTERLEAVE=BIP'

cmd = '%s %s -ot %s -of %s -co %s %s' % (gdalCall,inputFileName,outputDataType,outputDataFormat,outputInterleave,outputFileName+'.noflip')

if outputWinFirstPix is not None:
    #                                         xoff                  yoff                  xsize = width                     ysize= height
    #                                        1 --> 0               1 --> 0
    cmd = cmd + (' -srcwin %s %s %s %s' % (outputWinFirstPix-1,outputWinFirstLine-1,outputWinLastPix-outputWinFirstPix+1,outputWinLastLine-outputWinFirstLine+1))
    print(cmd)

failure = os.system(cmd)
if failure:
    print('%s: running %s failed' % (sys.argv[0],cmd))
    sys.exit(1)
#else:
#    os.rename(os.path.splitext(outputFileName)[0]+'.j00',outputFileName)



print('\nINFO     : Start flipping image to radar coordinate for a %s orbit.' % PASS)

# Flip using cpxfiddle
# get keyword asc or desc from .res
#  if resFile is not None:
#     inStream = open(resFile,'r')
#     find ascending tag then
# -m Y for ascending
# -m X for descending orbit

cpxfcall='cpxfiddle -q normal -o short '              # cpxfiddle -w 3788 -f ci2 -q normal -o short -m Y rs2_slc.raw.noflip  > rs2_slc.raw 

if PASS == 'Ascending':
    cmd = cpxfcall + ( '-w %s -f %s -m %s %s > %s ' % (outputWinLastPix-outputWinFirstPix+1,outputDataType2,'Y',outputFileName+'.noflip',outputFileName)) 
    print('INFO     : %s ' % cmd) 
else:
    cmd = cpxfcall + ( '-w %s -f %s -m %s %s > %s ' % (outputWinLastPix-outputWinFirstPix+1,outputDataType2,'X',outputFileName+'.noflip',outputFileName)) 
    print('INFO     : %s ' % cmd) 

failure = os.system(cmd)
if failure:
    print('%s: running %s failed' % (sys.argv[0],cmd))
    sys.exit(1)
else:
    os.remove(outputFileName+'.noflip')



# Manual call to update Doris res file
# check whether the file exist!!!
if resFile is not None:

    print(resFile)

    # load header
    #headerFileStream = open(os.path.splitext(outputFileName)[0]+'.hdr','r')
    headerFileStream = open(outputFileName+'.hdr','r')
    for line in headerFileStream:
        pair = line.split()
        if len(pair) > 1:
            vars()[pair[0]] = pair[2]   # set IMAGE_LINES and LINE_SAMPLES
#            print(vars()[pair[0]])

    # check whether the file exist
    outStream = open(resFile,'a')

    outStream.write('\n')
    outStream.write('*******************************************************************\n')
    outStream.write('*_Start_crop:			geotiff\n')
    outStream.write('*******************************************************************\n')
    outStream.write('Data_output_file:                       %s\n' % outputFileName)
    outStream.write('Data_output_format: 			complex_short\n') # HARDCODED

    if outputWinFirstPix is not None:
        outStream.write('First_line (w.r.t. original_image): 	%s\n'   % in_FirstLine)   
        outStream.write('Last_line (w.r.t. original_image): 	%s\n'   % in_LastLine)
        outStream.write('First_pixel (w.r.t. original_image): 	%s\n' % in_FirstPix)
        outStream.write('Last_pixel (w.r.t. original_image): 	%s\n'   % in_LastPix)

    elif outputWinFirstPix is None:
        outStream.write('First_line (w.r.t. original_image): 	1\n')
        outStream.write('Last_line (w.r.t. original_image): 	%s\n'   % lines )
        outStream.write('First_pixel (w.r.t. original_image): 	1\n')
        outStream.write('Last_pixel (w.r.t. original_image): 	%s\n'   % samples )

    outStream.write('*******************************************************************\n')
    outStream.write('* End_crop:_NORMAL\n')
    outStream.write('*******************************************************************\n')

    outStream.write('\n')
    outStream.write('    Current time: %s\n' % time.asctime())
    outStream.write('\n')

    # close output file
    outStream.close()

    # replace crop tag in result file
    sourceText   = "crop:			0"
    replaceText  = "crop:			1"
    inputStream  = open(resFile,'r')
    textStream   = inputStream.read()
    inputStream.close()
    outputStream = open(resFile, "w")
    outputStream.write(textStream.replace(sourceText, replaceText))
    outputStream.close()


#EOF
