#!/usr/bin/env python
import os,sys,time

def usage():
    print '\nUsage: python tsx_dump_data.py tsx_COSAR_product outputfile [l0 lN p0 pN] -res RESFILE'
    print '  where tsx_COSAR_product is the input filename'
    print '        outputfile        is the output filename'
    print '        l0                is the first azimuth line (starting at 1)'
    print '        lN                is the last azimuth line'
    print '        p0                is the first range pixel (starting at 1)'
    print '        pN                is the last range pixel'
    print '        RESFILE           DORIS result file that is to be updated for crop metadata'

try:
    inputFileName  = sys.argv[1]; 
    outputFileName = sys.argv[2]
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)

inputFileName    = sys.argv[1]
outputFileName   = sys.argv[2]

resFile = None
for element in range(len(sys.argv)):
    option = sys.argv[element];
    if option == '-res':
        resFile = str(sys.argv[element+1])
#        print resFile
        del sys.argv[element+1]
        del sys.argv[element]
        break

if len(sys.argv) == 3:
    outputWinFirstLine = None
    outputWinLastLine  = None
    outputWinFirstPix = None
    outputWinLastPix  = None
elif len(sys.argv) > 3 and len(sys.argv) < 9:
    outputWinFirstLine = int(sys.argv[3])-1    # gdal srcwin starting at 0
    outputWinLastLine  = int(sys.argv[4])      # Lastline --> yoff  (later)
    outputWinFirstPix  = int(sys.argv[5])-1    # gdal srcwin starting at 0
    outputWinLastPix   = int(sys.argv[6])      # Lastpix  --> yoff  (later)
elif len(sys.argv) > 3 and len(sys.argv) < 7:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
else:
    outputWinFirstLine = None
    outputWinLastLine  = None
    outputWinFirstPix = None
    outputWinLastPix  = None

# system parameters :: check whether its there
gdalCall = 'gdal_translate'

# Output data parameters
outputDataFormat = 'MFF'
outputDataType   = 'CInt16'

cmd = '%s %s -ot %s -of %s %s' % (gdalCall,inputFileName,outputDataType,outputDataFormat,outputFileName)

if outputWinFirstPix is not None:
    cmd = cmd + (' -srcwin %s %s %s %s' % (outputWinFirstPix,outputWinFirstLine,outputWinLastPix-outputWinFirstPix,outputWinLastLine-outputWinFirstLine))
    #print cmd 

failure = os.system(cmd)
if failure:
    print '%s: running %s failed' % (sys.argv[0],cmd)
    sys.exit(1)
else:
    os.rename(os.path.splitext(outputFileName)[0]+'.j00',outputFileName)
# call to cpxfiddle

# check whether the file exist!!!
if resFile is not None:

    print resFile

    # load header
    headerFileStream = open(os.path.splitext(outputFileName)[0]+'.hdr','r')
    for line in headerFileStream:
        pair = line.split()
        if len(pair) > 1:
            vars()[pair[0]] = pair[2]   # set IMAGE_LINES and LINE_SAMPLES
#            print vars()[pair[0]]

    # check whether the file exist
    outStream = open(resFile,'a')

    outStream.write('\n')
    outStream.write('*******************************************************************\n')
    outStream.write('*_Start_crop:			cosar\n')
    outStream.write('*******************************************************************\n')
    outStream.write('Data_output_file: 	%s\n' % outputFileName)
    outStream.write('Data_output_format: 			complex_short\n') # hardcoded

    if outputWinFirstPix is not None:
        outStream.write('First_line (w.r.t. original_image): 	%s\n' % (outputWinFirstLine+1))   # back to Doris convention
        outStream.write('Last_line (w.r.t. original_image): 	%s\n' % outputWinLastLine)
        outStream.write('First_pixel (w.r.t. original_image): 	%s\n' % (outputWinFirstPix+1))
        outStream.write('Last_pixel (w.r.t. original_image): 	%s\n' % outputWinLastPix)

    elif outputWinFirstPix is None:
        outStream.write('First_line (w.r.t. original_image): 	1\n')
        outStream.write('Last_line (w.r.t. original_image): 	%s\n' % IMAGE_LINES)
        outStream.write('First_pixel (w.r.t. original_image): 	1\n')
        outStream.write('Last_pixel (w.r.t. original_image): 	%s\n' % LINE_SAMPLES)

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
