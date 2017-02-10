#!/usr/bin/env python
import os,sys,time

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

import sentinel_1.main_code.resdata as resdata


def dump_data(input_file,res_file, output_file='', coordinates=[]):
    # This function dumps a .raw file from the original .tif sentinel data. The input_file is the .tif file and the
    # res_file is the .res file that corresponds with the output file. Coordinates is an optional variable which can be
    # called if the lines and pixels are not yet defined in the .res file.

    res_vars = resdata.ResData(filename=res_file)
    res_vars.res_read()

    # Check if information about crop is available

    if not coordinates:
        if res_vars.process_control['crop'] == '0':
            print 'There is no information available about how to crop this file!'
            return
        else:
            outputWinFirstPix = int(res_vars.processes['crop']['First_pixel (w.r.t. original_image)'])
            outputWinLastPix = int(res_vars.processes['crop']['Last_pixel (w.r.t. original_image)'])
            outputWinFirstLine = int(res_vars.processes['crop']['First_line (w.r.t. tiff_image)'])
            outputWinLastLine = int(res_vars.processes['crop']['Last_line (w.r.t. tiff_image)'])

    else:
        outputWinFirstPix = coordinates[0]
        outputWinLastPix = coordinates[1]
        outputWinFirstLine = coordinates[2]
        outputWinLastLine = coordinates[3]

    if not output_file:
        if res_vars.process_control['crop'] == '1':
            if 'Data_output_file' in res_vars.processes['crop'].keys():
                output_file = os.path.join(os.path.basename(res_file), res_vars.processes['crop']['Data_output_file'])
        if not output_file:
            output_file = res_file.split(".")[0] + '.raw'

    # system parameters :: check whether its there

    gdalCall = 'gdal_translate'

    # Output data parameters
    outputDataFormat = 'MFF'
    outputDataType   = 'CInt16'

    cmd = '%s %s -ot %s -of %s %s' % (gdalCall,input_file,outputDataType,outputDataFormat,output_file)

    if outputWinFirstPix is not None:
        cmd = cmd + (' -srcwin %s %s %s %s' % (outputWinFirstPix,outputWinFirstLine,outputWinLastPix-outputWinFirstPix+1,outputWinLastLine-outputWinFirstLine+1))

    failure = os.system(cmd)
    if failure:
        print '%s: running %s failed' % (sys.argv[0],cmd)
        sys.exit(1)
    else:
        os.rename(os.path.splitext(output_file)[0]+'.j00',output_file)
        os.remove(os.path.splitext(output_file)[0]+'.hdr')
        os.remove(os.path.splitext(output_file)[0]+'.hdr.aux.xml')


# Actually execute the code to unzip one data file.
if __name__ == "__main__":

    input_file = sys.argv[1]
    res_file = sys.argv[2]

    dump_data(input_file, res_file, output_file='', coordinates=[])


