#!/usr/bin/env python
import sys
import os
import numpy as np
from numpy import *
import gdal
from gdalconst import *
# from resdata import ResData
from get_ramp import get_ramp, freadbk, get_parameter

def usage():
    print '\nUsage: python  do_reramp_SLC.py dataFilename resFilename resampled'
    print '  where dataFilename        is the name of burst you want to deramp'
    print '        resFilename         is the .res file of burst              ' 
    print '        plotFlag            is a boolean var, to plot only         '
    print '                            default of doPlot is false             '
    print ' This python applies the inverse phase ramp to the burst pointed by DATAFILENAME (slc)'
    print ' and RESFILENAME (res) that was deramped by deramp_SLC.m. The phase screen'
    print ' must account for the new resampled grids PIXRGGRID and PIXAZGRID    '
    print ' [Nlines_mst x Nsamples_mst] that contain the time coordinates of the'
    print ' resampled image into the master grid:                               '
    print '  for example                                                        '
    print ' python   do_reramp_SLC.py slave_rsmp.raw slave.res False            '
    print ' matlab: Lorenzo Iannini'
    print ' Python: Wu Wenhao   Wuhan university   QQ:460249274'
try:
    dataFilename         = sys.argv[1] 
    resFilename          = sys.argv[2]
     
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
if len(sys.argv) == 3:
    resampled  = True
elif len(sys.argv) == 4:
    resampled  = sys.argv[3]
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1) 

##################################################################################

Link_ORIG_FILE=dataFilename+'.orig'

if (os.path.isfile(Link_ORIG_FILE)):
    os.rename(Link_ORIG_FILE,dataFilename)
os.rename(dataFilename,Link_ORIG_FILE)

# Read information
################################################################################

rsmpDataFormat = get_parameter('Data_output_format',resFilename,2,'*_Start_resample','* End_resample:_NORMAL')
print "rsmpDataFormat =",rsmpDataFormat

if rsmpDataFormat=='complex_real4':
    dataFormat = 'cpxfloat32'
    print 'dataFormat =',dataFormat
else:
    dataFormat = 'cpxint16'

# Number of lines
lNum = int(get_parameter('Number_of_lines_original',resFilename,1))

l0 = int(get_parameter('First_line (w.r.t. original_master)',resFilename,2,'*_Start_resample','* End_resample:_NORMAL'))
lN = int(get_parameter('Last_line (w.r.t. original_master)',resFilename,2,'*_Start_resample','* End_resample:_NORMAL'))
p0 = int(get_parameter('First_pixel (w.r.t. original_master)',resFilename,2,'*_Start_resample','* End_resample:_NORMAL'))
pN = int(get_parameter('Last_pixel (w.r.t. original_master)',resFilename,2,'*_Start_resample','* End_resample:_NORMAL'))

# Get resampled Slv size
Naz_res = lN-l0+1
Nrg_res = pN-p0+1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Obtain chirp %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if resampled == False:
    ChirpFilt = get_ramp(resFilename, resampled=0, type='chirp')
if resampled == True:
    ChirpFilt = get_ramp(resFilename, resampled=1, type='chirp')

################################################################################
# Read data

Path_MFF_HDR   =Link_ORIG_FILE[:-5] + '.hdr'

if dataFormat=='cpxfloat32':
    Link_DATA      =Link_ORIG_FILE[:-5] + '.x00'
else:
    Link_DATA      =Link_ORIG_FILE[:-5] + '.j00'

if (os.path.isfile(Path_MFF_HDR)):
    os.remove(Path_MFF_HDR)
if (os.path.isfile(Link_DATA)):
    os.remove(Link_DATA)

RAW_DATA_ABSOLUTE_PATH=os.path.abspath(Link_ORIG_FILE)
print "RAW_DATA_ABSOLUTE_PATH=", RAW_DATA_ABSOLUTE_PATH
os.symlink(RAW_DATA_ABSOLUTE_PATH,Link_DATA)

outStream      = open(Path_MFF_HDR,'w')
outStream.write('IMAGE_FILE_FORMAT = MFF\n')
outStream.write('FILE_TYPE = IMAGE\n')
outStream.write('IMAGE_LINES = %d\n' % int(Naz_res))
outStream.write('LINE_SAMPLES = %d\n'% int(Nrg_res))
outStream.write('BYTE_ORDER = LSB\n')
outStream.write('END\n')
outStream.close()

slc       = freadbk(Path_MFF_HDR,1, 1,int(Naz_res),int(Nrg_res))
slc       = slc.astype(complex128)

if (os.path.isfile(Path_MFF_HDR)):
    os.remove(Path_MFF_HDR)
if (os.path.isfile(Link_DATA)):
    os.remove(Link_DATA)
#######################################################################################

#%% Apply reramping
slc_reramped = slc*ChirpFilt
slc_reramped = slc_reramped.astype(np.complex64)

print "\nStart Write (save) reramped SLC !!!\n"
print "The data type is ",slc_reramped.dtype

cols = slc_reramped.shape[1]
rows = slc_reramped.shape[0]
print "cols = ",cols
print "rows = ",rows

fid = open(dataFilename,'wb')
if rsmpDataFormat=='complex_real4':
    for i_temp in arange(0,rows):   
        fid.write(slc_reramped[i_temp,:])
else:
     save_row_sl_deramped=zeros((2*cols,1),dtype='int16')
     for i_temp in arange(0,rows):        
         save_row_sl_deramped[0::2]= (slc_reramped[i_temp,:].real.astype(np.int16)).reshape(-1,1,order='F').copy()
         save_row_sl_deramped[1::2]= (slc_reramped[i_temp,:].imag.astype(np.int16)).reshape(-1,1,order='F').copy()
         fid.write(save_row_sl_deramped)
     if save_row_sl_deramped.dtype=='int16':
         print '\n Good ! data format is right!!!!\n'
     else:
         print 'Warning !The data format may be wrong !'
     del save_row_sl_deramped

fid.close()
print "\nThe reramp operation completed!\n"
