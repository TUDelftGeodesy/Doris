#!/usr/bin/env python
import os,sys,time
import numpy as np
from numpy import *

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

from sentinel_1.functions.get_ramp import get_ramp, freadbk, get_parameter


def usage():
    print '\nUsage: python do_deramp_SLC_nom.py dataFilename  resFilename plotFlag'
    print '  where dataFilename     is the name of burst you want to deramp'
    print '        resFilename      is the .res file of burst'
    print '        plotFlag         is a boolean var, to plot only'
    print '                         default of doPlot is false'
    print ' This function removes the phase ramp (Doppler centroid variations) from single burst of' 
    print ' RS2 or S1 TOPS acquisition. The original binary image at path '
    print " DATAFILENAME is saved in 'DATAFILENAME'.orig, whereas the new instance " 
    print ' will be characterized by baseband spectrum. The function also requires '
    print ' the .res file RESFILENAME.                                             '     
    print '  for example                                                           '
    print ' python  do_deramp_SLC.py   20140821_iw_2_burst_1.raw slave.res False   '
    print ' matlab : Lorenzo Iannini tudelft                                       '
    print ' Python: Wu Wenhao   Wuhan university   QQ:460249274                    '
try:
    dataFilename         = sys.argv[1] 
    resFilename          = sys.argv[2] 
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
if len(sys.argv) == 3:
    plotFlag  = False
elif len(sys.argv) == 4:
    plotFlag  = sys.argv[3]
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1) 

##################################################################################

#%% Make a backup of the original burst image
#% This is done in order not to overwrite the input data. DATAFILENAME at 
#% the end of the routine will in fact contain the new deramped burst

# DATAFILENAME is saved in 'DATAFILENAME'.orig, whereas the new instance 
# will be characterized by baseband spectrum. The function also requires 
# the .res file RESFILENAME.

Link_ORIG_FILE=dataFilename+'.orig'

if (os.path.isfile(Link_ORIG_FILE)):
    os.rename(Link_ORIG_FILE,dataFilename)
os.rename(dataFilename,Link_ORIG_FILE)

#*****************************************************************************#
# Calculate chirp for deramping
ChirpFilt = get_ramp(resFilename, resampled=0, type='chirp')

# Image size properties
if get_parameter('First_line (w.r.t. ovs_image)',resFilename,1): #oversampled data
    l0 = int(get_parameter('First_line (w.r.t. ovs_image)',resFilename,1))  
    lN = int(get_parameter('Last_line (w.r.t. ovs_image)',resFilename,1))  
    p0 = int(get_parameter('First_pixel (w.r.t. ovs_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. ovs_image)',resFilename,1))
    dataFormat = 'cpxfloat32'
else: #original data  
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resFilename,1))
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resFilename,1))
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resFilename,1))
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resFilename,1))
    dataFormat = 'cpxint16'

# Number of lines
lNum = int(get_parameter('Number_of_lines_original',resFilename,1))

# Image size
Naz = lN-l0+1
Nrg = pN-p0+1

#**********************************************************************#

Path_MFF_HDR   =dataFilename.split('.')[0]+'.hdr'

if dataFormat == 'cpxfloat32':
    Link_DATA  =dataFilename.split('.')[0]+'.x00'
elif dataFormat == 'cpxint16':
    Link_DATA  =dataFilename.split('.')[0]+'.j00'
else:
    print 'Unrecognized input'
    usage()
    sys.exit(1)
outStream      = open(Path_MFF_HDR,'w')
outStream.write('IMAGE_FILE_FORMAT = MFF\n')
outStream.write('FILE_TYPE = IMAGE\n')
outStream.write('IMAGE_LINES = %d\n' % int(Naz))
outStream.write('LINE_SAMPLES = %d\n'% int(Nrg ))
outStream.write('BYTE_ORDER = LSB\n')
outStream.write('END\n')
outStream.close()

if (os.path.isfile(Link_DATA)):
    os.remove(Link_DATA)
RAW_DATA_ABSOLUTE_PATH=os.path.abspath(Link_ORIG_FILE)
print "RAW_DATA_ABSOLUTE_PATH=", RAW_DATA_ABSOLUTE_PATH
os.symlink(RAW_DATA_ABSOLUTE_PATH,Link_DATA)

slc=freadbk(Path_MFF_HDR,1, 1,int(Naz),int(Nrg))
slc_deramped_nom = conj(ChirpFilt)*slc

##************************************************#
if (os.path.isfile(Link_DATA)):
    os.remove(Link_DATA)
if (os.path.isfile(Path_MFF_HDR)):
    os.remove(Path_MFF_HDR)
#******************************************************#

cols = slc_deramped_nom.shape[1]
rows = slc_deramped_nom.shape[0]

fid = open(dataFilename,'wb')
if dataFormat=='cpxfloat32':
    slc_deramped=slc_deramped.astype(np.complex64)  #It's not useful ,may be a bug
    if slc_deramped_nom.dtype=='complex64':
        print 'data format is right\n!!!!'
    else:
        print "Warning !The data format may be wrong !"
    for i_temp in arange(0,rows):
        fid.write(slc_deramped_nom[i_temp,:])
else:
    save_row_sl_deramped=zeros((2*cols,1),dtype='int16')
    for i_temp in arange(0,rows):
        save_row_sl_deramped[0::2]= (slc_deramped_nom[i_temp,:].real.astype(np.int16)).reshape(-1,1,order='F')
        save_row_sl_deramped[1::2]= (slc_deramped_nom[i_temp,:].imag.astype(np.int16)).reshape(-1,1,order='F')
        fid.write(save_row_sl_deramped)
    if save_row_sl_deramped.dtype=='int16':
        print '\n Good ! data format is right!!!!\n'
    else:
        print 'Warning !The data format may be wrong !'
    del save_row_sl_deramped
fid.close()