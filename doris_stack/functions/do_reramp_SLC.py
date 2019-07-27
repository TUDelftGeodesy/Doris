#!/usr/bin/env python
import numpy as np
from numpy import *
from doris.doris_stack.functions.get_ramp import get_ramp
from doris.doris_stack.functions.ESD_functions import freadbk
from doris.doris_stack.main_code.resdata import ResData
import sys

def usage():
    print('\nUsage: python  do_reramp_SLC.py dataFilename resFilename resampled')
    print('  where dataFilename        is the name of burst you want to deramp')
    print('        resFilename         is the .res file of burst              ')
    print(' This python applies the inverse phase ramp to the burst pointed by DATAFILENAME (slc)')
    print(' and RESFILENAME (res) that was deramped by deramp_SLC.m. The phase screen')
    print(' must account for the new resampled grids PIXRGGRID and PIXAZGRID    ')
    print(' [Nlines_mst x Nsamples_mst] that contain the time coordinates of the')
    print(' resampled image into the master grid:                               ')
    print('  for example                                                        ')
    print(' python   do_reramp_SLC.py slave_rsmp.raw slave.res False            ')
    print(' created by Gert Mulder')
    print(' Part of code adapted from Lorenzo Iannini and Wu Wenhao')
try:
    dataFilename = sys.argv[1]
    resFilename = sys.argv[2]
     
except:
    print('Unrecognized input')
    usage()
    sys.exit(1)
if len(sys.argv) == 3:
    resampled = True
elif len(sys.argv) == 4:
    resampled = sys.argv[3]
else:
    print('Unrecognized input')
    usage()
    sys.exit(1)

# Read information
################################################################################

res = ResData(resFilename)

# Image size properties
if res.process_control['resample'] == '1': #oversampled data
    l0 = int(res.processes['resample']['First_line (w.r.t. original_master)'])
    lN = int(res.processes['resample']['Last_line (w.r.t. original_master)'])
    p0 = int(res.processes['resample']['First_pixel (w.r.t. original_master)'])
    pN = int(res.processes['resample']['Last_pixel (w.r.t. original_master)'])
    dataFormat = 'complex64'
    resampled = True
else: # original data
    l0 = int(res.processes['crop']['First_line (w.r.t. original_image)'])
    lN = int(res.processes['crop']['Last_line (w.r.t. original_image)'])
    p0 = int(res.processes['crop']['First_pixel (w.r.t. original_image)'])
    pN = int(res.processes['crop']['Last_pixel (w.r.t. original_image)'])
    dataFormat = 'cpxint16'
    resampled = False

# Get resampled Slv size
Naz_res = lN-l0+1
Nrg_res = pN-p0+1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Obtain chirp %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not resampled:
    ChirpFilt = get_ramp(resFilename, resampled=0, type='chirp')
elif resampled:
    ChirpFilt = get_ramp(resFilename, resampled=1, type='chirp')
else:
    print('Resampled should be True or False')

################################################################################
# Read data

slc = freadbk(dataFilename, 1, 1, int(Naz_res), int(Nrg_res), dataFormat, int(Naz_res), int(Nrg_res))

#######################################################################################

newFilename = dataFilename[:-4] + '_reramped.raw'
fid = open(newFilename, 'wb')
slc_reramped = slc * ChirpFilt
del ChirpFilt
del slc

#%% Apply reramping
if dataFormat == 'complex64':
    slc_dat = slc_reramped.astype(np.complex64)
else:  # cpxint16
    slc_dat = np.zeros(shape=(int(Naz_res), int(Nrg_res) * 2)).astype('int16')
    slc_dat[:, 0::2] = np.real(slc_reramped)
    slc_dat[:, 1::2] = np.imag(slc_reramped)

fid.write(slc_dat)
fid.close()

print("\nReramp operation completed\n")
