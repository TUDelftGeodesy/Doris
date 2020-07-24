#!/usr/bin/env python
import numpy as np
from numpy import *
from get_ramp import get_ramp
from doris.doris_stack.functions.ESD_functions import freadbk
from doris.doris_stack.main_code.resdata import ResData
import sys


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
    print ' python  do_deramp_SLC.py   20140821_iw_2_burst_1.raw subordinate.res False   '
    print ' created by Gert Mulder'
    print ' Part of code adapted from Lorenzo Iannini and Wu Wenhao'
try:
    dataFilename = sys.argv[1]
    resFilename = sys.argv[2]
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)

#*****************************************************************************#
# Calculate chirp for deramping
ChirpFilt = get_ramp(resFilename, resampled=0, type='chirp')

res = ResData(resFilename)

# Image size properties
if res.process_control['oversample'] == '1': #oversampled data
    l0 = int(res.processes['oversample']['First_line (w.r.t. ovs_image)'])
    lN = int(res.processes['oversample']['Last_line (w.r.t. ovs_image)'])
    p0 = int(res.processes['oversample']['First_pixel (w.r.t. ovs_image)'])
    pN = int(res.processes['oversample']['Last_pixel (w.r.t. ovs_image)'])
    dataFormat = 'cpxfloat32'
else: # original data
    l0 = int(res.processes['crop']['First_line (w.r.t. original_image)'])
    lN = int(res.processes['crop']['Last_line (w.r.t. original_image)'])
    p0 = int(res.processes['crop']['First_pixel (w.r.t. original_image)'])
    pN = int(res.processes['crop']['Last_pixel (w.r.t. original_image)'])
    dataFormat = 'cpxint16'

# Image size
Naz_res = lN-l0+1
Nrg_res = pN-p0+1

################################################################################
# Read data

slc = freadbk(dataFilename, 1, 1, int(Naz_res), int(Nrg_res), dataFormat, int(Naz_res), int(Nrg_res))

#######################################################################################

newFilename = dataFilename[:-4] + '_deramped.raw'
fid = open(newFilename, 'wb')
slc_deramped = conj(ChirpFilt)*slc
del ChirpFilt
del slc

# %% Apply reramping
if dataFormat == 'complex64':
    slc_dat = slc_deramped.astype(np.complex64)
else:  # cpxint16
    slc_dat = np.zeros(shape=(int(Naz_res), int(Nrg_res) * 2)).astype('int16')
    slc_dat[:, 0::2] = np.real(slc_deramped)
    slc_dat[:, 1::2] = np.imag(slc_deramped)

fid.write(slc_dat)
fid.close()

print "\nDeramp operation completed\n"
