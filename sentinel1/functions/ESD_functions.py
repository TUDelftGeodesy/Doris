import os,sys,time
import math
import numpy as np
from numpy import *
import scipy as Sci
import scipy.linalg
import scipy.io as sio
from  scipy  import ndimage
import scipy.optimize as optimization
# import matplotlib.pyplot as plt
import struct,shutil
import gdal,gdalconst
from gdalconst import *
from scipy import stats
########################################################################################################################
# Function to get parameters from files
# Parameter = get_parameter(First_param,file_name,format_flag=1,Second_param=None,Third_param=None)
def get_parameter(First_param,file_name,format_flag=1,Second_param=None,Third_param=None):
    Read_contine_flag=0
    class set_class(object):
        pass
    orbit_info = set_class()
    time_temp = []
    x_temp = []
    y_temp = []
    z_temp = []
    value=None

    for line in open(file_name):
        if format_flag==1:
            if not (line.find(First_param)):
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                return value

        if format_flag==2:
            if  not (line.find(Second_param)):
                Read_contine_flag=1
            if (Read_contine_flag==1) and (not (line.find(First_param))):  ##Be careful
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
            if  (not (line.find(Third_param))):  ##Be careful
                Read_contine_flag=0
                return value

        if format_flag==3:
            if not (line.find(First_param)):
                index=line.find(':')
                pixel_time=(line[(index+1):].strip(' \n\t')).split(' ')[1].split(':')
                return pixel_time

        if format_flag==4:
            if not (line.find(First_param)):
                index=line.find(':')
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1
                continue
            if (Read_contine_flag>=1) :
                new_line = line.strip('\n').split()
                time_temp.append(float(new_line[0]))
                x_temp.append(float(new_line[1]))
                y_temp.append(float(new_line[2]))
                z_temp.append(float(new_line[3]))
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+1)):
                    setattr(orbit_info,'x',x_temp)
                    setattr(orbit_info,'y',y_temp)
                    setattr(orbit_info,'z',z_temp)
                    setattr(orbit_info,'time',time_temp)
                    return orbit_info

########################################################################################################################
# Definition to extract data
# In Matlab:
# thisBurstData = freadbk(['burst' num2str(nBurst)   '/cint.raw'],nofLines1,formatData1, line1, nofLines1,1,nofPixels1);
# In Python:
# thisBurstData = freadbk(path_file,line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt=np.dtype(np.float32), lines=0, pixels=0):
def freadbk(path_file,line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32', lines=0, pixels=0):

    # First use memmap to get a memory map of the full file, than extract the requested part.
    if dt == 'cpxint16':
        file_dat = np.memmap(path_file, dtype=np.dtype(np.int16), mode='r', shape=(lines, pixels * 2))
        data = 1j * file_dat[:,0::2].astype('float32',subok=False)
        data += file_dat[:,1::2].astype('float32',subok=False)
        data = data[line_start-1:line_start+nofLines-1,pixel_start-1:pixel_start+nofPixels-1]
    elif dt == 'cpxshort':
        file_dat = np.memmap(path_file, dtype=np.dtype(np.float16), mode='r', shape=(lines, pixels * 2))
        data = 1j * file_dat[:,1::2].astype('float32',subok=False)
        data += file_dat[:,0::2].astype('float32',subok=False)
        data = data[line_start-1:line_start+nofLines-1,pixel_start-1:pixel_start+nofPixels-1]
    else:
        dt = np.dtype(dt)
        file_dat = np.memmap(path_file, dtype=dt, mode='r', shape=(lines, pixels))
        data = file_dat[line_start-1:line_start+nofLines-1,pixel_start-1:pixel_start+nofPixels-1].astype(dt,subok=False)

    return data
########################################################################################################################
# Function to get the Burst Overlap of the two bursts
## Textread might be wrong
def get_BOL_lines(burstList):
#   Example: BOL_length, BOL_lines = get_BOL_lines([1 2])
#   burstList is an array with the two bursts of which you want to have the
#   burst overlap indexes of.
#   BOL_lines is a structure with BOL_lines.this and
#   BOL_lines.next, where BOL_lines.this contains an array with indexes of
#   the burstoverlap area of the first burst and BOL_lines.next contains an
#   array with indexes of the burstoverlap area of the second burst

    pr_dir = []

    if len(burstList) > 2:
        print 'burstList should be 2x1 array, indicating the two bursts to be processed.'
        return

    # Loop over bursts
    for nBurst in burstList:
        # Read metadata
        # Master res and ifg res files for current burst
        masterRes1 = 'burst-' + str(nBurst) + '/' + 'master.res'
        ifgRes1 = 'burst-' + str(nBurst) + '/' + 'ifgs.res'

        # Tries to read the Azimuth_time_interval otherwise points out a possible error
        try:
            PRF1 = get_parameter('Pulse_Repetition_Frequency (computed, Hz)',masterRes1,1)
            PRF1 = float(PRF1)
        except:
            print "Something is wrong reading Azimuth_time_interval, check get_parameter... "

        # Master res and ifg res files for next bursts
        masterRes2 = 'burst-' + str(nBurst+1) + '/' + 'master.res'
        ifgRes2 = 'burst-' + str(nBurst+1) + '/' + 'ifgs.res'

        PRF2 = get_parameter('Pulse_Repetition_Frequency (computed, Hz)',masterRes2,1)
        PRF2 = float(PRF2)

        # Lines and data format
        Number_lines=get_parameter('First_line (w.r.t. original_master)',ifgRes1,1)
        l1_1     = float(Number_lines)

        Number_lines=get_parameter('Last_line (w.r.t. original_master)',ifgRes1,1)
        l1_last     = float(Number_lines)
        burst_Length = l1_last-l1_1+1

        Number_lines=get_parameter('First_line (w.r.t. original_master)',ifgRes2,1)
        l1_2     = float(Number_lines)

        # Reads time and transform it to seconds
        #TODO for Taz_start near midnight this code will not work, because days are not used
        Taz_start = get_parameter('First_pixel_azimuth_time (UTC)',masterRes1,3)

        Taz_start1 = float(Taz_start[0])*3600 + float(Taz_start[1])*60 + float(Taz_start[2])
        Taz_start1 = Taz_start1 + ((l1_1-1)/PRF1)

        Taz_start = get_parameter('First_pixel_azimuth_time (UTC)',masterRes2,3)

        Taz_start2 = float(Taz_start[0])*3600 + float(Taz_start[1])*60 + float(Taz_start[2])
        Taz_start2 = Taz_start2 + ((l1_2-1)/PRF2)

        # Overlapping part goes for nextburstdata from 1 to line2
        BOL_Length = round((Taz_start1-Taz_start2)*PRF2)+burst_Length

        class BOL_lines():
            pass

        # Azimuth coordinate
        az_first     = int(get_parameter('First_line (w.r.t. output_image)',masterRes2,1))
        BOL_centre = az_first + int(round(BOL_Length / 2))

        BOL_lines.this = [(burst_Length-BOL_Length+1),burst_Length]
        BOL_lines.next = [1,BOL_Length]


        return BOL_Length, BOL_lines, BOL_centre

########################################################################################################################
# Function to calculate normalized Doppler Centroid frequency
# Df_DC = get_f_DC_difference(f_DC_1, f_DC_2, BOL_Length, BOL_lines, PRF, normalize)
def get_f_DC_difference(nBurst, BOL_lines, normalize = True):

    burst1 = 'burst-' + str(nBurst) + '/'
    burst2 = 'burst-' + str(nBurst+1)+ '/'
    this_m_resData = burst1 + 'master.res'
    this_s_resData = burst1 + 'slave.res'
    next_m_resData = burst2 + 'master.res'
    next_s_resData = burst2 + 'slave.res'

    # Get variables from first burst
    this_line_first     = int(get_parameter('First_line (w.r.t. output_image)',this_m_resData,1))
    this_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',this_m_resData,1))
    this_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',this_m_resData,1))
    this_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',this_m_resData,1))

    # Get variables from second burst
    next_line_first     = int(get_parameter('First_line (w.r.t. output_image)',next_m_resData,1))
    next_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',next_m_resData,1))
    next_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',next_m_resData,1))
    next_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',next_m_resData,1))

    # Read only the Burstoverlap
    if this_pixel_first < next_pixel_first:
        first_pixel = next_pixel_first
    elif this_pixel_first >= next_pixel_first:
        first_pixel = this_pixel_first

    if this_pixel_last > next_pixel_last:
        last_pixel = next_pixel_last
    elif this_pixel_last <= next_pixel_last:
        last_pixel = this_pixel_last



    master_fdc_parameters_file='burst-'+str(nBurst)+'/master_fdc_parameters.mat'
    f_DC_1     =scipy.io.loadmat(master_fdc_parameters_file)
    f_DC_1     =f_DC_1['f_DC']

    slave_fdc_parameters_file ='burst-'+str(nBurst+1)+'/rsmp_slave_deramp_parameters.mat'
    f_DC_2      =scipy.io.loadmat(slave_fdc_parameters_file)
    f_DC_2      =f_DC_2['f_DC']

    if normalize == True or normalize == 'True':
        f_DC_1 = f_DC_1 - np.nanmean(f_DC_1)
        f_DC_2 = f_DC_2 - np.nanmean(f_DC_2)

    this_1 = first_pixel - this_pixel_first
    next_1 = first_pixel - next_pixel_first
    this_2 = last_pixel - this_pixel_first
    next_2 = last_pixel - next_pixel_first
    Df_DC = f_DC_1[BOL_lines.this[0]-1:BOL_lines.this[1], this_1:this_2+1] - f_DC_2[BOL_lines.next[0]-1:BOL_lines.next[1], next_1:next_2+1]

    return Df_DC

########################################################################################################################
# Function to calculate the interferograms of both bursts
# thisburstdata, nextburstdata, diffBursts, PRF = get_interfero(nBurst, BOL_lines, BOL_Length)
def get_interfero(nBurst, BOL_lines, BOL_Length, Df_DC):

    burst1 = 'burst-' + str(nBurst) + '/'
    burst2 = 'burst-' + str(nBurst+1)+ '/'
    this_m_resData = burst1 + 'master.res'
    this_s_resData = burst1 + 'slave.res'
    next_m_resData = burst2 + 'master.res'
    next_s_resData = burst2 + 'slave.res'

    # cpxint16 and cpxfloat32
    dataFormat_m = 'cpxint16'
    dataFormat_s = 'complex64'

    # Get variables from first burst
    this_line_first     = int(get_parameter('First_line (w.r.t. output_image)',this_m_resData,1))
    this_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',this_m_resData,1))
    this_nr_oflines     = int(this_line_last) - int(this_line_first) +1
    this_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',this_m_resData,1))
    this_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',this_m_resData,1))
    this_nr_ofpixels    = int(this_pixel_last) - int(this_pixel_first) +1
    mastername_1        = get_parameter('Data_output_file',this_m_resData,1)
    slavename_1         = get_parameter('Data_output_file',this_s_resData,2,'*_Start_resample','* End_resample:_NORMAL')
    PRF_1               = float(get_parameter('Pulse_Repetition_Frequency (computed, Hz)',this_m_resData,1))
    this_pixel_length   = int(get_parameter('Number_of_pixels_original', this_m_resData,1))

    # Get variables from second burst
    next_line_first     = int(get_parameter('First_line (w.r.t. output_image)',next_m_resData,1))
    next_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',next_m_resData,1))
    next_nr_oflines     = int(next_line_last) - int(next_line_first) +1
    next_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',next_m_resData,1))
    next_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',next_m_resData,1))
    next_nr_ofpixels    = int(next_pixel_last) - int(next_pixel_first) +1
    mastername_2        = get_parameter('Data_output_file',next_m_resData,1)
    slavename_2         = get_parameter('Data_output_file',next_s_resData,2,'*_Start_resample','* End_resample:_NORMAL')
    PRF_2               = float(get_parameter('Pulse_Repetition_Frequency (computed, Hz)',next_m_resData,1))
    next_pixel_length   = int(get_parameter('Number_of_pixels_original', next_m_resData,1))

    PRF = PRF_1

    # Read only the Burstoverlap
    if this_pixel_first < next_pixel_first:
        first_pixel = next_pixel_first
    elif this_pixel_first >= next_pixel_first:
        first_pixel = this_pixel_first

    if this_pixel_last > next_pixel_last:
        pixel_length = next_pixel_last - first_pixel + 1
    elif this_pixel_last <= next_pixel_last:
        pixel_length = this_pixel_last - first_pixel + 1

    # Mas_1 = freadbk(burst1 + mastername_1, min(BOL_lines.this), first_pixel - this_pixel_first + 1, BOL_Length, pixel_length , dataFormat_m,  this_nr_oflines, this_nr_ofpixels)
    # Mas_2 = freadbk(burst2 + mastername_2, min(BOL_lines.next), first_pixel - next_pixel_first + 1, BOL_Length, pixel_length , dataFormat_m,  next_nr_oflines, next_nr_ofpixels)
    # Sla_1  = freadbk(burst1 + slavename_1, min(BOL_lines.this), first_pixel - this_pixel_first + 1, BOL_Length, pixel_length , dataFormat_s,  this_nr_oflines, this_nr_ofpixels)
    # Sla_2  = freadbk(burst2 + slavename_2, min(BOL_lines.next), first_pixel - next_pixel_first + 1, BOL_Length, pixel_length , dataFormat_s,  next_nr_oflines, next_nr_ofpixels)

    ifgs_1  = freadbk(burst1 + 'cint.raw', min(BOL_lines.this), first_pixel - this_pixel_first + 1, BOL_Length, pixel_length , dataFormat_s,  this_nr_oflines, this_nr_ofpixels)
    ifgs_2  = freadbk(burst2 + 'cint.raw', min(BOL_lines.next), first_pixel - next_pixel_first + 1, BOL_Length, pixel_length , dataFormat_s,  next_nr_oflines, next_nr_ofpixels)

    # Remove invalid data both in range and azimuth
    valid_range = []
    valid_azimuth = []
    for i in range(0,len(ifgs_1[0,:])):
        if np.nanmean(abs(ifgs_1[:,i])) != 0 and np.nanmean(abs(ifgs_2[:,i])) != 0:
            valid_range.append(i)

    for i in range(0,len(ifgs_1[:,0])):
        if np.nanmean(abs(ifgs_1[i,:])) != 0 and np.nanmean(abs(ifgs_2[i,:])) != 0:
            valid_azimuth.append(i)

    # Mas_1 = Mas_1[:, valid_range[:]]
    # Mas_2 = Mas_2[:, valid_range[:]]
    # Sla_1 = Sla_1[:, valid_range[:]]
    # Sla_2 = Sla_2[:, valid_range[:]]
    ifgs_1 = ifgs_1[:, valid_range[:]]
    ifgs_2 = ifgs_2[:, valid_range[:]]

    Df_DC = Df_DC[:, valid_range[:]]

    # Mas_1 = Mas_1[valid_azimuth[:], :]
    # Mas_2 = Mas_2[valid_azimuth[:], :]
    # Sla_1 = Sla_1[valid_azimuth[:], :]
    # Sla_2 = Sla_2[valid_azimuth[:], :]
    ifgs_1 = ifgs_1[valid_azimuth[:], :]
    ifgs_2 = ifgs_2[valid_azimuth[:], :]

    Df_DC = Df_DC[valid_azimuth[:], :]

    # thisBurstData = Mas_1*Sla_1.conj()
    # nextBurstData = Mas_2*Sla_2.conj()

    thisBurstData = ifgs_1
    nextBurstData = ifgs_2

    diffBursts = nextBurstData*thisBurstData.conj()

    real_diffBursts = ndimage.gaussian_filter(np.real(diffBursts),sigma=5)
    imag_diffBursts = ndimage.gaussian_filter(np.imag(diffBursts),sigma=5)
    diffBursts      = real_diffBursts +1j*imag_diffBursts

    abs_nextBurstData=ndimage.gaussian_filter(np.abs(nextBurstData)*np.abs(nextBurstData),sigma=5)+0.001
    abs_thisBurstData=ndimage.gaussian_filter(np.abs(thisBurstData)*np.abs(thisBurstData),sigma=5)+0.001
    ESD_coh          = np.abs(diffBursts)/sqrt(abs_nextBurstData*abs_thisBurstData)  #

    W = float(np.sum(ESD_coh > 0.5)) / float(np.sum(ESD_coh > 0))
    diffBursts = diffBursts[ESD_coh > 0.5]
    thisBurstData = thisBurstData[ESD_coh > 0.5]
    nextBurstData = nextBurstData[ESD_coh > 0.5]
    Df_DC = Df_DC[ESD_coh > 0.5]

    return thisBurstData, nextBurstData, diffBursts, PRF, Df_DC, W

########################################################################################################################
# Function to calculate pixel offset for each burst, according to Nida's proposed method. Threshold can be assigned by
# user.
# offset = apply_ESD_Nida(diffBursts, Df_DC, PRF, threshold)
def apply_ESD_Nida(diffBursts, Df_DC, PRF, threshold = 0.0001):

    # Determine phase of interferogram
    ph_esd = np.angle(diffBursts)

    # Do an intitial estimation based on the peak of the histogram
    N,X = np.histogram(ph_esd[:],50)
    idx = np.argmax(N)
    D_az_init = X[idx]*(PRF/(2*np.pi*np.nanmean(Df_DC[:])))

    # Create linspace from -0.1 to 0.1 in 7 steps
    D_az_span = -0.3*np.pi*(PRF/(2*np.pi*np.nanmean(Df_DC[:])))
    D_azs = np.linspace(D_az_init-D_az_span, D_az_init+D_az_span,num=7)
    del D_az_span

    print str(D_az_init)
    print str(D_azs[0])

    c = -1
    D_az_min = []
    # Keeps refinining until error with earlier estimation is lower than 0.0001
    while True:
        c += 1
        ph_test = np.ones(len(D_azs))
        # Calculate estimated phase, residual phase and test phase
        for k in range(0,len(D_azs)):
            D_az = D_azs[k]
            ph_est = (2*np.pi*Df_DC*D_az)/PRF
            ph_res = ph_esd - ph_est

            ph_test[k] = np.nanmean(np.angle(exp(1j * ph_res[:]))) # should be ph_test(k) = np.nanmean(exp(1i*ph_res[:]))
            print ph_test

        ind = np.argmin(abs(ph_test))
        D_az_min.append(D_azs[ind])

        print 'next iteration'
        # Break iteration when pixel shift is sufficient
        if c > 2 and abs(D_az_min[c]-D_az_min[c-1]) < threshold and abs(D_az_min[c-1]-D_az_min[c-2]) < threshold:
            ph_est_opt = np.nanmean(ph_est[:])
            offset = -D_az_min[c]
            break

        # Use a smaller difference for next iteration
        D_az_span = D_azs[1] - D_azs[0]
        D_azs = np.linspace(D_azs[ind]-D_az_span, D_azs[ind]+D_az_span,num=7)
        del ph_test

    print 'amount of loops in iteration ' + str(c)

    pix_offset = offset / (PRF/(2*np.pi*np.nanmean(Df_DC[:])))

    return offset, pix_offset
########################################################################################################################
# Function to calculate pixel offset for each burst, according to Wu's proposed method. Threshold can be assigned by
# user.
# pixel_shift, Total_shift, Iteration_ESD, azimuth_shift, Weight_shift, Valid_ESD_index,
# = apply_ESD_Wu(diffBursts, Df_DC, azimFreq1, Total_shift, Iteration_ESD, threshold, max_iterate)
def apply_ESD_Wu(nBurst, diffBursts, Df_DC, Total_shift, Iteration_ESD, azimuth_shift, Weight_shift, Valid_ESD_index, threshold = 0.0001, max_iterate = 25):

    masterRes = 'burst-' + str(nBurst) + '/' + 'master.res'
    Azimuth_time_interval = get_parameter('Azimuth_time_interval (s)',masterRes,1)
    azimFreq    = 1.0/float(Azimuth_time_interval)

     # Do initial estimation based on peak of histogram
    diffBursts = np.angle(diffBursts)/Df_DC/2/np.pi*azimFreq
    print "Data size is: ",diffBursts.shape
    hist,bins = np.histogram(diffBursts,500)
    print 'Currect variable value:Histogram_shift=',bins[np.argmax(hist)]

    width =0.7*(bins[1]-bins[0])
    bins=0.5*bins[1:]+0.5*bins[:-1]
    ##plt.subplot(121)
    ##plt.plot(bins,hist)
    ##plt.bar(bins, hist, align='center', width=width)
    ##plt.show()
    ##plt.close()

    # Iterate until fit is close enough to peak
    nloc,nscale=stats.norm.fit(diffBursts)
    print 'Currect variable value: nloc,nscale=',nloc,nscale

    iteration_temp=0

    print "\nStart iteration operation:\n",
    while np.abs(nloc - bins[argmax(hist)])> threshold:
        iteration_temp= iteration_temp+1
        #print "\niterative time:",iteration_temp
        if nloc > bins[argmax(hist)]:
            diffBursts_optimization_fit = diffBursts[where(diffBursts < (nloc+ nscale))]
            diffBursts_optimization_fit = diffBursts_optimization_fit[where(diffBursts_optimization_fit > (bins[argmax(hist)] - nscale))]
        else:
            diffBursts_optimization_fit = diffBursts[where(diffBursts < (bins[argmax(hist)] + nscale))]
            diffBursts_optimization_fit = diffBursts_optimization_fit[where(diffBursts_optimization_fit > (nloc - nscale))]

        nloc, nscale = stats.norm.fit(diffBursts_optimization_fit)
        # print 'Currect variable value:nloc,nscale=',nloc,nscale

        stat_val,p_val=stats.kstest(diffBursts_optimization_fit,'norm',(nloc,nscale))
        # print 'KS-statistic D = %6.3f p-value = %6.4f' % (stat_val, p_val)

        hist,bins = np.histogram(diffBursts_optimization_fit,500)
        # print 'Currect variable value: Histogram_shift=',bins[argmax(hist)]

        info = stats.describe(diffBursts_optimization_fit)
        # print "Data size is : "+ str(info[0])
        # print "Minimum value is: " + str(info[1][0])
        # print "Maximum value is: " + str(info[1][1])
        # print "Arithmetic mean is: " + str(info[2])
        # print "Unbiased variance is: " + str(info[3])
        # print "Biased skewness is: " + str(info[4])
        # print "Biased kurtosis is: " + str(info[5])+'\n'

        # Break if iteration takes too long
        if iteration_temp > max_iterate:
            break

    # Determine pixel shift
    if np.abs(nloc) < np.abs(bins[argmax(hist)]):
        pixel_shift=-nloc
    else:
        pixel_shift=-bins[argmax(hist)]

    # Calculate total pixel shift of all bursts and save variables
    if math.fabs(pixel_shift) < 0.025:
        Total_shift += pixel_shift
        print 'pixel_shift=',pixel_shift
        print '\n'
        azimuth_shift.append(pixel_shift)
        Weight_shift.append(nscale)
        Valid_ESD_index.append(nBurst+0.5)

    # Give warnings if values seem wrong
    if np.abs(nloc- bins[argmax(hist)]) > 0.0035 and np.abs(nloc) < 0.025 and np.abs(bins[argmax(hist)]) < 0.025:
        print "The difference between normfit and histogram",np.abs(nloc- bins[argmax(hist)])
        print "Warning !!!, please Check this result ! Maybe here is low coherence !!!"
        print "Warning !!!, please Check this result ! Maybe here is low coherence !!!"
        print "Warning !!!, please Check this result ! Maybe here is low coherence !!!"
        print "Warning !!!, please Check this result ! Maybe here is low coherence !!!"
        print "Warning !!!, please Check this result ! Maybe here is low coherence !!!"
    if  np.abs(nloc) > 0.025 and np.abs(bins[argmax(hist)]) > 0.025:
        print "Warning !!!, please Check this result ! Maybe caused by Geophysical signal !!!"
        print "Warning !!!, please Check this result ! Maybe caused by Geophysical signal !!!"
        print "Warning !!!, please Check this result ! Maybe caused by Geophysical signal !!!"
        print "Warning !!!, please Check this result ! Maybe caused by Geophysical signal !!!"
        print "Warning !!!, please Check this result ! Maybe caused by Geophysical signal !!!"

    return pixel_shift, Total_shift, Iteration_ESD, azimuth_shift, Weight_shift, Valid_ESD_index

