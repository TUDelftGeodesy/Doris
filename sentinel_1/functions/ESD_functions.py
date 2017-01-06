import numpy as np
from numpy import *
from sentinel_1.functions.get_ramp import get_ramp, get_parameter
import os, sys


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
def freadbk(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32', lines=0, pixels=0):
    # First use memmap to get a memory map of the full file, than extract the requested part.

    if dt == 'cpxint16':
        dtype = np.dtype([('re', np.int16), ('im', np.int16)])
        file_dat = np.memmap(path_file, dtype=dtype, mode='r', shape=(lines, pixels)).view(np.int16).astype(np.float32).view(np.complex64)
        data = file_dat[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1].astype(
            'complex64', subok=False)
    elif dt == 'cpxshort':

        file_dat = np.memmap(path_file, dtype=np.dtype(np.float16), mode='r', shape=(lines, pixels * 2))
        data = 1j * file_dat[:, 1::2].astype('float32', subok=False)
        data += file_dat[:, 0::2].astype('float32', subok=False)
        data = data[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1]

    else:
        dt = np.dtype(dt)
        file_dat = np.memmap(path_file, dtype=dt, mode='r', shape=(lines, pixels))
        data = file_dat[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1].astype(
            dt, subok=False)

    return data


########################################################################################################################
# Function to calculate normalized Doppler Centroid frequency
# Df_DC = get_f_DC_difference(f_DC_1, f_DC_2, BOL_Length, BOL_lines, PRF, normalize)
def get_f_DC_difference(nBurst):

    burst1 = 'burst_' + str(nBurst) + '/'
    burst2 = 'burst_' + str(nBurst + 1) + '/'

    this_m_resData = burst1 + 'master.res'
    next_m_resData = burst2 + 'master.res'

    os.chdir(os.getcwd() + '/' + burst1)
    f_DC_1 = get_ramp(os.path.basename(this_m_resData), resampled=0, type='DC')
    os.chdir(os.path.dirname(os.getcwd()))

    os.chdir(os.getcwd() + '/' + burst1)
    f_DC_2 = get_ramp(os.path.basename(next_m_resData), resampled=0, type='DC')
    os.chdir(os.path.dirname(os.getcwd()))

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    Df_DC = f_DC_1[line_start - 1:line_start + line_length - 1, first_pixel_this - 1:first_pixel_this + pixel_length - 1] - \
            f_DC_2[0:line_length, first_pixel_next - 1: first_pixel_next + pixel_length - 1]

    return Df_DC

def get_coordinates(nBurst):

    burst1 = 'burst_' + str(nBurst) + '/'
    burst2 = 'burst_' + str(nBurst+1) + '/'
    this_m_resData = burst1 + 'master.res'
    next_m_resData = burst2 + 'master.res'

    # Get variables from first burst
    this_line_first     = int(get_parameter('First_line (w.r.t. output_image)',this_m_resData,1))
    this_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',this_m_resData,1))
    this_nr_oflines     = int(this_line_last) - int(this_line_first) +1
    this_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',this_m_resData,1))
    this_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',this_m_resData,1))
    this_nr_ofpixels    = int(this_pixel_last) - int(this_pixel_first) +1
    PRF_1               = float(get_parameter('Pulse_Repetition_Frequency (computed, Hz)',this_m_resData,1))

    # Get variables from second burst
    next_line_first     = int(get_parameter('First_line (w.r.t. output_image)',next_m_resData,1))
    next_line_last      = int(get_parameter('Last_line (w.r.t. output_image)',next_m_resData,1))
    next_nr_oflines     = int(next_line_last) - int(next_line_first) +1
    next_pixel_first    = int(get_parameter('First_pixel (w.r.t. output_image)',next_m_resData,1))
    next_pixel_last     = int(get_parameter('Last_pixel (w.r.t. output_image)',next_m_resData,1))
    next_nr_ofpixels    = int(next_pixel_last) - int(next_pixel_first) +1

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

    first_pixel_this = first_pixel - this_pixel_first + 1
    first_pixel_next = first_pixel - next_pixel_first + 1

    line_length = this_line_last - next_line_first + 1
    line_start = this_nr_oflines - line_length + 1

    return line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, this_nr_ofpixels,\
           next_nr_oflines, next_nr_ofpixels, PRF


########################################################################################################################
# Function to calculate the interferograms of both bursts
# thisburstdata, nextburstdata, diffBursts, PRF = get_interfero(nBurst, BOL_lines, BOL_Length)
def get_offset(nBurst, Df_DC, coh_treshold=0.3):

    burst1 = 'burst_' + str(nBurst) + '/'
    burst2 = 'burst_' + str(nBurst+1) + '/'

    # cpxint16 and cpxfloat32
    dataFormat_s = 'complex64'

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    ifgs_1  = freadbk(burst1 + 'cint.raw.old', line_start, first_pixel_this, line_length, pixel_length , dataFormat_s,  this_nr_oflines, this_nr_ofpixels)
    ESD_coh_1  = freadbk(burst1 + 'coherence.raw', line_start, first_pixel_this, line_length, pixel_length , 'float32',  this_nr_oflines, this_nr_ofpixels)
    ifgs_2  = freadbk(burst2 + 'cint.raw.old', 1, first_pixel_next, line_length, pixel_length, dataFormat_s, next_nr_oflines, next_nr_ofpixels)
    ESD_coh_2 = freadbk(burst2 + 'coherence.raw', 1, first_pixel_next, line_length,
                     pixel_length, 'float32', next_nr_oflines, next_nr_ofpixels)
    ESD_coh = (ESD_coh_1 + ESD_coh_2) / 2

    #ifgs_1_total = freadbk(burst1 + 'cint.raw.old', 1, 1, this_nr_oflines, this_nr_ofpixels, dataFormat_s,  this_nr_oflines, this_nr_ofpixels)
    #ifgs_2_total = freadbk(burst2 + 'cint.raw.old', 1, 1, next_nr_oflines, next_nr_ofpixels, dataFormat_s,  next_nr_oflines, next_nr_ofpixels)

    # Remove invalid data both in range and azimuth
    valid_range = []
    valid_azimuth = []
    for i in range(0,len(ifgs_1[0,:])):
        if np.nanmean(abs(ifgs_1[:,i])) != 0 and np.nanmean(abs(ifgs_2[:,i])) != 0:
            valid_range.append(i)

    for i in range(0,len(ifgs_1[:,0])):
        if np.nanmean(abs(ifgs_1[i,:])) != 0 and np.nanmean(abs(ifgs_2[i,:])) != 0:
            valid_azimuth.append(i)

    if valid_range and valid_azimuth:
        ifgs_1 = ifgs_1[:, valid_range[:]]
        ifgs_2 = ifgs_2[:, valid_range[:]]
        ESD_coh = ESD_coh[:, valid_range[:]]

        ifgs_1 = ifgs_1[valid_azimuth[:], :]
        ifgs_2 = ifgs_2[valid_azimuth[:], :]
        ESD_coh = ESD_coh[valid_azimuth[:], :]

        Df_DC = Df_DC[:, valid_range[:]]
        Df_DC = Df_DC[valid_azimuth[:], :]

    # First downsample 2 * 10
    Nra = 10
    Naz = 2
    new_ra = ESD_coh.shape[1] / 10
    new_az = ESD_coh.shape[0] / 2

    ESD_coh = ESD_coh[0:new_az*Naz-1:Naz, 0:new_ra*Nra-1:Nra]
    ifgs_1_multilook = ifgs_1[:new_az*2, :new_ra*10].reshape([new_az, Naz, new_ra, Nra]).mean(3).mean(1)
    ifgs_2_multilook = ifgs_2[:new_az*2, :new_ra*10].reshape([new_az, Naz, new_ra, Nra]).mean(3).mean(1)
    Df_DC_multilook = Df_DC[:new_az * 2, :new_ra * 10].reshape([new_az, Naz, new_ra, Nra]).mean(3).mean(1)

    # Double difference and calculate weights according to Cramer Rao bound
    diffBursts = ifgs_1_multilook * ifgs_2_multilook.conj()
    weights = 2 * ESD_coh*ESD_coh / (1 - ESD_coh*ESD_coh)

    W = np.sum(weights[ESD_coh > coh_treshold])
    offset = np.angle(np.sum(diffBursts[ESD_coh > coh_treshold] * weights[ESD_coh > coh_treshold]) / W) * \
             (PRF / (2 * np.pi * np.nanmean(Df_DC_multilook[ESD_coh > coh_treshold] * weights[ESD_coh > coh_treshold] /
                                            np.mean(weights[ESD_coh > coh_treshold]))))

    return offset, W

########################################################################################################################
# Function to calculate pixel offset for each burst, according to Nida's proposed method. Threshold can be assigned by
# user.
# offset = apply_ESD_Nida(diffBursts, Df_DC, PRF, threshold)
def apply_ESD_Nida(diffBursts, Df_DC, PRF, threshold = 0.0001):

    # Determine phase of interferogram
    ph_esd = np.angle(diffBursts)

    # Do an intitial estimation based on the peak of the histogram
    N,X = np.histogram(ph_esd[:], 50)
    idx = np.argmax(N)
    D_az_init = X[idx]*(PRF/(2*np.pi*np.nanmean(Df_DC[:])))

    # Create linspace from -0.1 to 0.1 in 7 steps
    D_az_span = -0.3*np.pi*(PRF/(2*np.pi*np.nanmean(Df_DC[:])))
    D_azs = np.linspace(D_az_init-D_az_span, D_az_init+D_az_span,num=7)
    del D_az_span

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
            #print ph_test

        ind = np.argmin(abs(ph_test))
        D_az_min.append(D_azs[ind])

        # Break iteration when pixel shift is sufficient
        if c > 2 and abs(D_az_min[c]-D_az_min[c-1]) < threshold and abs(D_az_min[c-1]-D_az_min[c-2]) < threshold:
            ph_est_opt = np.nanmean(ph_est[:])
            offset = -D_az_min[c]
            break

        # Use a smaller difference for next iteration
        D_az_span = D_azs[1] - D_azs[0]
        D_azs = np.linspace(D_azs[ind]-D_az_span, D_azs[ind]+D_az_span,num=7)
        del ph_test

    #print 'amount of loops in iteration ' + str(c)

    pix_offset = offset / (PRF/(2*np.pi*np.nanmean(Df_DC[:])))

    return offset, pix_offset
