import numpy as np
import os, sys
from get_ramp import get_ramp
from datetime import datetime

def save_overlapping(stack_path, master_date, dates, overlap):

    nBurst, burst, next_burst = get_burst(overlap)

    esd_folder = os.path.join(stack_path, 'esd')
    if not os.path.exists(esd_folder):
        os.mkdir(esd_folder)
    overlap_path = os.path.join(stack_path, 'esd', overlap)
    if not os.path.exists(overlap_path):
        os.mkdir(overlap_path)

    path = swath_path(stack_folder, dates[0], master_date, burst)
    os.chdir(path)

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    burst1 = 'burst_' + str(nBurst) + '/'
    burst2 = 'burst_' + str(nBurst + 1) + '/'

    # First get the data for the master and df_dc
    master_path = os.path.join(overlap_path, master_date)
    df_dc_path = os.path.join(overlap_path, 'df_dc')

    if not os.path.exists(master_path + '_1') or not os.path.exists(master_path + '_2'):
        master_path_1 = dat_file(stack_path, burst, date='master')
        master_path_2 = dat_file(stack_path, next_burst, date='master')
        master_1 = freadbk(burst1 + master_path_1, line_start, first_pixel_this, line_length,
                          pixel_length, 'cpxint16', this_nr_oflines, this_nr_ofpixels)
        master_2 = freadbk(burst2 + master_path_2, 1, first_pixel_next, line_length, pixel_length,
                          'cpxint16', next_nr_oflines, next_nr_ofpixels)
        master_1_file = np.memmap(master_path + '_1', 'complex64', shape=master_1.shape, mode='w+')
        master_2_file = np.memmap(master_path + '_2', 'complex64', shape=master_2.shape, mode='w+')
        master_1_file[:, :] = master_1
        master_2_file[:, :] = master_2

    if not os.path.exists(df_dc_path):
        df_dc = get_f_DC_difference(nBurst)
        df_dc_file = np.memmap(df_dc_path, 'float32', shape=df_dc.shape, mode='w+')
        df_dc_file[:,:] = df_dc[:,:]

    # Then loop over the slaves
    for date in dates:
        path = swath_path(stack_folder, date, master_date, burst)
        os.chdir(path)

        data_path = os.path.join(overlap_path, date)

        burst1 = 'burst_' + str(nBurst) + '/'
        burst2 = 'burst_' + str(nBurst + 1) + '/'

        if not os.path.exists(data_path + '_1') or not os.path.exists(data_path + '_2'):
            slave_1 = freadbk(burst1 + 'slave_rsmp.raw.old', line_start, first_pixel_this, line_length, pixel_length , 'complex64',  this_nr_oflines, this_nr_ofpixels)
            slave_2 = freadbk(burst2 + 'slave_rsmp.raw.old', 1, first_pixel_next, line_length, pixel_length, 'complex64', next_nr_oflines, next_nr_ofpixels)
            slave_1_file = np.memmap(data_path + '_1', 'complex64', shape=slave_1.shape, mode='w+')
            slave_2_file = np.memmap(data_path + '_2', 'complex64', shape=slave_2.shape, mode='w+')
            slave_1_file[:,:] = slave_1
            slave_2_file[:,:] = slave_2


def find_ps_overlapping(stack_path, master_date, overlap):
    # This is used to find the ps point in overlapping areas

    nBurst, burst, next_burst = get_burst(overlap)

    esd_folder = os.path.join(stack_path, 'esd')
    if not os.path.exists(esd_folder):
        os.mkdir(esd_folder)
    files = os.listdir(os.path.join(esd_folder, overlap))
    dates = sorted([f[:-2] for f in files if f.endswith('_1')])

    path = swath_path(stack_folder, dates[0], master_date, burst)
    os.chdir(path)

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    # First calculate the ps point for first overlap
    files = os.listdir(os.path.join(esd_folder, overlap))
    first_files = [os.path.join(esd_folder, overlap, f) for f in files if f.endswith('_1')]
    first_name = os.path.join(esd_folder, overlap, 'first')
    first = np.memmap(first_name, 'float32', shape=(line_length, pixel_length, len(first_files)), mode='w+')

    for f, n in zip(first_files, range(len(first_files))):
        first_dat = np.memmap(f, 'complex64', mode='r+', shape=(line_length, pixel_length))
        first[:,:,n] = np.abs(first_dat[:,:])

    mean = np.mean(first, axis=2)
    std = np.std(first, axis=2)
    ps1 = (std/mean) < 0.3
    ps1_file = os.path.join(esd_folder, overlap, 'ps_1')
    ps1_dat = np.memmap(ps1_file, 'bool', mode= 'w+', shape=first_dat.shape)
    ps1_dat[:,:] = ps1[:,:]

    # Then calculate the ps point for second overlap
    second_files = [os.path.join(esd_folder, overlap, f) for f in files if f.endswith('_2')]
    second_name = os.path.join(esd_folder, overlap, 'second')
    second = np.memmap(second_name, 'float32', shape=(line_length, pixel_length, len(second_files)), mode='w+')

    for f, n in zip(second_files, range(len(second_files))):
        second_dat = np.memmap(f, 'complex64', mode='r+', shape=(line_length, pixel_length))
        second[:, :, n] = np.abs(second_dat[:, :])

    mean = np.mean(second, axis=2)
    std = np.std(second, axis=2)
    ps2 = (std/mean) < 0.3
    ps2_file = os.path.join(esd_folder, overlap, 'ps_2')
    ps2_dat = np.memmap(ps2_file, 'bool', mode= 'w+', shape=second_dat.shape)
    ps2_dat[:,:] = ps2[:,:]

    ps_file = os.path.join(esd_folder, overlap, 'ps')
    ps_dat = np.memmap(ps_file, 'bool', mode='w+', shape=second_dat.shape)
    ps_dat[:, :] = ((ps1_dat * ps2_dat) == 1)


def network_esd_ps(stack_folder, master_date, overlap, max_t_baseline=500):
    # Based on ps point esd is calculated using a network approach

    esd_folder = os.path.join(stack_folder, 'esd')
    if not os.path.exists(esd_folder):
        os.mkdir(esd_folder)

    files = os.listdir(os.path.join(esd_folder, overlap))
    dates = sorted([f[:-2] for f in files if (f.endswith('_1') and len(f) > 10)])

    diff_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    std_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    weight = 0

    path = swath_path(stack_folder, dates[0], master_date, burst)
    os.chdir(path)

    nBurst = int(overlap.split('_')[4])
    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    ps_file = os.path.join(esd_folder, overlap, 'ps')
    ps_dat = np.memmap(ps_file, 'bool', mode='r', shape=(line_length, pixel_length))

    ps_id = np.where(ps_dat == 1)
    if not ps_id:  # If there are no ps points
        return diff_matrix, std_matrix, weight
    else:
        weight = np.array([len(ps_id[0])])

    df_dc_file = np.memmap(os.path.join(esd_folder, overlap, 'df_dc'), 'float32', mode='r+', shape=(line_length, pixel_length))
    df_dc_ps = df_dc_file[ps_id]

    for date, n in zip(dates, range(len(dates))):
        for date_2, num in zip(dates[n+1:], range(n+1, len(dates))):
            time_diff = np.abs(
                (datetime.strptime(date, '%Y-%m-%d') - datetime.strptime(date_2, '%Y-%m-%d')).days)

            if time_diff < max_t_baseline:
                first_master = np.memmap(os.path.join(esd_folder, overlap, date + '_1'), 'complex64', mode='r+',
                                         shape=(line_length, pixel_length))
                first_slave = np.memmap(os.path.join(esd_folder,overlap, date_2 + '_1'), 'complex64', mode='r+',
                                        shape=(line_length, pixel_length))
                second_master = np.memmap(os.path.join(esd_folder,overlap, date + '_2'), 'complex64', mode='r+',
                                          shape=(line_length, pixel_length))
                second_slave = np.memmap(os.path.join(esd_folder,overlap, date_2 + '_2'), 'complex64', mode='r+',
                                         shape=(line_length, pixel_length))

                double_diff = (first_slave[ps_id] * first_master[ps_id].conj()) * \
                              (second_slave[ps_id] * second_master[ps_id].conj()).conj()
                pixel_std = np.std(np.angle(double_diff) * (PRF/(2*np.pi*df_dc_ps)))
                pixel_diff = np.angle(np.sum(double_diff)) * (PRF/(2*np.pi*np.nanmean(df_dc_ps)))

                std_matrix[0, n, num] = pixel_std
                diff_matrix[0, n, num] = pixel_diff

    return diff_matrix, std_matrix, weight, dates

def get_burst(overlap):
    s = overlap.split('_')
    burst = s[0] + '_' + s[1] + '_' + s[2] + '_' + s[3] + '_' + s[4]
    next_burst = s[5] + '_' + s[6] + '_' + s[7] + '_' + s[8] + '_' + s[9]

    nBurst = int(s[4])

    return nBurst, burst, next_burst

def swath_path(stack_folder, date, master_date, key):
    date_folder = master_date[:4] + master_date[5:7] + master_date[8:10] + '_' + date[:4] + date[5:7] + date[8:10]
    swath_burst = key.split('_')
    file_path = os.path.join(stack_folder, date_folder, swath_burst[0] + '_' + swath_burst[1])

    return file_path


def dat_file(stack_path, key, date, full_path=False, swath=False):
    # TODO refactor string lengt code
    # This function converts combinations of dates and keys to a datafile name

    string = '_iw_' + key[6] + '_burst_' + key[17:]

    if date == 'master' or date == 'slave':
        string = date + string + '.raw'
    else:  # In case it is a date
        string = date[:4] + date[5:7] + date[8:10] + string + '.raw'

    if full_path is True:
        date_folder = date[:4] + date[5:7] + date[8:10]
        swath_folder = key[:10]
        burst_folder = key[11:]
        if swath is False:
            string = os.path.join(stack_path, date_folder, swath_folder, burst_folder, string)
        elif swath is True:
            string = os.path.join(stack_path, date_folder, swath_folder, string)

    return string

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


# Actually execute the code...
if __name__ == "__main__":

    if len(sys.argv) == 4:
        stack_folder  = sys.argv[1]
        burst         = sys.argv[2]
        ps_select     = sys.argv[3]  # Use 1 if needed, use 0 if not
    else:
        sys.exit('usage: stack_folder burst')

    print('stack folder is ' + stack_folder)
    print('burst is ' + burst)
    print('ps select is ' + ps_select)

    # first get all the dates from the stack:
    ifgs = [f for f in os.listdir(stack_folder) if len(f) == 17]
    dates = [f[9:13] + '-' + f[13:15] + '-' + f[15:] for f in ifgs]
    master_date = ifgs[0][0:4] + '-' + ifgs[0][4:6] + '-' + ifgs[0][6:8]

    # Then run the overlap cutout / ps selection and
    save_overlapping(stack_folder, master_date, dates, burst)

    # If we want to select ps points
    if ps_select == '1':
        find_ps_overlapping(stack_folder, master_date, burst)

    # Get the esd results for the overlapping areas
    diff_matrix, std_matrix, weight, dates = network_esd_ps(stack_folder, master_date, burst)

    # And save them in the corresponding folder:
    folder = os.path.join(stack_folder, 'esd', burst)

    np.save(os.path.join(folder, 'diff_matrix'), diff_matrix)
    np.save(os.path.join(folder, 'std_matrix'), std_matrix)
    np.save(os.path.join(folder, 'weight'), weight)
    np.save(os.path.join(folder, 'dates'), dates)


