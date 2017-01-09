import numpy as np
import os, sys
from datetime import datetime

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

from sentinel_1.functions.get_ramp import freadbk
from sentinel_1.functions.ESD_functions import get_f_DC_difference, get_coordinates


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
    to_angle_matrix = np.zeros(shape=(1, len(dates), len(dates)))
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
                # Phase ramp per pixel
                to_angle_matrix[0, n, num] = (PRF/(2*np.pi*np.nanmean(df_dc_ps))) * (line_start - 1)

    return diff_matrix, std_matrix, to_angle_matrix, weight, dates

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
    diff_matrix, std_matrix, to_angle_matrix, weight, dates = network_esd_ps(stack_folder, master_date, burst)

    # And save them in the corresponding folder:
    folder = os.path.join(stack_folder, 'esd', burst)

    np.save(os.path.join(folder, 'diff_matrix'), diff_matrix)
    np.save(os.path.join(folder, 'std_matrix'), std_matrix)
    np.save(os.path.join(folder, 'to_angle_matrix'), to_angle_matrix)
    np.save(os.path.join(folder, 'weight'), weight)
    np.save(os.path.join(folder, 'dates'), dates)


