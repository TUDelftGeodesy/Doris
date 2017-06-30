import numpy as np
import os, sys
from datetime import datetime
from datetime import timedelta
from doris.doris_stack.functions.ESD_functions import get_f_DC_difference, get_coordinates, freadbk


def save_overlapping(stack_folder, master_date, dates, overlap):

    nBurst, burst, next_burst = get_burst(overlap)

    esd_folder = os.path.join(stack_folder, 'esd')
    if not os.path.exists(esd_folder):
        print('ESD folder does not exist')
        return
    overlap_path = os.path.join(stack_folder, 'esd', overlap)
    if not os.path.exists(overlap_path):
        os.mkdir(overlap_path)

    path = swath_path(stack_folder, dates[0], burst)
    print('reading metadata from ' + path)
    os.chdir(path)

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    burst1 = 'burst_' + str(nBurst) + '/'
    burst2 = 'burst_' + str(nBurst + 1) + '/'

    # First get the data for the master and df_dc
    master_path = os.path.join(overlap_path, master_date)
    df_dc_path = os.path.join(overlap_path, 'df_dc')

    master_1 = master_file(burst)
    master_2 = master_file(next_burst)

    if not os.path.exists(master_path + '_1') or not os.path.exists(master_path + '_2'):
        master_1 = freadbk(burst1 + master_1, line_start, first_pixel_this, line_length,
                          pixel_length, 'cpxint16', this_nr_oflines, this_nr_ofpixels)
        master_2 = freadbk(burst2 + master_2, 1, first_pixel_next, line_length, pixel_length,
                          'cpxint16', next_nr_oflines, next_nr_ofpixels)
        master_1_file = np.memmap(master_path + '_1', 'complex64', shape=master_1.shape, mode='w+')
        master_2_file = np.memmap(master_path + '_2', 'complex64', shape=master_2.shape, mode='w+')
        master_1_file[:] = master_1
        master_2_file[:] = master_2
        master_1_file.flush()
        master_2_file.flush()

    if not os.path.exists(df_dc_path):
        df_dc = get_f_DC_difference(nBurst)
        df_dc_file = np.memmap(df_dc_path, 'float32', shape=df_dc.shape, mode='w+')
        df_dc_file[:,:] = df_dc[:,:]

    # Then loop over the slaves
    for date in dates:
        if date == master_date:
            continue

        path = swath_path(stack_folder, date, burst)
        os.chdir(path)

        print(path)
        data_path = os.path.join(overlap_path, date)

        burst1 = 'burst_' + str(nBurst) + '/'
        burst2 = 'burst_' + str(nBurst + 1) + '/'

        if not os.path.exists(data_path + '_1') or not os.path.exists(data_path + '_2'):
            slave_1 = freadbk(burst1 + 'slave_rsmp_reramped.raw', line_start, first_pixel_this, line_length, pixel_length , 'complex64',  this_nr_oflines, this_nr_ofpixels)
            slave_2 = freadbk(burst2 + 'slave_rsmp_reramped.raw', 1, first_pixel_next, line_length, pixel_length, 'complex64', next_nr_oflines, next_nr_ofpixels)
            slave_1_file = np.memmap(data_path + '_1', 'complex64', shape=slave_1.shape, mode='w+')
            slave_2_file = np.memmap(data_path + '_2', 'complex64', shape=slave_2.shape, mode='w+')
            slave_1_file[:] = slave_1
            slave_2_file[:] = slave_2
            slave_1_file.flush()
            slave_2_file.flush()


def find_ps_overlapping(stack_folder, overlap):
    # This is used to find the ps point in overlapping areas

    nBurst, burst, next_burst = get_burst(overlap)
    esd_folder = os.path.join(stack_folder, 'esd')
    overlap_path = os.path.join(esd_folder, overlap)
    files = os.listdir(os.path.join(overlap_path))
    dates = sorted([f[:-2] for f in files if (f.endswith('_1') and len(f) > 10)])
    folder = dates[0][0:4] + dates[0][5:7] + dates[0][8:10]
    os.chdir(os.path.join(stack_folder, folder, overlap[0:7]))

    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    # Remove earlier generated files for ps
    files = next(os.walk(overlap_path))[2]
    files = [os.path.join(overlap_path, name) for name in files if (name.endswith('ps') or name.startswith('ps'))]
    for filename in files:
        os.remove(filename)

    # Gather data in one matrix
    first, second = gather_stack(overlap_path, line_length, pixel_length)

    # First calculate the ps point for first overlap
    mean = np.mean(first, axis=2)
    std = np.std(first, axis=2)

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(std, mean)
        c[~ np.isfinite(c)] = 10000
    ps1 = (c < 0.3)
    ps1_file = os.path.join(overlap_path, 'ps_1')
    ps1_dat = np.memmap(ps1_file, 'bool', mode= 'w+', shape=(line_length, pixel_length))
    ps1_dat[:,:] = ps1[:,:]

    # Then calculate the ps point for second overlap
    mean = np.mean(second, axis=2)
    std = np.std(second, axis=2)

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(std, mean)
        c[~ np.isfinite(c)] = 10000
    ps2 = (c < 0.3)
    ps2_file = os.path.join(overlap_path, 'ps_2')
    ps2_dat = np.memmap(ps2_file, 'bool', mode= 'w+', shape=(line_length, pixel_length))
    ps2_dat[:,:] = ps2[:,:]

    ps_file = os.path.join(overlap_path, 'ps')
    ps_dat = np.memmap(ps_file, 'bool', mode='w+', shape=(line_length, pixel_length))
    ps_dat[:, :] = ((ps1_dat * ps2_dat) == 1)


def select_ps_data(stack_folder, overlap):
    # This function creates seperate files for the values of the ps points only.

    esd_folder = os.path.join(stack_folder, 'esd')
    overlap_path = os.path.join(esd_folder, overlap)
    files = os.listdir(os.path.join(overlap_path))
    dates = sorted([f[:-2] for f in files if (f.endswith('_1') and len(f) > 10)])

    folder = dates[0][0:4] + dates[0][5:7] + dates[0][8:10]
    os.chdir(os.path.join(stack_folder, folder, overlap[0:7]))

    nBurst = int(overlap.split('_')[3])
    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    ps_file = os.path.join(overlap_path, 'ps')
    ps_dat = np.memmap(ps_file, 'bool', mode='r', shape=(line_length, pixel_length))
    ps_num = np.sum(ps_dat)

    # Save only the ps points to file.
    for date in dates:
        slave_ps_name = os.path.join(overlap_path, date + '_1_ps')
        master_ps_name = os.path.join(overlap_path, date + '_2_ps')

        if not os.path.exists(slave_ps_name) or not os.path.exists(master_ps_name):
            slave_ps = np.memmap(slave_ps_name, 'complex64', mode='w+', shape=(ps_num))
            master_ps = np.memmap(master_ps_name, 'complex64', mode='w+', shape=(ps_num))
            slave = np.memmap(os.path.join(overlap_path, date + '_1'), 'complex64', mode='r',
                                     shape=(line_length, pixel_length))
            master = np.memmap(os.path.join(overlap_path, date + '_2'), 'complex64', mode='r',
                                    shape=(line_length, pixel_length))
            if ps_num > 0:
                slave_ps[:] = slave[ps_dat]
                master_ps[:] = master[ps_dat]
    # Do the same for the df_dc file
    if not os.path.exists(os.path.join(overlap_path, 'df_dc_ps')):
        df_dc_ps = np.memmap(os.path.join(overlap_path, 'df_dc_ps'), 'float32', mode='w+', shape=(ps_num))
        df_dc = np.memmap(os.path.join(overlap_path, 'df_dc'), 'float32', mode='r',
                          shape=(line_length, pixel_length))
        if ps_num > 0:
            df_dc_ps[:] = df_dc[ps_dat]


def network_esd_ps(stack_folder, overlap, master_date, max_baseline, max_offset=0.02):
    # Based on ps point esd is calculated using a network approach

    dates, overlap_path, diff_matrix, var_matrix, to_angle_matrix, weight_matrix, processing = prepare_esd(stack_folder, overlap)
    folder = dates[0][0:4] + dates[0][5:7] + dates[0][8:10]
    os.chdir(os.path.join(stack_folder, folder, overlap[0:7]))

    nBurst = int(overlap.split('_')[3])
    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    ps_file = os.path.join(overlap_path, 'ps')
    ps_dat = np.memmap(ps_file, 'bool', mode='r', shape=(line_length, pixel_length))

    ps_id = np.where(ps_dat == 1)
    if not ps_id:  # If there are no ps points
        return diff_matrix, var_matrix, weight_matrix, to_angle_matrix
    else:
        ps_num = len(ps_id[0])

    df_dc_ps = np.memmap(os.path.join(overlap_path, 'df_dc_ps'), 'float32', mode='r+', shape=ps_num)[:]

    for date, n in zip(dates, range(len(dates))):
        for date_2, num in zip(dates, range(len(dates))):
            # Only calculate the upper triangle, as the others will be the same
            if processing[n, num] == 1:
                continue

            timediff = datetime.strptime(date_2, '%Y-%m-%d') - datetime.strptime(date, '%Y-%m-%d')
            if timediff > timedelta(minutes=1) and timediff < timedelta(days=max_baseline):

                first_master = np.memmap(os.path.join(overlap_path, date + '_1_ps'), 'complex64', mode='r',
                                         shape=ps_num)
                first_slave = np.memmap(os.path.join(overlap_path, date_2 + '_1_ps'), 'complex64', mode='r',
                                        shape=ps_num)
                second_master = np.memmap(os.path.join(overlap_path, date + '_2_ps'), 'complex64', mode='r',
                                          shape=ps_num)
                second_slave = np.memmap(os.path.join(overlap_path, date_2 + '_2_ps'), 'complex64', mode='r',
                                         shape=ps_num)

                double_diff = (first_master * first_slave.conj()) * (second_master * second_slave.conj()).conj()

                # Now select all pixels with an offset of less than x milipixel
                double_diff[np.isnan(double_diff)] = 0.050
                val = (np.abs(np.angle(double_diff)) < max_offset)
                w = np.sum(val)

                if w > 0:
                    pixel_diff = np.angle(np.sum(double_diff[val])) * (PRF / (2 * np.pi * np.nanmean(df_dc_ps[val])))
                    pixel_var = np.var(np.angle(double_diff[val]) * (PRF/(2*np.pi*df_dc_ps[val])))
                    temp_baseline_w = np.exp(-(float(timediff.days) / 100))
                    weight_matrix[0, n, num] = temp_baseline_w * w
                    var_matrix[0, n, num] = pixel_var * temp_baseline_w
                    diff_matrix[0, n, num] = pixel_diff
                    # Phase ramp per pixel
                    to_angle_matrix[0, n, num] = (PRF/(2*np.pi*np.nanmean(df_dc_ps[val]))) * (line_start - 1)

    return diff_matrix, var_matrix, to_angle_matrix, weight_matrix, dates


def network_esd_coh(stack_folder, overlap, master_date, max_baseline, ra=10, az=2):

    dates, overlap_path, diff_matrix, var_matrix, to_angle_matrix, weight_matrix, processed = prepare_esd(stack_folder, overlap)
    folder = master_date[0:4] + master_date[5:7] + master_date[8:10] + "_" + dates[0][0:4] + dates[0][5:7] + dates[0][8:10]
    os.chdir(os.path.join(stack_folder, folder, overlap[0:7]))

    nBurst = int(overlap.split('_')[3])
    line_start, line_length, first_pixel_this, first_pixel_next, pixel_length, this_nr_oflines, \
    this_nr_ofpixels, next_nr_oflines, next_nr_ofpixels, PRF = get_coordinates(nBurst)

    # Gather data in one matrix
    first, second = gather_stack(overlap_path, line_length, pixel_length, abs=False)

    # Remove the empty rows / columns
    columns = np.where((np.min(np.abs(np.sum(first, axis=0)), axis=1) != 0) *
                       (np.min(np.abs(np.sum(second, axis=0)), axis=1) != 0) == True)[0]
    rows = np.where((np.min(np.abs(np.sum(first, axis=1)), axis=1) != 0) *
                       (np.min(np.abs(np.sum(second, axis=1)), axis=1) != 0) == True)[0]
    first = first[rows[0]:rows[-1]+1, columns[0]:columns[-1]+1, :]
    second = second[rows[0]:rows[-1]+1, columns[0]:columns[-1]+1, :]

    # Multilook the df_dc
    df_dc = np.memmap(os.path.join(overlap_path, 'df_dc'), 'float32', mode='r+', shape=(line_length, pixel_length))
    df_dc_ml = multilook(df_dc[rows[0]:rows[-1]+1, columns[0]:columns[-1]+1], az, ra)

    for date, n in zip(dates, range(len(dates))):

        # First select the dates we want to compare with
        c_dates = []
        nums = []
        for date_2, num in zip(dates, range(len(dates))):
            # Only calculate the upper triangle, as the others will be the same
            timediff = datetime.strptime(date_2, '%Y-%m-%d') - datetime.strptime(date, '%Y-%m-%d')

            if timediff > timedelta(minutes=1) and timediff < timedelta(days=max_baseline) and processed[n, num] != 1:
                c_dates.append(date_2)
                nums.append(num)

        if len(c_dates) != 0:
            # Then create ifgs of first and second
            shape_ifg = (first.shape[0], first.shape[1], len(nums))
            first_ifg = np.memmap(os.path.join(overlap_path, 'first_ifg'), 'complex64', shape=shape_ifg, mode='w+')
            first_ifg[:] = first[: ,: ,n][:, :, None] * first[:, :, nums].conj()
            second_ifg = np.memmap(os.path.join(overlap_path, 'second_ifg'), 'complex64', shape=shape_ifg, mode='w+')
            second_ifg[:] = second[:, :, n][:, :, None] * second[:, :, nums].conj()

            # And the double difference
            double_diff = np.memmap(os.path.join(overlap_path, 'double_diff'), 'complex64', shape=shape_ifg, mode='w+')
            double_diff[:] = first_ifg * second_ifg.conj()
            double_diff = multilook(double_diff, az, ra, summation=True)
            diff_phase = np.angle(double_diff)
            diff_amp = np.abs(double_diff)

            # Calculate coherence
            amp_sq = np.memmap(os.path.join(overlap_path, 'amp_sq'), 'float32', shape=shape_ifg, mode='w+')
            amp_sq[:] = (np.real(first_ifg) ** 2 + np.imag(first_ifg) ** 2) * (np.real(second_ifg) ** 2 + np.imag(second_ifg) ** 2)
            coh_amp_sq = multilook(amp_sq, az, ra, summation=True)
            coh = diff_amp / np.sqrt(coh_amp_sq)

            # Calculate weights
            coh[coh == 1] = 0
            weight = 2 * az * ra * coh**2 / (1 - coh**2)
            tot_weight = np.sum(np.sum(weight, axis=0), axis=0)

            to_angle = PRF / (2 * np.pi * df_dc_ml)
            shift_pix = diff_phase * to_angle[:, :, None]
            shift = np.sum(np.sum(shift_pix * weight, axis=0), axis=0) / tot_weight
            to_angle_weighted = np.sum(np.sum(to_angle[:, :, None] * weight, axis=0), axis=0) / tot_weight
            var = np.sum(np.sum((shift[None, None, :] - shift_pix)**2 * weight, axis=0), axis=0) / tot_weight

            var_matrix[0, n, nums] = var
            diff_matrix[0, n, nums] = shift
            to_angle_matrix[0, n, nums] = to_angle_weighted
            weight_matrix[0, n, nums] = tot_weight

    return diff_matrix, var_matrix, to_angle_matrix, weight_matrix, dates


def prepare_esd(stack_folder, overlap, esd_type='ps', load_existing=False):
    # Get some standard variables for the esd processing, used in both methods.

    esd_folder = os.path.join(stack_folder, 'esd')
    overlap_path = os.path.join(esd_folder, overlap)
    files = os.listdir(os.path.join(overlap_path))
    dates = sorted([f[:-2] for f in files if (f.endswith('_1') and len(f) > 10)])

    diff_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    var_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    to_angle_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    weight_matrix = np.zeros(shape=(1, len(dates), len(dates)))
    processed = np.zeros(shape=(len(dates), len(dates)))

    # Now load the existing results if they are available.
    if load_existing == True:
        diff_m = np.load(os.path.join(overlap_path, esd_type + '_diff_matrix.npy'))
        var_m = np.load(os.path.join(overlap_path, esd_type + '_var_matrix.npy'))
        to_angle_m = np.load(os.path.join(overlap_path, esd_type + '_to_angle_matrix.npy'))
        weight_m = np.load(os.path.join(overlap_path, esd_type + '_weight_matrix.npy'))
        old_dates = np.load(os.path.join(overlap_path, esd_type + '_dates.npy'))

        dates_overlap = [id for id, date in zip(range(len(dates)), old_dates) if date in dates]
        diff_matrix[np.ix_(1, dates_overlap, dates_overlap)] = diff_m
        var_matrix[np.ix_(1, dates_overlap, dates_overlap)] = var_m
        to_angle_matrix[np.ix_(1, dates_overlap, dates_overlap)] = to_angle_m
        weight_matrix[np.ix_(1, dates_overlap, dates_overlap)] = weight_m

        processed[np.ix_(dates_overlap, dates_overlap)] = 1

    return dates, overlap_path, diff_matrix, var_matrix, to_angle_matrix, weight_matrix, processed


def gather_stack(overlap_path, line_length, pixel_length, abs=True):
    # This function gathers all data from different date into 2 3D matrices

    files = os.listdir(overlap_path)
    first_files = [os.path.join(overlap_path, f) for f in files if f.endswith('_1') and len(f) > 10]
    first_name = os.path.join(overlap_path, 'first')
    if os.path.exists(first_name):
        os.remove(first_name)
    if abs:
        first = np.memmap(first_name, 'float32', shape=(line_length, pixel_length, len(first_files)), mode='w+')
    else:
        first = np.memmap(first_name, 'complex64', shape=(line_length, pixel_length, len(first_files)), mode='w+')

    for f, n in zip(first_files, range(len(first_files))):
        first_dat = np.memmap(f, 'complex64', mode='r', shape=(line_length, pixel_length))
        if abs:
            first[:, :, n] = np.abs(first_dat[:, :])
        else:
            first[:, :, n] = first_dat[:, :]

    second_files = [os.path.join(overlap_path, f) for f in files if f.endswith('_2') and len(f) > 10]
    second_name = os.path.join(overlap_path, 'second')
    if os.path.exists(second_name):
        os.remove(second_name)
    if abs:
        second = np.memmap(second_name, 'float32', shape=(line_length, pixel_length, len(second_files)), mode='w+')
    else:
        second = np.memmap(second_name, 'complex64', shape=(line_length, pixel_length, len(second_files)), mode='w+')

    for f, n in zip(second_files, range(len(second_files))):
        second_dat = np.memmap(f, 'complex64', mode='r', shape=(line_length, pixel_length))
        if abs:
            second[:, :, n] = np.abs(second_dat[:, :])
        else:
            second[:, :, n] = second_dat[:, :]

    return first, second


def multilook(matrix, az, ra, summation=False):
    # This function multilooks a matrix, which is either in 2D or 3D. In the case of 3D the third dimension is
    # considered the time dimension. If summation is True we do not average but sum the values in the multilooking area.
    # Multilooking always starts at the first range/azimuth pixel.

    # First downsample 2 * 10
    new_ra = matrix.shape[1] / ra
    new_az = matrix.shape[0] / az

    size = matrix.shape
    if len(size) == 3 and summation == False:
        matrix_multilook = matrix[:new_az * az, :new_ra * ra, :].reshape([new_az, az, new_ra, ra, size[2]]).mean(3).mean(1)
    elif len(size) == 2 and summation == False:
        matrix_multilook = matrix[:new_az * az, :new_ra * ra].reshape([new_az, az, new_ra, ra]).mean(3).mean(1)
    elif len(size) == 3 and summation == True:
        matrix_multilook = matrix[:new_az * az, :new_ra * ra, :].reshape([new_az, az, new_ra, ra, size[2]]).mean(3).mean(1)
    elif len(size) == 2 and summation == True:
        matrix_multilook = matrix[:new_az * az, :new_ra * ra].reshape([new_az, az, new_ra, ra]).mean(3).mean(1)
    else:
        print('matrix does not have the right size')
        return []

    matrix_multilook = matrix_multilook.astype(matrix_multilook.dtype, subok=False)

    return matrix_multilook


def get_burst(overlap):
    s = overlap.split('_')
    burst = s[0] + '_' + s[1] + '_' + s[2] + '_' + s[3]
    next_burst = s[4] + '_' + s[5] + '_' + s[6] + '_' + s[7]

    nBurst = int(s[3])

    return nBurst, burst, next_burst


def swath_path(stack_folder, date, key):
    date_folder = date[:4] + date[5:7] + date[8:10]
    swath_burst = key.split('_')
    file_path = os.path.join(stack_folder, date_folder, swath_burst[0] + '_' + swath_burst[1])

    return file_path


def master_file(key):
    # This function converts combinations of dates and keys to a datafile name
    string = '_iw_' + key[6] + '_burst_' + key[14:]
    string = 'master' + string + '.raw'

    return string

# Actually execute the code...
if __name__ == "__main__":

    ps_select = '1'
    if len(sys.argv) == 7:
        stack_folder  = sys.argv[1]
        overlap       = sys.argv[2]
        esd_type      = sys.argv[3]
        max_baseline  = sys.argv[4]
        master_date   = sys.argv[5]
        ps_select     = sys.argv[6]  # Use 1 if needed, use 0 if not
    elif len(sys.argv) == 6:
        stack_folder  = sys.argv[1]
        overlap       = sys.argv[2]
        esd_type      = sys.argv[3]
        max_baseline  = sys.argv[4]
        master_date   = sys.argv[5]
    else:
        sys.exit('usage: stack_folder type burst')

    print('stack folder is ' + stack_folder)
    print('burst is ' + overlap)
    print('ps select is ' + ps_select)
    print('max baseline is ' + max_baseline)
    print('master date is ' + master_date)
    print('type ESD is ' + esd_type)

    # first get all the dates from the stack:
    master_key = master_date[:4] + master_date[5:7] + master_date[8:]

    ifgs = [f for f in os.listdir(stack_folder) if len(f) == 8]
    dates = [f[:4] + '-' + f[4:6] + '-' + f[6:8] for f in ifgs if f != master_key]

    # Then run the overlap cutout / ps selection and
    save_overlapping(stack_folder, master_date, dates, overlap)

    # If we want to select ps points
    if ps_select == '1':
        find_ps_overlapping(stack_folder, overlap)

    # Get the esd results for the overlapping areas either based on ps or coherence
    max_baseline = int(max_baseline)
    if esd_type == 'ps':
        select_ps_data(stack_folder, overlap)
        diff_matrix, var_matrix, to_angle_matrix, weight, dates = network_esd_ps(stack_folder, overlap, master_date, max_baseline)
    elif esd_type == 'coh':
        diff_matrix, var_matrix, to_angle_matrix, weight, dates = network_esd_coh(stack_folder, overlap, master_date, max_baseline)
    else:
        sys.exit('Type should either be coh or ps')

    # And save them in the corresponding folder:
    folder = os.path.join(stack_folder, 'esd', overlap)

    np.save(os.path.join(folder, esd_type + '_diff_matrix'), diff_matrix)
    np.save(os.path.join(folder, esd_type + '_var_matrix'), var_matrix)
    np.save(os.path.join(folder, esd_type + '_to_angle_matrix'), to_angle_matrix)
    np.save(os.path.join(folder, esd_type + '_weight_matrix'), weight)
    np.save(os.path.join(folder, esd_type + '_dates'), dates)
