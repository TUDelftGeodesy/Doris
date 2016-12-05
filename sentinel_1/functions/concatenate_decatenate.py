from sentinel_1.main_code.resdata import ResData
import numpy as np
import os, sys


def decatenate(date_folder, image_file, burst_file):
    # Decatenate data.
    master = os.path.join(date_folder, image_file)

    # Load .res files
    image_res, burst_res = read_res(date_folder)

    # Read image size
    bursts = burst_res.keys()
    no_lines = int(burst_res[bursts[0]][type].processes['readfiles']['Number_of_lines_output_image'])
    no_pixels = int(burst_res[bursts[0]][type].processes['readfiles']['Number_of_pixels_output_image'])

    # First use memmap to get a memory map of the full file.
    full_image = np.memmap(master, dtype=dt, mode='r', shape=(no_lines, no_pixels))

    for burst in bursts:
        # Finally write all data from individual bursts to master file. We assume a simple 20 pixel offset from
        # the side to prevent copying data without information.

        if burst_file == 'master':
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '.raw'
        elif burst_file == 'master_deramped':
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '.raw.orig'
        else:
            string = burst_file

        burst_dat = os.path.join(date_folder, burst[0:7], burst[9:], string)

        line_0 = int(burst_res[burst][type].processes['readfiles']['First_line (w.r.t. output_image)'])
        line_1 = int(burst_res[burst][type].processes['readfiles']['Last_line (w.r.t. output_image)'])
        pix_0 = int(burst_res[burst][type].processes['readfiles']['First_pixel (w.r.t. output_image)'])
        pix_1 = int(burst_res[burst][type].processes['readfiles']['Last_pixel (w.r.t. output_image)'])

        burst_pix = pix_1 - pix_0 + 1
        burst_line = line_1 - line_0 + 1

        # Cut out data with border of 20 px and write to file.
        burst_image = np.memmap(burst_dat, dtype=dt, mode='w+', shape=(burst_line,burst_pix))
        burst_image[:,:] = full_image[line_0-1:line_1,pix_0-1:pix_1]
        burst_image.flush()


def concatenate(date_folder, image_file, burst_file):
    # Concatenate data.
    master = os.path.join(date_folder, image_file)

    # Load .res files
    image_res, burst_res = read_res(date_folder)

    # Read image size
    bursts = burst_res.keys()
    no_lines = int(burst_res[bursts[0]][type].processes['readfiles']['Number_of_lines_output_image'])
    no_pixels = int(burst_res[bursts[0]][type].processes['readfiles']['Number_of_pixels_output_image'])

    # First use memmap to get a memory map of the full file.
    full_image = np.memmap(master, dtype=dt, mode='w+', shape=(no_lines, no_pixels))

    for burst in bursts:
        # Finally write all data from individual bursts to master file. We assume a simple 20 pixel offset from
        # the side to prevent copying data without information.

        if burst_file == 'master':
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '.raw'
        elif burst_file == 'master_deramped':
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '.raw.orig'
        else:
            string = burst_file

        burst_dat = os.path.join(date_folder, burst[0:7], burst[9:], string)

        line_0 = int(burst_res[burst][type].processes['readfiles']['First_line (w.r.t. output_image)'])
        line_1 = int(burst_res[burst][type].processes['readfiles']['Last_line (w.r.t. output_image)'])
        pix_0 = int(burst_res[burst][type].processes['readfiles']['First_pixel (w.r.t. output_image)'])
        pix_1 = int(burst_res[burst][type].processes['readfiles']['Last_pixel (w.r.t. output_image)'])

        burst_pix = pix_1 - pix_0 + 1
        burst_line = line_1 - line_0 + 1

        # Cut out data with border of 20 px and write to file.
        burst_image = np.memmap(burst_dat, dtype=dt, mode='r', shape=(burst_line,burst_pix))
        full_image[(line_0+19):(line_1-20),(pix_0+19):(pix_1-20)] = burst_image[20:-20,20:-20]


def read_res(date_folder):
    # Read .res data to the burst objects. Generally done after a processing step.

    swaths = next(os.walk(os.path.join(date_folder, fold)))[1]
    swaths = [fol for fol in swaths if len(fol) == 10]

    res_burst = dict()
    res_image = dict()

    for swath in swaths:

        bursts = next(os.walk(os.path.join(self.stack_folder, fold, swath)))[1]
        bursts = [burst for burst in bursts if burst.startswith('burst')]

        for burst in bursts:
            slave_res = os.path.join(date_folder, swath, burst, 'slave.res')
            master_res = os.path.join(date_folder, swath, burst, 'master.res')
            ifgs_res = os.path.join(date_folder, swath, burst, 'ifgs.res')

            burst_name = swath + '_' + burst
            res_burst[burst_name] = dict()

            if os.path.exists(slave_res):
                res_burst[burst_name]['slave'] = ResData(filename=slave_res)
            if os.path.exists(master_res):
                res_burst[burst_name]['master'] = ResData(filename=master_res)
            if os.path.exists(ifgs_res):
                res_burst[burst_name]['ifgs'] = ResData(filename=ifgs_res)

    slave_res = os.path.join(date_folder, 'slave.res')
    master_res = os.path.join(date_folder, 'master.res')
    ifgs_res = os.path.join(date_folder, 'ifgs.res')

    if os.path.exists(slave_res):
        res_image['slave'] = ResData(filename=slave_res)
    if os.path.exists(master_res):
        res_image['master'] = ResData(filename=master_res)
    if os.path.exists(ifgs_res):
        res_image['ifgs'] = ResData(filename=ifgs_res)

    return res_image, res_burst

# Actually execute the code...
if __name__ == "__main__":

    if len(sys.argv) == 4:
        date_folder         = sys.argv[1]
        type                = sys.argv[2]
        image_file          = sys.argv[3]
        burst_file          = sys.argv[4]
    else:
        sys.exit('usage: date_folder type image_file burst_file')


    if type == 'decatenate':
        decatenate(date_folder, image_file, burst_file)
    elif type == 'concatenate':
        concatenate(date_folder, image_file, burst_file)
    else:
        sys.exit('type should either be decatenate or concatenate')
