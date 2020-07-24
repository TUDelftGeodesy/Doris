from doris.doris_stack.main_code.resdata import ResData
import numpy as np
import os, sys


def decatenate(date_folder, image_file, burst_file, datatype, multilooked='none', res_type='main'):
    # Decatenate data.
    # Multilooks can be none > not concatenated, only > only multilooked images are concatenated or >
    # all > both original and multilooked images are concatenated.

    if len(multilooked) == 7:
        main = os.path.join(date_folder, image_file[:-4] + '_' + multilooked + '.raw')
    else:
        main = os.path.join(date_folder, image_file)

    # Load .res files
    image_res, burst_res = read_res(date_folder, type=res_type)

    # Read image size
    bursts = burst_res.keys()
    if multilooked != 'none':
        try:
            no_lines = int(burst_res[bursts[0]].processes['readfiles']['Number_of_ml_lines_output_image'])
            no_pixels = int(burst_res[bursts[0]].processes['readfiles']['Number_of_ml_pixels_output_image'])
        except:
            print('Not able to load multilooking parameters for ' + image_file)
            return
    else:
        no_lines = int(burst_res[bursts[0]].processes['readfiles']['Number_of_lines_output_image'])
        no_pixels = int(burst_res[bursts[0]].processes['readfiles']['Number_of_pixels_output_image'])

    # First use memmap to get a memory map of the full file.
    full_image = np.memmap(main, dtype=datatype, mode='r', shape=(no_lines, no_pixels))

    for burst in bursts:
        # Finally write all data from individual bursts to main file. We assume a simple 20 pixel offset from
        # the side to prevent copying data without information.

        burst_dat, line_0, line_1, pix_0, pix_1, burst_pix, burst_line, az_offset, ra_offset = \
            burst_info(burst, burst_file, burst_res, multilooked)

        # Cut out data with border of 20 px and write to file.
        burst_image = np.memmap(burst_dat, dtype=datatype, mode='w+', shape=(burst_line, burst_pix))
        burst_image[:, :] = full_image[line_0-1:line_1, pix_0-1:pix_1]
        burst_image.flush()


def concatenate(date_folder, image_file, burst_file, datatype, multilooked='none', res_type='main'):
    # Concatenate data.

    if len(multilooked) == 7:
        main = os.path.join(date_folder, image_file[:-4] + '_' + multilooked + '.raw')
    else:
        main = os.path.join(date_folder, image_file)

    # Load .res files
    image_res, burst_res = read_res(date_folder, type=res_type)

    # Read image size
    bursts = burst_res.keys()
    if multilooked != 'none':
        try:
            no_lines = int(burst_res[bursts[0]].processes['readfiles']['Number_of_ml_lines_output_image'])
            no_pixels = int(burst_res[bursts[0]].processes['readfiles']['Number_of_ml_pixels_output_image'])
        except:
            print('Not able to load multilooking parameters for ' + image_file)
            return
    else:
        no_lines = int(burst_res[bursts[0]].processes['readfiles']['Number_of_lines_output_image'])
        no_pixels = int(burst_res[bursts[0]].processes['readfiles']['Number_of_pixels_output_image'])

    # First use memmap to get a memory map of the full file.
    full_image = np.memmap(main, dtype=datatype, mode='w+', shape=(no_lines, no_pixels))

    for burst in bursts:
        # Finally write all data from individual bursts to main file. We assume a simple 20 pixel offset from
        # the side to prevent copying data without information. (corrected with multilooking factor)

        burst_dat, line_0, line_1, pix_0, pix_1, burst_pix, burst_line, daz, dra = \
            burst_info(burst, burst_file, burst_res, multilooked)

        # Cut out data with border of 20 px and write to file.
        burst_image = np.memmap(burst_dat, dtype=datatype, mode='r', shape=(burst_line,burst_pix))
        full_image[(line_0+(daz-1)):(line_1-daz), (pix_0+(dra-1)):(pix_1-dra)] = burst_image[daz:-daz, dra:-dra]


def burst_info(burst, burst_file, burst_res, multilooked='none'):
    # Information about this specific burst

    if burst_file == 'main.raw':
        if multilooked:
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '_' + multilooked + '.raw'
        else:
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '.raw'
    elif burst_file == 'main_deramped.raw':
        if multilooked:
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '_deramped' + '_' + multilooked + '.raw'
        else:
            string = '_iw_' + burst[6] + '_burst_' + burst[17:] + '_deramped.raw'
    else:
        string = burst_file

    if len(multilooked) == 7:
        burst_dat = os.path.join(date_folder, burst[0:7], burst[8:], string[:-4] + '_' + multilooked + '.raw')
        line_0 = int(burst_res[burst].processes['readfiles']['First_line_' + multilooked])
        line_1 = int(burst_res[burst].processes['readfiles']['Last_line_ml_' + multilooked])
        pix_0 = int(burst_res[burst].processes['readfiles']['First_pixel_ml_' + multilooked])
        pix_1 = int(burst_res[burst].processes['readfiles']['Last_pixel_ml_' + multilooked])
        ra = int(burst_res[burst].processes['readfiles']['Multilook_range_' + multilooked])
        az = int(burst_res[burst].processes['readfiles']['Multilook_azimuth_' + multilooked])
    else:
        burst_dat = os.path.join(date_folder, burst[0:7], burst[8:], string)
        line_0 = int(burst_res[burst].processes['readfiles']['First_line (w.r.t. output_image)'])
        line_1 = int(burst_res[burst].processes['readfiles']['Last_line (w.r.t. output_image)'])
        pix_0 = int(burst_res[burst].processes['readfiles']['First_pixel (w.r.t. output_image)'])
        pix_1 = int(burst_res[burst].processes['readfiles']['Last_pixel (w.r.t. output_image)'])
        ra = 1
        az = 1

    az_offset = int(np.ceil(20 / float(az)))
    ra_offset = int(np.ceil(100 / float(ra)))

    burst_pix = pix_1 - pix_0 + 1
    burst_line = line_1 - line_0 + 1

    return burst_dat, line_0, line_1, pix_0, pix_1, burst_pix, burst_line, az_offset, ra_offset


def read_res(date_folder, type='main'):
    # Read .res data to the burst objects. Generally done after a processing step.

    swaths = next(os.walk(date_folder))[1]
    swaths = [fol for fol in swaths if len(fol) == 7]

    res_burst = dict()

    for swath in swaths:

        bursts = next(os.walk(os.path.join(date_folder, swath)))[1]
        bursts = [burst for burst in bursts if burst.startswith('burst')]

        for burst in bursts:
            subordinate_res = os.path.join(date_folder, swath, burst, 'subordinate.res')
            main_res = os.path.join(date_folder, swath, burst, 'main.res')

            burst_name = swath + '_' + burst

            if type == 'main' and os.path.exists(main_res):
                res_burst[burst_name] = ResData(filename=main_res)
            elif type == 'subordinate' and os.path.exists(subordinate_res):
                res_burst[burst_name]= ResData(filename=subordinate_res)
            elif os.path.exists(main_res):
                res_burst[burst_name] = ResData(filename=main_res)
            elif os.path.exists(subordinate_res):
                res_burst[burst_name]= ResData(filename=subordinate_res)
            else:
                print('No burst main or subordinate image available')
                return

    subordinate_res = os.path.join(date_folder, 'subordinate.res')
    main_res = os.path.join(date_folder, 'main.res')

    if type == 'main' and os.path.exists(main_res):
        res_image = ResData(filename=main_res)
    elif type == 'subordinate' and os.path.exists(subordinate_res):
        res_image = ResData(filename=subordinate_res)
    elif os.path.exists(main_res):
        res_image = ResData(filename=main_res)
    elif os.path.exists(subordinate_res):
        res_image = ResData(filename=subordinate_res)
    else:
        print('No image main or subordinate image available')
        return

    return res_image, res_burst

# Actually execute the code...
if __name__ == "__main__":

    date_folder         = sys.argv[1]
    type                = sys.argv[2]
    burst_file          = sys.argv[3]
    datatype            = sys.argv[4]
    if len(sys.argv) > 5:
        multilooked         = sys.argv[5]
    else:
        multilooked = 'False'
    if len(sys.argv) > 6:
        res_type = sys.argv[6]
    else:
        res_type = 'main'

    print('concatenate folder is ' + date_folder)
    print('burst_file is ' + burst_file)
    print('datatype is ' + datatype)
    print('concatenating multilooked image is ' + multilooked)

    if len(multilooked) == 7:
        # Multilooked should be given in a 7 length string azimuth multilook _ range multiloop
        # example '004_020'
        # The script detects whether this multilooking is available. Otherwise it will produce an error.
        multilooked = multilooked
    else:
        multilooked = 'none'
    image_file = burst_file

    if datatype == 'cpxint16' or datatype == 'complex_short':
        datatype = np.dtype([('re', np.int16), ('im', np.int16)])

    if type == 'decatenate':
        decatenate(date_folder, image_file, burst_file, datatype, multilooked, res_type)
    elif type == 'concatenate':
        concatenate(date_folder, image_file, burst_file, datatype, multilooked, res_type)
    else:
        sys.exit('type should either be decatenate or concatenate')
