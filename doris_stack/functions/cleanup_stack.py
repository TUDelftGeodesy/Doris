# This script can be used to cleanup a datastack
# In this script we use the following rational
# 1. Data from the master date and ESD are always kept (can change when these scripts evolve)
# 2. The last processing step is saved (or if different processing steps are connected and the last one is not finished yet
# 3. All .res and .ras files are kept.
# 4. Following steps are never deleted > resampled slave / subt_refdem / filtphase (multilooked) / coherence (multilooked) / unwrap

import os
import sys
from doris.doris_stack.main_code.resdata import ResData

def res_file_selection(path, dat_type='burst'):
    # Give an oversight of all files connected to different steps:

    del_files = []
    res_dat = dict()

    # First check if slave / master / ifgs res files are there
    for res_type in ['master', 'slave', 'ifgs']:
        if os.path.exists(os.path.join(path, res_type + '.res')):
            res_dat[res_type] = ResData(os.path.join(path, res_type + '.res'), type=res_type)
        else:
            print('No data found in ' + path)
            return []

    if res_dat['slave'].process_control['readfiles'] == '1' and res_dat['slave'].process_control['crop'] == '1':
        if res_dat['slave'].processes['readfiles']['deramp'] == '1':
            # Remove the ramped image data.
            del_files.append(res_dat['slave'].processes['crop']['Data_output_file'][:-12] + '.raw')
    if res_dat['ifgs'].process_control['dem_assist'] == '1':
        # Remove the temporary files for dem_assist
        del_files.extend(['dac_delta_demline.temp', 'dac_delta_dempixel.temp', 'dac_m_demline.temp', 'dac_m_dempixel.temp'])
    if res_dat['slave'].process_control['resample'] == '1' and res_dat['ifgs'].process_control['dem_assist'] == '1':
        del_files.extend(['dac_delta_line.raw', 'dac_delta_pixel.raw'])
    if res_dat['slave'].process_control['resample'] == '1':
        # After resampling the deramped file is not needed anymore. Also processing data from resampling is not needed.
        del_files.append(res_dat['slave'].processes['crop']['Data_output_file'])
        del_files.extend(['rsmp_orig_slave_line.raw', 'rsmp_orig_slave_pixel.raw'])

    # Now the resampled slave stays.
    if res_dat['ifgs'].process_control['subtr_refphase'] == '1':
        # If the reference phase is subtracted both the reramped slave and interferogram can be removed
        del_files.extend(['cint.raw', 'slave_rsmp_reramped.raw'])
    if res_dat['ifgs'].process_control['subtr_refdem'] == '1':
        # When the dem phase is removed, the interferogram with subtracted reference phase can be removed.
        del_files.extend(['cint.raw', 'cint_srp.raw', 'demcrop.raw', 'dem_radar.raw', 'master_slave.crd'])
    if res_dat['ifgs'].process_control['filtphase'] == '1':
        # When after the the removal of the reference dem the filtphase step is done, the subtrefdem ifgs is removed.
        del_files.extend(['cint.raw', 'cint_srp.raw', 'cint_srd.raw', 'demcrop.raw', 'dem_radar.raw', 'master_slave.crd'])

    # Finally if it is about a burst image, we can check whether the filtered and coherence files are already
    # concatenated, which means they can be removed.
    image_ifgs = os.path.join(os.path.dirname(os.path.dirname(path)), 'ifgs.res')
    if dat_type == 'burst' and os.path.exists(image_ifgs):
        # If there is a ifgs.res file for the full image.

        image_res = ResData(image_ifgs)
        if image_res.process_control['coherence'] == '1':
            del_files.append('coherence.raw')
        if image_res.process_control['filtphase'] == '1':
            del_files.append('cint.0.2filtered')
        if image_res.process_control['filtphase'] == '1':
            del_files.append('cint.0.2filtered')

    elif dat_type == 'image':
        # If we are looking at the full image, is it important to remove the non-multilooked coherence or filtered image
        if res_dat['ifgs'].process_control['filtphase'] == '1':
            if res_dat['ifgs'].processes['filtphase']['Data_output_file'].endswith('ml.raw'):
                del_files.append('cint.0.2filtered')
        if res_dat['ifgs'].process_control['coherence'] == '1':
            if res_dat['ifgs'].processes['coherence']['Data_output_file'].endswith('ml.raw'):
                del_files.append('coherence.raw')

    del_files = set(del_files)
    del_files = [os.path.join(path, del_file) for del_file in del_files]
    del_files = [del_file for del_file in del_files if os.path.exists(del_file)]

    return del_files

def cleanup_stack(path, master_key):
    # This is the main cleanup function.

    folders = next(os.walk(path))[1]
    if not master_key in folders:
        print('master folder not found in path')
        return
    else:
        folders.remove(master_key)

    del_files = []

    for folder in folders:
        del_dat = res_file_selection(os.path.join(path, folder), dat_type='image')
        del_files.extend(del_dat)

        burst_folders = find_burst_folders(os.path.join(path, folder))
        for burst_folder in burst_folders:
            del_dat = res_file_selection(os.path.join(path, burst_folder), dat_type='burst')
            del_files.extend(del_dat)

    for filename in del_files:
        print('removing: ' + filename)
        os.remove(filename)
    print('Succes! path ' + path + ' is cleaned from temporary files.')


def find_burst_folders(folder):
    # This function finds all the burst folders in an image folder
    folders = []
    swaths = next(os.walk(folder))[1]
    for swath in swaths:
        bursts = next(os.walk(os.path.join(folder, swath)))[1]
        for burst in bursts:
            folders.append(os.path.join(folder, swath, burst))
    return folders

# Actually execute the code...
if __name__ == "__main__":
    path = sys.argv[1]
    master_key = sys.argv[2]

    print('path to be cleaned ' + path)
    print('master key is ' + master_key)

    cleanup_stack(path, master_key)
