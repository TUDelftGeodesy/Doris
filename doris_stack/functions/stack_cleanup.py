"""
This script is meant to cleanup datastacks and remove unnecessary intermediate results. The most general options are the
settings for a PS approach and a distributed scatterers approach, which is also used for atmospheric observations.
"""

import os
import shutil

# Test:
stack_folder = '/media/gert/Data/datastacks/netherlands/asc_t88/stack'
cleanup_ps = False
cleanup_ds = True
full_swath_rm = []
full_swath_keep = []
burst_rm = []
burst_keep = []
remove = False
files = cleanup(stack_folder, cleanup_ps, cleanup_ds, full_swath_rm, full_swath_keep, burst_rm, burst_keep, remove)

def cleanup(stack_folder, cleanup_ps=True, cleanup_ds=False, full_swath_rm=[], full_swath_keep=[], burst_rm=[],
            burst_keep=[], remove=True):
    """
    This is the main script to decide what should be removed and what not for an individual full swath.
    Basically a list is made of the files that should be removed from the full swath and burst folders based on the
    ending of the file name.
    The function returns a dictionary where are deleted files are listed. If you only want to list all the files that
    will be deleted use the remove=False option.

    Default options are:

    remove                          abbreviation            PS          DS
    interferogram                   ifg                     yes         yes
    ifg earth phase corrected       ifg_srp                 yes         yes
    ifg dem phase corrected         ifg_srd                 no          yes
    ifg phase filtered              ifg_filt                NA (yes)    no
    ifg coherence                   ifg_coh                 yes         no
    ifg unwrapped                   ifg_unw                 NA (yes)    no
    dem pixel shift (pixel)         dac_delta_p             yes         yes
    dem pixel shift (line)          dac_delta_l             yes         yes
    subordinate image                     s_ramp                  yes         yes
    subordinate image deramped            s_deramp                no          yes
    main image                    m_ramp                  yes         yes
    main image deramped           m_deramp                no (*)      yes
    dem phase                       dem                     no          no
    latitude                        phi                     no (*)      no
    longitude                       lam                     no (*)      no

    burst folder                    b_folder                yes         yes
        burst raw files             b_raw                   yes         yes
        burst ras files             b_ras                   yes         yes
        burst res files             b_res                   yes         yes

    * Only one needed in a single main stack. Is not implemented yet.
    """

    # First check what should be removed.
    if cleanup_ps:
        swath_clean = {'ifg': True, 'ifg_srp': True, 'ifg_srd': False, 'ifg_filt': True, 'ifg_coh': True,
                       'ifg_unw': True, 's_ramp': True, 's_deramp': False, 'm_ramp': True, 'm_deramp': False,
                       'dem': False, 'phi': False, 'lam': False, 'dac_delta_p': True, 'dac_delta_l': True}
        burst_clean = {'b_folder': True, 'b_raw': True, 'b_ras': True, 'b_res': True}
    elif cleanup_ds:
        swath_clean = {'ifg': True, 'ifg_srp': True, 'ifg_srd': True, 'ifg_filt': False, 'ifg_coh': False,
                       'ifg_unw': False, 's_ramp': True, 's_deramp': True, 'm_ramp': True, 'm_deramp': True,
                       'dem': False, 'phi': False, 'lam': False, 'dac_delta_p': True, 'dac_delta_l': True}
        burst_clean = {'b_folder': True, 'b_raw': True, 'b_ras': True, 'b_res': True}
    else:  # Otherwise nothing is removed unless indicated
        swath_clean = {'ifg': False, 'ifg_srp': False, 'ifg_srd': False, 'ifg_filt': False, 'ifg_coh': False,
                       'ifg_unw': False, 's_ramp': False, 's_deramp': False, 'm_ramp': False, 'm_deramp': False,
                       'dem': False, 'phi': False, 'lam': False, 'dac_delta_p': False, 'dac_delta_l': False}
        burst_clean = {'b_folder': False, 'b_raw': False, 'b_ras': False, 'b_res': False}

    # Add additional options.
    for dat in full_swath_rm:
        swath_clean[dat] = True
    for dat in full_swath_keep:
        swath_clean[dat] = False
    for dat in burst_rm:
        burst_clean[dat] = True
    for dat in burst_keep:
        burst_clean[dat] = False

    # Then create the strings with which these parts end
    swath_endings = {'ifg': 'int.raw', 'ifg_srp': 'srp.raw', 'ifg_srd': 'srd.raw', 'ifg_filt': 'filt.raw',
                     'ifg_coh': 'rence.raw', 'ifg_unw': 'unwrapped.raw', 's_ramp': 'rsmp.raw',
                     's_deramp': 'rsmp_deramped.raw', 'm_ramp': 'ster.raw', 'm_deramp': 'ster_deramped.raw',
                     'dem': 'dem_radar.raw', 'phi': 'phi.raw', 'lam': 'lam.raw',
                     'dac_delta_p': 'delta_pixel.raw', 'dac_delta_l': 'delta_line.raw'}
    burst_endings = {'b_folder': '', 'b_raw': '.raw', 'b_ras': '.ras', 'b_res': '.res'}

    # Finally, make a list of which endings should be deleted
    swath_remove = [dat for key, dat in swath_endings.iteritems() if swath_clean[key]]
    burst_remove = [dat for key, dat in burst_endings.iteritems() if burst_clean[key]]

    # Check the total ifgs in the stack
    swath_folders = scan_stack(stack_folder)

    # Go over these ifgs and remove intermediate steps in full swath and bursts.
    deleted = dict()
    for swath_folder in swath_folders:
        deleted[swath_folder] = dict()

        if burst_clean['b_folder']:
            folders = remove_burst_folders(swath_folder, remove)
            deleted[swath_folder]['folders'] = folders
        else:
            filenames = remove_burst_files(swath_folder, burst_remove, remove)
            deleted[swath_folder]['burst_files'] = filenames

        filenames = remove_file(swath_folder, swath_remove, remove)
        deleted[swath_folder]['swath_files'] = filenames

    return deleted


def scan_stack(stack_folder):
    # This function enters the children directories and checks whether a main.res, ifgs.res, swath folder and ifg
    # exist.

    swath_folders = []
    root, dirs, files = os.walk(stack_folder).next()

    for folder in dirs:
        r, folders, files = os.walk(os.path.join(root, folder)).next()

        if 'swath_1' in folders and 'main.res' in files and 'ifgs.res' in files and 'cint.raw' in files:
            swath_folders.append(os.path.join(root, folder))

    return swath_folders


def remove_burst_folders(swath_folder, remove):
    # Remove all burst folders from swath folder

    folder_names = []
    root, dirs, files = os.walk(swath_folder).next()

    for folder in dirs:
        if folder.startswith('swath'):
            if remove:
                shutil.rmtree(os.path.join(root, folder))
            folder_names.append(os.path.join(root, folder))

    return folder_names


def remove_file(swath_folder, file_endings, remove):
    # Remove the files in the main folder.

    file_names = []
    root, dirs, files = os.walk(swath_folder).next()

    for filename in files:
        for end in file_endings:
            if filename.endswith(end):
                if remove:
                    os.remove(os.path.join(root, filename))
                file_names.append(os.path.join(root, filename))

    return file_names


def remove_burst_files(swath_folder, file_endings, remove):
    # Remove the files in the burst folders.

    file_names = []
    root1, swaths, files = os.walk(swath_folder).next()

    if len(swaths) == 0:
        'Files seems to be deleted already!'
        return file_names

    for swath in swaths:
        root2, bursts, files = os.walk(os.path.join(root1, swath)).next()
        for burst in bursts:
            root3, burst_fold, files = os.walk(os.path.join(root2, burst)).next()
            for filename in files:
                for end in file_endings:
                    if filename.endswith(end) and remove:
                        if remove:
                            os.remove(os.path.join(root3, filename))
                        file_names.append(os.path.join(root3, filename))

    return file_names
