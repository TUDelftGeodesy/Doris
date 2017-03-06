import numpy as np
import os, sys
from ESD_functions import freadbk
from resdata import ResData


def correct_slaves(burst_key, stack_folder, master_date, cor_phase_ifg = 'cint_srd.raw', cor_phase_s = 'slave_corrected.raw'):
    # This function corrects a slave image using the earth and dem phase corrected ifg.

    master_folder = master_date[0:4] + master_date[5:7] + master_date[8:10] + '_'

    folders = next(os.walk(stack_folder))[1]
    folders = [fold for fold in folders if fold.startswith(master_folder) and len(fold) == 17]

    swath = burst_key[:10]
    burst = burst_key[11:]
    burst_folders = [os.path.join(f, swath, burst) for f in folders]

    # Load the master file:
    m_file = os.path.join(burst_folders[0], 'master_iw' + swath[5] + '_' + burst)
    m_res = os.path.join(burst_folders[0], 'master.res')
    res_dat = ResData(m_res)

    pixels = int(res_dat.processes['crop']['Last_line (w.r.t. original_image)']) - \
             int(res_dat.processes['crop']['First_line (w.r.t. original_image)']) + 1
    lines = int(res_dat.processes['crop']['Last_pixel (w.r.t. original_image)']) - \
            int(res_dat.processes['crop']['First_pixel (w.r.t. original_image)']) + 1
    m_data = freadbk(m_file, lines=lines, pixels=pixels, dt='cpxint')

    # Load the phase corrected files and divide by the master.
    for fold in burst_folders:
        c_ifg = freadbk(os.path.join(fold, cor_phase_ifg), dt='complex64', lines=lines, pixels=pixels)

        # Calculate the corrected slave and save to disk.
        c_slave_file = os.path.join(fold, cor_phase_s)
        c_slave_dat = np.memmap(c_slave_file, dtype='complex64', shape=(lines, pixels))
        c_slave_dat[:, :] = (c_ifg / m_data).conj()
        c_slave_file.flush()


# Actually execute the code...
if __name__ == "__main__":

    if len(sys.argv) == 4:
        stack_folder  = sys.argv[1]
        burst_key     = sys.argv[2]
        master_date   = sys.argv[3]
    else:
        sys.exit('usage: stack_folder type burst')

    print('phase corrected slave will be calculated')
    print('stack folder is ' + stack_folder)
    print('burst is ' + burst_key)

    # Run the script
    correct_slaves(burst_key, stack_folder, master_date)
