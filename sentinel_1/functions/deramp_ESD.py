# This script is used to remove the residual ramp in the ifgs of individual bursts based on ESD estimates.

import numpy as np
import os, sys

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

from sentinel_1.functions.get_ramp import freadbk
from sentinel_1.functions.resdata import ResData

def remove_ramp(burst_path, file, angle_pixel):
    # Remove ramp from burst

    res_file = os.path.join(burst_path, 'master.res')
    res_dat = ResData(res_file, 'master')
    crop = res_dat.processes['crop']
    lines = int(crop['Last_line (w.r.t. original_image)']) - int(crop['First_line (w.r.t. original_image)']) + 1
    pixels = int(crop['Last_pixel (w.r.t. original_image)']) - int(crop['First_pixel (w.r.t. original_image)']) + 1

    n = np.arange(lines)
    phase_diff = n * angle_pixel
    complex_diff = np.cos(phase_diff).astype('complex64') + 1j * np.sin(phase_diff).astype('complex64')

    dat_file = np.memmap(os.path.join(burst_path, file), dtype='complex64', mode='r+', shape=(lines, pixels))
    dat_file[:, :] = complex_diff * dat_file.conj()
    dat_file.flush()

if __name__ == "main":
    # If calling script directly we run the remove_ramp function.

    if len(sys.argv) == 4:
        burst_path      = sys.argv[1]
        file            = sys.argv[2]
        angle_pixel     = sys.argv[3]  # Use 1 if needed, use 0 if not
    else:
        sys.exit('usage: burst_folder, file, angle_per_pixel')

    print('burst folder is ' + burst_path)
    print('file we will deramp is ' + file)
    print('the angle per pixel is ' + angle_pixel)

    angle_pixel = float(angle_pixel)
    remove_ramp(burst_path, file, angle_pixel)

