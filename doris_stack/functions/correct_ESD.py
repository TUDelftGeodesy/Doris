# This script is used to remove the residual ramp in the ifgs of individual bursts based on ESD estimates.

import numpy as np
import os, sys

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

from .resdata import ResData


def remove_ramp(file, angle_pixel):
    # Remove ramp from burst

    res_file = 'master.res'
    res_dat = ResData(res_file, 'master')
    crop = res_dat.processes['crop']
    lines = int(crop['Last_line (w.r.t. original_image)']) - int(crop['First_line (w.r.t. original_image)']) + 1
    pixels = int(crop['Last_pixel (w.r.t. original_image)']) - int(crop['First_pixel (w.r.t. original_image)']) + 1

    n = np.arange(lines)
    phase_diff = n * angle_pixel
    complex_diff = np.cos(phase_diff).astype('complex64') + 1j * np.sin(phase_diff).astype('complex64')

    dat_file = np.memmap(file, dtype='complex64', mode='r+', shape=(lines, pixels))
    p_before = np.nanmean(np.angle(dat_file))
    print('Average phase before is ' + str(p_before))

    dat_file[:, :] = dat_file * complex_diff.conj()[:, None]
    p_after = np.nanmean(np.angle(dat_file))
    print('Average phase after is ' + str(p_after))
    dat_file.flush()

if __name__ == "__main__":
    # If calling script directly we run the remove_ramp function.

    if len(sys.argv) == 3:
        file            = sys.argv[1]
        angle_pixel     = sys.argv[2]  # Use 1 if needed, use 0 if not
    else:
        sys.exit('usage: burst_folder, file, angle_per_pixel')

    print('file we will deramp is ' + file)
    print('the angle per pixel is ' + angle_pixel)

    angle_pixel = float(angle_pixel)
    remove_ramp(file, angle_pixel)

