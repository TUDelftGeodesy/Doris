# Based on the orbit of the swath the orbits of the individual burst is calculated.

from .orbit_coordinates import lph2xyz, xyz2ell, intrp_orbit
import os
import numpy as np
import collections
from datetime import datetime
from doris.doris_stack.functions.resdata import ResData
from shapely.geometry import Polygon


def burst_header(resID):

    meta = collections.OrderedDict()

    meta['row_1'] = ['===============================================\n']
    meta['MASTER RESULTFILE:'] = resID
    meta['Created by'] = 'G.Mulder TU Delft'
    meta['row_2'] = 'Doris (Delft o-o Radar Interferometric Software)'
    meta['Version'] = 'Version (2015) (For TOPSAR)'
    meta['FFTW library'] = 'used'
    meta['VECLIB library'] = 'not used'
    meta['LAPACK library'] = 'not used'
    meta['Compiled at'] = 'XXXXXXXX'
    meta['By GUN gcc'] = 'XXXXXXXX'
    meta['row_3'] = ['===============================================\n']

    return meta


def burst_readfiles(meta, burst_num, burst_center, burst_border, swath_data):
    # First copy swath metadata for burst and create a georef dict which stores information about the geo reference of
    # the burst.
    meta['Burst_number_index'] = str(burst_num)
    aux = meta['aux']
    aux['azimuthPRF'] = [meta['Pulse_Repetition_Frequency (computed, Hz)']]
    meta.pop('aux')

    # First find coordinates of center and optionally the corners
    meta['Scene_centre_longitude'] = str(burst_center[0])
    meta['Scene_centre_latitude'] = str(burst_center[1])
    meta['Scene_ul_corner_latitude'] = str(burst_border[0][1])
    meta['Scene_ur_corner_latitude'] = str(burst_border[1][1])
    meta['Scene_lr_corner_latitude'] = str(burst_border[2][1])
    meta['Scene_ll_corner_latitude'] = str(burst_border[3][1])
    meta['Scene_ul_corner_longitude'] = str(burst_border[0][0])
    meta['Scene_ur_corner_longitude'] = str(burst_border[1][0])
    meta['Scene_lr_corner_longitude'] = str(burst_border[2][0])
    meta['Scene_ll_corner_longitude'] = str(burst_border[3][0])

    # Find doppler centroid frequency and azimuth reference time
    doppler_times = [np.datetime64(aux['doppler_azimuth_Time'][i] + '-00') for i in range(len(aux['doppler_azimuth_Time']))]
    frequency_times = [np.datetime64(aux['doppler_azimuth_Time'][i] + '-00') for i in range(len(aux['doppler_azimuth_Time']))]
    burst_start_time = np.datetime64(aux['azimuthTimeStart'][burst_num-1] + '-00')
    meta['First_pixel_azimuth_time (UTC)'] = burst_start_time.astype(datetime).strftime('%Y-%b-%d %H:%M:%S.%f')

    # First index after start burst for doppler and azimuth
    doppler_id = np.where(doppler_times > burst_start_time)
    frequency_id = np.where(frequency_times > burst_start_time)

    # Assign DC values to metadata
    if len(doppler_id[0]) > 0:
        doppler_id = np.min(doppler_id)
        parameter = aux['dopplerCoeff'][doppler_id].split()
        meta['DC_reference_azimuth_time'] = np.datetime64(aux['doppler_azimuth_Time'][doppler_id] + '-00').astype(datetime).strftime('%Y-%b-%d %H:%M:%S.%f')
        meta['DC_reference_range_time'] = aux['doppler_range_Time'][doppler_id]
    else:
        parameter = ['0.0','0.0','0.0']
        meta['DC_reference_azimuth_time'] = '0.0 0.0'
        meta['DC_reference_range_time'] = aux['doppler_range_Time'][-1]
    # Assign parameters
    meta['Xtrack_f_DC_constant (Hz, early edge)'] = parameter[0]
    meta['Xtrack_f_DC_linear (Hz/s, early edge)'] = parameter[1]
    meta['Xtrack_f_DC_quadratic (Hz/s/s, early edge)'] = parameter[2]

    # Assign FM values to metadata
    if len(frequency_id[0]) > 0:
        frequency_id = np.min(frequency_id)
        if aux['azimuthFmRate_c0']:
            parameter = [aux['azimuthFmRate_c0'][frequency_id],aux['azimuthFmRate_c1'][frequency_id],aux['azimuthFmRate_c2'][frequency_id]]
        else:
            parameter = aux['azimuthFmRatePolynomial'][frequency_id].split()
        meta['FM_reference_azimuth_time'] = np.datetime64(aux['azimuthFmRate_reference_Azimuth_time'][frequency_id] + '-00').astype(datetime).strftime('%Y-%b-%d %H:%M:%S.%f')
        meta['FM_reference_range_time'] = aux['azimuthFmRate_reference_Range_time'][frequency_id]
    else:
        parameter = ['0.0','0.0','0.0']
        meta['doppler_azimuth_Time'] = '0.0 0.0'
        meta['FM_reference_range_time'] = aux['azimuthFmRate_reference_Range_time'][-1]
    # Assign parameters
    meta['FM_polynomial_constant_coeff (Hz, early edge)'] = parameter[0]
    meta['FM_polynomial_linear_coeff (Hz/s, early edge)'] = parameter[1]
    meta['FM_polynomial_quadratic_coeff (Hz/s/s, early edge)'] = parameter[2]

    # Add information about swath
    meta['row_1'] = ['******************************************************************']
    meta['Datafile'] = os.path.basename(swath_data)
    meta['Dataformat'] = 'tiff'
    meta['Number_of_lines_original'] = aux['imageLines'][0]
    meta['Number_of_pixels_original'] = aux['imagePixels'][0]
    meta['deramp'] = '0'
    meta['reramp'] = '0'
    meta['ESD_correct'] = '0'

    return meta


def burst_crop(meta,burst_num,swath_data,new_burst_num):
    # This function returns a description of the crop files part, which defines how the burst is cropped out of the
    # swath data. This function is generally called when the .res and raw data is written to the datastack and uses one
    # of the outputs from the burst_readfiles

    crop = collections.OrderedDict()

    last_sample =  [int(x) for x in meta['aux']['lastValidSample'][burst_num-1].split()]
    first_sample = [int(x) for x in meta['aux']['firstValidSample'][burst_num-1].split()]

    swath_data = os.path.basename(swath_data)
    crop['Data_output_file'] = 'slave_iw_' + swath_data[6] + '_burst_' + str(new_burst_num) + '.raw'
    crop['Data_output_format'] = 'complex_short'

    # Start line of this burst in total swath product
    lines = meta['aux']['imageLines'][0]
    extra_lines = int(lines) * (burst_num-1)

    # All coordinates are one based (e.g. we start at pixel 1)
    halfway = int(len(last_sample)/2)

    crop['First_line (w.r.t. tiff_image)'] = str(1 + last_sample[:halfway].count(-1) + extra_lines)
    crop['Last_line (w.r.t. tiff_image)'] = str(len(last_sample) - last_sample[halfway:].count(-1) + extra_lines)
    crop['First_line (w.r.t. original_image)'] = str(1 + last_sample[:halfway].count(-1))
    crop['Last_line (w.r.t. original_image)'] = str(len(last_sample) - last_sample[halfway:].count(-1))
    crop['First_pixel (w.r.t. original_image)'] = str(max(first_sample))
    crop['Last_pixel (w.r.t. original_image)'] = str(max(last_sample))

    return crop


def center_shape_from_res(resfile):
    # This function reads the shape and center of a burst from a .res file.

    res = ResData(resfile)
    res.res_read()
    meta = res.processes['readfiles']

    center = (float(meta['Scene_centre_longitude']), float(meta['Scene_centre_latitude']))
    ul = (float(meta['Scene_ul_corner_longitude']), float(meta['Scene_ul_corner_latitude']))
    ur = (float(meta['Scene_ur_corner_longitude']), float(meta['Scene_ur_corner_latitude']))
    lr = (float(meta['Scene_lr_corner_longitude']), float(meta['Scene_lr_corner_latitude']))
    ll = (float(meta['Scene_ll_corner_longitude']), float(meta['Scene_ll_corner_latitude']))

    coverage = Polygon([ul, ur, lr, ll])

    return center, coverage