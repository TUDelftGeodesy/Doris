# Based on the orbit of the swath the orbits of the individual burst is calculated.
from sentinel_1.functions.orbit_coordinates import lph2xyz, xyz2ell, intrp_orbit
import os
import numpy as np
import collections
from datetime import datetime
import time
from shapely.geometry import Polygon
from sentinel_1.functions.precise_read import interpolate_orbit, orbit_read


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


def pixel_location(orbit_info,line,pixel,burst_num):
    # This function uses orbit + lat/lon/height information of a known location, to calculate the lat/lon/height of
    # another location. Info of the first is stored in the orbit_info. Pixel location is given by variables line and pixel.

    norm_orbit, norm_orbit_line = intrp_orbit(line,orbit_info,burst_num)
    coord_xyz = lph2xyz(line, pixel, orbit_info, norm_orbit_line,
                        orbit_info['centroid_lon'], orbit_info['centroid_lat'], orbit_info['height'])
    phi_lam_height = xyz2ell(coord_xyz)

    return phi_lam_height


def burst_location(aux,burst_num,corners=False):

    # Create orbit dict to interpolate burst lat/lon/height locations
    orbit_info = dict()
    orbit_info['orbitTime'] = aux['orbitTime']
    orbit_info['orbitX'] = aux['orbitX']
    orbit_info['orbitY'] = aux['orbitY']
    orbit_info['orbitZ'] = aux['orbitZ']
    orbit_info['azimuthPRF'] = aux['azimuthPRF']
    orbit_info['azimuthTimeStart'] = aux['azimuthTimeStart']
    orbit_info['rangeTimePix'] = aux['rangeTimePix']
    orbit_info['rangeRSR'] = aux['rangeRSR']

    # Calculate centre of burst based on swath orbit information.
    # Image centroid parameters
    line_location_index_down    = burst_num * np.array(aux['imageLines'][0], 'int') # Last line
    line_location_index_up      = (burst_num - 1) * np.array(aux['imageLines'][0], 'int') # First line may be not accurate
    cen_line                    = (line_location_index_down + line_location_index_up) / 2 - line_location_index_up
    cen_pixel                   = np.array(aux['imagePixels'][0], 'int') / 2

    # Find the lon/lat information for this burst and calculate image centroid lat/lon/height.
    line_nums   = np.array(aux['sceneCenLine_number'],'int')
    pixel_nums  = np.array(aux['sceneCenPixel_number'],'int')
    distances   = (line_nums - cen_line)**2 + (cen_pixel - pixel_nums)**2
    id          = distances.argmin()
    orbit_info['centroid_lat'] = float(aux['sceneCenLat'][id])
    orbit_info['centroid_lon'] = float(aux['sceneCenLon'][id])
    orbit_info['height'] = float(aux['height'][id])

    # Finally calculate the new coordinates
    phi_lam_height = pixel_location(orbit_info,cen_line,cen_pixel,burst_num-1)
    lon = phi_lam_height[1]
    lat = phi_lam_height[0]
    coverage = list()

    # If we also want to calculate the coordinates of the corners in our script.
    if corners:
        # First update lon/lat/height because corners are calculated with reference to the centre of the burst.
        orbit_info['centroid_lat'] =  phi_lam_height[0]
        orbit_info['centroid_lon'] =  phi_lam_height[1]
        orbit_info['height'] = phi_lam_height[2]

        # Now calculate the coordinates of the corners.
        line        = ['', '']
        line[0]     = int(aux['lastValidSample'][burst_num-1].split(' ')[0:int(aux['imageLines'][0])/2].count('-1'))
        line[1]     = int(aux['imageLines'][0]) - aux['lastValidSample'][burst_num-1].split(' ')[:int(aux['imageLines'][0])/2:-1].count('-1')
        firstValidSample_num        =" ".join(aux['firstValidSample']).split(' ')
        firstValidSample_num_new    =list(set(firstValidSample_num))
        firstValidSample_num_new.sort()
        firstValidSample_num_new.remove('-1')
        pixel       = ['', '']
        pixel[0]    = int(firstValidSample_num_new[0])+1
        pixel[1]    = int(max("".join(aux['lastValidSample']).split(' ')))+1

        lin_pix = [[line[0],pixel[0]],[line[1],pixel[0]],[line[1],pixel[1]],[line[0],pixel[1]]]

        for xy in lin_pix:
            lat_lon = pixel_location(orbit_info, xy[0], xy[1], burst_num-1)
            coverage.append([lat_lon[1],lat_lon[0]])

        coverage = Polygon(coverage)

    return lat, lon, coverage


def burst_readfiles(meta,burst_num,swath_data,corners = False):
    # First copy swath metadata for burst and create a georef dict which stores information about the geo reference of
    # the burst.
    meta['Burst_number_index'] = str(burst_num)
    aux = meta['aux']
    aux['azimuthPRF'] = [meta['Pulse_Repetition_Frequency (computed, Hz)']]
    meta.pop('aux')

    # First find coordinates of center and optionally the corners
    lat, lon, coverage = burst_location(aux,burst_num=burst_num,corners=corners)
    meta['Scene_centre_longitude'] = lon
    meta['Scene_centre_latitude'] = lat

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

    return meta, coverage


def burst_datapoints(meta,burst_num):
    # First check which datapoints should be included.
    burst_start_time = np.datetime64(meta['aux']['azimuthTimeStart'][burst_num-1] + '-00')

    swath_start = burst_start_time - np.timedelta64(100,'s')
    swath_stop = burst_start_time + np.timedelta64(100,'s')

    datapoints = collections.OrderedDict()

    datapoints['row_1'] = ['t(s)','X(m)','Y(m)','Z(m)']
    t = meta['aux']['orbitTime']
    x = meta['aux']['orbitX']
    y = meta['aux']['orbitY']
    z = meta['aux']['orbitZ']
    datapoints['NUMBER_OF_DATAPOINTS'] = str(len(t))

    i = 0
    for n in range(len(t)):
        t_point = np.datetime64(t[n])

        if t_point > swath_start and t_point < swath_stop:
            t_s = datetime.strptime(t[n],'%Y-%m-%dT%H:%M:%S.%f')
            t_s = float(t_s.hour * 3600 + t_s.minute * 60 + t_s.second) + float(t_s.microsecond) / 1000000
            t_s = "%.6f" % t_s
            datapoints['row_' + str(i + 2)] = [t_s, x[n], y[n], z[n]]
            i += 1

    datapoints['NUMBER_OF_DATAPOINTS'] = str(i)

    return datapoints


def burst_precise(meta,burst_num,precise_folder,type='POE'):
    # This function utilizes the orbit_read script to read precise orbit files and export them to the resfile format.
    # Additionally it removes the burst_datapoints part, as it is not needed anymore.

    # First check whether the precise orbit file exists and load data if that is the case.

    t_s = datetime.strptime(meta['aux']['azimuthTimeStart'][burst_num-1], '%Y-%m-%dT%H:%M:%S.%f')
    date_orbit = int(t_s.hour * 3600 + t_s.minute * 60 + t_s.second)

    input_time = np.arange(date_orbit-100,date_orbit+100)
    date = time.mktime(time.strptime(meta['aux']['azimuthTimeStart'][burst_num-1], '%Y-%m-%dT%H:%M:%S.%f'))

    X,Y,Z=interpolate_orbit(precise_folder, input_time, date, type, 'spline', satellite =meta['Product type specifier'])

    if len(X) == 0 and type == 'POE':
        X, Y, Z = interpolate_orbit(precise_folder, input_time, date, 'RES', 'spline', satellite =meta['Product type specifier'])
        print 'There is no precise orbit file available, we try the restituted files'

    if len(X) == 0:
        print 'There is no precise or restituted orbit file available'
        return

    datapoints = collections.OrderedDict()
    datapoints['row_1'] = ['t(s)','X(m)','Y(m)','Z(m)']
    datapoints['NUMBER_OF_DATAPOINTS'] = str(len(input_time))

    for n in range(len(input_time)):
        datapoints['row_' + str(n + 2)] = [str(input_time[n]), str(X[n]), str(Y[n]), str(Z[n])]

    return datapoints


def burst_crop(meta,burst_num,swath_data,new_burst_num):
    # This function returns a description of the crop files part, which defines how the burst is cropped out of the
    # swath data. This function is generally called when the .res and raw data is written to the datastack and uses one
    # of the outputs from the burst_readfiles

    crop = collections.OrderedDict()

    last_sample =  [int(x) for x in meta['aux']['lastValidSample'][burst_num-1].split()]
    first_sample = [int(x) for x in meta['aux']['firstValidSample'][burst_num-1].split()]

    swath_data = os.path.basename(swath_data)
    crop['Data_output_file'] = swath_data[15:23] + '_iw_' + swath_data[6] + '_burst_' + str(new_burst_num) + '.raw'
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

