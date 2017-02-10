from datetime import datetime
import collections
from sentinel_1.functions.precise_read import interpolate_orbit
from shapely.geometry import Polygon
import numpy as np
from scipy.interpolate import RectBivariateSpline
import time


def swath_datapoints(meta):
    # First check which datapoints should be included.
    swath_start_time = np.datetime64(meta['aux']['azimuthTimeStart'][0])

    swath_start = swath_start_time - np.timedelta64(100, 's')
    swath_stop = swath_start_time + np.timedelta64(100, 's')

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


def swath_precise(meta, precise_folder, dat_type='POE'):
    # This function utilizes the orbit_read script to read precise orbit files and export them to the resfile format.
    # Additionally it removes the burst_datapoints part, as it is not needed anymore.

    # First check whether the precise orbit file exists and load data if that is the case.

    t_s = datetime.strptime(meta['aux']['azimuthTimeStart'][0], '%Y-%m-%dT%H:%M:%S.%f')
    date_orbit = int(t_s.hour * 3600 + t_s.minute * 60 + t_s.second)

    input_time = np.arange(date_orbit-100,date_orbit+100)
    date = time.mktime(time.strptime(meta['aux']['azimuthTimeStart'][0], '%Y-%m-%dT%H:%M:%S.%f'))

    if dat_type == 'POE' or dat_type == 'RES':
        X, Y, Z = interpolate_orbit(precise_folder, input_time, date, dat_type, 'spline', satellite =meta['Product type specifier'])

        if len(X) == 0 and dat_type == 'POE':
            dat_type = 'RES'
            X, Y, Z = interpolate_orbit(precise_folder, input_time, date, 'RES', 'spline', satellite =meta['Product type specifier'])
            print('There is no precise orbit file available, we try the restituted files')

    if len(X) == 0 or dat_type == 'XML':
        print('There is no precise or restituted orbit file available we use the datapoints from the .xml file')
        datapoints = swath_datapoints(meta)
        return datapoints

    datapoints = collections.OrderedDict()
    datapoints['row_1'] = ['t(s)','X(m)','Y(m)','Z(m)']
    datapoints['NUMBER_OF_DATAPOINTS'] = str(len(input_time))

    for n in range(len(input_time)):
        datapoints['row_' + str(n + 2)] = [str(input_time[n]), str(X[n]), str(Y[n]), str(Z[n])]

    return datapoints


def swath_pixel_line_coordinate(meta):
    # This function converts the given datapoints from xml to line/pixel coordinates with respect to the upper left 
    # pixel of the swath. This information can be used to make a first guess of the corners and center locations of 
    # different burst.

    aux = meta['aux']
    lats = np.array([float(l) for l in aux['sceneCenLat']])
    lons = np.array([float(l) for l in aux['sceneCenLon']])
    lines = [int(l) for l in aux['sceneCenLine_number']]
    pixels = [int(l) for l in aux['sceneCenPixel_number']]

    az_time = [np.datetime64(t) - np.datetime64(aux['azimuthTime'][0]) for t in aux['azimuthTime']]
    az_step = np.timedelta64(int(float(meta['Azimuth_time_interval (s)']) * 1000000000000), 'ps')
    new_lines = [int(round(t/az_step)) for t in az_time]

    lines, count_l = np.unique(new_lines, return_counts=True)
    pixels, count_p = np.unique(pixels, return_counts=True)
    lats = lats.reshape((count_p[0], count_l[0]))
    lons = lons.reshape((count_p[0], count_l[0]))

    lat_interp = RectBivariateSpline(lines, pixels, lats)
    lon_interp = RectBivariateSpline(lines, pixels, lons)
    # These functions can be called using: lat_interp.ev(li, pi). Where li and pi are lists of pixels and lines.

    return lat_interp, lon_interp


def burst_coverage(meta, corners=True, shape=True):
    # This function returns the lat, lon of the corners of all bursts in this swath. If polygon is True also the poly
    # gons are generated.

    # First get the interpolation from pix/line to lat/lon
    lat_interp, lon_interp = swath_pixel_line_coordinate(meta)

    # Now calculate the centre pixels of individual bursts.
    l_b = int(meta['aux']['imageLines'][0])
    p_b = int(meta['aux']['imagePixels'][0])

    # Calculate first lines
    start_times = meta['aux']['azimuthTimeStart']
    az_time = [np.datetime64(t) - np.datetime64(start_times[0]) for t in start_times]
    az_step = np.timedelta64(int(float(meta['Azimuth_time_interval (s)']) * 1000000000000), 'ps')
    start_lines = [int(round(t/az_step)) for t in az_time]

    burst_center = []
    burst_corners = []
    burst_shapes = []

    # Now loop over all the bursts and calc the center pixel / corners / polygon
    for l in start_lines:
        center = [(lat_interp(l+np.floor(l_b/2), np.floor(p_b/2))[0][0], lon_interp(l+np.floor(l_b/2), np.floor(p_b/2))[0][0])]
        burst_center.append(center)
        if corners == True or shape == True:
            ul = (lat_interp(l, 0)[0][0], lon_interp(l , 0)[0][0])
            ur = (lat_interp(l, p_b-1)[0][0], lon_interp(l , p_b-1)[0][0])
            lr = (lat_interp(l+l_b-1, p_b-1)[0][0], lon_interp(l+l_b-1, p_b-1)[0][0])
            ll = (lat_interp(l+l_b-1, 0)[0][0], lon_interp(l+l_b-1, 0)[0][0])
            burst_corners.append([ul, ur, lr, ll])

            if shape == True:
                burst_shapes.append(Polygon([ul, ur, lr, ll]))

    return burst_center, burst_corners, burst_shapes


def swath_coverage(meta):
    # This function calculates the total coverage of the swath.

    aux = meta['aux']
    lats = np.array([float(l) for l in aux['sceneCenLat']])
    lons = np.array([float(l) for l in aux['sceneCenLon']])
    lines = [int(l) for l in aux['sceneCenLine_number']]
    pixels = [int(l) for l in aux['sceneCenPixel_number']]

    lines, count_l = np.unique(lines, return_counts=True)
    pixels, count_p = np.unique(pixels, return_counts=True)
    lats = lats.reshape((count_p[0], count_l[0]))
    lons = lons.reshape((count_p[0], count_l[0]))

    # ul, ur, lr, ll
    swath_corners = [(lats[0,0],lons[0,0]),(lats[0,-1],lons[0,-1]),(lats[-1,-1],lons[-1,-1]),(lats[-1,0],lons[-1,0])]
    swath_shapes = Polygon(swath_corners)

    return swath_corners, swath_shapes
