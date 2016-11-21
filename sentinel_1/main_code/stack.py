# This script gathers available sentinel files from the database and checks the coverage in space and time for this data
# stack.

import os, sys
import numpy as np
import warnings
from datetime import datetime, timedelta
from shapely.geometry import Polygon, shape, mapping, box
from shapely.ops import cascaded_union
import fiona
import image
import subprocess
from copy import deepcopy
from collections import Counter, OrderedDict
from sentinel_1.functions.sentinel_dump_data_function import dump_data

# TODO Make a function that describes an iteration over the burst in the datastack
# TODO Split this function in initialization (single products) and processing part (ifgs)

class StackData(object):
    # This function holds information for a full datastack of sentinel data and is used to select relevant images and
    # bursts based on a given area of interest.

    def __init__(self, track_dir, shape_dat, buffer=0.02, start_date='2014-04-01',end_date='', polarisation=['vh'], path='', db_type=1, precise_dir=''):
        # Initialize variables:

        # Datastack folder where all data from the datastack is stored (database,shapes,burst data, res files and results)
        # You should have write acces to this folder!
        self.path = []

        # Search path, shape, buffer and polarisation for this datastack. Currently only one polarisation is implemented. If
        # needed this can be extended later.
        self.search_path = []
        self.start_date = []
        self.end_date = []
        self.master_date = []
        self.shape = []
        self.shape_filename = []
        self.buffer = []
        self.polarisation = ''
        self.precise_orbits = ''

        # All images included in this datastack. Including lists of shapes and acquistion dates
        self.images = []
        self.image_files = []
        self.image_shapes = []
        self.image_dates = []

        # The shapes and names of the swaths
        self.swath_names = list()
        self.swath_shapes = list()

        # The resulting dates, with underlying bursts. This information can be used to create res files and images
        # of individual bursts. Variable data includes a structure using (dates > swaths > bursts).
        self.dates = list()
        self.datastack = dict()
        self.concatenated = dict()
        self.coordinates = dict()
        self.burst_names = list()
        self.burst_shapes = list()
        self.burst_centers = list()

        # Some important characteristics of the dataset of bursts. (matrices with names on y and dates on x axis)
        self.burst_availability = []
        self.burst_lon = []
        self.burst_lat = []
        self.burst_baselines = []

        # Temporary variable to store images which are not yet checked for shape or dates.
        self.image_dump = []

        ####################################################################

        # This function initializes the datastack using a search path, start/end dates and a buffer shape. The start and
        # end date should have the format yyyy-mm-dd and the shape should either be a shapefile or a list of coordinate
        # pairs. [[lat,lon][lat,lon] enz. ]. Minimum input is a folder which contains the .SAFE folders (a track folder)

        if not track_dir or not os.path.exists(track_dir):
            warnings.warn('This function needs an existing path as input!')
            return

        self.search_files(track_dir,db_type=db_type)

        if shape:
            self.create_shape(shape_dat,buffer)

        if end_date:
            self.end_date = np.datetime64(end_date).astype('datetime64[s]') + np.timedelta64(1, 'D').astype('timedelta64[s]')
        else:
            self.end_date = np.datetime64('now').astype('datetime64[s]')
        self.start_date = np.datetime64(start_date).astype('datetime64[s]')

        if isinstance(polarisation, basestring):
            polarisation = [polarisation]
            for i in polarisation:
                if not i in ['hh','vv','hv','vh']:
                    warnings.warn('This polarisation does not exist for sentinel data.')
                    return
        self.polarisation = polarisation
        self.search_path = track_dir

        if not path:
            warnings.warn('You did not specify an output path. Please do so later on using the add_path function')
        else:
            self.add_path(path)

        if precise_dir:
            if os.path.exists(precise_dir):
                self.precise_orbits = precise_dir
            else:
                print 'Precise orbit path does not exist'

    def read_existing(self,path):
        # This function reads an existing datastack. If needed also the interferograms which are already created.

        print('In development')

    def add_path(self,path):
        # This function adds the output path.

        if os.path.isdir(path):
            self.path = path
        elif os.path.isdir(os.path.dirname(path[:-1])):
            os.mkdir(path)
            self.path = path
        else:
            warnings.warn('Neither the directory itself nor the parents directory exists. Choose another path.')

    def search_files(self,track_dir,db_type=1):
        # This function searches for files within a certain folder. This folder should contain images from the same
        # track. These images are added to the variable image dump.
        images = list()

        if db_type == 1:
            image_dirs=next(os.walk(track_dir))[1]
            for i in image_dirs:
                if i.endswith('.SAFE'):
                    images.append(os.path.join(track_dir, i))

        if db_type == 2:
            for dates in next(os.walk(track_dir))[1]:
                print(dates)
                date = os.path.join(track_dir,dates)

                for files in next(os.walk(date))[2]:

                    data = os.path.join(date,files)
                    if data.endswith('.SAFE.zip'):
                        images.append(data)

        if images:
            self.image_dump = sorted(images)
        else:
            warnings.warn('No images found! Please choose another data folder.')

    def create_shape(self,shape_dat,buffer=0.02):
        # This function creates a shape to make a selection of usable bursts later on. Buffer around shape is in
        # degrees.

        if not shape_dat:
            warnings.warn('Please provide a shapefile or coordinates.')

        try:
            if isinstance(shape_dat,list):
                self.shape = Polygon(shape_dat)

            else: # It should be a shape file. We always select the first shape.
                sh = fiona.open(shape_dat).next()
                self.shape = shape(sh['geometry'])

            # Now we have the shape we add a buffer and simplify first to save computation time.
            self.shape = self.shape.simplify(buffer / 2)
            self.shape = self.shape.buffer(buffer)

        except:
            warnings.warn('Unrecognized shape')
            return

        self.buffer = buffer

    def unpack_image(self):
        # This program unpacks the images which are needed for processing.
        for image in self.images:
            image.unzip()
            image.init_unzipped()

    def del_unpacked_image(self):
        # This program unpacks the images which are needed for processing.
        for image in self.images:
            image.del_unzip()

    def select_image(self,start_date='',end_date=''):
        # This function selects usable images based on .kml files and dates

        if not self.shape:
            warnings.warn('There is no shape loaded to select images. Please use the create_shape function to do so.')

        if start_date:
            self.start_date = np.datetime64(start_date).astype('datetime64[s]')
        if end_date:
            self.end_date = np.datetime64(end_date).astype('datetime64[s]') + np.timedelta64(1, 'D').astype('timedelta64[s]')

        # First select images based on dates and check if polygons intersect.
        for i in self.image_dump:
            d = os.path.basename(i)[17:32]
            acq_time = np.datetime64(d[0:4] + '-' + d[4:6] + '-' + d[6:11] + ':' + d[11:13] + ':' + d[13:] + '-0000')
            if acq_time >= self.start_date and acq_time <= self.end_date:
                im = image.ImageMeta(path=i)
                im.read_kml()
                im_shape = deepcopy(im.coverage)

                if im_shape.intersects(self.shape):
                    self.images.append(im)
                    self.image_shapes.append(im_shape)
                    self.image_dates.append(acq_time)
                    self.image_files.append(os.path.basename(i))
        self.dates = sorted(list(set([d.astype('datetime64[D]') for d in self.image_dates])))

    def select_burst(self,date=''):
        # This function selects the usefull bursts at one epoch (user defined or automatically selected) and searches
        # usefull burst at other dates. This function uses the extend_burst function, which is intended to search for
        # bursts at other dates. This function can be run later on to update the datastack.

        image_dates = [im.astype('datetime64[D]') for im in self.image_dates]
        # First select which date will be the master
        if date:
            date = np.datenum64(date).astype('datetime64[D]')
            date = deepcopy(image_dates[min(abs(self.image.dates-date))])
        else:  # if no date is specified
            date = Counter(image_dates).most_common(1)[0][0]

        # Load the metadata for this date for the bursts and swaths
        image_id = np.where(image_dates == date)[0]
        for i in image_id:
            self.images[i].meta_swath()
            self.images[i].meta_burst()

        # Order the selected images by acquisition time.
        image_id = [x for (y,x) in sorted(zip([self.image_dates[i] for i in image_id],image_id))]
        date = date.astype(datetime).strftime('%Y-%m-%d')
        self.master_date = date
        self.datastack[date] = dict()

        for p in self.polarisation:
            for swath in ['1','2','3']:
                swath_id = 'iw' + swath + '-slc-' + p
                swath_name = 'swath_' + swath + '_' + p
                burst_no = 1

                for i in image_id:
                    # Check which swath should be selected.
                    swath_names = [os.path.basename(xml) for xml in self.images[i].swaths_xml]

                    swath_no = [no for no in range(len(swath_names)) if swath_id in swath_names[no]][0]
                    if not swath_no and swath_no is not 0: # If there is no data for this swath
                        break

                    for burst in self.images[i].swaths[swath_no].bursts:
                        if burst.burst_coverage.intersects(self.shape):

                            # If there are bursts in this swath, which intersect, add burst to list.
                            if swath_name not in self.datastack[date].keys(): # If there are burst and no burst list exists
                                self.datastack[date][swath_name] = dict()
                                if swath_name not in self.swath_names:
                                    self.swath_names.append(swath_name)

                            # Assign burst to data stack
                            burst_name = 'burst_' + str(burst_no)
                            burst.new_burst_num = burst_no
                            burst.res_burst(precise_folder=self.precise_orbits)
                            self.datastack[date][swath_name][burst_name] = burst

                            # Assign coverage, center coordinates and burst name.
                            self.burst_shapes.append(burst.burst_coverage)
                            self.burst_centers.append(burst.burst_center)
                            self.burst_names.append(swath_name + '_' + burst_name)
                            burst_no += 1

    def swath_coverage(self):
        # Create a convex hull for the different swaths.

        for swath in self.swath_names:
                # Assign coverage of swath using convex hull.
                burst_id = [i for i in range(len(self.burst_names)) if swath in self.burst_names[i]]
                swath_shape = cascaded_union([self.burst_shapes[i] for i in burst_id])
                self.swath_shapes.append(swath_shape)

    def extend_burst(self):
        # Searches for burst at dates other than the dates that are already available. This means that if there are
        # images added for dates which are already indexed, this data will not be used!

        # first read all metadata if not done already.
        processed_dates = [np.datetime64(d) for d in self.datastack.keys()]
        image_set = [d for d in self.dates if d not in processed_dates]
        image_dates = [d.astype('datetime64[D]') for d in self.image_dates]

        for date in image_set:
            # Append data to datastack variable
            date_str = date.astype(datetime).strftime('%Y-%m-%d')
            self.datastack[date_str] = dict()

            # Load the metadata for this date for the bursts and swaths
            image_id = np.where(image_dates == date)[0]
            for i in image_id:
                self.images[i].meta_swath()
                self.images[i].meta_burst(corners=False)
                # TODO Change precise orbits script to image wise approach (will be the same for all bursts in image)

            for swath in self.swath_names:
                # Add swath to datastack
                self.datastack[date_str][swath] = OrderedDict()

                for i in image_id:
                    # Select correct swath in image
                    swath_files = [os.path.basename(xml) for xml in self.images[i].swaths_xml]
                    identifier = swath[-4] + '-slc-' + swath[-2:]
                    swath_no = [no for no in range(len(swath_files)) if identifier in swath_files[no]][0]

                    no = 0
                    for burst in self.images[i].swaths[swath_no].bursts:
                        x_dist = np.array([xy[0] - burst.burst_center[0] for xy in self.burst_centers])
                        y_dist = np.array([xy[1] - burst.burst_center[1] for xy in self.burst_centers])

                        dist = np.sqrt(x_dist**2 + y_dist**2)
                        burst_id = np.argmin(dist)
                        no = no + 1
                        print str(dist) + ' ' + date_str + ' ' + swath + ' ' + str(no)

                        if dist[burst_id] < 0.1:
                            # Assign burst to data stack
                            burst.new_burst_num = int(self.burst_names[burst_id][17:])
                            burst.res_burst(precise_folder=self.precise_orbits)
                            self.datastack[date_str][swath][self.burst_names[burst_id][11:]] = burst

    def define_burst_coordinates(self,slaves=False):
        # This function defines the exact coordinates in pixels of every burst based on the lower left corner of the first
        # burst image. In this way the total overlap of these bursts can easily be monitored. Results are written to the
        # coordinates variable

        if slaves is True:
            dates = self.datastack.keys()
        else:
            dates = [self.master_date]

        self.coordinates = OrderedDict()

        for date in dates:
            self.coordinates[date] = OrderedDict()
            self.coordinates[date]['shapes'] = []
            ref = False

            min_line = 1; min_pixel = 1
            max_line = 1; max_pixel = 1
            for swath in self.datastack[date].keys():

                self.coordinates[date][swath] = OrderedDict()
                self.coordinates[date][swath]['corners'] = np.zeros([len(self.datastack[date][swath].keys()), 4, 2],dtype='int')

                b = 0
                for burst in sorted(self.datastack[date][swath].keys(), key = lambda x: int(x[6:])):

                    if ref is False:
                        self.coordinates[date]['ref_az_time'] = self.datastack[date][swath][burst].processes['readfiles']['First_pixel_azimuth_time (UTC)']
                        self.coordinates[date]['ref_range_time'] = self.datastack[date][swath][burst].processes['readfiles']['Range_time_to_first_pixel (2way) (ms)']
                        ref = True

                    az_first = self.datastack[date][swath][burst].processes['readfiles']['First_pixel_azimuth_time (UTC)']
                    az_samp = self.datastack[date][swath][burst].processes['readfiles']['Pulse_Repetition_Frequency (computed, Hz)']
                    first_line = int(self.datastack[date][swath][burst].processes['crop']['First_line (w.r.t. original_image)'])
                    last_line = int(self.datastack[date][swath][burst].processes['crop']['Last_line (w.r.t. original_image)'])

                    range_first = self.datastack[date][swath][burst].processes['readfiles']['Range_time_to_first_pixel (2way) (ms)']
                    range_samp = self.datastack[date][swath][burst].processes['readfiles']['Range_sampling_rate (computed, MHz)']
                    first_pixel = int(self.datastack[date][swath][burst].processes['crop']['First_pixel (w.r.t. original_image)'])
                    last_pixel = int(self.datastack[date][swath][burst].processes['crop']['Last_pixel (w.r.t. original_image)'])

                    no_lines = int(self.datastack[date][swath][burst].processes['readfiles']['Number_of_lines_original'])
                    no_pixels = int(self.datastack[date][swath][burst].processes['readfiles']['Number_of_pixels_original'])

                    # Calculate difference w.r.t. reference point.
                    range_time_diff = (float(range_first) - float(self.coordinates[date]['ref_range_time']))
                    pixel_diff = int(round((float(range_samp) * 1e3) * range_time_diff))

                    az_time1 = datetime.strptime(az_first, '%Y-%b-%d %H:%M:%S.%f')
                    az_time2 = datetime.strptime(self.coordinates[date]['ref_az_time'], '%Y-%b-%d %H:%M:%S.%f')
                    az_time_diff = (az_time1 - az_time2).total_seconds()
                    line_diff = int(round(float(az_samp) * az_time_diff))

                    # Calculate final corner coordinates.
                    ll = np.array([line_diff + first_line, pixel_diff + first_pixel], ndmin=2)
                    ul = np.array([line_diff + last_line, pixel_diff + first_pixel], ndmin=2)
                    ur = np.array([line_diff + last_line, pixel_diff + last_pixel], ndmin=2)
                    lr = np.array([line_diff + first_line, pixel_diff + last_pixel], ndmin=2)

                    self.coordinates[date][swath]['corners'][b,:,:] = np.vstack([ll,ul,ur,lr])
                    b += 1

                    # Check for max/min line/pixel to prevent negative pixel numbers.
                    min_line = min(1 + line_diff, min_line)
                    max_line = max(line_diff + no_lines, max_line)
                    min_pixel = min(1 + pixel_diff, min_pixel)
                    max_pixel = max(pixel_diff + no_pixels, max_pixel)

            if min_line < 1:
                max_line = max_line + (1 - min_line)
            if min_pixel < 1:
                max_pixel = max_pixel + (1 - min_pixel)

            for swath in self.datastack[date].keys():
                # If one of the lines or pixels is lower then 1, correct all coordinates to have positive coordinates for all
                # bursts.
                if min_line < 1:
                    self.coordinates[date][swath]['corners'][:,:,0] = self.coordinates[date][swath]['corners'][:,:,0] + (1 - min_line)
                if min_pixel < 1:
                    self.coordinates[date][swath]['corners'][:,:,1] = self.coordinates[date][swath]['corners'][:,:,1] + (1 - min_pixel)

                for burst in range(len(self.datastack[date][swath].keys())):
                    shape_c = self.coordinates[date][swath]['corners'][burst,:,:]
                    shape = box(shape_c[0,1],shape_c[0,0],shape_c[2,1],shape_c[2,0])
                    self.coordinates[date]['shapes'].append(shape)

                    # Finally add information to the .res file if already loaded
                    if self.datastack[date][swath]['burst_' + str(burst+1)].processes['readfiles']:

                        read = self.datastack[date][swath]['burst_' + str(burst+1)].processes['readfiles']
                        read['First_line (w.r.t. output_image)'] = str(shape_c[0,0])
                        read['Last_line (w.r.t. output_image)'] = str(shape_c[2,0])
                        read['First_pixel (w.r.t. output_image)'] = str(shape_c[0,1])
                        read['Last_pixel (w.r.t. output_image)'] = str(shape_c[2,1])
                        read['Number_of_pixels_output_image'] = str(max_pixel)
                        read['Number_of_lines_output_image'] = str(max_line)
                        self.datastack[date][swath]['burst_' + str(burst+1)].processes['readfiles'] = read
                    else:
                        print 'No resfile available, so information is not added to resfile'

    def correct_jitter(self):
        # This function corrects for jitter errors in the data. This will be based on first line azimuth time and PRF.
        # Likely this error is resolved after June 2015, so for later products it is not necessary.

        dates = self.datastack.keys()

        for date in dates:
            ref = False

            for swath in self.datastack[date].keys():
                for burst in sorted(self.datastack[date][swath].keys(), key = lambda x: int(x[6:])):
                    # First select azimuth time and reference time

                    if ref is False:
                        ref_az_time = self.datastack[date][swath][burst].processes['readfiles']['First_pixel_azimuth_time (UTC)']
                        ref_az_time = datetime.strptime(ref_az_time, '%Y-%b-%d %H:%M:%S.%f')
                        prf = self.datastack[date][swath][burst].processes['readfiles']['Pulse_Repetition_Frequency (computed, Hz)']
                        ref = True

                    az_time = self.datastack[date][swath][burst].processes['readfiles']['First_pixel_azimuth_time (UTC)']
                    az_time = datetime.strptime(az_time, '%Y-%b-%d %H:%M:%S.%f')

                    # Then round the time values to values that correspond with full pixels.
                    time_diff = (az_time - ref_az_time).total_seconds()
                    jitter_diff = (round(float(prf) * time_diff) / float(prf)) - time_diff
                    az_time = (az_time + timedelta(seconds=jitter_diff)).strftime('%Y-%b-%d %H:%M:%S.%f')
                    self.datastack[date][swath][burst].processes['readfiles']['First_pixel_azimuth_time (UTC)'] = az_time

    def lat_lon_availability(self):
        # This function searches for the lat/lon and availability of all bursts in the dataset. Every time the dataset
        # is extended this function should be updated.

        date_no = len(self.dates)
        burst_no = len(self.burst_names)

        self.burst_availability = np.zeros((date_no,burst_no),dtype=bool)
        self.burst_lon = np.zeros((date_no,burst_no))
        self.burst_lat = np.zeros((date_no,burst_no))

        d = -1
        for date in self.dates:
            date = date.astype(datetime).strftime('%Y-%m-%d')
            d += 1
            b = -1
            for name in self.burst_names:
                swath = name[:10]
                burst = name[11:]
                b += 1
                if burst in self.datastack[date][swath].keys():
                    self.burst_availability[d,b] = True
                    self.burst_lon[d,b] = self.datastack[date][swath][burst].burst_center[0]
                    self.burst_lat[d,b] = self.datastack[date][swath][burst].burst_center[1]

    def baseline(self):
        # This function gives the baseline based on doris function. The selected burst will be given as reference.

        if not self.burst_lat:
            self.lat_lon_availability()

        date_no = len(self.dates)
        burst_no = len(self.burst_names)

        print('In progress')

    def make_network(self):
        # This function will create a network of interferograms

        print('In progress')

    def write_stack(self,write_path='',no_data=False):
        # This function writes the full datastack to a given folder using the dates / swaths / bursts setup. This
        # also generates the res readfiles data.
        if write_path and os.path.exists(write_path):
            self.path = write_path
        if (not write_path or not os.path.exists(write_path)) and not self.path:
            warnings.warn('Please specify a path that exists to write the data')
            return

        for date in self.datastack.keys():

            date_basic = date.translate(None,'-')
            date_path = os.path.join(self.path,date_basic)
            if not os.path.exists(date_path):
                os.mkdir(date_path)

            for swath in self.datastack[date].keys():

                swath_path = os.path.join(date_path,swath)
                if not os.path.exists(swath_path):
                    os.mkdir(swath_path)

                for burst in self.datastack[date][swath].keys():
                    # Finally write the bursts with their res files and precise orbits
                    xml = self.datastack[date][swath][burst].swath_xml
                    data = self.datastack[date][swath][burst].swath_data
                    image_no = str(self.datastack[date][swath][burst].burst_num)
                    stack_no = burst[6:]
                    xml_base = os.path.basename(xml)
                    res_name = os.path.join(swath_path,xml_base[15:23] + '_iw_' + xml_base[6] + '_burst_' + stack_no + '.res')
                    swath_res = os.path.join(swath_path,xml_base[15:23] + '_iw_' + xml_base[6] + '.res')
                    outdata =  os.path.join(swath_path,xml_base[15:23] + '_iw_' + xml_base[6] + '_burst_' + stack_no + '.raw')

                    # First create res file

                    # Old method
                    #create_res = ('python sentinel_dump_header2doris_single.py ' + xml + ' ' + stack_no + ' ' + image_no
                    #              + ' ' + swath_path)
                    #os.system(create_res)

                    #if not no_data:
                    #    # Then add precise orbits
                    #    if precise_orbits:
                    #        try:
                    #            add_precise = 'python orbit_read.py ' + res_name + ' ' + precise_orbits + ' cubic POE'
                    #            os.system(add_precise)
                    #        except:
                    #            warnings.warn('File will be written without precise orbits')

                    #    # Finally write data
                    #    crop = 'python get_parameter_coord_single.py ' + swath_res + ' ' + stack_no
                    #    crop = subprocess.Popen(crop , shell=True, stdout=subprocess.PIPE,).communicate()[0][:-1]
                    #    write_data = 'python sentinel_dump_data.py ' + data + ' ' + outdata + ' ' + crop + ' -res ' + res_name
                    #    os.system(write_data)

                    # Write res file
                    self.datastack[date][swath][burst].write(res_name)

                    # Save data
                    # crop = 'python get_parameter_coord_single.py ' + swath_res + ' ' + stack_no
                    # crop = subprocess.Popen(crop , shell=True, stdout=subprocess.PIPE,).communicate()[0][:-1]
                    # write_data = 'python sentinel_dump_data.py ' + data + ' ' + outdata + ' ' + crop + ' -res ' + res_name
                    # os.system(write_data)

                    dump_data(data,res_name,output_file=outdata)

    def create_concatenate_stack(self,slaves=True):
        # This function creates a folder and res files for the concatenated files

        if self.coordinates and self.datastack:
           if slaves is True and not len(self.coordinates.keys()) is len(self.datastack.keys()):
                print 'If you want to create the base for the concatenated data, please run the define_burst_coordinates with slaves=True!'
                return
           # Otherwise it is fine...
        else:
            print 'Please run functions to create datastack and coordinates first!'

    def setup_single_master(self,start_date='',end_date='',master=''):
        # This function makes a setup for the single master of sentinel data.

        print('in progress')

    def do_coreg(self,dem=True):
        # This function does the coregistration for this datastack
        print('in progress')

    def create_coreg_concatenate(self):
        # This function reads
        print('in progress')

    def write_shapes(self,coverage=True,images=True,swaths=True,bursts=True,coordinates=True):
        # This function writes shapefiles of the area of interest, images and bursts.

        if not coverage and not images and not bursts:
            warnings.warn('Select at least one of shape types')
            return

        shapes = list(); shape_names = list(); shape_files = list()
        if coverage and self.shape:
            shape_files.append(self.path + 'area_of_interest.shp')
            if self.shape.type == 'MultiPolygon':
                shapes.append([sh for sh in self.shape])
                shape_names.append([('coverage_with_buffer_of_' + str(self.buffer) + '_degrees_' + str(i)) for i in range(len(shapes[0]))])
            elif self.shape.type == 'Polygon':
                shape_names.append(['coverage_with_buffer_of_' + str(self.buffer) + '_degrees'])
                shapes.append([self.shape])
        if images and self.image_shapes:
            shape_files.append(self.path + 'image_coverage.shp')
            shape_names.append([date.astype(datetime).strftime('%Y-%m-%dT%H:%M:%S') for date in self.image_dates])
            shapes.append(self.image_shapes)
        if swaths and self.swath_shapes:
            shape_files.append(self.path + 'swath_coverage.shp')
            shape_names.append(self.swath_names)
            shapes.append(self.swath_shapes)
        if bursts and self.burst_shapes:
            shape_files.append(self.path + 'burst_coverage.shp')
            shape_names.append(self.burst_names)
            shapes.append(self.burst_shapes)
        if coordinates and self.coordinates:
            shape_files.append(self.path + 'burst_pixels.shp')
            shape_names.append(self.burst_names)
            shapes.append(self.coordinates[self.master_date]['shapes'])

        shape_setup = {
            'geometry': 'Polygon',
            'properties': {'name': 'str'},
        }

        for i in range(len(shape_files)):
            with fiona.open(shape_files[i], 'w', 'ESRI Shapefile', shape_setup) as sh:
                for n in range(len(shapes[i])):
                    sh.write({
                        'geometry': mapping(shapes[i][n]),
                        'properties': {'name': shape_names[i][n]},
                    })
