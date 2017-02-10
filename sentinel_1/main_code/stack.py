# This script gathers available sentinel files from the database and checks the coverage in space and time for this data
# stack.

import os
import warnings
from collections import Counter, OrderedDict
from datetime import datetime

import fiona
import numpy as np
from shapely.geometry import shape, mapping, box

import sentinel_1.main_code.image as image
from jobs import Jobs
from sentinel_1.functions.load_shape_unzip import extract_kml_preview, shape_im_kml
from sentinel_1.functions.load_shape_unzip import load_shape
from sentinel_1.main_code.dorisparameters import DorisParameters
from sentinel_1.functions.burst_metadata import center_shape_from_res


class StackData(object):
    # This function holds information for a full datastack of sentinel data and is used to select relevant images and
    # bursts based on a given area of interest.

    def __init__(self, track_dir, shape_dat, buffer=0.02, start_date='2014-04-01', end_date='', polarisation='vh', path='', db_type=1, precise_dir=''):
        # Initialize variables:

        # Datastack folder where all data from the datastack is stored (database,shapes,burst data, res files and results)
        # You should have write acces to this folder!
        self.path = path

        # Search path, shape, buffer and polarisation for this datastack. Currently only one polarisation is implemented. If
        # needed this can be extended later.
        self.search_path = []
        self.start_date = []
        self.end_date = []
        self.master_date = []
        self.shape = []
        self.shape_filename = shape_dat
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
        self.burst_no = 0
        self.burst_availability = []
        self.burst_lon = []
        self.burst_lat = []
        self.burst_baselines = []

        # Temporary variable to store images which are not yet checked for shape or dates.
        self.image_dump = []

        # parallel computing:
        doris_parameters = DorisParameters(os.path.dirname(os.path.dirname(self.path)))
        self.doris_parameters = doris_parameters
        self.nr_of_jobs = doris_parameters.nr_of_jobs
        self.parallel = doris_parameters.parallel
        self.function_path = doris_parameters.function_path

        ####################################################################

        # This function initializes the datastack using a search path, start/end dates and a buffer shape. The start and
        # end date should have the format yyyy-mm-dd and the shape should either be a shapefile or a list of coordinate
        # pairs. [[lat,lon][lat,lon] enz. ]. Minimum input is a folder which contains the .SAFE folders (a track folder)

        if not track_dir or not os.path.exists(track_dir):
            warnings.warn('This function needs an existing path as input!')
            return

        self.search_files(track_dir)

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

    def add_path(self,path):
        # This function adds the output path.

        if os.path.isdir(path):
            self.path = path
        elif os.path.isdir(os.path.dirname(path[:-1])):
            os.mkdir(path)
            self.path = path
        else:
            warnings.warn('Neither the directory itself nor the parents directory exists. Choose another path.')

    def search_files(self,track_dir):
        # This function searches for files within a certain folder. This folder should contain images from the same
        # track. These images are added to the variable image dump.
        images = list()

        top_dir = next(os.walk(track_dir))

        for data in top_dir[2]:
            if data.endswith('.zip'):
                images.append(os.path.join(track_dir, data))
        for data in top_dir[1]:
            if data.endswith('.SAFE'):
                images.append(os.path.join(track_dir, data))
            else:  # Likely images are stored in a folder.
                sec_dir = next(os.walk(os.path.join(track_dir, data)))

                for dat in sec_dir[1]:
                    if dat.endswith('.SAFE'):
                        images.append(os.path.join(track_dir, data, dat))
                for dat in sec_dir[2]:
                    if dat.endswith('.zip'):
                        images.append(os.path.join(track_dir, data, dat))

        if images:
            images = sorted(images)
        else:
            warnings.warn('No images found! Please choose another data folder. Track_dir = ' + str(track_dir))

        base = []  # Remove double hits because data is already unzipped.
        for i in images: # Drop all .zip files which are unpacked already.
            if i.endswith('.SAFE.zip'):
                base.append(os.path.basename(i[:-9]))
            elif i.endswith('.zip'):
                base.append(os.path.basename(i[:-4]))
            elif i.endswith('.SAFE'):
                base.append(os.path.basename(i[:-5]))
        b, id = np.unique(base, return_index=True)

        rem = []
        for i in range(len(base)):
            if i in id:
                self.image_dump.append(images[i])
            else:
                rem.append(images[i])

        if rem:
            print('removed the following zip files from stack:')
            for r in rem:
                print(r)
            print('It is advised to work with zipfiles instead of unpacked data. This saves diskspace and will ')

    def create_shape(self,shape_dat,buffer=0.02):
        # This function creates a shape to make a selection of usable bursts later on. Buffer around shape is in
        # degrees.

        self.shape = load_shape(shape_dat, buffer)
        self.buffer = buffer

    def check_new_images(self, master):
        # This function checks which images are already processed, and which are not. If certain dates are already
        # processed they are removed from the list. You have to specify the master date, otherwise the script will not
        # know how many burst are expected per date.

        # Which dates are available?
        image_dates = [im.astype('datetime64[D]') for im in self.image_dates]
        # What is the master date
        date = np.datetime64(master).astype('datetime64[D]')

        date_folders = [d for d in next(os.walk(self.path))[1] if len(d) == 8]
        rm_id = []

        if date_folders:
            dates = [np.datetime64(d[0:4] + '-' + d[4:6] + '-' + d[6:8]) for d in date_folders]

            if date in dates:
                date_folder = date_folders[np.where(dates == date)[0][0]]

                # Check existing files in master folder
                swaths = dict()
                swath_folders = next(os.walk(os.path.join(self.path, date_folder)))[1]
                if len(swath_folders) == 0:
                    print('No swaths in master folder')
                    #return

                for swath in swath_folders:
                    swaths[swath] = sorted([b[14:] for b in next(os.walk(os.path.join(self.path, date_folder, swath)))[2]])

                # Now check if the burst also in slave folders exist....
                for folder, d in zip(date_folders, dates):
                    # Check existing files in master folder
                    try:
                        swath_folders = next(os.walk(os.path.join(self.path, folder)))[1]
                        if not set(swath_folders) == set(swaths.keys()):
                            raise LookupError('Amount of swaths is not the same for ' + folder)

                        for swath in swath_folders:
                            bursts = sorted([b[14:] for b in next(os.walk(os.path.join(self.path, folder, swath)))[2]])
                            if not set(bursts) == set(swaths[swath]):
                                raise LookupError('Amount of bursts is not the same for ' + folder)

                            if d == date:
                                # If the master is already processed we have to create the list of center and coverage of
                                # bursts.
                                # TODO make this robust for the case no seperate input data folders are created.
                                res_files = sorted([b for b in next(os.walk(os.path.join(self.path, folder, swath)))[2] if b.endswith('res')])
                                for res_file in res_files:
                                    res_file = os.path.join(self.path, folder, swath, res_file)
                                    center, coverage = center_shape_from_res(resfile=res_file)
                                    # Assign coverage, center coordinates and burst name.
                                    self.burst_shapes.append(coverage)
                                    self.burst_centers.append(center)
                                    self.burst_names.append(res_file[14:])
                                    self.burst_no += 1

                        # If all bursts are the same these files are not processed.
                        for id in np.where(image_dates == d)[0][::-1]:
                            del self.image_dates[id]
                            del image_dates[id]
                            del self.images[id]
                            del self.image_files[id]

                    except LookupError as error:
                        print(error)

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
                kml, png = extract_kml_preview(i, png=False)
                succes = shape_im_kml(self.shape, kml)

                if succes:
                    self.images.append(im)
                    self.image_dates.append(acq_time)
                    self.image_files.append(os.path.basename(i))
        self.dates = sorted(list(set([d.astype('datetime64[D]') for d in self.image_dates])))

    def select_burst(self,date=''):
        # This function selects the usefull bursts at one epoch (user defined or automatically selected) and searches
        # usefull burst at other dates. This function uses the extend_burst function, which is intended to search for
        # bursts at other dates. This function can be run later on to update the datastack.

        if self.burst_names:
            print('Master data is already loaded from earlier processed data.')
            return

        image_dates = [im.astype('datetime64[D]') for im in self.image_dates]
        # First select which date will be the master
        if date:
            date = np.datetime64(date).astype('datetime64[D]')
            # date = deepcopy(image_dates[min(abs(self.image_dates-date))])
        else:  # if no date is specified
            date = Counter(image_dates).most_common(1)[0][0]

        # Load the metadata for this date for the bursts and swaths
        image_id = np.where(image_dates == date)[0]
        for i in image_id:
            self.images[i].meta_swath()

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
                data = []

                for i in image_id:
                    # Check which swath should be selected.
                    swath_names = [os.path.basename(xml) for xml in self.images[i].swaths_xml]

                    swath_no = [no for no in range(len(swath_names)) if swath_id in swath_names[no]][0]
                    if not swath_no and swath_no is not 0:  # If there is no data for this swath
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
                            self.datastack[date][swath_name][burst_name] = burst

                            # Assign coverage, center coordinates and burst name.
                            self.burst_shapes.append(burst.burst_coverage)
                            self.burst_centers.append(burst.burst_center)
                            self.burst_names.append(swath_name + '_' + burst_name)
                            burst_no += 1

                            # Finally add also to the number of bursts from the image
                            self.images[i].burst_no += 1
                            self.burst_no += 1

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
            data = []

            # Load the metadata for this date for the bursts and swaths
            image_id = np.where(image_dates == date)[0]
            for i in image_id:
                self.images[i].meta_swath()

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

                        if dist[burst_id] < 0.1:
                            # Assign burst to data stack
                            burst.new_burst_num = int(self.burst_names[burst_id][17:])
                            if data:
                                data = burst.res_burst(precise_folder=self.precise_orbits, data=data)
                            else:
                                data = burst.res_burst(precise_folder=self.precise_orbits)
                            self.datastack[date_str][swath][self.burst_names[burst_id][11:]] = burst
                        else:
                            print('No corresponding burst found! Closest is ' + str(dist[burst_id]) + ' ' +
                                  date_str + ' ' + swath + ' ' + self.burst_names[burst_id])

    def remove_incomplete_images(self):
        # This function removes all the images with less than maximum bursts. This will make a stack more consistent.

        for key in self.datastack.keys():
            burst_no = 0
            for key_swath in self.datastack[key].keys():
                burst_no += len(self.datastack[key][key_swath])
            if burst_no != self.burst_no:
                self.datastack.pop(key)
                print('Number of burst for ' + key + ' is ' + str(burst_no) + ' instead of ' + str(self.burst_no) +
                      ' and is removed from the datastack.')

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

    def write_stack(self,write_path='',no_data=False):
        # This function writes the full datastack to a given folder using the dates / swaths / bursts setup. This
        # also generates the res readfiles data.
        if write_path and os.path.exists(write_path):
            self.path = write_path
        if (not write_path or not os.path.exists(write_path)) and not self.path:
            warnings.warn('Please specify a path that exists to write the data')
            return

        jobs = []
        burst_num = []

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
                    outdata =  os.path.join(swath_path,xml_base[15:23] + '_iw_' + xml_base[6] + '_burst_' + stack_no + '.raw')

                    self.datastack[date][swath][burst].write(res_name)
                    if not os.path.exists(res_name) or not os.path.exists(outdata):

                        jobs.append('python ' + self.function_path + 'sentinel_dump_data_function.py ' + data + ' ' + res_name)
                        burst_num.append(stack_no + '_' + xml_base[6] + '_' + xml_base[15:23])

        # Burst are sorted in such a way that mainly read from different data files sorted by burst then swath then date.
        ids = sorted(range(len(burst_num)), key=lambda x: burst_num[x])
        jobs = jobs(ids)

        jobList1 = []
        for job in jobs:
            jobList1.append([self.path, job])
            if not self.parallel:
                os.chdir(self.path)
                # Resample
                os.system(job)
        if self.parallel:
            jobs = Jobs(self.nr_of_jobs, self.doris_parameters)
            jobs.run(jobList1)

    def unpack_image(self, dest_folder=''):

        if not dest_folder:
            dest_folder = os.path.join(self.path, 'slc_data_files')
        if not os.path.exists(dest_folder):
            os.mkdir(dest_folder)

        jobList1 = []

        # This program unpacks the images which are needed for processing. If unpacking fails, they are removed..
        for imagefile in self.images:

            zipped_folder = imagefile.zip_path
            if zipped_folder.endswith('.SAFE.zip'):
                imagefile.unzip_path = os.path.join(dest_folder, os.path.basename(zipped_folder[:-9] + '.SAFE'))
            elif zipped_folder.endswith('.zip'):
                imagefile.unzip_path = os.path.join(dest_folder, os.path.basename(zipped_folder[:-4] + '.SAFE'))
            shapefile = self.shape_filename
            pol = self.polarisation[0]
            overwrite = False
            command1 = ('python ' + self.function_path + 'load_shape_unzip.py ' + zipped_folder + ' ' + dest_folder +
                        ' ' + shapefile + ' ' + pol + ' ' + str(overwrite))
            jobList1.append([self.path, command1])
            if not self.parallel:
                os.chdir(self.path)
                # Resample
                os.system(command1)

        if self.parallel:
            jobs = Jobs(self.nr_of_jobs, self.doris_parameters)
            jobs.run(jobList1)

    def del_unpacked_image(self):
        # This program unpacks the images which are needed for processing.
        for image in self.images:
            os.remove(image.unzip_path)

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
