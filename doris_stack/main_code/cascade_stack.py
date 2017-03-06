import os
import numpy as np
from datetime import datetime
from collections import OrderedDict
from baselines import baselines
import copy
from copy import deepcopy
import shutil
from resdata import ResData
from ESD_functions import get_f_DC_difference, get_offset
from scipy import linalg
import resdata as resdata
from dorisparameters import DorisParameters
import collections

from jobs import Jobs

from single_master_stack import SingleMaster

class Cascade(SingleMaster):

    def __init__(self,datastack=[],stack_read=False,start_date='',end_date='',stack_folder='',processing_folder='',
                 input_files='',doris_path='', cpxfiddle_folder=''):
        # This function loads in a datastack to create a single master stack. Optional are the start date, end date and
        # master date. If they are not defined all dates will be loaded. The master date can be loaded later using the
        # master function. If you want a random master value, choose master_date='random'
        # Dates should be given as 'yyyy-mm-dd'. If stack_read is True information is read from the stack_folder. Otherwise
        # the datastack from an StackData object should be given as input.

        doris_parameters = DorisParameters(os.path.dirname(os.path.dirname(processing_folder))) # (assuming it is stored in the stackfolder
        self.doris_parameters = doris_parameters
        
        if not start_date:
            self.start_date = doris_parameters.start_date_default
        else:
            self.start_date = start_date
        if not end_date:
            self.end_date = doris_parameters.end_date_default
        else:
            self.end_date = end_date

        self.start_date = datetime.strptime(self.start_date,'%Y-%m-%d')
        self.end_date = datetime.strptime(self.end_date, '%Y-%m-%d')

        self.nr_of_jobs = doris_parameters.nr_of_jobs
        self.parallel = doris_parameters.parallel

        self.stack = dict()
        self.full_swath = dict()
        self.folder = processing_folder
        self.stack_folder = stack_folder

        self.doris_path = doris_parameters.doris_path
        self.cpxfiddle = doris_parameters.cpxfiddle_path  # '/...../cpxfiddle'
        self.function_path = doris_parameters.function_path
        self.input_files = input_files
        self.ESD_shift = dict()
        self.pi_shift = dict()
        self.swath = dict()
        self.overlapping = []

        if stack_folder:
            self.datastack = dict()
            self.stack_load()
            self.stack_read()

        if datastack:
            self.datastack = datastack
            self.stack_read()

        if processing_folder and not self.datastack:
            self.processing_read()

    def stack_load(self):
        # This function reads a datastack based on the data in a stack folder

        folders = next(os.walk(self.stack_folder))[1]
        folders = [fold for fold in folders if len(fold) == 8]

        for fold in folders:
            date = fold[:4] + '-' + fold[4:6] + '-' + fold[6:8]
            date_im = datetime.strptime(date,'%Y-%m-%d')

            if self.start_date <= date_im <= self.end_date:
                self.datastack[date] = dict()
                swaths = next(os.walk(os.path.join(self.stack_folder, fold)))[1]
                swaths = [fol for fol in swaths if len(fol) == 10]

                for swath in swaths:
                    self.datastack[date][swath] = dict()

                    bursts = next(os.walk(os.path.join(self.stack_folder, fold, swath)))[2]
                    bursts = [burst for burst in bursts if burst.endswith('.res')]

                    for burst in bursts:
                        slave_res = os.path.join(self.stack_folder, fold, swath, burst)
                        burst_name = 'burst_' + burst.split('.')[0][20:]
                        self.datastack[date][swath][burst_name] = ResData(slave_res)

    def stack_read(self):
        dates = self.datastack.keys()
        dates = dates.sort()

        for date_m, date_s in zip(dates[0:-1], dates[1:]):
            date_im_m = datetime.strptime(date_m,'%Y-%m-%d')
            date_im_s = datetime.strptime(date_s, '%Y-%m-%d')

            if (self.start_date <= date_m <= self.end_date) and (self.start_date <= date_s <= self.end_date):
                date = date_m + '_' + date_s
                self.stack[date] = OrderedDict()

                for swath in self.datastack[date].keys():
                    for burst in self.datastack[date][swath].keys():
                        new_name = swath + '_' + burst

                        self.stack[date][new_name] = dict()
                        self.stack[date][new_name]['slave'] = self.datastack[date][swath][burst]
                self.full_swath[date] = dict()
        if self.folder:
            print 'Because a datastack is defined using stack_folder or the datastack variable the processing folder is not read'

    def processing_read(self):
        # This function reads a processing datastack based on the processing folder

        folders = next(os.walk(self.folder))[1]
        folders = [fold for fold in folders if len(fold) == 17]

        self.stack = dict()

        for fold in folders:
            m_date = fold[:4] + '-' + fold[4:6] + '-' + fold[6:8]
            s_date = fold[9:13] + '-' + fold[13:15] + '-' + fold[15:17]
            m_datetime = datetime.strptime(m_date,'%Y-%m-%d')
            s_datetime = datetime.strptime(s_date,'%Y-%m-%d')

            if (self.start_date <= m_datetime <= self.end_date) and (self.start_date <= s_datetime <= self.end_date):
                self.stack[m_date + '_' + s_date] = dict()
                self.full_swath[m_date + '_' + s_date] = dict()

                swaths = next(os.walk(os.path.join(self.folder, fold)))[1]
                swaths = [fol for fol in swaths if len(fol) == 10]

                for swath in swaths:
                    bursts = next(os.walk(os.path.join(self.folder, fold, swath)))[1]

                    for burst in bursts:
                        burst_name = swath + '_' + burst
                        self.stack[s_date][burst_name] = dict()

        self.read_res()

    def remove_finished(self, step='filtphase'):
        # Checks which processing folders are already finished based on the final step.
        # possible steps are ('dem_assist', 'comp_coregpm', 'interfero', 'coherence', 'subtr_refphase', 'subtr_refdem',
        # 'filtphase', 'unwrap', 'geocoding'
        self.read_res()

        for date in self.stack.keys():
            if 'ifgs' in self.full_swath[date].keys():
                if self.full_swath[date]['ifgs'].process_control[step] == '1':
                    self.stack.pop(date)
                    self.full_swath.pop(date)

    def initialize(self,path='',cascade=False):
        # If path is not defined the stack will be initialized in the same folder as the datastack. This function does:
        # - copy .res and .raw files

        os.chdir(self.folder)

        for date in self.stack.keys():
            m_date = date[0:10]
            s_date = date[11:]

            date_folder = self.image_path(date)
            if not os.path.exists(date_folder):
                os.mkdir(date_folder)

            for burst in self.stack[date].keys():
                swath_folder = self.swath_path(date,burst)
                if not os.path.exists(swath_folder):
                    os.mkdir(swath_folder)

                burst_folder = self.burst_path(date,burst)
                if not os.path.exists(burst_folder):
                    os.mkdir(burst_folder)
                os.chdir(burst_folder)

                # Copy data files
                m_source = self.dat_file(burst,m_date,full_path=True,swath=True)
                m_dest = self.dat_file(burst,date='master')
                s_source = self.dat_file(burst,s_date,full_path=True,swath=True)
                s_dest = self.dat_file(burst,date='slave')
                if not os.path.exists(m_dest):
                    shutil.copy(m_source,m_dest)
                if not os.path.exists(s_dest):
                    shutil.copy(s_source,s_dest)

                # Write res files
                new_filename = os.path.join(burst_folder,'slave.res')
                if not os.path.exists(new_filename):
                    res = deepcopy(self.stack[s_date][burst]['slave'])
                    res.processes['crop']['Data_output_file'] = 'slave' + res.processes['crop']['Data_output_file'][8:]
                    res.write(new_filename=new_filename)

                new_filename = os.path.join(burst_folder,'master.res')
                if not os.path.exists(new_filename):
                    res = deepcopy(self.stack[m_date][burst]['slave'])
                    res.processes['crop']['Data_output_file'] = 'master' + res.processes['crop']['Data_output_file'][8:]
                    res.write(new_filename=new_filename)

        self.create_full_swath()

        del self.stack[self.master_date]

    def create_full_swath(self):
        # Create folders with full swath for individual interferogram.

        dates = self.stack.keys()

        # Change the res files.
        for date in dates:
            bursts = self.stack[date].keys()
            self.full_swath[date] = copy.deepcopy(self.stack[date][bursts[0]])


            # Information slave images
            az_time = self.full_swath[date]['slave'].processes['readfiles']['First_pixel_azimuth_time (UTC)']
            az_time = datetime.strptime(az_time,'%Y-%b-%d %H:%M:%S.%f')
            range_time = float(self.full_swath[date]['slave'].processes['readfiles']['Range_time_to_first_pixel (2way) (ms)'])

            # First adjust pixel and range times for pixel (1,1)
            for burst in bursts:
                az_burst = self.stack[date][burst]['slave'].processes['readfiles']['First_pixel_azimuth_time (UTC)']
                az_burst = datetime.strptime(az_burst,'%Y-%b-%d %H:%M:%S.%f')
                range_burst = float(self.stack[date][burst]['slave'].processes['readfiles']['Range_time_to_first_pixel (2way) (ms)'])

                if az_burst < az_time:
                    az_time = az_burst
                if range_burst < range_time:
                    range_time = range_burst

            az_time = az_time.strftime('%Y-%b-%d %H:%M:%S.%f')
            range_time = "{0:.15f}".format(range_time)
            self.full_swath[date]['slave'].processes['readfiles']['First_pixel_azimuth_time (UTC)'] = az_time
            self.full_swath[date]['slave'].processes['readfiles']['Range_time_to_first_pixel (2way) (ms)'] = range_time

            # Then change information on image size and crop.
            no_lines = self.full_swath[date]['slave'].processes['readfiles']['Number_of_lines_output_image']
            no_pixels = self.full_swath[date]['slave'].processes['readfiles']['Number_of_pixels_output_image']

            self.full_swath[date]['slave'].processes['readfiles']['Number_of_lines_original'] = no_lines
            self.full_swath[date]['slave'].processes['readfiles']['Number_of_pixels_original'] = no_pixels

            # Change readfiles
            self.full_swath[date]['slave'].processes['readfiles']['First_line (w.r.t. output_image)'] = '1'
            self.full_swath[date]['slave'].processes['readfiles']['Last_line (w.r.t. output_image)'] = no_lines
            self.full_swath[date]['slave'].processes['readfiles']['First_pixel (w.r.t. output_image)'] = '1'
            self.full_swath[date]['slave'].processes['readfiles']['Last_pixel (w.r.t. output_image)'] = no_pixels

            # Change in crop
            self.full_swath[date]['slave'].processes['crop']['First_line (w.r.t. original_image)'] = '1'
            self.full_swath[date]['slave'].processes['crop']['Last_line (w.r.t. original_image)'] = no_lines
            self.full_swath[date]['slave'].processes['crop']['First_pixel (w.r.t. original_image)'] = '1'
            self.full_swath[date]['slave'].processes['crop']['Last_pixel (w.r.t. original_image)'] = no_pixels
            self.full_swath[date]['slave'].processes['crop']['Data_output_file'] = date + '.raw'
            self.full_swath[date]['slave'].processes['crop'].pop('First_line (w.r.t. tiff_image)')
            self.full_swath[date]['slave'].processes['crop'].pop('Last_line (w.r.t. tiff_image)')

        # Write data to folder
        for date in dates:
            folder = self.image_path(date)
            os.chdir(folder)

            m_date = date[0:10]
            s_date = date[11:]

            # Create master.res
            master_path = self.image_path(date,'master.res')
            self.full_swath[m_date]['slave'].write(new_filename=master_path)

            # Create slave.res
            slave_path = self.image_path(date,'slave.res')
            self.full_swath[s_date]['slave'].write(new_filename=slave_path)

