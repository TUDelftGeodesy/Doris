# This function creates networks of interferograms based on a single master stack. Basically the following steps are
# done here:
# 1. Read in of a single master stack
# 2. Definition of a network. We want to implement the following options:
#       - Cascade stack > only subsequent images are processed
#       - Full stack > all interferograms are created (This can create a huge amount of interferograms!!)
#       - Small baseline > only small baselines are processed (treshold temporal, perpendicular or added in meter/days)
#       - Baseline model > this functionality creates a baseline model based on calculated coherence and a coherence
#           threshold input: model based on perpendicular/temporal baseline and coherence treshold.
# 3. Creation of network of interferograms (based on single master stack with earth/topographic phase removed)
# 4. Phase filtering
# 5. Multilooking
# 6. Unwrapping

import numpy as np
import os, sys
from single_master_stack import SingleMaster
from datetime import datetime
from resdata import ResData
from dorisparameters import DorisParameters
from network_functions import


class NetworkStack(SingleMaster):

    def __init__(self, stack_folder='', start_date='2014-01-01',end_date='2020-10-20',master_date='',processing_folder='',
                 input_files='',doris_path='', cpxfiddle='', start_file='cint_srd.raw', start_step='comp_refdem'):

        if master_date:
            self.master_date = master_date
        else:
            print('Give the master date for the singel master stack as input.')

        doris_parameters = DorisParameters(
            os.path.dirname(os.path.dirname(processing_folder)))  # (assuming it is stored in the stackfolder)
        self.doris_parameters = doris_parameters

        self.start_date = datetime.strptime(start_date, '%Y-%m-%d')
        self.end_date = datetime.strptime(end_date, '%Y-%m-%d')

        self.nr_of_jobs = doris_parameters.nr_of_jobs
        self.parallel = doris_parameters.parallel

        self.stack = dict()
        self.datastack = dict()
        self.full_swath = dict()
        self.folder = processing_folder
        self.stack_folder = stack_folder

        if not doris_path:
            self.doris_path = doris_parameters.doris_path
        else:
            self.doris_path = doris_path
        if not cpxfiddle:
            self.cpxfiddle = doris_parameters.cpxfiddle_path
        else:
            self.cpxfiddle = cpxfiddle

        self.function_path = doris_parameters.function_path

        # Search for single master stack data in stack folder
        self.read_single_master(start_file=start_file)

    def read_single_master(self, start_file):
        # This function reads the data from our single master stack

        master_folder = self.master_date[0:4] + self.master_date[5:7] + self.master_date[8:10] + '_'

        folders = next(os.walk(self.stack_folder))[1]
        folders = [fold for fold in folders if fold.startswith(master_folder)]

        for fold in folders:
            date = fold[9:13] + '-' + fold[14:16] + '-' + fold[17:19]
            date_im = datetime.strptime(date,'%Y-%m-%d')

            if self.start_date <= date_im <= self.end_date:
                self.datastack[date] = dict()
                swaths = next(os.walk(os.path.join(self.stack_folder, fold)))[1]
                swaths = [fol for fol in swaths if len(fol) == 10]

                for swath in swaths:
                    self.datastack[date][swath] = dict()

                    bursts = next(os.walk(os.path.join(self.stack_folder, fold, swath)))[1]

                    for burst in bursts:
                        ifgs_res = os.path.join(self.stack_folder, fold, swath, burst, 'ifgs.res')
                        self.datastack[date][swath][burst] = dict()
                        self.datastack[date][swath][burst]['ifgs'] = ResData(ifgs_res)
                        self.datastack[date][swath][burst]['start_file'] = \
                            os.path.join(self.stack_folder, fold, swath, burst, start_file)


    def define_network(self, function_handle, per_baseline_max, temp_baseline_max, model_coherence_min):
        # Based on the given dates and baselines a network is created.



    def create_network(self):


    def read_res_network(self, coreg_dates=False, dates=[]):


    def update_res_network(self, dates=[]):

