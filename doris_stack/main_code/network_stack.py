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
from doris_parameters import DorisParameters
from define_network import cascade_network, threshold_network
import copy
from jobs import Jobs


class NetworkStack(SingleMaster):

    def __init__(self, stack_folder='', start_date='2014-01-01',end_date='2020-10-20',master_date='',processing_folder='',
                 input_files='',doris_path='', cpxfiddle='', start_file='cint_srd.raw'):

        if master_date:
            self.master_date = master_date
        else:
            print('Give the master date for the singel master stack as input.')

        doris_parameters = DorisParameters(
            os.path.dirname(os.path.dirname(processing_folder)))  # (assuming it is stored in the stackfolder)
        self.doris_parameters = doris_parameters

        Jobs.id = 0

        self.nr_of_jobs = doris_parameters.nr_of_jobs
        self.parallel = doris_parameters.parallel
        self.doris_path = doris_parameters.doris_path
        self.cpxfiddle = doris_parameters.cpxfiddle_path  # '/...../cpxfiddle'
        self.function_path = doris_parameters.function_path

        self.start_date = datetime.strptime(start_date, '%Y-%m-%d')
        self.end_date = datetime.strptime(end_date, '%Y-%m-%d')

        self.nr_of_jobs = doris_parameters.nr_of_jobs
        self.parallel = doris_parameters.parallel

        self.stack = dict()
        self.datastack = dict()
        self.full_swath_datastack = dict()
        self.full_swath = dict()
        self.folder = processing_folder
        self.stack_folder = stack_folder
        self.ifgs_list = []
        self.ifgs_keys = []

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
                    bursts = next(os.walk(os.path.join(self.stack_folder, fold, swath)))[1]

                    for burst in bursts:
                        ifgs_res = os.path.join(self.stack_folder, fold, swath, burst, 'ifgs.res')
                        burst_key = swath + '_' + burst
                        self.datastack[date][burst_key] = dict()
                        self.datastack[date][burst_key]['ifgs'] = ResData(ifgs_res)
                        self.datastack[date][burst_key]['start_file'] = \
                            os.path.join(self.stack_folder, fold, swath, burst, start_file)

            # Read in the data for the full image.
            ifgs_res = os.path.join(self.stack_folder, fold, 'ifgs.res')
            self.full_swath_datastack[date] = dict()
            self.full_swath_datastack[date]['ifgs'] = ResData(ifgs_res)
            self.full_swath_datastack[date]['start_file'] = \
                os.path.join(self.stack_folder, fold, start_file)


    def read_network(self):
        # This function is used to read an already existing network of ifgs.

        folders = next(os.walk(self.stack_folder))[1]
        folders = [fold for fold in folders if fold.endswith('_ifg')]

        for fold in folders:
            date_1 = fold[9:13] + '-' + fold[14:16] + '-' + fold[17:19]
            date_2 = fold[0:4] + '-' + fold[5:7] + '-' + fold[8:10]
            date_im1 = datetime.strptime(date_1,'%Y-%m-%d')
            date_im2 = datetime.strptime(date_2, '%Y-%m-%d')

            # Check if it falls between start and end date.
            if self.start_date <= date_im1 <= self.end_date and self.start_date <= date_im2 <= self.end_date:
                im_key = date_1 + '_' + date_2
                self.stack[im_key] = dict()
                swaths = next(os.walk(os.path.join(self.stack_folder, fold)))[1]
                swaths = [fol for fol in swaths if len(fol) == 7]

                for swath in swaths:
                    bursts = next(os.walk(os.path.join(self.stack_folder, fold, swath)))[1]

                    # Load the .res files for the Bursts
                    for burst in bursts:
                        ifgs_res = os.path.join(self.stack_folder, fold, swath, burst, 'ifgs.res')
                        master_res = os.path.join(self.stack_folder, fold, swath, burst, 'master.res')
                        slave_res = os.path.join(self.stack_folder, fold, swath, burst, 'slave.res')
                        burst_key = swath + '_' + burst
                        self.stack[im_key][burst_key] = dict()
                        self.stack[im_key][burst_key]['ifgs'] = ResData(ifgs_res)
                        self.stack[im_key][burst_key]['master'] = ResData(master_res)
                        self.stack[im_key][burst_key]['slave'] = ResData(slave_res)

                # Load the .res files for the images
                ifgs_res = os.path.join(self.stack_folder, fold, 'ifgs.res')
                if os.path.exists(ifgs_res):
                    self.full_swath[im_key] = dict()
                    self.full_swath[im_key] = ResData(ifgs_res)


    def define_network(self, per_baseline_max, temp_baseline_max, model_coherence_min, network_type='cascade'):
        # Based on the given dates and baselines a network is created.

        dates = self.datastack.keys()
        swath = self.datastack[dates[0]].keys()[0]
        baselines = [float(self.datastack[d][swath]['burst_1']['ifgs'].processes['coarse_orbits']['Bperp']) for d in dates]

        if network_type == 'cascade':
            self.ifgs_list, self.ifgs_keys = cascade_network(dates)
        elif network_type == 'threshold':
            self.ifgs_list, self.ifgs_keys = threshold_network(dates, baselines, per_baseline_max, temp_baseline_max)
        else:
            print('Other networks than cascade and temporal/perpendicular threshold are not implemented yet!')


    def create_network(self):
        # This function creates the network interferograms for the single master stack. This includes:
        # > creation of new data folders
        # > copying of the ifgs.res to the right folders
        # > making ifgs from former ifgs to create the new ifgs. (This is a current hack, which should be improved
        #       later on.)

        # First create folder structure and master/slave result files
        for ifg_key, ifg_pair in zip(self.ifgs_keys, self.ifgs_list):
            # Loop over all the different ifgs pairs.
            folder = os.path.join(self.stack_folder, 'ifg_key' + '_ifg')
            if not os.path.exists(folder):
                os.mkdir(folder)
            self.stack[ifg_key] = dict()

            date_1 = ifg_key[:10]
            date_2 = ifg_key[11:21]
            for burst_key in self.datastack[date_1].keys():
                swath = burst_key[:7]
                burst = burst_key[8:]

                swath_folder = os.path.join(self.stack_folder, ifg_key + '_ifg', swath)
                if not os.path.exists(swath_folder):
                    os.mkdir(swath_folder)

                burst_folder = os.path.join(self.stack_folder, ifg_key + '_ifg', swath, burst)
                if not os.path.exists(burst_folder):
                    os.mkdir(burst_folder)
                ifgs_path = os.path.join(self.stack_folder, ifg_key + '_ifg', swath, burst, 'ifgs.res')
                master_path = os.path.join(self.stack_folder, ifg_key + '_ifg', 'master.res')
                slave_path = os.path.join(self.stack_folder, ifg_key + '_ifg', 'slave.res')
                self.stack[ifg_key][burst_key] = dict()

                if date_1 == self.master_date:
                    self.stack[ifg_key][burst_key]['ifgs'] = copy.deepcopy(self.datastack[date_2][burst_key]['ifgs'])
                else:
                    self.stack[ifg_key][burst_key]['ifgs'] = copy.deepcopy(self.datastack[date_1][burst_key]['ifgs'])

                # For the master and slave files always the slave files are used except the single_master master date
                # is used.
                if not date_1 == self.master_date:
                    self.stack[ifg_key][burst_key]['master'] = copy.deepcopy(self.datastack[date_1][burst_key]['slave'])
                else:
                    self.stack[ifg_key][burst_key]['master'] = copy.deepcopy(self.datastack[date_2][burst_key]['master'])
                if not date_2 == self.master_date:
                    self.stack[ifg_key][burst_key]['slave'] = copy.deepcopy(self.datastack[date_2][burst_key]['slave'])
                else:
                    self.stack[ifg_key][burst_key]['slave'] = copy.deepcopy(self.datastack[date_1][burst_key]['master'])

                self.stack[ifg_key][burst_key]['ifgs'].res_path = ifgs_path
                self.stack[ifg_key][burst_key]['ifgs'].res_path = master_path
                self.stack[ifg_key][burst_key]['ifgs'].res_path = slave_path

            # Create the full_swath resfiles.
            ifgs_path = os.path.join(self.stack_folder, ifg_key + '_ifg', 'ifgs.res')
            self.full_swath[ifg_key] = dict()
            self.full_swath[ifg_key]['ifgs'] = copy.deepcopy(self.full_swath_datastack[date_1]['ifgs'])
            self.full_swath[ifg_key]['ifgs'].res_path = ifgs_path

        # Remove the interferogram, coherence, filtering and unwrapping steps if needed.
        self.del_process('interfero', type='ifgs')
        self.del_process('coherence', type='ifgs')
        self.del_process('filtphase', type='ifgs')
        self.del_process('unwrap', type='ifgs')
        # Now write the .res files to disk.
        self.update_res()

        # Create the corrected slave images (earth phase / dem phase)
        # This step is done per burst and paralellized, as we will load the
        date = self.stack.keys()[0]
        job_list = []
        for burst in self.stack[date].keys():
            path = self.stack_folder
            command = self.function_path + ' corrected_slave.py ' + self.stack_folder + ' ' + burst + ' ' + self.master_date
            job_list.append({"path": path, "command": command})
            if (not(self.parallel)):
                os.chdir(path)
                os.system(command)
        if (self.parallel):
            jobs = Jobs(self.nr_of_jobs, self.doris_parameters)
            jobs.run(job_list)

        # Create links to the different images for coherence calculations.
        for ifg_key, ifg_pair in zip(self.ifgs_keys, self.ifgs_list):
            # Loop over all the different ifgs pairs.
            date_1 = ifg_key[:10]
            date_2 = ifg_key[11:21]
            master_key = self.master_date[:4] + self.master_date[5:7] + self.master_date[8:]

            for burst_key in self.datastack[ifg_key].keys():
                swath = burst_key[:7]
                burst = burst_key[8:]
                new_slave = os.path.join(self.stack_folder, ifg_key + '_ifg', swath, burst, 'slave_corrected.raw')
                new_master = os.path.join(self.stack_folder, ifg_key + '_ifg', swath, burst, 'master_corrected.raw')

                # Make a link to the right data file depending on whether we have a slave or master image.
                if not date_1 == self.master_date:
                    old_master = os.path.join(self.stack_folder, master_key + '_' + date_1, swath, burst,
                                self.stack[ifg_key][burst_key]['master'].processes['resample']['Data_output_file'])
                    self.stack[ifg_key][burst_key]['master'].processes['resample']['Data_output_file'] = 'master_corrected.raw'
                else:
                    old_master = os.path.join(self.stack_folder, date_1 + '_' + date_2, swath, burst,
                                 self.stack[ifg_key][burst_key]['master'].processes['crop']['Data_output_file'])
                    self.stack[ifg_key][burst_key]['master'].processes['crop']['Data_output_file'] = 'master_corrected.raw'
                if not date_2 == self.master_date:
                    old_slave = os.path.join(self.stack_folder, master_key + '_' + date_2, swath, burst,
                                self.stack[ifg_key][burst_key]['slave'].processes['resample']['Data_output_file'])
                    self.stack[ifg_key][burst_key]['slave'].processes['resample']['Data_output_file'] = 'slave_corrected.raw'
                else:
                    old_slave = os.path.join(self.stack_folder, date_2 + '_' + date_1 , swath, burst,
                                self.stack[ifg_key][burst_key]['slave'].processes['crop']['Data_output_file'])
                    self.stack[ifg_key][burst_key]['slave'].processes['crop']['Data_output_file'] = 'slave_corrected.raw'

                os.link(old_master, new_master)
                os.link(old_slave, new_slave)

        self.update_res()


    def skip_earth_dem_phase(self):
        # This function is to skip the generation earth phase and dem phase images, because that is done already.
        # Basically it renames the created ifgs to the output for the ifgs where reference phases are subtracted.

        self.read_res()

        for ifg_key in self.stack.keys():
            for burst_key in self.stack[ifg_key][burst_key]:

                if (self.stack[ifg_key][burst_key]['ifgs'].process_control['interfero'] == '1' and
                       self.stack[ifg_key][burst_key]['ifgs'].process_control['subt_refdem'] == '1'):

                    path = os.path.join(self.stack_folder, ifg_key + '_ifg', burst_key[:7], burst_key[8:])
                    os.rename(os.path.join(path, 'cint.raw'), os.path.join(path, 'cint_srd.raw'))

                else:
                    print('First create interferograms from the corrected slaves.')




class APS_network(NetworkStack):
    # This class is specifically used to solve the ambiguity in APS per epoch based on a network approach. This is based
    # on a former defined network of interferograms, which are already unwrapped. Coherence is used to improve the
    # quality of the final product. To solve the network several approaches are implemented:
    # 1. A solution referring to the average values of the different interferograms.
    # 2. A solution which includes a set of user defined model outcomes. Generally, these will be based on good model
    #       results during clear winter days.
    # 3. A solution where all model solutions are included


    def include_Harmonie(self):
        # This function includes Harmonie input data.
        print('Working on this 2')

    def create_A_matrix(self):
        # Here the A matrix is made
        print('Working on this!')


