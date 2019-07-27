from datetime import datetime
import os
from .doris_config import DorisConfig
import xml.etree.ElementTree as ET


class DorisParameters():

    """This class contains all parameters that are used in the execution of test_dat_EDS_5.py.
       path parameters are checked for existance,
       all parameters are printed to stdout
    """

    def __init__(self, stack_path):

        grs_config = DorisConfig()

        self.doris_path = grs_config.doris_path
        self.cpxfiddle_path = grs_config.cpxfiddle_path
        self.job_handler_script = grs_config.job_handler_script
        self.function_path = grs_config.function_path
        self.source_path = grs_config.source_path

        self.verbose = True

        tree = ET.parse(os.path.join(stack_path, 'doris_input.xml'))
        self.settings = tree.getroot()

        project_path = self._settings_get('.datastack_folder')
        self.project_path = project_path
        data_path = self._settings_get('.sar_data_folder')
        self.data_path = data_path
        #
        # used in single_master.py
        #
        #
        # used in test_dat_ESD
        #
        # TODO DLE FIX shape path
        self.shape_dat = self._settings_get('.shape_file_path')
        self.track_dir = data_path
        self.stack_path = project_path + '/stack/'
        self.precise_orbits = self._settings_get('.orbits_folder')
        # Start data of datastack. If end date not given it searches till current.
        self.input_files = project_path + '/input_files/'

        self.parallel = self._settings_compare('.parallel', 'yes')
        self.nr_of_jobs = int(self._settings_get('.cores'))
        self.initialize_flag = self._settings_compare('.initialize_flag', 'yes')

        self.profile_log = project_path + '/profile_log'
        self.doris_parallel_flag_dir = project_path + '/.Doris_parallel'
        self.between_sleep_time = 1
        self.end_sleep_time = 1

        self.do_coarse_orbits = self._settings_compare('.do_coarse_orbits', 'yes')
        self.do_deramp = self._settings_compare('.do_deramp', 'yes')
        self.do_fake_fine_coreg_bursts = self._settings_compare('.do_fake_fine_coreg_bursts', 'yes')
        self.do_dac_bursts = self._settings_compare('.do_dac_bursts', 'yes')
        self.do_fake_coreg_bursts = self._settings_compare('.do_fake_coreg_bursts', 'yes')
        self.do_resample = self._settings_compare('.do_resample', 'yes')
        self.do_reramp = self._settings_compare('.do_reramp', 'yes')
        self.do_interferogram = self._settings_compare('.do_interferogram', 'yes')
        self.do_compref_phase = self._settings_compare('.do_compref_phase', 'yes')
        self.do_compref_dem = self._settings_compare('.do_compref_dem', 'yes')
        self.do_coherence = self._settings_compare('.do_coherence', 'yes')
        self.do_esd = self._settings_compare('.do_esd', 'yes')
        self.do_network_esd = self._settings_compare('.do_network_esd', 'yes')
        self.do_ESD_correct = self._settings_compare('.do_ESD_correct', 'yes')
        self.do_ref_phase = self._settings_compare('.do_ref_phase', 'yes')
        self.do_ref_dem = self._settings_compare('.do_ref_dem', 'yes')
        self.do_phasefilt = self._settings_compare('.do_phasefilt', 'yes')
        self.do_calc_coordinates = self._settings_compare('.do_calc_coordinates', 'yes')
        self.do_multilooking = self._settings_compare('.do_multilooking', 'yes')
        self.do_unwrap = self._settings_compare('.do_unwrap', 'yes')
        #
        # used in Jobs
        #
        # self.job_handler_script = source_path + "/sentinel1/main_code/jobHandlerScript"
        #
        # Print parameters, check if paths exist
        #

        print('self.shape_dat: ' + self.shape_dat)
        # self._check_path_exists(self.shape_dat)
        print('self.track_dir:	' + self.track_dir)
        self._check_path_exists(self.track_dir)
        print('self.stack_path:	' + self.stack_path)
        # self._check_path_exists(self.stack_path)
        print('self.precise_orbits:	' + self.precise_orbits)
        self._check_path_exists(self.precise_orbits)
        print('self.input_files:	' + self.input_files)
        # self._check_path_exists(self.input_files)
#        print('self.main_code_folder:	' + self.main_code_folder)
#        self._check_path_exists(self.main_code_folder)
#        print('self.script_folder:	' + self.script_folder)
#        self._check_path_exists(self.script_folder)
        print('self.nr_of_jobs:	' + str(self.nr_of_jobs))
        print('self.initialize_flag:	' + str(self.initialize_flag))
        print('self.jobHandlerScript:	' + self.job_handler_script)
        self._check_path_exists(self.job_handler_script)

    def _check_path_exists(self, path):
        if not(os.path.exists(path)):
            print('Error Doris_Parameters: path ' + path + ' does not exist')
            
    def _settings_get(self, string):
        return self.settings.find('*/' + string).text


    def _settings_compare(self, string, comp_string):
        if (self.settings.find('*/' + string).text.lower()==comp_string.lower()):
            return True
        return False
