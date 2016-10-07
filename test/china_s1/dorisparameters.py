from datetime import datetime
import os
from grs_config import GrsConfig

class DorisParameters():

    """This class contains all parameters that are used in the execution of test_dat_EDS_5.py.
       path parameters are checked for existance,
       all parameters are printed to stdout
    """

    def __init__(self):

        grs_config = GrsConfig()
        self.doris_path = grs_config.doris_path
        self.cpxfiddle_path = grs_config.cpxfiddle_path
        self.job_handler_script = grs_config.job_handler_script

        self.verbose = True

        self.start_date = datetime.strptime('2016-01-13','%Y-%m-%d')
        self.master_date = '2016-01-13'
        self.end_date = datetime.strptime('2016-02-06','%Y-%m-%d')

    	project_path = '/data2/Projects/datastacks/china_s1'
        data_path = '/data2/Projects/datastacks/china_s1_dat'
        #
        # used in single_master.py
        #
        #
        # used in test_dat_ESD
        #
        self.shape_dat = project_path + '/shape/AOI.shp'
        self.track_dir = data_path
        self.stack_path = project_path + '/stack/'
        self.precise_orbits = project_path + '/orbits/'
        # Start data of datastack. If end date not given it searches till current.
        self.input_files = project_path + '/input_files/'


#        self.main_code_folder = source_path + '/sentinel1/main_code/'
#        self.script_folder = source_path + '/sentinel1/functions/'

        self.parallel = True
        self.nr_of_jobs = 1
        self.between_sleep_time = 1
        self.end_sleep_time = 5
        self.initialize_flag = True
        self.profile_log = project_path + '/profile_log'
        self.doris_parallel_flag_dir = project_path + '/.Doris_parallel'
        #
        # used in Jobs
        #
#        self.job_handler_script = source_path + "/sentinel1/main_code/jobHandlerScript"
        #
        # Print parameters, check if paths exist
        #

        print 'self.shape_dat: ' + self.shape_dat
        self._check_path_exists(self.shape_dat)
        print 'self.track_dir:	' + self.track_dir
        self._check_path_exists(self.track_dir)
        print 'self.stack_path:	' + self.stack_path
        self._check_path_exists(self.stack_path)
        print 'self.precise_orbits:	' + self.precise_orbits
        self._check_path_exists(self.precise_orbits)
        print 'self.start_date:	' + str(self.start_date)
        print 'self.master_date:	' + self.master_date
        print 'self.end_date:	' + str(self.end_date)
        print 'self.input_files:	' + self.input_files
        self._check_path_exists(self.input_files)
#        print 'self.main_code_folder:	' + self.main_code_folder
#        self._check_path_exists(self.main_code_folder)
#        print 'self.script_folder:	' + self.script_folder
#        self._check_path_exists(self.script_folder)
        print 'self.nr_of_jobs:	' + str(self.nr_of_jobs)
        print 'self.initialize_flag:	' + str(self.initialize_flag)
        print 'self.jobHandlerScript:	' + self.job_handler_script
        self._check_path_exists(self.job_handler_script)

    def _check_path_exists(self, path):
        if not(os.path.exists(path)):
            print 'Error Doris_Parameters: path ' + path + ' does not exist'
            
            
