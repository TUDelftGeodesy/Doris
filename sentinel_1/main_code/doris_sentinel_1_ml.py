# This file is created by:
# Gert Mulder
# TU Delft
# 04-07-2016
##############################

import sys
import os

class DorisSentinel1(object):


    def run(self, doris_parameters_path, start_date, end_date, master_date):

        print 'start sentinel 1 processing'

        print(sys.path)

        from sentinel_1.main_code.create_datastack import prepare_datastack
        from sentinel_1.main_code.dorisparameters import DorisParameters
        from sentinel_1.main_code.grs_profile import GRS_Profile

        #Set your input variables here. You should use absolute paths.

        global dorisParameters
        dorisParameters = DorisParameters(doris_parameters_path)
        sys.path.extend([dorisParameters.function_path])

        print sys.path

        if not os.path.exists(dorisParameters.project_path):
            prepare_datastack(dorisParameters.project_path, dorisParameters.source_shapefile,
                              dorisParameters.dem_source_path, dorisParameters.satellite, doris_parameters_path)

        profile = GRS_Profile(dorisParameters.profile_log + '_' + str(dorisParameters.nr_of_jobs), dorisParameters.verbose)
        # doris executable
        doris_path = dorisParameters.doris_path

        # cpxfiddle executable
        cpxfiddle_folder = dorisParameters.cpxfiddle_path  #'/...../cpxfiddle'

        # function

        # The shapefile to select the area of interest. You can easily find a shapefile countries or regions on the internet.
        # For example via diva-gis.org. Shapefile for the Netherlands can be found in the same folder under shapes.
        shape_dat = dorisParameters.shape_dat # '/...../test.shp'

        # The folder where SLC images from a specific track are stored. This data will be used as input for the script
        track_dir = dorisParameters.track_dir # '/......'

        # This is the output folder.
        stack_path = dorisParameters.stack_path  #'/.....'

        # Folder where the precise or restituted orbits are stored. Precise orbits can be found via the following link:
        # 'https://qc.sentinel1.eo.esa.int/aux_poeorb/'. The script will assume that there are precise orbits if this folder is
        # defined, but falls back to other algorithms if needed.
        precise_orbits =  dorisParameters.precise_orbits #'/......'

        # Here the doris inputfiles are stored. (Comes with the python scripts)
        input_files = dorisParameters.input_files  #'/......'

        # Specify here the path to the functions folder with seperate python functions.
#        script_folder = dorisParameters.script_folder  #'/....../doris_v5.0.0_beta/sentinel1/functions'

        #################################

        #sys.path.append(main_code_folder)
        #sys.path.append(script_folder)

#        sys.path.append(script_folder)
        from stack import StackData
        import single_master_stack

        profile.log_time_stamp('start')

        # Now we import the script to create a single master interferogram
        processing = single_master_stack.SingleMaster(master_date=master_date, start_date=start_date,
                                                      end_date=end_date, stack_folder=stack_path,
                                                      input_files=input_files, processing_folder=stack_path,
                                                      doris_path=doris_path, cpxfiddle_folder=cpxfiddle_folder)

        # These lines can be used if you want to skip the initialize step because a some calculation steps are already performed....
        del processing.stack[master_date]
        del processing.full_swath[master_date]
        processing.read_res()

        #processing.del_process('coherence', type='ifgs')
        #processing.coherence()
        #profile.log_time_stamp('coherence')

        #processing.phasefilt()
        #profile.log_time_stamp('phasefilt')
        # Multilook filtered image and coherence image
        #processing.multilook(step='coherence', ra=100, az=25)
        processing.multilook(step='filtphase', ra=100, az=25)
        #profile.log_time_stamp('multilooking')
        # Unwrap image
        processing.del_process('unwrap', type='ifgs', images=True)
        processing.unwrap()
        profile.log_time_stamp('unwrapping')

        #profile.log_time_stamp('end')

    print 'end sentinel 1 processing'

