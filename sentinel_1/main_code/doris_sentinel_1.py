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

        from sentinel_1.main_code.prepare_processing_folder import prepare_datastack
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
        # Create a datastack using the stack function
        stack = StackData(track_dir=track_dir,shape_dat=shape_dat,start_date=start_date,end_date=end_date,polarisation='vv',path=stack_path,db_type=2,precise_dir=precise_orbits)
        profile.log_time_stamp('StackData')
        # All images which correspond with the start and end date are selected
        stack.select_image()
        # Then these images are unzipped
        stack.unpack_image()
        # Based on the shape file bursts are selected for one date
        stack.select_burst()
        # And also for the other dates the needed bursts are selected
        stack.extend_burst()
        # Now the exact coordinates of the different burst in the concatenated image is calculated
        stack.define_burst_coordinates(slaves=True)
        # Write the datastack to the stack_path folder
        stack.write_stack(write_path=stack_path,no_data=False)
        # A few auxiliary functions which are not strictly necessary.
        # Calculate the coverage of the different sub-swaths
        stack.swath_coverage()
        # Calculate the centre point of the different burst and create a list with all available burst. (To see whether all data
        # files are available for these dates.)
        stack.lat_lon_availability()
        # Write the shapes from the bursts and swaths to a shapefile to check in a GIS program like Qgis.
        stack.write_shapes()
        profile.log_time_stamp('stack preparation finished')
        # Finally delete unzipped images
        stack.del_unpacked_image()

        # Now we import the script to create a single master interferogram
        processing = single_master_stack.SingleMaster(master_date=master_date, start_date=start_date,
                                                      end_date=end_date, stack_folder=stack_path,
                                                      input_files=input_files, processing_folder=stack_path,
                                                      doris_path=doris_path, cpxfiddle_folder=cpxfiddle_folder)

        # These lines can be used if you want to skip the initialize step because a some calculation steps are already performed....
        #del processing.stack[master_date.strftime('%Y-%m-%d')]
        #del processing.full_swath[master_date.strftime('%Y-%m-%d')]
        #processing.read_res()

        # Copy the necessary files to start processing
        processing.initialize()
        profile.log_time_stamp('initialize')
        # Calculate the coarse orbits of individual bursts
        processing.coarse_orbits()
        profile.log_time_stamp('coarse_orbits')
        # Calculate the coarse correlation of individual bursts
        processing.coarse_correlation()
        profile.log_time_stamp('coarse_correlation')
        # Gather all results from former step and calculate average shifts for the full swath
        processing.correct_coarse_correlation()
        profile.log_time_stamp('correct_coarse_correlation')
        # Deramp the data of both slave and master
        processing.deramp()
        profile.log_time_stamp('deramp')
        # Perform icc coregistration per burst
        processing.icc_burst()
        profile.log_time_stamp('icc_burst')
        # Combine coregistration windows from all bursts for full swath
        processing.coreg_full_swath()
        profile.log_time_stamp('coreg_full_swath')
        # Perform DEM coregistration for individual bursts
        processing.dac_bursts()
        profile.log_time_stamp('dac_bursts')
        # Combine results from burst DEM coregistration to full swath
        processing.dac_full_swath()
        profile.log_time_stamp('dac_full_swath')
        # Calculate polynomial for residuals icc coregistration and DEM coregistration for full swath and write to individual bursts
        processing.coreg_bursts()
        profile.log_time_stamp('coreg_bursts')
        # Resample individual bursts
        processing.resample()
        profile.log_time_stamp('resample')
        # Reramp burst
        processing.reramp(master=True)
        profile.log_time_stamp('reramp')
        # Make interferograms for individual bursts
        processing.interferogram(concatenate=True)
        profile.log_time_stamp('interferogram')

        # Perform enhanced spectral diversity for full swath
        processing.ESD()
        profile.log_time_stamp('ESD')
        # Correct for ESD shift
        processing.ESD_correct()
        profile.log_time_stamp('ESD_correct')
        # Remove resample and interferogram steps
        processing.del_process('resample' ,type='slave')
        processing.del_process('interfero' , type='ifgs')
        # Resample again with additional ESD shift
        processing.resample()
        profile.log_time_stamp('resample')
        # Reramp data based on last resampling
        processing.reramp()
        profile.log_time_stamp('reramp')
        # Create interferogram and combine to full swath
        processing.interferogram()
        profile.log_time_stamp('interferogram')
        # Combine all slave bursts to full swath
        processing.combine_slave()
        profile.log_time_stamp('combine_slave')
        # Combine all master bursts to full swath
        processing.combine_master()
        profile.log_time_stamp('combine_master')

        # Remove earth reference phase from interferograms and combine for full swath
        processing.ref_phase()
        profile.log_time_stamp('ref_phase')
        # Remove height effects from interferograms and combine for full swath
        processing.ref_dem()
        profile.log_time_stamp('ref_dem')
        # Compute coherence
        processing.coherence()
        profile.log_time_stamp('coherence')
        # Apply phase filtering
        processing.phasefilt()
        profile.log_time_stamp('phasefilt')
        # Geocode data
        processing.calc_coordinates()
        profile.log_time_stamp('calc_coordinates')
        # Multilook filtered image and coherence image
        processing.multilook(step='coherence')
        processing.multilook(step='filtphase')
        profile.log_time_stamp('multilooking')
        # Unwrap image
        processing.unwrap()
        profile.log_time_stamp('unwrapping')

        profile.log_time_stamp('end')

    print 'end sentinel 1 processing'

