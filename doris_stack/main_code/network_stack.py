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
from sentinel_1.main_code.single_master_stack import SingleMaster

class NetworkStack(SingleMaster):

    def __init__(self, stack_folder='', start_date='',end_date='',master_date='',processing_folder='',
                 input_files='',doris_path='', cpxfiddle_folder=''):

        # Search for single master stack data in stack folder

        #

    def define_network(self, function_handle, per_baseline_max, temp_baseline_max, model_coherence_min):


    def create_network(self):


    def read_res_network(self, coreg_dates=False, dates=[]):


    def update_res_network(self, dates=[]):

