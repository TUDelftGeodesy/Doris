# This function creates the bash script in the data folder to run the final datastack.
# please note that you still have to define your start and end dates!

import os

def create_bash(stack_folder, python_source_path, nodes):

    f = open(os.path.join(stack_folder, 'run.sh'), 'w')

    doris_run_script = os.path.join(python_source_path, 'main_code', 'doris_main.py')
    processing = stack_folder

    f.write('#!/bin/bash \n')
    f.write('\n')
    f.write('#PBS -l nodes=1:ppn=' + nodes + ' \n')
    f.write('\n')
    f.write('export PYTHONPATH=' + python_source_path + '/main_code/:' + python_source_path + '/functions/:$PYTHONPATH \n')
    f.write('export PATH=' + python_source_path + '/main_code/:' + python_source_path + '/functions/:$PATH \n')
    f.write('python ' + doris_run_script + ' -p ' + processing + ' -s ' + 'yyyy-mm-dd' + ' -e ' + 'yyyy-mm-dd' + ' -m ' + 'yyyy-mm-dd \n')

    f.close()