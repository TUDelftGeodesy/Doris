# This function creates the bash script in the data folder to run the final datastack.
# please note that you still have to define your start and end dates!

import os

class CreateBash(object):


    def __init__(self):
        return

    def create(self, stack_folder, root_folder, nodes):
        file_path=os.path.join(stack_folder, 'run.sh')

        f = open(file_path, 'w')

        doris_run_script = os.path.join(root_folder, 'doris_stack', 'main_code', 'doris_main.py')
        processing = stack_folder

        f.write('#!/bin/bash \n')
        f.write('\n')
        f.write('#PBS -l nodes=1:ppn=' + nodes + ' \n')
        f.write('\n')
        f.write('source_path=/data/src/Doris_s1_git/\n')
        f.write('export PYTHONPATH=$source_path/doris_stack/main_code/:$source_path/doris_stack/functions/:$PYTHONPATH \n')
        f.write('export PATH=$source_path/bin:$source_path/sar_tools:$PATH \n')
        f.write('python ' + doris_run_script + ' -p ' + processing + ' -s ' + 'yyyy-mm-dd' + ' -e ' + 'yyyy-mm-dd' + ' -m ' + 'yyyy-mm-dd \n')

        f.close()
        # make sure the file is executable
        os.chmod(file_path, 0744)
