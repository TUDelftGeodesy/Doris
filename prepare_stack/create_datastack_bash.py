# This function creates the bash script in the data folder to run the final datastack.
# please note that you still have to define your start and end dates!

import os
import xml.etree.ElementTree as ET

class CreateBash(object):


    def __init__(self):
        return

    def create(self, stack_folder, root_folder, nodes):

        xml_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                'install/doris_config.xml')
        tree = ET.parse(xml_file)
        settings = tree.getroot()

        source_path = settings.find('.source_path').text
        doris_folder = os.path.dirname(settings.find('.doris_path').text)
        cpxfiddle_folder = os.path.dirname(settings.find('.cpxfiddle_path').text)
        snaphu_folder = os.path.dirname(settings.find('.snaphu_path').text)

        file_path=os.path.join(stack_folder, 'doris_stack.sh')

        f = open(file_path, 'w')

        doris_run_script = os.path.join(source_path, 'doris_stack', 'main_code', 'doris_main.py')
        processing = stack_folder

        f.write('#!/bin/bash \n')
        f.write('\n')
        f.write('#PBS -l nodes=1:ppn=' + nodes + ' \n')
        f.write('\n')
        f.write('source_path=' + source_path + '\n')
        f.write('export PYTHONPATH=$source_path/doris_stack/main_code/:$source_path/doris_stack/functions/:$PYTHONPATH \n')
        f.write('export PATH=' + doris_folder + ':' + cpxfiddle_folder + ':' + snaphu_folder + ':' + '$PATH \n')
        f.write('python ' + doris_run_script + ' -p ' + processing + ' -s ' + 'yyyy-mm-dd' + ' -e ' + 'yyyy-mm-dd' + ' -m ' + 'yyyy-mm-dd \n')

        f.close()
        # make sure the file is executable
        os.chmod(file_path, 0744)
