'''
	GrsConfig defines paths that are local to the source tree.
	They are copied into DorisParameters for use in Doris python scripts
'''

import xml.etree.ElementTree as ET
import sys, os

class GrsConfig(object):

    def __init__(self):

        xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'grs_config.xml')
        tree = ET.parse(xml_file)
        settings = tree.getroot()

        self.source_path = settings.find('.source_path').text
        self.doris_path = settings.find('.doris_path').text
        self.cpxfiddle_path = settings.find('.cpxfiddle_path').text

        self.job_handler_script = self.source_path + "/sentinel_1/main_code/jobHandlerScript"
        self.function_path = self.source_path + "/sentinel_1/functions/"
        self.main_code_path = self.source_path + "/sentinel_1/main_code/"

        # Extend path
        sys.path.extend([self.function_path])
        sys.path.extend([self.main_code_path])
