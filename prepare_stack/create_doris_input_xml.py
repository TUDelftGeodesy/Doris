import xml.etree.ElementTree as ET
import os

class CreateDorisInputXml(object):

    def __init__(self, input_file):
        self.input_file_dict={}
        if(len(input_file)==0):
            self._create_xml()
        else:
            self._read_xml(input_file)


    def _create_xml(self):
        # This will first create the framework with data folders the stackfolder should contain
        # a dorisparameters file.


        input = False
        while input == False:
            user_input = raw_input("Enter the path to the folder of the sar data: ")
            if os.path.exists(user_input):
                self.input_file_dict['sar_data_folder'] = user_input
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            self.input_file_dict['datastack_folder'] = raw_input("Enter the path to the folder of new datastack: ")
            if os.path.exists(self.input_file_dict['datastack_folder']):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            self.input_file_dict['shape_file_path'] = raw_input("Enter full path to the shapefile: ")
            if os.path.exists(self.input_file_dict['shape_file_path']) and self.input_file_dict['shape_file_path'].endswith('.shp'):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = raw_input("Enter the path to the folder of the orbit files: ")
            if os.path.exists(user_input):
                self.input_file_dict['orbits_folder'] = user_input
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = raw_input("Do you want to generate the DEM file automaticly (Yes/No): ").lower()
            if user_input == 'yes' or user_input == 'no':
                self.input_file_dict['generate_dem'] = user_input
                input = True
            else:
                print('You should use either yes or no')

        input = False
        while input == False:
            self.input_file_dict['dem_folder'] = raw_input("Enter path to the dem folder: ")
            if os.path.exists(self.input_file_dict['dem_folder']):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = raw_input("Do you want to use parallel computing (Yes/No): ").lower()
            if user_input == 'yes' or user_input == 'no':
                self.input_file_dict['parallel'] = user_input
                input = True
            else:
                print('You should use either yes of no')

        if user_input == 'yes':
            nodes = raw_input("How many cores do you want to use: ")
            self.input_file_dict['cores'] = nodes

        xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'doris_input_template.xml')
        tree = ET.parse(xml_file)
        self.settings = tree.getroot()

        for key in self.input_file_dict.keys():
            self.settings.find('*/' + key).text = self.input_file_dict.get(key)

        tree.write(os.path.join(self.input_file_dict['datastack_folder'], 'doris_input.xml'))

        return self.input_file_dict

    def _read_xml(self, input_file):
        xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), input_file)
        tree = ET.parse(xml_file)
        self.settings = tree.getroot()

    def get_value(self, key):
        return self.settings.find('*/' + key).text
