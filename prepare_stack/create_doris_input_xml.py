import xml.etree.ElementTree as ET
import os
from datetime import datetime

class CreateDorisInputXml(object):

    def __init__(self, input_file):
        self.input_file_dict={}
        self.dem_folder=''
        if(len(input_file)==0):
            self._create_xml()
        else:
            self._read_xml(input_file)


    def _create_xml(self):
        # This will first create the framework with data folders the stackfolder should contain
        # a dorisparameters file.

        input = False
        while input == False:
            user_input = input("Enter the path to the archive data folder: ")
            if os.path.exists(user_input):
                self.input_file_dict['sar_data_folder'] = user_input
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = input("Which polarisation do you want to use (vv,hh,vh,hv): ")
            if user_input in ['vv', 'hh', 'vh', 'hv']:
                self.input_file_dict['polarisation'] = user_input
                input = True
            else:
                print('This polarisation does not exist')

        input = False
        while input == False:
            user_input = input("Which track do you want to work with? (explore on https://scihub.copernicus.eu/dhus/) : ")
            try:
                input = str(int(user_input)).zfill(3)
                self.input_file_dict['track'] = user_input
                input = True
            except:
                print('This track does not exist')

        input = False
        while input == False:
            user_input = input("Is this track ascending or descending? (asc/dsc) : ")
            if user_input in ['asc', 'dsc']:
                self.input_file_dict['direction'] = user_input
                input = True
            else:
                print('Input should either be asc or dsc')

        input = False
        while input == False:
            self.input_file_dict['datastack_folder'] = input("Enter the path to the folder of new datastack: ")
            if os.path.exists(self.input_file_dict['datastack_folder']):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            self.input_file_dict['shape_file_path'] = input("Enter full path to the shapefile: ")
            if os.path.exists(self.input_file_dict['shape_file_path']) and self.input_file_dict['shape_file_path'].endswith('.shp'):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = input("Enter the path to the folder of the orbit files: ")
            if os.path.exists(user_input):
                self.input_file_dict['orbits_folder'] = user_input
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = input("Do you want to generate the DEM file automaticly (Yes/No): ").lower()
            if user_input == 'yes' or user_input == 'no':
                self.input_file_dict['generate_dem'] = user_input
                input = True
            else:
                print('You should use either yes or no')

        input = False
        while input == False:
            self.input_file_dict['dem_processing_folder'] = input("Enter path to the dem folder: ")
            self.input_file_dict['dem_folder'] = os.path.join(self.input_file_dict['datastack_folder'], 'dem')
            if os.path.exists(self.input_file_dict['dem_processing_folder']):
                input = True
            else:
                print('The path is incorrect, use another path')

        input = False
        while input == False:
            user_input = input("Do you want to use parallel computing (Yes/No): ").lower()
            if user_input == 'yes' or user_input == 'no':
                self.input_file_dict['parallel'] = user_input
                input = True
            else:
                print('You should use either yes of no')

        if user_input == 'yes':
            nodes = input("How many cores do you want to use: ")
            self.input_file_dict['cores'] = nodes

        input = False
        while input == False:
            user_input = input("What is the start date of your stack in yyyy-mm-dd (can be changed later): ").lower()
            try:
                date = datetime.strptime(user_input, '%Y-%m-%d')
                self.input_file_dict['start_date'] = user_input
                input = True
            except:
                print('Format not recognized, 01-01-2014 chosen')
                self.input_file_dict['start_date'] = user_input

        input = False
        while input == False:
            user_input = input("What is the end date of your stack in yyyy-mm-dd (can be changed later): ").lower()
            try:
                date = datetime.strptime(user_input, '%Y-%m-%d')
                self.input_file_dict['end_date'] = user_input
                input = True
            except:
                print('Format not recognized, 01-01-2050 chosen')
                self.input_file_dict['end_date'] = user_input

        input = False
        while input == False:
            user_input = input("What is the master date of your stack in yyyy-mm-dd (can be changed later): ").lower()
            try:
                date = datetime.strptime(user_input, '%Y-%m-%d')
                self.input_file_dict['master_date'] = user_input
                input = True
            except:
                print('Format not recognized, 01-01-2016 chosen. Check https://scihub.copernicus.eu/dhus/#/home for valid date')
                self.input_file_dict['master_date'] = user_input

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
