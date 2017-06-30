# Function created by Gert Mulder
# Institute TU Delft
# Date 9-11-2016
# Part of Doris 5.0

# This function makes a setup of the processing folder, based on a single shapefile.
# Inputs are
# - the satellite sensor used
# - the shapefile
# - the processing folder we want to create
# - the dem source folder where the intermediate DEM data is stored

import os
from create_dem import CreateDem
from create_inputfiles import CreateInputFiles
from create_doris_input_xml import CreateDorisInputXml
from create_datastack_bash import CreateBash
import xml.etree.ElementTree as ET

class PrepareDatastack(object):

    def __init__(self):
        return

    def prepare(self, inputfile):

        # This will first create the framework with data folders the stackfolder should contain
        # a dorisparameters file.

        doris_input_xml = CreateDorisInputXml(inputfile)

        folders = ['input_files', 'stack', 'dem']
        for foldername in folders:
            if not os.path.exists(os.path.join(doris_input_xml.get_value('datastack_folder'), foldername)):
                os.mkdir(os.path.join(doris_input_xml.get_value('datastack_folder'), foldername))

        # Then create the dem file
        password_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'install', 'doris_config.xml')
        tree_doris = ET.parse(password_file)
        settings_doris = tree_doris.getroot()
        srtm_username = settings_doris.find('.usgs_username').text
        srtm_password = settings_doris.find('.usgs_password').text
        dem_out = os.path.join(doris_input_xml.get_value('dem_folder'), 'dem.raw')
        dem_var = dem_out + '.var'

        dem = CreateDem()
        if doris_input_xml.get_value('generate_dem'.lower())=='yes':
            dem.create(doris_input_xml.get_value('shape_file_path'), dem_out, dem_var, resample=None,
                                    doris_input=True, rounding=0.1, border=1.5,
                                    data_folder=doris_input_xml.get_value('dem_processing_folder'), quality='SRTM3',
                                    password=srtm_password, username=srtm_username)

        ## Then create the inputfiles
        inputfiles_folder = os.path.join(doris_input_xml.get_value('datastack_folder'), 'input_files')
        xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inputfile_template.xml')
        satellite = 'sentinel-1'
        CreateInputFiles(dem_var, xml_file, satellite).create(inputfiles_folder)

        root_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        CreateBash().create(doris_input_xml.get_value('datastack_folder'), root_folder, doris_input_xml.get_value('cores'))

