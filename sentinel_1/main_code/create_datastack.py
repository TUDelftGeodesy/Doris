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

import os, sys
import shutil

if __name__ == "__main__":
    # If calling script directly we have to load the package first to our python path
    folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(folder)
    sys.path.extend([folder])

from sentinel_1.main_code.create_dem import create_binary
from sentinel_1.main_code.create_inputfiles import CreateInputFiles
from sentinel_1.main_code.init_datastack_xml import init_datastack_xml
from sentinel_1.main_code.dorisparameters import DorisParameters
from sentinel_1.main_code.init_datastack_bash import create_bash


def prepare_datastack():
    # This will first create the framework with data folders the stackfolder should contain
    # a dorisparameters file.

    input = False
    while input == False:
        folder = raw_input("Enter path of new datastack: ")
        if os.path.exists(folder):
            input = True
        else:
            print('The path is incorrect, use another path')

    input = False
    while input == False:
        dem_folder = raw_input("Enter path to folder with raw dem files: ")
        if os.path.exists(dem_folder):
            input = True
        else:
            print('The path is incorrect, use another path')

    input = False
    while input == False:
        shapefile = raw_input("Enter path to shapefile: ")
        if os.path.exists(shapefile) and shapefile.endswith('.shp'):
            input = True
        else:
            print('The path is incorrect, use another path')

    satellite = 'sentinel-1'

    # input = False
    # while input == False:
    #    shapefile = raw_input("Enter sensor type:  (default is sentinel-1)")
    #    if os.path.exists(shapefile) and shapefile.endswith('.shp'):
    #        input = True
    #    else:
    #        print('The path is incorrect, use another path')

    if not os.path.exists(folder):
        os.mkdir(folder)

    # Now initialize the first .xml file with basic information
    nodes = init_datastack_xml(folder)
    # And create the .bash script
    python_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    create_bash(folder, python_folder, nodes)

    dorisparameters = DorisParameters(folder)

    folders = ['dem', 'input_files', 'shape', 'stack']
    for foldername in folders:
        if not os.path.exists(os.path.join(folder, foldername)):
            os.mkdir(os.path.join(folder, foldername))

    # move the shapefile
    suffixes = ['.shp', '.dbf', '.prj', '.qpj', '.shx']

    for suffix in suffixes:
        shape_dest = os.path.join(folder, 'shape', 'AOI' + suffix)
        shutil.copyfile(shapefile[:-4] + suffix, shape_dest)

    # Then create the dem file
    dem_out = os.path.join(folder, 'dem', 'dem.raw')
    dem_out, dem_var, dem_inputfile = create_binary(os.path.join(folder, 'shape', 'AOI.shp'), dem_out, resample=None,
                                                    doris_input=True, rounding=0.1, border=1.5,
                                                    data_folder=dem_folder, quality='SRTM1')

    # Then create the inputfiles
    inputfiles_folder = os.path.join(folder, 'input_files')
    xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inputfile.xml')
    CreateInputFiles(dem_var, inputfiles_folder, xml_file, satellite)


# Actually execute the code...
if __name__ == "__main__":

    prepare_datastack()
