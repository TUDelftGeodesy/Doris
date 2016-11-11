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

from create_dem import create_binary
from create_inputfiles import CreateInputFiles
import os
import shutil

#folder = '/media/gert/Data/datastacks/zuid_holland_t88'
#shapefile = '/media/gert/Data/shapes/netherlands/zuid-holland.shp'
#satellite = 'sentinel-1'
#dem_source_folder = '/media/gert/Data/dem'


def prepare_datastack(folder, shapefile, dem_source_folder, satellite, dorisparameters):
    # This will first create the framework with data folders

    if not os.path.exists(folder):
        os.mkdir(folder)

    folders = ['dem', 'input_files', 'shape', 'stack']
    for foldername in folders:
        if not os.path.exists(os.path.join(folder, foldername)):
            os.mkdir(os.path.join(folder, foldername))

    # move the shapefile
    suffixes = ['.shp', '.cpg', '.dbf', '.prj', '.qpj', '.shx']

    for suffix in suffixes:
        shape_dest = os.path.join(folder, 'shape', 'AOI' + suffix)
        shutil.copyfile(shapefile[:-4] + suffix, shape_dest)

    # Then create the dem file
    dem_out = os.path.join(folder, 'dem', 'dem.raw')
    dem_out, dem_var, dem_inputfile = create_binary(os.path.join(folder, 'shape', 'AOI.shp'), dem_out, resample=None,
                                                    doris_input=True, rounding=0.1, border=1.5,
                                                    data_folder=dem_source_folder, quality='SRTM1')

    # Then create the inputfiles
    inputfiles_folder = os.path.join(folder, 'input_files')
    CreateInputFiles(dem_var, inputfiles_folder, 'inputfile.xml', satellite)

    # Finally copy dorisparameters file to stack folder.
    dorisparemeters_new = os.path.join(folder, 'dorisparameters.py')
    shutil.copy(dorisparameters, dorisparemeters_new)
