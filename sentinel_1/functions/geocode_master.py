import os
import shutil
import resdata
import collections

# Example for the Netherlands
#dem_inputfile = '/home/gertmulder/software/sentinel/python/inputfiles/netherlands/input.createdem'
#geocode_inputfile = '/home/gertmulder/software/sentinel/python/inputfiles/netherlands/input.geocode'
#folder = '/home/gertmulder/datastacks/s1a_t88_netherlands/20151030'
#doris_dir = '/home/gertmulder/bin/doris/current_version/'

def geocode_master(folder, geocode_inputfile, dem_inputfile, doris_dir=''):
    # this function geocode the bursts of a master file, based on a DEM
    # Run the create DEM command

    os.chdir(folder)
    shutil.copy('ifgs.res','create_dem.res')
    resultfile = resdata.ResData(filename='create_dem.res')

    if not dem_inputfile == False:
        if not doris_dir:
            os.system('doris ' + dem_inputfile)
        else:
            os.system(doris_dir + ' ' + dem_inputfile)

    # Read the .res file
    resultfile.res_read()

    # Add the slant2height information. This is meant to fake the doris script
    sl2h_dat = collections.OrderedDict()
    sl2h_dat['Method'] = 'schwabisch'
    sl2h_dat['Data_output_file'] = 'dem_radar.raw'
    sl2h_dat['Data_output_format'] = 'real4'
    sl2h_dat['First_line (w.r.t. original_master)'] = resultfile.processes['comp_refdem']['First_line (w.r.t. original_master)']
    sl2h_dat['Last_line (w.r.t. original_master)'] = resultfile.processes['comp_refdem']['Last_line (w.r.t. original_master)']
    sl2h_dat['First_pixel (w.r.t. original_master)'] = resultfile.processes['comp_refdem']['First_pixel (w.r.t. original_master)']
    sl2h_dat['Last_pixel (w.r.t. original_master)'] = resultfile.processes['comp_refdem']['Last_pixel (w.r.t. original_master)']
    sl2h_dat['Multilookfactor_azimuth_direction'] = resultfile.processes['comp_refdem']['Multilookfactor_azimuth_direction']
    sl2h_dat['Multilookfactor_range_direction'] = resultfile.processes['comp_refdem']['Multilookfactor_range_direction']
    sl2h_dat['Ellipsoid (name,a,b)'] = 'WGS84 6.37814e+06 6.35675e+06'

    resultfile.insert(sl2h_dat, 'slant2h')
    resultfile.write()

    # Generate lat / lon files
    if not doris_dir:
        os.system('doris ' + geocode_inputfile)
    else:
        os.system(doris_dir + ' ' + geocode_inputfile)




