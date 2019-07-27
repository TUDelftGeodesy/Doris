# Function created by Gert Mulder
# Institute TU Delft
# Date 9-11-2016
# Part of Doris 5.0

# This function creates a dem based on either a shape/kml file or a given bounding box. If a shape/kml file is given a
# minimum offset of about 0.1 degrees is used.
# All grids are based on the WGS84 projection. 
# Downloaded data is based on SRTM void filled data:
# Documentation: https://lpdaac.usgs.gov/sites/default/files/public/measures/docs/NASA_SRTM_V3.pdf

# Description srtm data: https://lpdaac.usgs.gov/dataset_discovery/measures/measures_products_table/SRTMGL1_v003
# Description srtm q data: https://lpdaac.usgs.gov/node/505

import numpy as np
import gdal
import gdalconst
import osr
from HTMLParser import HTMLParser
import pickle
import requests
import os
import zipfile
import fastkml
import shutil
from fiona import collection
import xml.etree.ElementTree as ET
import sys

# This are some example values...
# shape_filename = '/media/gert/Data/shapes/netherlands/netherland.shp'
# data_folder = '/media/gert/Data/dem/'
# dem_path = '/media/gert/Data/dem/test/'
# resample = 'regular_grid'
# doris_input = True
# quality = 'SRTM1'
# out_file = dem_path + 'output.dem'
# lats = np.arange(50.1, 53, 0.0005)
# lons = np.arange(2, 8, 0.0005)
# rounding = 1
# border = 0.1
# create_binary('', out_file, resample, doris_input, lats, lons, rounding, border, data_folder, quality)


class CreateDem:

    def __init__(self):
        return

    def create(self, shape_filename='', out_file='' ,var_file='', resample='regular_grid', doris_input=True, lats=[], lons=[],
                      rounding=1, border=0.1, data_folder='', quality='SRTM1', password='', username='', password_file=''):
        # This function stitches the different files together. If no data is available values will be zero. Which is
        # generally true because it is above sealevel.
        # The resampling is either
        # - none
        # - regular_grid (based on vectors of lats/lons)
        # - irregular_grid (based on lats/lons)

        if not out_file and not resample == 'irregular_grid':
            print('Please specify an output filename!')

        if shape_filename:
            latlim, lonlim = self._kml_shp_2_bb(shape_filename)
        else:
            try:
                latlim = [min(lats), max(lats)]
                lonlim = [min(lons), max(lons)]
            except:
                print('Not possible to create DEM grid.')
                return

        # Now add the rounding and borders to add the sides of our image
        # Please use rounding as a n/60 part of a degree (so 1/15 , 1/10 or 1/20 of a degree for example..)
        latlim = [np.floor((latlim[0] - border) / rounding) * rounding, np.ceil((latlim[1] + border) / rounding) * rounding]
        lonlim = [np.floor((lonlim[0] - border) / rounding) * rounding, np.ceil((lonlim[1] + border) / rounding) * rounding]

        # First download needed .hgt files. Quality is either 1, 3 or 30. If possible the files are downloaded. Otherwise
        # we fall back to lower quality. This is done using the elevation package
        tiles, q_tiles, tiles_30 = self._download_dem_files(latlim, lonlim, quality, data_folder, password=password, username=username)

        # Then create the final grid. This depends on the needed data type and possible resampling...
        if quality == 'SRTM1':
            pixel_degree = 3600
        elif quality == 'SRTM3':
            pixel_degree = 1200
        elif quality == 'SRTM30':
            pixel_degree = 120
        else:
            print('quality should be either SRTM1, SRTM3 or SRTM30!')
            return

        lat_size = int((latlim[1] - latlim[0]) * pixel_degree) + 1
        lon_size = int((lonlim[1] - lonlim[0]) * pixel_degree) + 1
        print('Bounding box is:')
        print('from ' + str(latlim[0]) + ' latitude to ' + str(latlim[1]))
        print('from ' + str(lonlim[0]) + ' longitude to ' + str(lonlim[1]))

        # Create final grid and add tiles.
        if quality == 'SRTM1' or quality == 'SRTM3':
            # Create file for DEM data
            bin_data = np.memmap(out_file, dtype=np.int16, shape=(lat_size, lon_size), mode='w+')
            bin_data = self._add_tiles(bin_data, tiles, quality, latlim, lonlim)

            # Create file for quality data
            bin_q_data = np.memmap(out_file + '.q', dtype=np.uint8, shape=(lat_size, lon_size), mode='w+')
            bin_q_data = self._add_tiles(bin_q_data, q_tiles, quality, latlim, lonlim, quality_file=True)

            # Save quality file to geotiff
            temp_q = os.path.join(data_folder, 'temp_q.tiff')
            n_latlim = [latlim[0] - 0.5 / pixel_degree, latlim[1] + 0.5 / pixel_degree]
            n_lonlim = [lonlim[0] - 0.5 / pixel_degree, lonlim[1] + 0.5 / pixel_degree]
            q_tiff = self._create_georeference(n_latlim, n_lonlim, 1.0 / pixel_degree, 1.0 / pixel_degree, 'uint8', temp_q)
            q_tiff.GetRasterBand(1).WriteArray(bin_q_data)
            q_tiff = None

        else:  # We are creating a 30 seconds map
            # Create file for DEM data
            bin_data = np.memmap(out_file, dtype=np.uint16, shape=(lat_size, lon_size))
            bin_data = self._add_tiles(bin_data, quality, tiles_30, latlim, lonlim)

            n_latlim = latlim
            n_lonlim = lonlim

        # Save the binary data as a geotiff
        print('Save data to geotiff')
        dem_tiff = os.path.join(data_folder, 'temp_dem.tiff')
        dem_data = self._create_georeference(n_latlim, n_lonlim, 1.0 / pixel_degree, 1.0 / pixel_degree, 'int16', dem_tiff)
        dem_data.GetRasterBand(1).WriteArray(bin_data)
        dem_data = None

        # First remove the egm96 from this raster.
        # Pre assign the egm tiff file first.
        print('Calculate geoid correction for SRTM data')
        egm_tiff = os.path.join(data_folder, 'egm96_resample.tiff')
        egm_data = self._create_georeference(n_latlim, n_lonlim, 1.0 / pixel_degree, 1.0 / pixel_degree, 'float32', egm_tiff)
        egm_data = None

        print('Correct DEM for geoid')
        self._add_egm96(dem_tiff, egm_tiff, data_folder)

        if resample == 'regular_grid':
            # If regular grid is used, we convert using gdal.
            # First create a geotiff file, then resample geotiff file. We always use cubic interpolation.
            print('Resampling to new regular grid')
            dlat = lats[1] - lats[0]
            dlon = lons[1] - lons[0]
            dem_tiff_final = os.path.join(data_folder, 'dem.tiff')
            dem_data_final = self._create_georeference(latlim, lonlim, dlat, dlon, 'float32', dem_tiff_final)
            dem_data_final = None
            dem_data_final = gdal.Open(dem_tiff_final, gdal.GA_Update)
            dem_data = gdal.Open(dem_tiff, gdal.GA_Update)

            gdal.ReprojectImage(dem_data, dem_data_final, None, None, gdal.GRA_Cubic)
            dem_data_final = None
            dem_tiff = dem_tiff_final

        elif resample == 'irregular_grid':
            # Use a simple bilinear approach to find values for specific points.
            print('Resampling to new irregular grid')
            heights = self._simple_bilinear(lats, lons, dem_tiff)

            return heights

        if doris_input == True:
            # Create a binary output file
            command = 'gdal_translate -of MFF ' + dem_tiff + ' ' + dem_tiff[:-5] + '.raw'
            os.system(command)

            if not os.path.exists(os.path.dirname(out_file)):
                os.makedirs(os.path.dirname(out_file))
            shutil.move(dem_tiff[:-5] + '.r00', out_file)

            # And create the scripts that go with them.
            self._output_doris_inputfiles(dem_tiff, out_file, var_file)

            # Finally also move .tiff file
            shutil.move(dem_tiff, out_file[:-4] + '.tiff')

        return


    def _add_tiles(self, outputdata, tiles, quality, latlim, lonlim, quality_file=False):
        # This function adds tiles to np.memmap file

        if quality == 'SRTM1':
            shape = (3601, 3601)
            s_size = 1.0 / 3600.0
            step_lat = 1
            step_lon = 1
        elif quality == 'SRTM3':
            shape = (1201, 1201)
            s_size = 1.0 / 1200.0
            step_lat = 1
            step_lon = 1
        elif quality == 'SRTM30':
            shape = (6000, 4800)
            s_size = 1.0 / 120.0
            step_lat = 50.0 - s_size
            step_lon = 40.0 - s_size
        else:
            print('quality should be either SRTM1, SRTM3 or SRTM30!')
            return

        print('total file size is ' + str(outputdata.shape[0]) + ' in latitude and ' + str(outputdata.shape[1]) + ' in longitude')

        for tile in tiles:
            if quality_file:
                image = np.fromfile(tile, dtype='>u1').reshape(shape)
            else:
                image = np.fromfile(tile, dtype='>i2').reshape(shape)

            if os.path.basename(tile)[7] == 'N':
                lat = float(os.path.basename(tile)[8:10])
            else:
                lat = - float(os.path.basename(tile)[8:10])
            if os.path.basename(tile)[10] == 'E':
                lon = float(os.path.basename(tile)[11:14])
            else:
                lon = - float(os.path.basename(tile)[11:14])
            if quality == 'SRTM30':
                lat = lat - 50 + (s_size / 2)
                lon += (s_size / 2)

            print('adding ' + tile)

            # Find the coordinates of the part of the tile that should be written to the output data.
            t_latlim = [max(lat, latlim[0]), min(lat + step_lat, latlim[1])]
            t_lonlim = [max(lon, lonlim[0]), min(lon + step_lon, lonlim[1])]
            t_latid = [shape[0] - int(round((t_latlim[0] - lat) / s_size)), shape[0] - (int(round((t_latlim[1] - lat) / s_size)) + 1)]
            t_lonid = [int(round((t_lonlim[0] - lon) / s_size)), int(round((t_lonlim[1] - lon) / s_size)) + 1]
            latsize = int(round((latlim[1] - latlim[0]) / s_size)) + 1
            latid = [latsize - int(round((t_latlim[0] - latlim[0]) / s_size)), latsize - (int(round((t_latlim[1] - latlim[0]) / s_size)) + 1)]
            lonid = [int(round((t_lonlim[0] - lonlim[0]) / s_size)), int(round((t_lonlim[1] - lonlim[0]) / s_size)) + 1]

            print('Adding tile lat ' + str(t_latid[1] + 1) + ' to ' + str(t_latid[0]) + ' into dem file ' +
                  str(latid[1] + 1) + ' to ' + str(latid[0]))
            print('Adding tile lon ' + str(t_lonid[0] + 1) + ' to ' + str(t_lonid[1]) + ' into dem file ' +
                  str(lonid[0] + 1) + ' to ' + str(lonid[1]))

            # Assign values from tiles to outputdata
            if quality == 'SRTM30':
                outputdata[latid[1]: latid[0]-2, lonid[0]: lonid[1]-2] = image[t_latid[1]: t_latid[0]-2, t_lonid[0]: t_lonid[1]-2]
            else:
                outputdata[latid[1]: latid[0]-1, lonid[0]: lonid[1]-1] = image[t_latid[1]: t_latid[0]-1, t_lonid[0]: t_lonid[1]-1]

        return outputdata


    def _download_dem_files(self, latlim, lonlim, quality, data_folder, username='', password='', password_xml=''):
        # This function downloads data either in 1,3 or 30 arc seconds. When you choose either 1 or 3 seconds also 30
        # seconds is downloaded to fill the voids.
        # In this operation the needed files are also extracted...

        # Check the username and password and load from config file if needed.
        if not username or not password:
            if os.path.exists(password_xml):
                tree = ET.parse(password_xml)
                settings = tree.getroot()

                username = settings.find('.usgs_username').text
                password = settings.find('.usgs_password').text
            else:
                print('You should specify a username and password to download SRTM data. ')
                return

        filelist = self._srtm_listing(data_folder, username=username, password=password)
        outfiles = []
        q_files = []
        outfiles_30 = []

        lats = np.arange(np.floor(latlim[0]), np.ceil(latlim[1]))
        lons = np.arange(np.floor(lonlim[0]), np.ceil(lonlim[1]))

        if quality == 'SRTM1' or quality == 'SRTM3':
            for lat in lats:
                for lon in lons:

                    lat = int(lat)
                    lon = int(lon)

                    if lat < 0:
                        latstr = 'S' + str(abs(lat)).zfill(2)
                    else:
                        latstr = 'N' + str(lat).zfill(2)
                    if lon < 0:
                        lonstr = 'W' + str(abs(lon)).zfill(3)
                    else:
                        lonstr = 'E' + str(lon).zfill(3)

                    # Check if file exists in filelist
                    if str(lat) not in filelist[quality]:
                        continue
                    elif str(lon) not in filelist[quality][str(lat)]:
                        continue

                    filename = os.path.join(data_folder, latstr + lonstr + 'SRTMGL3.hgt.zip')
                    q_file = os.path.join(data_folder, latstr + lonstr + 'SRTMGL3.q.zip')
                    extracted_file = os.path.join(data_folder, quality + '__' + latstr + lonstr + '.hgt')
                    q_extracted = os.path.join(data_folder, quality + '__' + latstr + lonstr + '.q')

                    if not os.path.exists(extracted_file) or not os.path.exists(q_extracted):
                        # Download and unzip
                        download_dem = filelist[quality][str(lat)][str(lon)]
                        download_q = download_dem[:-7] + 'num.zip'

                        command = 'wget ' + download_dem + ' --user ' + username + ' --password ' + password + ' -O ' + filename
                        q_command = 'wget ' + download_q + ' --user ' + username + ' --password ' + password + ' -O ' + q_file
                        try:
                            os.system(command)
                            zip_data = zipfile.ZipFile(filename)
                            source = zip_data.open(zip_data.namelist()[0])
                            target = open(extracted_file, 'wb')
                            shutil.copyfileobj(source, target, length=-1)
                            target.close()
                            outfiles.append(extracted_file)
                            os.remove(filename)

                            os.system(q_command)
                            zip_data = zipfile.ZipFile(q_file)
                            source = zip_data.open(zip_data.namelist()[0])
                            target = open(q_extracted, 'wb')
                            shutil.copyfileobj(source, target, length=-1)
                            target.close()
                            q_files.append(q_extracted)
                            os.remove(q_file)
                        except:
                            print('Failed to download or process ' + filename)

                    else:
                        outfiles.append(extracted_file)
                        q_files.append(q_extracted)

        for lat in lats:
            for lon in lons:
                # Now check with which 30 seconds file this corresponds
                lat50 = int(np.floor(float(lat + 10) / 50)) * 50 + 40
                lon40 = int(np.floor(float(lon + 20) / 40)) * 40 - 20

                if lat50 < 0:
                    latstr = 'S' + str(abs(lat50)).zfill(2)
                else:
                    latstr = 'N' + str(lat50).zfill(2)
                if lon40 < 0:
                    lonstr = 'W' + str(abs(lon40)).zfill(3)
                else:
                    lonstr = 'E' + str(lon40).zfill(3)

                # Does this file exist..
                if str(lat50) not in filelist['SRTM30']:
                    continue
                elif str(lon40) not in filelist['SRTM30'][str(lat50)]:
                    continue

                filename = os.path.join(data_folder, latstr + lonstr + 'SRTMGL3.hgt.zip')
                extracted_file = os.path.join(data_folder, 'SRTM30_' + latstr + lonstr + '.hgt')

                # Download and unzip file if possible and add to list
                if not os.path.exists(extracted_file):
                    # Download and unzip
                    user = 'USERNAME'
                    password = 'PASSWORD'

                    command = 'wget ' + filelist['SRTM30'][str(lat50)][str(lon40)] + ' --user ' + user + ' --password ' + password + \
                              ' -O ' + filename
                    try:
                        os.system(command)
                        zip_data = zipfile.ZipFile(filename)
                        source = zip_data.open(zip_data.namelist()[0])
                        target = open(extracted_file, 'wb')
                        shutil.copyfileobj(source, target, length=-1)
                        target.close()
                        outfiles_30.append(extracted_file)

                    except:
                        print('Failed to download or process ' + filename)
                elif extracted_file not in outfiles_30:
                    outfiles_30.append(extracted_file)

        return outfiles, q_files, outfiles_30


    def _simple_bilinear(self, lats, lons, dem_tiff, data_folder):
        # To resample to datapoints, we convert the tiff file to a binary file
        bin_file = os.path.join(data_folder, 'dem.raw')
        command = 'gdal_translate -of MFF ' + dem_tiff + ' ' + bin_file
        os.system(command)
        shutil.move(bin_file[:-4] + '.r00', bin_file)

        # Read it as a memmap function
        dem_data = gdal.Open(dem_tiff, gdal.GA_Update)
        size = (dem_data.RasterYSize, dem_data.RasterXSize)
        data = np.memmap(bin_file, shape=size, dtype=np.dtype('float32'))
        r = dem_data.GetGeoTransform()

        # And find the values of corresponding points using bilinear interpolation.
        x_id = np.floor((lons - r[0]) / r[1]).astype('int32')
        x_diff = (((lons - r[0]) / r[1]) - x_id)
        y_id = np.floor((lats - r[3]) / r[5]).astype('int32')
        y_diff = (((lats - r[3]) / r[5]) - y_id)

        # Calculate final results
        ll_cont = data[y_id, x_id] * (1-x_diff) * (1-y_diff)
        ul_cont = data[y_id + 1, x_id] * (1-x_diff) * y_diff
        ur_cont = data[y_id + 1, x_id + 1] * x_diff * y_diff
        lr_cont = data[y_id, x_id + 1] * x_diff * (1-y_diff)
        heights = ll_cont + ul_cont + ur_cont + lr_cont

        return heights


    def _kml_shp_2_bb(self, filename):
        # from kml and shape file to a bounding box. We will always use a bounding box to create the final product.

        if filename.endswith('.shp'):
            with collection(filename, "r") as inputshape:

                shapes = [shape for shape in inputshape]
                # only first shape
                dat = shapes[0]['geometry']['coordinates']

                lon = [l[0] for l in dat[0]]
                lat = [l[1] for l in dat[0]]

                latlim = [min(lat), max(lat)]
                lonlim = [min(lon), max(lon)]

        elif filename.endswith('.kml'):
            doc = file(filename).read()
            k = fastkml.KML()
            k.from_string(doc)
            dat = list(list(k.features())[0].features())[0].geometry[0].exterior.coords[:]

            lon = [l[0] for l in dat[0]]
            lat = [l[1] for l in dat[0]]

            latlim = [min(lat), max(lat)]
            lonlim = [min(lon), max(lon)]
        else:
            print('format not recognized! Pleas creat either a .kml or .shp file.')
            return []

        return latlim, lonlim


    def _add_egm96(self, dem_tiff, egm_tiff, data_folder):
        # This function adds the geoid from egm96

        # Load egm96 grid and resample to input grid using gdal.
        # (For this purpose the grid is downloaded from:
        # http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/binary/binarygeoid.html
        # In principle this is converted to geotiff here,

        filename = os.path.join(data_folder, 'EGM96_15min.dat')
        egm_source_tiff = os.path.join(data_folder, 'EGM96_15min.tiff')

        if not os.path.exists(egm_source_tiff):
            if not os.path.exists(filename):
                # Download egm96 file
                command = 'wget http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/binary/WW15MGH.DAC -O ' + filename
                os.system(command)

            # Get georeference
            latlim = [-90.125, 90.125]
            lonlim = [-0.125, 359.875]
            dlat = 0.25
            dlon = 0.25
            egm_data = self._create_georeference(latlim, lonlim, dlat, dlon, 'float32', egm_source_tiff)

            # Load data
            egm96 = np.fromfile(filename, dtype='>i2').reshape((721, 1440)).astype('float32')

            # Save as geotiff
            egm_data.GetRasterBand(1).WriteArray(egm96 / 100)
            egm_data = None

        egm_source = gdal.Open(egm_source_tiff, gdal.GA_Update)
        egm_data = gdal.Open(egm_tiff, gdal.GA_Update)

        # Now resample to new dataset. (Copy information from gdal_dem)
        gdal.ReprojectImage(egm_source, egm_data, None, None, gdalconst.GRA_Bilinear)
        egm_data = None

        dem_new = dem_tiff + '.new'

        # Finally open original dataset and subtract
        command = 'gdal_calc.py -A ' + dem_tiff + ' -B ' + egm_tiff + ' --outfile=' + dem_new + ' --calc="A+B"'
        os.system(command)

        shutil.move(dem_new, dem_tiff)


    def _output_doris_inputfiles(self, dem_tiff, out_file, var_file):
        # This script outputs doris inputfiles which can be copied for further processing.

        dem = gdal.Open(dem_tiff)

        xsize = dem.RasterXSize
        ysize = dem.RasterYSize
        georef = dem.GetGeoTransform()
        dlat = georef[1]
        dlon = abs(georef[5])
        lonmin = georef[0] + (dlon * 0.5)
        latmax = georef[3] - (dlat * 0.5)

        output_txt = out_file + '.doris_inputfile'
        output_var = var_file
        output_var = open(output_var, 'w')
        txtfile = open(output_txt, 'w')

        dem_var = dict()
        dem_var['in_dem'] = out_file
        dem_var['in_format'] = 'r4'
        dem_var['in_size'] = str(ysize) + " " + str(xsize)
        dem_var['in_delta'] = str(dlat) + " " + str(dlon)
        dem_var['in_ul'] = str(latmax) + " " + str(lonmin)
        dem_var['in_nodata'] = '-32768'
        pickle.dump(dem_var, output_var)
        output_var.close()

        txtfile.write("# The processing cards generated by $(basename $0) script. \n")
        txtfile.write("# Using parameters: $@ \n")
        txtfile.write('# Copy the section(s) that is/are necessary to your processing setup. \n')
        txtfile.write("c         ___             ___ \n")
        txtfile.write("comment   ___SIM AMPLITUDE___ \n")
        txtfile.write("c                             \n")
        txtfile.write("SAM_IN_DEM     " + out_file + " \n")
        txtfile.write("SAM_IN_FORMAT   r4 \t\t\t // default is short integer \n")
        txtfile.write("SAM_IN_SIZE    " + str(ysize) + " " + str(xsize) + " \n")
        txtfile.write("SAM_IN_DELTA   " + str(dlat) + " " + str(dlon) + " \n")
        txtfile.write("SAM_IN_UL      " + str(latmax) + " " + str(lonmin) + " \n")
        txtfile.write("SAM_IN_NODATA  -32768 \n")
        txtfile.write("  \n")
        txtfile.write("  \n")
        txtfile.write("c         ___          ___ \n")
        txtfile.write("comment   ___DEM ASSIST___ \n")
        txtfile.write("c                             \n")
        txtfile.write("DAC_IN_DEM     $dempath/$outfile5 \n")
        txtfile.write("DAC_IN_FORMAT   r4 \t\t\t // default is short integer \n")
        txtfile.write("DAC_IN_SIZE    " + str(ysize) + " " + str(xsize) + " \n")
        txtfile.write("DAC_IN_DELTA   " + str(dlat) + " " + str(dlon) + " \n")
        txtfile.write("DAC_IN_UL      " + str(latmax) + " " + str(lonmin) + " \n")
        txtfile.write("DAC_IN_NODATA  -32768 \n")
        txtfile.write("  \n")
        txtfile.write("  \n")
        txtfile.write("c         ___             ___ \n")
        txtfile.write("comment   ___REFERENCE DEM___ \n")
        txtfile.write("c                             \n")
        txtfile.write("## CRD_METHOD   DEMINTRPMETHOD \n")
        txtfile.write("CRD_IN_DEM     $dempath/$outfile5 \n")
        txtfile.write("CRD_IN_FORMAT   r4 \t\t\t // default is short integer \n")
        txtfile.write("CRD_IN_SIZE    " + str(ysize) + " " + str(xsize) + " \n")
        txtfile.write("CRD_IN_DELTA   " + str(dlat) + " " + str(dlon) + " \n")
        txtfile.write("CRD_IN_UL      " + str(latmax) + " " + str(lonmin) + " \n")
        txtfile.write("CRD_IN_NODATA  -32768 \n")

        txtfile.close()


    def _srtm_listing(self, data_folder, username, password):
        # This script makes a list of all the available 1,3 and 30 arc second datafiles.
        # This makes it easier to detect whether files do or don't exist.

        data_file = os.path.join(data_folder, 'filelist')
        if os.path.exists(data_file):
            dat = open(data_file, 'r')
            filelist = pickle.load(dat)
            dat.close()
            return filelist

        server = "http://e4ftl01.cr.usgs.gov"

        folders = 'SRTM/SRTMGL1.003/2000.02.11/', 'SRTM/SRTMGL3.003/2000.02.11/', 'SRTM/SRTMGL30.002/2000.02.11/'
        keys = ['SRTM1', 'SRTM3', 'SRTM30']
        filelist = dict()
        filelist['SRTM1'] = dict()
        filelist['SRTM3'] = dict()
        filelist['SRTM30'] = dict()

        for folder, key_value in zip(folders, keys):

            conn = requests.get(server + '/' + folder, auth=(username, password))
            if conn.status_code == 200:
                print("status200 received ok")
            else:
                print("an error occurred during connection")

            data = conn.text
            parser = parseHTMLDirectoryListing()
            parser.feed(data)
            files = parser.getDirListing()

            if key_value == 'SRTM1' or key_value == 'SRTM3':
                files = [f for f in files if f.endswith('hgt.zip')]
                north = [int(filename[1:3]) for filename in files]
                east = [int(filename[4:7]) for filename in files]
                for i in [i for i, filename in enumerate(files) if filename[0] == 'S']:
                    north[i] = north[i] * -1
                for i in [i for i, filename in enumerate(files) if filename[3] == 'W']:
                    east[i] = east[i] * -1
            else:
                files = [f for f in files if f.endswith('dem.zip')]
                north = [int(filename[5:7]) for filename in files]
                east = [int(filename[1:4]) for filename in files]
                for i in [i for i, filename in enumerate(files) if filename[4] == 's']:
                    north[i] = north[i] * -1
                for i in [i for i, filename in enumerate(files) if filename[0] == 'w']:
                    east[i] = east[i] * -1

            for filename, n, e in zip(files, north, east):
                if not str(n) in filelist[key_value]:
                    filelist[key_value][str(n)] = dict()
                filelist[key_value][str(n)][str(e)] = server + '/' + folder + filename

        file_list = open(os.path.join(data_folder, 'filelist'), 'w')
        pickle.dump(filelist, file_list)
        file_list.close()

        return filelist


    def _create_georeference(self, latlim, lonlim, dlat, dlon, dtype='int16', filename=''):
        # This function also creates a dem file but creates an geotiff file instead. This functionality is used to allow
        # resampling in later steps using gdal.

        conversion = {
            "uint8": 1,
            "int8": 1,
            "uint16": 2,
            "int16": 3,
            "uint32": 4,
            "int32": 5,
            "float32": 6,
            "float64": 7,
            "complex64": 10,
            "complex128": 11,
        }

        if filename:
            driver = gdal.GetDriverByName('Gtiff')
            dataset = driver.Create(filename,
                                    int(np.round((lonlim[1] - lonlim[0]) / dlon)),
                                    int(np.round((latlim[1] - latlim[0]) / dlat)),
                                    1,
                                    conversion[dtype], ['COMPRESS=LZW', 'BIGTIFF=YES'])
        else:
            driver = gdal.GetDriverByName('mem')
            dataset = driver.Create('',
                                    int((lonlim[1] - lonlim[0]) / dlon),
                                    int((latlim[1] - latlim[0]) / dlat),
                                    1,
                                    conversion[dtype])

        dataset.SetGeoTransform((
            lonlim[0],  # 0
            dlat,  # 1
            0,  # 2
            latlim[1],  # 3
            0,  # 4
            -dlon))

        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(4326)

        dataset.SetProjection(spatial_ref.ExportToWkt())

        return dataset

    # This function is not used because data does not contain gaps....
    def _fill_voids(self, outputdata, output_30sec, quality):
        # This function fills the voids in our output data.

        if quality == 'SRTM1':
            s_size = 1.0 / 3600.0
        elif quality == 'SRTM3':
            s_size = 1.0 / 1200.0
        else:
            print('quality should be either SRTM1 or SRTM3!')
            return

        # This script assumes that the final dem is at least divided in 1 arc minute blocks so the two input fields are
        # perfectly alligned.
        id_void = np.argwhere(outputdata == -32767)

        if id:  # If there are any voids...
            id_30 = np.floor(id_void * 120.0 * s_size)
            outputdata[id_void] = output_30sec[id_30]

        return outputdata

# Following code is adapted from srtm-1.py > downloaded from
# https://svn.openstreetmap.org/applications/utils/import/srtm2wayinfo/python/srtm.py
class parseHTMLDirectoryListing(HTMLParser):
    def __init__(self):
        # print("parseHTMLDirectoryListing.__init__")
        HTMLParser.__init__(self)
        self.title = "Undefined"
        self.isDirListing = False
        self.dirList = []
        self.inTitle = False
        self.inHyperLink = False
        self.currAttrs = ""
        self.currHref = ""

    def handle_starttag(self, tag, attrs):
        # print("Encountered the beginning of a %s tag" % tag)
        if tag == "title":
            self.inTitle = True
        if tag == "a":
            self.inHyperLink = True
            self.currAttrs = attrs
            for attr in attrs:
                if attr[0] == 'href':
                    self.currHref = attr[1]

    def handle_endtag(self, tag):
        # print("Encountered the end of a %s tag" % tag)
        if tag == "title":
            self.inTitle = False
        if tag == "a":
            # This is to avoid us adding the parent directory to the list.
            if self.currHref != "":
                self.dirList.append(self.currHref)
            self.currAttrs = ""
            self.currHref = ""
            self.inHyperLink = False

    def handle_data(self, data):
        if self.inTitle:
            self.title = data
            print("title=%s" % data)
            if "Index of" in self.title:
                # print("it is an index!!!!")
                self.isDirListing = True
        if self.inHyperLink:
            # We do not include parent directory in listing.
            if "Parent Directory" in data:
                self.currHref = ""

    def getDirListing(self):
        return self.dirList

# Actually execute the code...
if __name__ == "__main__":

    stack_folder = sys.argv[1]
    if len(sys.argv) > 2:
        quality = sys.argv[2]
        if quality not in ['SRTM1', 'SRTM3']:
            quality = 'SRTM3'
    else:
        quality = 'SRTM3'

    xml_file = os.path.join(os.path.join(stack_folder, 'doris_input.xml'))
    print('reading xml file stack ' + xml_file)
    tree_stack = ET.parse(xml_file)
    settings_stack = tree_stack.getroot()[0]

    xml_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'install',
                                 'doris_config.xml')
    print('reading xml file settings doris ' + xml_file)
    tree_doris = ET.parse(xml_file)
    settings_doris = tree_doris.getroot()

    # Then create the dem file
    shapefile = settings_stack.find('.shape_file_path').text
    dem_calc_folder = settings_stack.find('.dem_processing_folder').text
    dem_out_folder = settings_stack.find('.dem_folder').text

    dem_out = os.path.join(dem_out_folder, 'dem.raw')
    dem_var = os.path.join(dem_out_folder, 'var.raw')

    srtm_username = settings_doris.find('.usgs_username').text
    srtm_password = settings_doris.find('.usgs_password').text

    dem = CreateDem()
    dem.create(shapefile, dem_out, dem_var, resample=None,
            doris_input=True, rounding=1, border=1,
            data_folder=dem_calc_folder, quality=quality,
            password=srtm_password, username=srtm_username)
