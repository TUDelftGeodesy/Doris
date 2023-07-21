import os, sys
import zipfile
import shutil
import warnings
import xml.etree.cElementTree as etree
from shapely.geometry import Polygon, shape
import fiona
import numpy as np


def unzip_folder(zipped_folder, dest_folder, shapefile='', pol='', data=True, swath='all', overwrite=False, check_valid=False):
    # This function unzips a data folder. The total amount of data to unpack can be reduced by selecting which data
    # files are unpacked. In first instance all files except the data files are extracted. Then the other files are
    # extracted based on the given input:
    # - pol > only this polarisation is unpacked ('vv', 'hh', 'hv', 'vh')
    # - shape > only the swaths which overlap with this shape are unpacked shape is a shapefile file
    # - swath > you can choose either 1/2/3
    # If you do not want to extract any data file use data=False. Finally, the script will skip files which are already
    # unpacked unless you set overwrite to True

    zipdat = zipfile.ZipFile(zipped_folder)

    # First check zipfile quality
    if check_valid == True:
        test = zipdat.testzip()
        if test:
            print('Some files in ' + zipped_folder + ' are corrupted.')
            return

    # First check wether the shape overlaps...
    kml_file, png_file = extract_kml_preview(zipped_folder, dir=dest_folder, overwrite=overwrite)

    if shapefile:
        shp = load_shape(shapefile, buffer=0.02)
        overlap = shape_im_kml(shp, kml_file)
        if not overlap:
            print('The image and kml_file do not overlap')
            return

    swaths = []
    for filename in zipdat.namelist():
        absname = os.path.join(dest_folder, filename)

        if not os.path.exists(absname) or overwrite == True:
            if not filename.endswith('.tiff'):
                zipdat.extract(filename, dest_folder)
            else:
                swaths.append(filename)

    if swath in ['1', '2', '3']:  # If only one of the swaths is extracted.
        swaths = [s for s in swaths if os.path.basename(s)[6] == str(swath)]
    if pol in ['vv','vh','hh','hv']:
        swaths = [s for s in swaths if os.path.basename(s)[12:14] == pol]
    if not swaths or data == False: # If there is nothing left, stop unpacking.
        return
    if shapefile:
        d_swath = []
        for s in swaths:
            xml_file = os.path.join(dest_folder, os.path.dirname(os.path.dirname(s)), 'annotation', os.path.basename(s)[:-4] + 'xml')
            if shape_swath_xml(shp, xml_file):
                d_swath.append(s)
        swaths = d_swath

    # Finally unpack the needed swaths
    for s in swaths:
        zipdat.extract(s, dest_folder)


def extract_kml_preview(zipped_folder, dir='', kml=True, png=True, overwrite=False):
    # Extracts quicklook and/or .kml files.
    #print(zipped_folder)
    zipdat = zipfile.ZipFile(zipped_folder)
    if not dir:
        dir = os.path.dirname(zipped_folder)

    png_name = ''
    kml_name = ''

    for filename in zipdat.namelist():  # Unzip and save .kml file
        if filename.endswith('map-overlay.kml') and kml == True:
            kml_name = os.path.join(dir, os.path.basename(zipped_folder)[:-4] + '.kml')
            zipped_kml = zipdat.open(filename)
            if not os.path.exists(kml_name) or overwrite == True:
                kml_file = file(kml_name, "wb")
                with zipped_kml, kml_file:
                    shutil.copyfileobj(zipped_kml, kml_file)

    for filename in zipdat.namelist():  # Unzip and save quicklook
        if filename.endswith('quick-look.png') and png == True:
            png_name = os.path.join(dir, os.path.basename(zipped_folder)[:-4] + '.png')
            zipped_png = zipdat.open(filename)
            if not os.path.exists(png_name) or overwrite == True:
                png_file = file(png_name, "wb")
                with zipped_png, png_file:
                    shutil.copyfileobj(zipped_png, png_file)

    return kml_name, png_name


def shape_im_kml(shp, kml_file):
    # This script extracts a Fiona/polygon shape of the footprint given in the .xml file of the image and checks whether
    # it overlaps

    # First check is .kml file exist
    if not os.path.exists(kml_file):
        warnings.warn('.kml file does not exist.')
        return False

    try:
        in_kml = etree.parse(kml_file)
        in_kml = in_kml.getroot()
        coor = in_kml[0][1][1][2][0].text
        coor = [i.split(',') for i in coor.split(' ')]
        coverage = Polygon([[float(i[0]),float(i[1])] for i in coor])
    except:
        warnings.warn('.kml file is corrupt')
        return False

    if coverage.intersects(shp):
        return True
    else:
        return False


def shape_swath_xml(shp, xml_file):
    # This script extracts a Fiona/polygon shape of the footprint given in the .xml file of the image and checks whether
    # it overlaps

    # First check is .xml file exist
    if not os.path.exists(xml_file):
        warnings.warn('.xml file does not exist.')
        return False

    try:
        in_xml = etree.parse(xml_file)
        in_xml = in_xml.getroot()
        coor = in_xml.find('.geolocationGrid').find('.geolocationGridPointList').findall('.geolocationGridPoint')
        lats = []
        lons = []
        line = []
        pixel = []
        for c in coor:
            lats.append(float(c.find('.latitude').text))
            lons.append(float(c.find('.longitude').text))
            pixel.append(int(c.find('.pixel').text))
            line.append(int(c.find('.line').text))
        maxpixel = np.max(pixel)
        coor = [[lat, lon] for lon, lat, p in zip(lats, lons, pixel) if p == 0]
        coor.extend([[lat, lon] for lon, lat, p in zip(lats[::-1], lons[::-1], pixel[::-1]) if p == maxpixel])
        coverage = Polygon(coor)

    except:
        warnings.warn('.xml file is corrupt')
        return False

    if coverage.intersects(shp):
        return True
    else:
        return False


def load_shape(shapefile, buffer=0.02):
    # This function creates a shape to make a selection of usable bursts later on. Buffer around shape is in
    # degrees.

    if not shapefile:
        warnings.warn('Please provide a shapefile or coordinates.')

    try:
        if isinstance(shapefile, list):  # If the coordinates are already loaded. (for example bounding box)
            shp = Polygon(shapefile)
        else:  # It should be a shape file. We always select the first shape.
            sh = next(iter(fiona.open(shapefile)))#fiona.open(shapefile).next()
            shp = shape(sh['geometry'])

        # Now we have the shape we add a buffer and simplify first to save computation time.
        shp = shp.simplify(buffer / 2)
        shp = shp.buffer(buffer)
    except:
        warnings.warn('Unrecognized shape')
        return

    return shp


# Testing ---------------------------------------------
# zipped_folder = '/media/gert/Data/radar_database/sentinel-1/s1_asc_t88/IW_SLC__1SDV_VVVH/20141116/S1A_IW_SLC__1SDV_20141116T172443_20141116T172510_003310_003D5C_E92F.SAFE.zip'
# dest_folder = '/media/gert/Data/radar_database/sentinel-1/s1_asc_t88/IW_SLC__1SDV_VVVH/20141116/test'
# shapefile = '/media/gert/Data/shapes/netherlands/netherland.shp'
# pol = ''
# data = True
# swath = 'all'
# overwrite = False
# check_valid = True
# ------------------------------------------------------

# Actually execute the code to unzip one data file.
if __name__ == "__main__":

    data = True
    check_valid = False
    swath = 'all'

    zipped_folder = sys.argv[1]
    dest_folder = sys.argv[2]
    shapefile = sys.argv[3]
    pol = sys.argv[4]
    overwrite = sys.argv[5]

    if overwrite == 'False':
        overwrite = False
    else:
        overwrite = True

    print('zipfile is ' + zipped_folder)
    print('shapefile is ' + shapefile)
    print('destination folder is ' + dest_folder)

    unzip_folder(zipped_folder, shapefile=shapefile, pol=pol, dest_folder=dest_folder, overwrite=overwrite, swath=swath, check_valid=check_valid, data=data)

