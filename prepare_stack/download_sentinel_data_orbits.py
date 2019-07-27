# This file contains a function to check which files for sentinel are available, which ones are downloaded and a quality
# check for the files which are downloaded.

import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import ssl
import re
import os, sys
import datetime
import base64
import subprocess
from fiona import collection
from fastkml import kml
from lxml import etree
import xml.etree.ElementTree as ET


def sentinel_available(start_day='', end_day='', sensor_mode='', product='', level='', track='', polarisation='', orbit_direction='', ROI='', user='', password=''):
    # All available sentinel 1 images are detected and printed on screen.
    # Following variables can be used to make a selection.
    # shape > defining shape file or .kml
    # start_day > first day for downloads (default one month before now) [yyyymmdd]
    # end_day > last day for downloads (default today)
    # track > the tracks we want to check (default all)
    # polarisation > which polarisation will be used. (default all)

    # string is the field we enter as url
    string = ''

    if sensor_mode:
        string = string + ' AND ' + 'sensoroperationalmode:' + sensor_mode
    if product:
        string = string + ' AND ' + 'producttype:' + product
    if level:
        string = string + ' AND ' + level
    if orbit_direction:
        string = string + ' AND ' + 'orbitdirection:' + orbit_direction
    if track:
        string = string + ' AND ' + 'relativeorbitnumber:' + track
    if start_day:
        start = datetime.datetime.strptime(start_day, '%Y-%m-%d').strftime('%Y-%m-%d')
    else:
        start = (datetime.datetime.now() - datetime.timedelta(days=350)).strftime('%Y-%m-%d')
    if end_day:
        end = datetime.datetime.strptime(end_day, '%Y-%m-%d').strftime('%Y-%m-%d')
    else:
        end = datetime.datetime.now().strftime('%Y-%m-%d')
    if polarisation:
        string = string + ' AND ' + 'polarisationmode:' + polarisation
    if ROI:
        shape_str = load_shape_info(ROI)
        string = string + ' AND footprint:"Intersects(POLYGON(' + shape_str + '))"'

    date_string = 'beginPosition:[' + start + 'T00:00:00.000Z TO ' + end + 'T23:59:59.999Z] AND endPosition:[' + start + 'T00:00:00.000Z TO ' + end + 'T23:59:59.999Z]'
    string = string + ' AND ' + date_string

    # Finally we do the query to get the search result.
    string = string[5:] + '&rows=1000'
    url = 'https://scihub.copernicus.eu/dhus/search?q=' + urllib.parse.quote_plus(string)
    print(url)

    print('Requesting available products: ' + url)
    request = urllib.request.Request(url)
    base64string = base64.b64encode('%s:%s' % (user, password))
    request.add_header("Authorization", "Basic %s" % base64string)

    # connect to server. Hopefully this works at once
    try:
        dat = urllib.request.urlopen(request)
    except:
        print('not possible to connect this time')
        return [], [], []

    html_dat = ''
    for line in dat:
        html_dat = html_dat + line

    parser = etree.HTMLParser()
    tree = etree.fromstring(html_dat, parser)
    products = [data for data in tree.iter(tag='entry')]
    links = [data.find('link').attrib for data in tree.iter(tag='entry')]
    dates = [data.findall('date')[1].text for data in tree.iter(tag='entry')]

    print('Following products will be downloaded')
    for link in links:
        print(link)


    return products, links, dates


def load_shape_info(shapefile):
    # This script converts .kml, .shp and .txt files to the right format. If multiple shapes are available the script
    # will select the first one.

    if shapefile.endswith('.shp'):
        with collection(shapefile, "r") as inputshape:
            for shape in inputshape:
                # only first shape
                dat = shape['geometry']['coordinates']

                st='('
                for p in dat[0]:
                    st = st + str(p[0]) + ' ' + str(p[1]) + ','
                st = st[:-1] + ')'

                break
    elif shapefile.endswith('.kml'):
        doc = file(shapefile).read()
        k = kml.KML()
        k.from_string(doc)
        shape = list(list(k.features())[0].features())[0].geometry.exterior.coords[:]
        st='('
        for p in shape:
            st = st + str(p[0]) + ' ' + str(p[1]) + ','
        st = st[:-1] + ')'
    else:
        print('format not recognized! Pleas creat either a .kml or .shp file.')
        return []

    return st


def sentinel_check_validity(products=[], destination_folder='', user='', password='', remove=True):
    # Check if the downloaded files are valid and remove if not

    valid_files = []
    invalid_files = []

    if not products:
        print('Nothing to check')
        return

    for product in products:
        date = str(product.findall('date')[1].text)
        date = datetime.datetime.strptime(date[:19], '%Y-%m-%dT%H:%M:%S')

        name = str(product.find('title').text)

        track = str(product.find('int[@name="relativeorbitnumber"]').text)
        data_type = str(product.find(".//str[@name='filename']").text)[4:16]
        pol = str(product.find(".//str[@name='polarisationmode']").text).replace(' ', '')
        direction = str(product.find(".//str[@name='orbitdirection']").text)
        if direction == 'ASCENDING':
            direction = 'asc'
        elif direction == 'DESCENDING':
            direction = 'dsc'

        trackfolder = os.path.join(destination_folder, 's1_' + direction + '_t' + track)
        typefolder = os.path.join(trackfolder, data_type + '_' + pol)
        datefolder = os.path.join(typefolder, date.strftime('%Y%m%d'))

        xml_dir = os.path.join(datefolder, name + '.xml')
        file_dir = os.path.join(datefolder, name + '.SAFE.zip')
        kml_dir = os.path.join(datefolder, name + '.kml')
        preview_dir = os.path.join(datefolder, name + '.jpg')

        # First check the file
        if os.path.exists(file_dir):
            uuid = product.find('id').text
            valid_dat = sentinel_quality_check(file_dir, uuid, user, password)
        else:
            valid_dat = False

        if not valid_dat:
            if os.path.exists(file_dir) and remove == True:
                os.system('rm ' + file_dir)
            if os.path.exists(xml_dir) and remove == True:
                os.system('rm ' + xml_dir)
            if os.path.exists(kml_dir) and remove == True:
                os.system('rm ' + kml_dir)
            if os.path.exists(preview_dir) and remove == True:
                os.system('rm ' + preview_dir)

            invalid_files.append(file_dir)
        else:
            valid_files.append(file_dir)

    return invalid_files, valid_files


def sentinel_download(products=[], xml_only=False,  destination_folder='', project_folder='', user='', password=''):
    # Download the files which are found by the sentinel_available script.

    if not products:
        print('No files to download')
        return

    wget_base = 'wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 --continue --tries=20 --no-check-certificate --user=' + user + ' --password=' + password + ' '

    for product in products:
        date = str(product.findall('date')[1].text)
        date = datetime.datetime.strptime(date[:19], '%Y-%m-%dT%H:%M:%S')

        url = str('"'+product.findall('link')[0].attrib['href'][:-6]+ urllib.parse.quote_plus('$value') +'"')
        name = str(product.find('title').text)

        track = str(product.find('int[@name="relativeorbitnumber"]').text)
        data_type = str(product.find(".//str[@name='filename']").text)[4:16]
        pol = str(product.find(".//str[@name='polarisationmode']").text).replace(' ', '')
        direction = str(product.find(".//str[@name='orbitdirection']").text)
        if direction == 'ASCENDING':
            direction = 'asc'
        elif direction == 'DESCENDING':
            direction = 'dsc'

        trackfolder = os.path.join(destination_folder, direction + '_t' + track.zfill(3))
        if not os.path.exists(trackfolder):
            os.mkdir(trackfolder)
        typefolder = os.path.join(trackfolder, data_type + '_' + pol)
        if not os.path.exists(typefolder):
            os.mkdir(typefolder)
        datefolder = os.path.join(typefolder, date.strftime('%Y%m%d'))
        if not os.path.exists(datefolder):
            os.mkdir(datefolder)

        xml_dir = os.path.join(datefolder, name + '.xml')
        file_dir = os.path.join(datefolder, name + '.SAFE.zip')
        kml_dir = os.path.join(datefolder, name + '.kml')
        preview_dir = os.path.join(datefolder, name + '.jpg')

        if project_folder:
            datefolder = os.path.join(project_folder, 's1', date.strftime('%Y%m%d') + '_t' + track)
            if not os.path.exists(datefolder):
                os.mkdir(datefolder)
            sentinel_folder = os.path.join(datefolder, 'sentinel_1')
            if not os.path.exists(sentinel_folder):
                os.mkdir(sentinel_folder)

            xml_project = os.path.join(datefolder, 'sentinel_1', name + '.xml')
            link_project = os.path.join(datefolder, 'sentinel_1', name + '.SAFE.zip')
            kml_project = os.path.join(datefolder, 'sentinel_1', name + '.kml')
            preview_project = os.path.join(datefolder, 'sentinel_1', name + '.jpg')

        # Save .xml files
        prod = etree.ElementTree(product)

        if not os.path.exists(xml_dir):
            prod.write(xml_dir, pretty_print = True)
        if project_folder:
            if not os.path.exists(xml_project):
                prod.write(xml_project, pretty_print = True)

        prev = "'preview'"
        png = "'quick-look.png'"
        kml = "'map-overlay.kml'"
        dat = "'" + name + ".SAFE'"

        preview_url = url[:-10] + '/Nodes(' + dat + ')/Nodes(' + prev + ')/Nodes(' + png + ')/' + urllib.parse.quote_plus('$value') + '"'
        kml_url = url[:-10] + '/Nodes(' + dat + ')/Nodes(' + prev + ')/Nodes(' + kml + ')/' + urllib.parse.quote_plus('$value') + '"'

        # Download data files and create symbolic link
        if xml_only == False: # So we also download the file
            if not os.path.exists(file_dir):
                wget_data = wget_base + url + ' -O ' + file_dir
                print('download url is:' + wget_data)
                os.system(wget_data)

                # Finally check whether the file is downloaded correctly. Otherwise delete file and wait for next round of
                # downloads.
                uuid = product.find('id').text
                valid = sentinel_quality_check(file_dir, uuid, user, password)
            else: # If the file already exists we assume it is valid.
                valid = True

            if valid == True:
                # First download additional files
                if not os.path.exists(preview_dir):
                    wget_preview = wget_base + preview_url + ' -O ' + preview_dir
                    os.system(wget_preview)
                if not os.path.exists(kml_dir):
                    wget_kml = wget_base + kml_url + ' -O ' + kml_dir
                    os.system(wget_kml)

                # Then copy to user folder and create links if project folder is used
                if project_folder:
                    if not os.path.exists(preview_project):
                        os.system('cp ' + preview_dir + ' ' + preview_project)
                    if not os.path.exists(kml_project):
                        os.system('cp ' + kml_dir + ' ' + kml_project)
                    if not os.path.exists(link_project):
                        os.system('ln -s ' + file_dir + ' ' + link_project)
            else:
                os.system('rm ' + file_dir)
                os.system('rm ' + xml_dir)
                os.system('rm ' + xml_project)

def sentinel_quality_check(filename, uuid, user, password):
    # Check whether the zip files can be unpacked or not. This is part of the download procedure.

    checksum_url = "https://scihub.copernicus.eu/dhus/odata/v1/Products('" + uuid + "')/Checksum/Value/" + urllib.parse.quote_plus('$value')
    request = urllib.request.Request(checksum_url)
    base64string = base64.b64encode('%s:%s' % (user, password))
    request.add_header("Authorization", "Basic %s" % base64string)

    # connect to server. Hopefully this works at once
    try:
        dat = urllib.request.urlopen(request)
    except:
        print('not possible to connect this time')
        return False

    html_dat = ''
    for line in dat:
        html_dat = html_dat + line

    # Check file on disk
    if sys.platform == 'darwin':
        md5 = subprocess.check_output('md5 ' + filename, shell=True)[-33:-1]
    elif sys.platform == 'linux2':
        md5 = subprocess.check_output('md5sum ' + filename, shell=True)[:32]
    else:
        'This function only works on mac or linux systems!'
        return False

    if md5 == html_dat.lower():
        return True
    else:
        return False

def download_orbits(start_date, end_date, pages=30, precise_folder='', restituted_folder =''):
    # This script downloads all orbits files from the precise orbits website, when pages is set to a very high number.
    # By default only the first page for the last two days (restituted) is checked.

    pages_res = min(pages, 60)
    pages_poe = pages # every day there are 8 restituted orbit files
    last_precise = '' # Last precise orbit file. Used to remove unused restituted orbit files.

    start_num = int(start_date[0:4] + start_date[5:7] + start_date[8:10])
    end_num = int(end_date[0:4] + end_date[5:7] + end_date[8:10])

    if precise_folder:
        for i in range(pages_poe):
            # First extract the orbitfiles from the page.

            url = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/?page=' + str(i + 1)
            gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
            try:
                page = urllib.request.urlopen(url, context=gcontext)
            except TypeError:
                page = urllib.request.urlopen(url)

            html = page.read().split('\n')
            orb_files = []

            for line in html:
                if re.search('<a .*href=.*>', line):
                    if re.search('EOF', line):
                        dat = re.search('<a href=.*>(.*)</a>', line)
                        orb_files.append(dat.group(1))

            if not last_precise:
                last_precise = orb_files[0]

            for orb in orb_files:
                # Download the orbit files
                filename = os.path.join(precise_folder, orb)

                if int(orb[42:50]) + 1 <= end_num and int(orb[42:50]) + 1 >= start_num:
                    url = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/' + orb
                    if not os.path.exists(filename):
                        try:
                            urllib.request.urlretrieve(url, filename, context=gcontext)
                        except TypeError:
                            urllib.request.urlretrieve(url, filename)
                        print(orb + ' downloaded')
                    else:
                        print(orb + ' already downloaded')
                else:
                    print(orb + ' is out of date range')

            if len(orb_files) > 0:
                if int(orb[42:50]) < start_num:
                    break

    if restituted_folder:
        now = datetime.datetime.now()
        last_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')
        diff = datetime.timedelta(days=25)

        print('Time difference to last date is ' + str((now - last_date).days))

        if now - last_date < diff:  # only run when precise orbits will not be available
            for i in range(pages_res):
                # First extract the orbitfiles from the page.

                url = 'https://qc.sentinel1.eo.esa.int/aux_resorb/?page=' + str(i + 1)
                gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
                try:
                    page = urllib.request.urlopen(url, context=gcontext)
                except TypeError:
                    page = urllib.request.urlopen(url)

                html = page.read().split('\n')
                orb_files = []

                for line in html:
                    if re.search('<a .*href=.*>', line):
                        if re.search('EOF', line):
                            dat = re.search('<a href=.*>(.*)</a>', line)
                            orb_files.append(dat.group(1))

                for orb in orb_files:
                    # Download the orbit files
                    filename = os.path.join(precise_folder, orb)

                    if int(orb[42:50]) + 1 <= end_num and int(orb[42:50]) + 1 >= start_num:
                        url = 'https://qc.sentinel1.eo.esa.int/aux_poeorb/' + orb
                        if not os.path.exists(filename):
                            try:
                                urllib.request.urlretrieve(url, filename, context=gcontext)
                            except TypeError:
                                urllib.request.urlretrieve(url, filename)
                            print(orb + ' downloaded')
                        else:
                            print(orb + ' already downloaded')
                    else:
                        print(orb + ' is out of date range')

                if len(orb_files) > 0:
                    if int(orb[42:50]) < start_num:
                        break


# Actually execute the code...
if __name__ == "__main__":

    stack_folder = sys.argv[1]

    xml_file = os.path.join(os.path.join(stack_folder, 'doris_input.xml'))
    tree = ET.parse(xml_file)
    settings = tree.getroot()[0]
    print('reading xml file stack ' + xml_file)

    ROI = settings.find('.shape_file_path').text
    polarisation = settings.find('.polarisation').text
    archive_folder = settings.find('.sar_data_folder').text
    track = settings.find('.track').text
    orbit_folder = settings.find('.orbits_folder').text

    start_date = settings.find('.start_date').text
    end_date = settings.find('.end_date').text

    # Standard settings
    level = 'L1'
    sensor_mode = 'IW'
    product = 'SLC'

    # user settings
    xml_name = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_xml_file = os.path.join(os.path.join(xml_name, 'install', 'doris_config.xml'))
    print('reading xml file settings doris ' + config_xml_file)
    tree = ET.parse(config_xml_file)
    settings = tree.getroot()
    user = settings.find('.scihub_username').text
    password = settings.find('.scihub_password').text

    products, links, dates = sentinel_available(start_day=start_date, end_day=end_date, ROI=ROI,
                                                polarisation=polarisation, sensor_mode=sensor_mode, track=track,
                                                orbit_direction='', level=level, product=product,user=user,
                                                password=password)

    sentinel_download(products, destination_folder=archive_folder, user=user, password=password)

    precise_folder = os.path.join(orbit_folder, 'precise')
    if not os.path.exists(precise_folder):
        os.makedirs(precise_folder)
    restituted_folder = os.path.join(orbit_folder, 'restituted')
    if not os.path.exists(restituted_folder):
        os.makedirs(restituted_folder)

    download_orbits(start_date, end_date, pages=100, precise_folder=precise_folder, restituted_folder=restituted_folder)
