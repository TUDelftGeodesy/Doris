# This files defines a class for metadata objects of sentinel images. Large part of the work depends on python readers
# from the tops toolbox.
import os
import warnings
import zipfile

from sentinel_1.main_code.swath import SwathMeta


class ImageMeta(object):
    # Object for image files for sentinel data

    def __init__(self, path='' , pol='all' ,swath_no=['1','2','3']):
        # Initialize function variables

        # This will contain a list of swath objects
        self.swaths = []
        self.pol = pol
        self.swath_no = swath_no

        # The following contain the path of xml and data files
        self.swaths_xml = ''
        self.swaths_data = ''
        self.image_kml = ''

        # This variable contains the convex hull of all swaths together
        self.metadata = []
        self.coverage = []

        # The following variables store data on which processing steps are performed and the results of these steps.
        self.steps = []
        self.steps_res = []

        # The following variable is to check the total number of bursts for this image. This is used to remove time
        # slots with less bursts.
        self.burst_no = 0

        # Check if the data is unzipped or not. If unzipped run the further initialization.
        self.zip_path = ''
        self.unzip_path = ''
        if path.endswith('.zip'):
            self.zip_path = path
        else:
            self.unzip_path = path

        # orbit information for this image
        self.orbit = ''

    ######################################################

    def init_unzipped(self, unzip_path=''):
        # This function creates an image object and searches for available data and xml files. It gives an error when
        # either the path does not exist, no data or xml files can be found or the data and xml files do not match.
        # It is possible to choose one of the available polarisations or swaths using the pol and swath variables.
        xml_dir = os.path.join(self.unzip_path, 'annotation')
        xml = [f for f in os.listdir(xml_dir) if os.path.isfile(os.path.join(xml_dir, f))]

        data_dir = os.path.join(self.unzip_path, 'measurement')
        data = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f))]

        # Select polarisation
        if any(s in self.pol for s in ('hh','vv','hv','vh')):
            xml = [x for x in xml if x[12:14] in self.pol]
            data = [x for x in data if x[12:14] in self.pol]
        elif self.pol != 'all':
            warnings.warn('Polarisation not recognized, using default (all)')

        # Select swaths
        xml = sorted([os.path.join(self.unzip_path,'annotation',x) for x in xml if x[6] in self.swath_no])
        data = sorted([os.path.join(self.unzip_path,'measurement',x) for x in data if x[6] in self.swath_no])

        # Initialize function values
        dat = [os.path.basename(d) for d in data]
        self.swaths_xml = [x for x in xml if os.path.basename(x)[:-4] + '.tiff' in dat]
        self.swaths_data = data

        # Check if the data is there and if the filenames coincide.
        if len(self.swaths_xml) == 0:
            warnings.warn('There are no xml files')


    def unzip(self, unzip_path=''):
        # This function unzips the corresponding image, based on some requirements.
        # Note that this is a backup function, while most unpacking is done by load_shape_unzip.py
        if not os.path.exists(self.path):
            try:
                zip = zipfile.ZipFile(self.path + '.zip')
                path = os.path.abspath(os.path.join(self.path, os.pardir))
                zip.extractall(path)
                return True
            except:
                print('Failed to unpack!')
                return False
        else:
            return True

    def meta_swath(self, precise_folder=''):
        # This function reads and stores metadata of different swaths in the swath objects.
        orbits = []

        if not self.swaths_xml:
            self.init_unzipped()

        if not self.swaths:
            for i in range(len(self.swaths_xml)):
                xml = self.swaths_xml[i]
                data = self.swaths_data[i]

                # Initialize swath and load data from xml file
                swath = SwathMeta(xml=xml, data=data)
                swath.meta_swath()

                # Calculate the orbits for this swath and reuse it for other swaths if it is already calculated
                if not orbits:
                    orbits = swath.orbits_swath(precise_folder=precise_folder)
                else:
                    swath.orbits = orbits

                # Define the resdata for the individual burst. And append the swath to the image object.
                swath.meta_burst()
                self.swaths.append(swath)
