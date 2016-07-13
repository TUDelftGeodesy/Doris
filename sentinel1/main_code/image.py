# This files defines a class for metadata objects of sentinel images. Large part of the work depends on python readers
# from the tops toolbox.
import os
from swath import SwathMeta
import warnings
import xml.etree.cElementTree as etree
from shapely.geometry import Polygon
import copy


class ImageMeta(object):
    # Object for image files for sentinel data

    def __init__(self, path='' , pol='all' ,swath_no=['1','2','3']):
        # Initialize function variables

        # This will contain a list of swath objects
        self.swaths = []

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

        ######################################################

        # This function creates an image object and searches for available data and xml files. It gives an error when
        # either the path does not exist, no data or xml files can be found or the data and xml files do not match.
        # It is possible to choose one of the available polarisations or swaths using the pol and swath variables.

        xml_dir = os.path.join(path, 'annotation')
        xml = [f for f in os.listdir(xml_dir) if os.path.isfile(os.path.join(xml_dir, f))]

        data_dir = os.path.join(path, 'measurement')
        data = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f))]

        # Select polarisation
        if any(s in pol for s in ('hh','vv','hv','vh')):
            xml = [x for x in xml if x[12:14] in pol]
            data = [x for x in data if x[12:14] in pol]
        elif pol != 'all':
            warnings.warn('Polarisation not recognized, using default (all)')

        # Select swaths
        xml = sorted([os.path.join(path,'annotation',x) for x in xml if x[6] in swath_no])
        data = sorted([os.path.join(path,'measurement',x) for x in data if x[6] in swath_no])

        # Check if the data is there and if the filenames coincide.
        if len(xml) == 0:
            warnings.warn('There are no xml files')
        else:
            for x in range(0,len(xml)):
                if os.path.basename(xml[x])[0:-3] != os.path.basename(data[x])[0:-4]:
                    warnings.warn('xml files and data files are not the same')

        # Initialize function values
        self.swaths_xml = xml
        self.swaths_data = data

    def meta_swath(self):
        # This function reads and stores metadata of different swaths in the swath objects.
        if not self.swaths:
            for i in range(len(self.swaths_xml)):
                xml = self.swaths_xml[i]
                data = self.swaths_data[i]

                swath = SwathMeta(xml=xml,data=data)
                swath.meta_swath()
                self.swaths.append(swath)

    def read_kml(self):
        # This function reads the kml file of this image, which can be used to select relevant images for a datastack

        # First check is .kml file exist
        if self.swaths_xml:
            self.image_kml = os.path.join(os.path.dirname(os.path.dirname(self.swaths_xml[0])), 'preview' , 'map-overlay.kml')
        else:
            warnings.warn('There is no xml file identified for this image.')
            return
        if not os.path.exists(self.image_kml):
            warnings.warn('.kml file does not exist.')
            return

        in_kml = etree.parse(self.image_kml)
        in_kml = in_kml.getroot()
        coor = in_kml[0][1][1][2][0].text
        coor = [i.split(',') for i in coor.split(' ')]
        self.coverage = Polygon([[float(i[0]),float(i[1])] for i in coor])

    def meta_burst(self,corners=True,precise_dir=''):
        # This function reads and stores metadata of different bursts in the bursts objects.
        if not self.swaths:
            self.meta_swath()

        print(precise_dir)

        for i in range(len(self.swaths)):
            if not self.swaths[0].bursts:
                for i in range(len(self.swaths)):
                    self.swaths[i].meta_burst(corners=corners)
            elif corners == True and not self.swaths[0].bursts[0].burst_coverage:
                for i in range(len(self.swaths)):
                    self.swaths[i].meta_burst(corners=corners)

    def write_res(self):
        # This function writes all metadata and performed processing steps to a .res file
        print 'In progress!'

    def write_xml(self):
        # This function writes all metadata and performed processing steps to a .xml file
        print 'In progress!'
