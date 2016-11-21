# This files defines a class for metadata objects of sentinel images. Large part of the work depends on python readers
# from the tops toolbox.

from sentinel_1.functions.xml_query import xml_query
import warnings
import os
from sentinel_1.main_code.burst import BurstMeta

class SwathMeta(object):
    # Class which stores and gathers information of a swath in a sentinel dataset

    def __init__(self, path='', swath_no='1', pol='vv', xml='', data=''):
        # Initialize object variables

        # This will contain a list of burst objects
        self.bursts = []

        # The following contain the path of xml and data file. Also swath_no and polarisation are given.
        self.swath_xml = ''
        self.swath_data = ''
        self.swath_no = ''
        self.swath_pol = ''

        # These variables contain the metadata and the convex hull of the burst
        self.metadata = []
        self.coverage = []

        # This function creates an swath object and searches for available data and xml files. It gives an error when
        # either the path does not exist, no data or xml files can be found or the data and xml files do not match.'
        if (not xml or not data) and not path:
            warnings.warn('Please provide either a product path or xml and data path')
            return

        if not xml or not data:
            warnings.warn('Data path is now used')

            xml_dir = os.path.join(path, 'annotation')
            xml = [f for f in os.listdir(xml_dir) if os.path.isfile(os.path.join(xml_dir, f))]

            data_dir = os.path.join(path, 'measurement')
            data = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f))]

            # Select polarisation
            if not any(s in pol for s in ('hh','vv','hv','vh')):
                warnings.warn('Polarisation not recognized, using default (vv)')
                pol = 'vv'
            if not swath_no in ('1','2','3'):
                warnings.warn('Swath number not recognized, using default (1)')

            xml = [os.path.join(path,'annotation',x) for x in xml if x[12:14] in pol and x[6] == swath_no]
            data = [os.path.join(path,'measurement',x) for x in data if x[12:14] in pol and x[6] == swath_no]

        # Check if the data is there and if the filenames coincide.
        # print xml + str(len(xml))
        # print data + str(len(data))

        if type(xml) is str:
            xml = [xml]
        if type(data) is str:
            data = [data]
        if (len(xml) != 1 and type(xml) is list) or len(data) != 1:
            warnings.warn('Total number of files should be one!')
        if not os.path.exists(xml[0]) or not os.path.exists(data[0]):
            warnings.warn('Either xml or data path does not exist')
        if xml[0:-3] != data[0:-4]:
            warnings.warn('xml and data file do not correspond.')

        self.swath_xml = xml[0]
        self.swath_data = data[0]

    def meta_swath(self):
        # This function reads and stores metadata of different swaths in the swath objects.
        self.metadata = xml_query(self.swath_xml)

    def meta_burst(self,corners=True,precise_folder=''):
        # This function reads and stores metadata of different bursts in the bursts objects.

        if not self.metadata:
            self.meta_swath()

        bursts_num = len(self.metadata['aux']['azimuthTimeStart'])
        if self.bursts:
            self.bursts = []

        for no in range(bursts_num):
            self.bursts.append(BurstMeta(path='',swath_no=self.swath_no, pol=self.swath_pol, burst_num=no + 1,
                                         xml=self.swath_xml, data=self.swath_data))
            self.bursts[no].meta_burst(swath_meta=self.metadata,corners=corners)

    def write_res(self):
        # This function writes all metadata and performed processing steps to a .res file
        print 'In progress!'

    def write_xml(self):
        # This function writes all metadata and performed processing steps to a .xml file
        print 'In progress!'