# This files defines a class for metadata objects of sentinel images. Large part of the work depends on python readers
# from the tops toolbox.

from sentinel_1.functions.burst_metadata import burst_header, burst_readfiles, burst_crop
from sentinel_1.main_code.resdata import ResData
import warnings
import os
import copy


class BurstMeta(ResData):
    # Class which holds and gathers information of a specific burst for sentinel 1 data.

    def __init__(self, path='', swath_no='1', pol='vv', burst_num=1, xml='', data=''):
        # Initialize variables

        # This indicates the burst number in the swath, acquisition date, centre location and the coverage of the burst (if available).
        self.burst_num = []
        self.new_burst_num = []
        self.swath_num = []
        self.burst_date = []
        self.burst_center = []
        self.burst_coverage = []
        self.burst_corners = []
        self.swath_meta = []

        # The following contain the path of xml and data file for swath and burst.
        self.swath_xml = ''
        self.swath_data = ''
        self.burst_res = ''
        self.burst_data = ''

        # orbits
        self.datapoints = []
        self.orbit_type = []

        #############################################################

        # This function creates an swath object and searches for available data and xml files. It gives an error when
        # either the path does not exist, no data or xml files can be found or the data and xml files do not match.'

        # Create resdata for this object
        super(BurstMeta, self).__init__(type='single')

        if not xml or not data:
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
        if not burst_num:
            warnings.warn('No burst number given')

        self.swath_xml = xml[0]
        self.swath_data = data[0]
        self.burst_num = burst_num
        self.swath_num = int(os.path.basename(xml[0])[6])

    def meta_burst(self, swath_meta=[], corners=True):
        # This function reads and stores metadata of the burst based on the swath xml file. If

        readfiles = burst_readfiles(copy.deepcopy(self.swath_meta), self.burst_num, self.burst_center, self.burst_corners, self.swath_data)
        crop = burst_crop(self.swath_meta, self.burst_num, self.swath_data, self.new_burst_num)

        # Read metadata from xml and inserts in resdata of burst
        # Insert the different steps (readfiles, orbits and crop)
        self.header = burst_header('master.res')
        self.insert(self.datapoints, process=self.orbit_type)
        self.insert(readfiles, process='readfiles')
        self.insert(crop, process='crop')



