# This files defines a class for metadata objects of sentinel images. Large part of the work depends on python readers
# from the tops toolbox.

from sentinel_1.functions.xml_query import xml_query
from sentinel_1.functions.burst_metadata import burst_header, burst_readfiles, burst_datapoints, burst_crop, burst_precise
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
        self.swath_meta = []

        # The following contain the path of xml and data file for swath and burst.
        self.swath_xml = ''
        self.swath_data = ''
        self.burst_res = ''
        self.burst_data = ''

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
        if not burst_num:
            warnings.warn('No burst number given')

        self.swath_xml = xml[0]
        self.swath_data = data[0]
        self.burst_num = burst_num
        self.swath_num = int(os.path.basename(xml[0])[6])

    def meta_burst(self,swath_meta=[],corners=True):
        # This function reads and stores metadata of the burst based on the xml file.

        if not swath_meta:
            if not self.swath_meta:
                self.swath_meta = xml_query(self.swath_xml)
        else:
            self.swath_meta = swath_meta

        readfiles, coverage = burst_readfiles(copy.deepcopy(self.swath_meta),self.burst_num,self.swath_data,corners=corners)

        # Read metadata from xml and inserts in resdata of burst
        self.insert(readfiles,process='readfiles')
        self.burst_coverage = coverage
        self.burst_center = [readfiles['Scene_centre_longitude'],readfiles['Scene_centre_latitude']]

    def res_burst(self,swath_meta=[],precise_folder=''):
        # Here the details for the res file are loaded.

        if not swath_meta:
            if not self.swath_meta:
                self.swath_meta = xml_query(self.swath_xml)
        else:
            self.swath_meta = swath_meta

        header = burst_header('master.res')
        crop = burst_crop(self.swath_meta,self.burst_num,self.swath_data,self.new_burst_num)

        # Read metadata from xml and inserts in resdata of burst
        self.header = header

        if not precise_folder:
            datapoints = burst_datapoints(self.swath_meta,self.burst_num)
            self.insert(datapoints,process='leader_datapoints')
        else:
            if os.path.exists(precise_folder):
                datapoints = burst_precise(self.swath_meta,self.burst_num,precise_folder,type='POE')
                if datapoints: # If it is empty, fall back to leader datapoints.
                    self.insert(datapoints,process='precise_orbits')
                else:
                    datapoints = burst_datapoints(self.swath_meta,self.burst_num)
                    self.insert(datapoints,process='leader_datapoints')
            else:
                datapoints = burst_datapoints(self.swath_meta,self.burst_num)
                self.insert(datapoints,process='leader_datapoints')

        # Insert crop after precise orbits..
        self.insert(crop,process='crop')

    def precise_orbits(self,swath_meta=[],precise_folder=''):
        # If the precise orbits where not added in an earlier step it can be done here.

        if not swath_meta:
            swath_meta = xml_query(self.swath_xml)

        if os.path.exists(precise_folder):
            datapoints = burst_precise(swath_meta,self.burst_num,precise_folder,type='POE')
            self.insert(datapoints,process='precise_orbits')
            self.del_process(process='leader_datapoints')
        else:
            print 'Precise orbit folder does not exist'

    def create_datafile(self,path='',filename=''):
        # This function creates a datafile based on metadata from the metadata file.
        if not self.metadata:
            warnings.warn('There is not metadata variable created yet! Use self.meta_burst to create one.')
        else:
            if not filename:
                filename = os.path.basename(self.swath_data)[0:-5] + '_burst_' + str(self.burst_num)
                burst_dir = os.path.join(path,filename)
                # create file (Maybe better to do this in a seperate folder)
