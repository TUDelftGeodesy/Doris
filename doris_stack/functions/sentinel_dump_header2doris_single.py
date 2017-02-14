#!/usr/bin/env python


# live code for parsing of XML Sentinel file into python data structures
# and from there into DORIS res file structure

# this is rather fair implementation and should be used as a structure
# for the future implementation of XML/GeoTIFF data readers

#[BO] EPD (Entought) running on Centos fails at from lxml import etree.
#     The exception tries another import. The cElementTree has c bindings.
#     I believe it is faster than the regular ElementTree. 

#from lxml import etree
try:
    import xml.etree.cElementTree as etree
except:
    try:
        from lxml import etree
    except:
        #import xml.etree.ElementTree as etree
        print 'Failed to load lxml.etree or xml.etree.cElementTree'
        sys.exit(1)

import string, time, sys, os
import math
import numpy
import numpy as np

def usage():
    print '\nUsage: python Sentinel_dump_header2doris.py Sentinel_XML_product outputfile'
    print '  where Sentinel_XML_product is the input filename'
    print '        outputfile      is the output DORIS resultfile'
    print 'For example:python Sentinel_dump_header2doris.py 20140801.xml 3 5 /home/gertmulder/datastack/2015_10_30/'
    print 'Adapted from Wu Wenhao QQ:460249274   Wuhan'
    print 'Changed by Gert Mulder, TU Delft'
try:
    xml_file        =sys.argv[1]
    stack_num       =sys.argv[2]
    image_num       =sys.argv[3]
    write_path      =sys.argv[4]
    # xml_file defines the pathname of the used xml source
    # stack_num is the burst number in the stack for this swath
    # image_num is the burst number in the swath data file of the product
    # write_path is the path where we write the data to

except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)


import time, sys, os
import math
from orbit_coordinates import xyz2ell, lph2xyz, intrp_orbit, hms2sec
import warnings

if not os.path.exists(xml_file):
    warnings.warn('xml file does not exist!')
if not os.path.exists(write_path):
    warnings.warn('path where file is written does not exist!')
if isinstance(stack_num,basestring):
    stack_num = int(stack_num)
if isinstance(image_num,basestring):
    image_num = int(image_num)

inTree = etree.parse(xml_file)

xml_file = os.path.basename(xml_file)

swath_num = xml_file[6]
name_ID = xml_file[15:23]

# query syntax for every field
queryList = {\
             # mission info
             'mission'      : './/adsHeader/missionId',\
             'image_mode'   : './/adsHeader/mode',\
             'polarisation' :  './/adsHeader/polarisation',\
             #Swath info
             'numberOfSamples_Swath': './/imageAnnotation/imageInformation/numberOfSamples',\
             'numberOfLines_Swath': './/imageAnnotation/imageInformation/numberOfLines',\
             'azimuthPixelSpacing':'.//imageAnnotation/imageInformation/azimuthPixelSpacing',\
             'rangePixelSpacing'  :'.//imageAnnotation/imageInformation/rangePixelSpacing',
             'Swath_startTime'  :'.//adsHeader/startTime',
             'Swath_stopTime'   :'.//adsHeader/stopTime',
             #Swath info
             #raw data info
             'Azimuth_steering_rate'  : './/generalAnnotation/productInformation/azimuthSteeringRate',\
             'PRF_raw_data'  : './/generalAnnotation/downlinkInformationList/downlinkInformation/prf',\
              #
             # Burst imageData file
            # 'imageData'  : './/productComponents/imageData/file/location/filename',\
             'imageLines' : './/swathTiming/linesPerBurst',\
             'imagePixels': './/swathTiming/samplesPerBurst',\
            # ValidSample
              'firstValidSample' : './/swathTiming/burstList/burst/firstValidSample',\
              'lastValidSample' : './/swathTiming/burstList/burst/lastValidSample',\
              # volume info
            # 'volFile' : './/productComponents/annotation/file/location/filename',\
             'volID'   : './/adsHeader/missionDataTakeId',\
            # 'volRef'  : './/generalHeader/referenceDocument',\
             # product info
             'radarFrequency' : './/generalAnnotation/productInformation/radarFrequency',\
             'productSpec'    : './/generalHeader/referenceDocument',\
             'productVolDate' : './/setup//IOCSAuxProductGenerationTimeUTC',\
             'productDate'    : './/generalHeader/generationTime',\
             'productFacility': './/productInfo/generationInfo/level1ProcessingFacility',\

             # scene info
             'scenePol'     : './/adsHeader/polarisation',\
             'sceneMode'    : './/adsHeader/mode',\
             'sceneSwath'    : './/adsHeader/swath',\
             'sceneCenLine_number'  : './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/line',\
             'sceneCenPixel_number'  : './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/pixel',\
             'sceneCenLat'  : './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/latitude',\
             'sceneCenLon'  : './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/longitude',\
             'height'   : './/geolocationGrid/geolocationGridPointList/geolocationGridPoint/height',\
             'sceneRecords' : './/imageDataInfo/imageRaster/numberOfRows',\
             # orbit info
             'orbitABS' : './/adsHeader/absoluteOrbitNumber',\
             'orbitDir' : './/generalAnnotation/productInformation/pass',\
             'orbitTime': './/generalAnnotation/orbitList/orbit/time',\
             'orbitX'   : './/generalAnnotation/orbitList/orbit/position/x',\
             'orbitY'   : './/generalAnnotation/orbitList/orbit/position/y',\
             'orbitZ'   : './/generalAnnotation/orbitList/orbit/position/z',\
             # range
             'rangeRSR'     :'.//generalAnnotation/productInformation/rangeSamplingRate',\
             'rangeBW'      :'.//imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/processingBandwidth',\
             'rangeWind'    :'.//imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/rangeProcessing/windowType',\
             'rangeTimePix' :'.//imageAnnotation/imageInformation/slantRangeTime',\
             # azimuth
             'azimuthPRF'       :'.//imageAnnotation/imageInformation/azimuthFrequency',\
             'azimuthTimeInterval' :'.//imageAnnotation/imageInformation/azimuthTimeInterval',\
             'azimuthBW'        :'.//imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/azimuthProcessing/totalBandwidth',\
             'azimuthWind'      :'.//imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/azimuthProcessing/windowType',\
             'azimuthTimeStart' : './/swathTiming/burstList/burst/azimuthTime',\
             'heading' : './/generalAnnotation/productInformation/platformHeading',\
             # doppler
             'doppler_azimuth_Time'  :'.//dopplerCentroid/dcEstimateList/dcEstimate/azimuthTime',\
             'doppler_range_Time'  :'.//dopplerCentroid/dcEstimateList/dcEstimate/t0',\
             'dopplerCoeff' :'.//dopplerCentroid/dcEstimateList/dcEstimate/dataDcPolynomial',\
             'azimuthFmRate_reference_Azimuth_time':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthTime',\
             'azimuthFmRate_reference_Range_time':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/t0',\
             # FMRate  old version
             'azimuthFmRate_c0':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/c0',\
             'azimuthFmRate_c1':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/c1',\
             'azimuthFmRate_c2':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/c2',\
             # FMRate  new version
             'azimuthFmRatePolynomial':'.//generalAnnotation/azimuthFmRateList/azimuthFmRate/azimuthFmRatePolynomial',\

             }

# temp variables and parameters
container     = {}
# containerTemp = {}
events        = ('end',)

for key in queryList.keys():

    try:
        vars()[key]
    except KeyError or NameError:
        vars()[key] = []

    for nodes in inTree.findall(queryList[key]):

        vars()[key].append(nodes.text)

    container[key] = vars()[key]


# ---------------------------------------------------------------------------------------------------------
spacing=55             # space key      format output
dummyVar = 'DUMMY'
Swath_number=str(container['sceneSwath'][0])

outputFile_Name_res = name_ID + '_iw_' + swath_num + '_burst_' + str(stack_num) + '.res'
outStream      = open(os.path.join(write_path,outputFile_Name_res),'w')
#Doppler Parameter
doppler_all=container['dopplerCoeff'][image_num-1]
doppler_parameter=doppler_all.split() # Doppler centroid parameter
# Image Centroid
container['Image_centroid_lat']=0
container['Image_centroid_lon']=0
Line_location_index_down    = image_num*float(container['imageLines'][0]) # Last line
Line_location_index_up      =(image_num-1)*float(container['imageLines'][0]) #First line  May be not accurate
cen_line                    = (Line_location_index_down +Line_location_index_up)/2
cen_pixel                   = int(container['imagePixels'][0])/2

Distance=float('-inf')
for i in range(len(container['sceneCenLat'])):
     lines_number=float(container['sceneCenLine_number'][i])
     pixel_number=float(container['sceneCenPixel_number'][i])

     if Distance<(lines_number-cen_line)**2+(cen_pixel-pixel_number)**2:
         Distance=(lines_number-cen_line)**2+(cen_pixel-pixel_number)**2
for i in range(len(container['sceneCenLat'])):
     lines_number=float(container['sceneCenLine_number'][i])
     pixel_number=float(container['sceneCenPixel_number'][i])

     if Distance==((lines_number-cen_line)**2+(cen_pixel-pixel_number)**2):
         norm_orbit,norm_orbit_line=intrp_orbit(int(math.fabs(cen_line-Line_location_index_up)),container,image_num-1)
         coord_xyz                =lph2xyz(int(math.fabs(cen_line-Line_location_index_up)),cen_pixel,container,norm_orbit_line,\
                          float(container['sceneCenLon'][i]),float(container['sceneCenLat'][i]),float(container['height'][i]))
         phi_lam_height=xyz2ell(coord_xyz)
         container['Image_centroid_lon']=phi_lam_height[1]
         container['Image_centroid_lat']=phi_lam_height[0]
         height            =phi_lam_height[2]

# ---------------------------------------------------------------------------------------------------------
outStream.write('===============================================\n')
outStream.write('MASTER RESULTFILE:                    %s\n' % outputFile_Name_res)
outStream.write('Created by:                           %s\n')
outStream.write('InSAR Processor: Doris (Delft o-o Radar Interferometric Software)\n')
outStream.write('Version:     Version   (2015) (For TOPSAR)\n')
outStream.write('FFTW library:                        used\n')
outStream.write('VECLIB library:                  not used\n')
outStream.write('LAPACK library:                  not used\n')
outStream.write('Compiled at:                     XXXXXXXX\n')
outStream.write('By GUN gcc                       XXXXXXXX\n')
outStream.write('===============================================\n')

outStream.write('\n')
outStream.write('\n')
outStream.write('Start_process_control\n')
outStream.write('readfiles:		1\n')
outStream.write('precise_orbits:         0\n')
outStream.write('crop:			0\n')
outStream.write('sim_amplitude:		0\n')
outStream.write('master_timing:		0\n')
outStream.write('oversample:		0\n')
outStream.write('resample:		0\n')
outStream.write('filt_azi:		0\n')
outStream.write('filt_range:		0\n')
outStream.write('NOT_USED:		0\n')
outStream.write('End_process_control\n')
outStream.write('\n')
# ---------------------------------------------------------------------------------------------------------
outStream.write('\nSentinel_dump_header2doris.py v1.0, doris software,QQ:460249274 (China)\n')
outStream.write('*******************************************************************\n')
outStream.write('*_Start_readfiles:\n')
outStream.write('*******************************************************************\n')
outStream.write(' '.join(['%s%s%s' %('Volume_file:'.ljust(spacing),           dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Volume_ID:'.ljust(spacing),             container['volID'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Volume_identifier:'.ljust(spacing),     dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Volume_set_identifier:'.ljust(spacing), dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Number of records in ref. file:'.ljust(spacing),dummyVar,'\n')]))

if cmp(str(container['mission'][0]),'RS2')==0:
     outStream.write('SAR_PROCESSOR:					RadarSAT-2\n')
elif cmp(str(container['mission'][0]),'S1A')==0:
     outStream.write('SAR_PROCESSOR:					Sentinel-1A\n')
else:
     outStream.write('SAR_PROCESSOR:					Unknown\n')
outStream.write(' '.join(['%s%s%s' %('SWATH:'.ljust(spacing),                container['sceneSwath'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('PASS:'.ljust(spacing),                 container['orbitDir'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('IMAGE_MODE:'.ljust(spacing),           container['image_mode'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('polarisation:'.ljust(spacing),         container['polarisation'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Product type specifier:'.ljust(spacing), container['mission'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Logical volume generating facility:'.ljust(spacing), dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Logical volume creation date:'.ljust(spacing), dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Location and date/time of product creation:'.ljust(spacing), dummyVar,'\n')]))
outStream.write(' '.join(['%s%s%s' %('Number_of_lines_Swath:'.ljust(spacing), container['numberOfLines_Swath'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Number_of_pixels_Swath:'.ljust(spacing), container['numberOfSamples_Swath'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('rangePixelSpacing:'.ljust(spacing),    container['rangePixelSpacing'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('azimuthPixelSpacing:'.ljust(spacing),  container['azimuthPixelSpacing'][0],'\n')]))
outStream.write(' '.join(['%s%d%s' %('total_Burst:'.ljust(spacing),          len(container['azimuthTimeStart']),'\n')]))
outStream.write(' '.join(['%s%d%s' %('Burst_number_index:'.ljust(spacing),   image_num,'\n')]))
outStream.write(' '.join(['%s%s%s' %('RADAR_FREQUENCY (HZ):'.ljust(spacing), container['radarFrequency'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Scene identification:'.ljust(spacing), ('Orbit: '+container['orbitABS'][0]),'\n')]))
outStream.write(' '.join(['%s%s%s' %('Scene location:'.ljust(spacing),('lat: '+ container['sceneCenLat'][0] +' lon:'+container['sceneCenLon'][0]),'\n')]))
outStream.write(' '.join(['%s%s%s' %('Sensor platform mission identifer:'.ljust(spacing), container['mission'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Scene_center_heading:'.ljust(spacing),  container['heading'][0],'\n')]))
outStream.write(' '.join(['%s%f%s' %('Scene_centre_latitude:'.ljust(spacing), container['Image_centroid_lat'],'\n')]))
outStream.write(' '.join(['%s%f%s' %('Scene_centre_longitude:'.ljust(spacing),container['Image_centroid_lon'],'\n')]))
outStream.write(' '.join(['%s%.9f%s' %('Radar_wavelength (m):'.ljust(spacing),(299792458.0/float(container['radarFrequency'][0])),'\n')]))
outStream.write(' '.join(['%s%s%s' %('Azimuth_steering_rate (deg/s):'.ljust(spacing), container['Azimuth_steering_rate'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Pulse_Repetition_Frequency_raw_data(TOPSAR):'.ljust(spacing), container['PRF_raw_data'][0],'\n')]))
outStream.write(' '.join(['%s%s%s%s' %('First_pixel_azimuth_time (UTC):'.ljust(spacing),time.strftime("%d-%b-%Y",time.strptime(container['azimuthTimeStart'][image_num-1].split('T')[0],"%Y-%m-%d")),' ' +container['azimuthTimeStart'][image_num-1].split('T')[1],'\n')]))

if (cmp(str(container['mission'][0]),'RS2')==0):
     outStream.write(' '.join(['%s%.9f%s' %('Pulse_Repetition_Frequency (computed, Hz):'.ljust(spacing),(1.0/float(container['azimuthTimeInterval'][0])),'\n')]))
elif cmp(str(container['mission'][0]),'S1A')==0:
     outStream.write(' '.join(['%s%s%s' %('Pulse_Repetition_Frequency (computed, Hz):'.ljust(spacing),container['azimuthPRF'][0],'\n')]))
else:
     outStream.write(' '.join(['%s%s%s' %('Pulse_Repetition_Frequency (computed, Hz):'.ljust(spacing),'Unknown','\n')]))
outStream.write(' '.join(['%s%s%s' %('Azimuth_time_interval (s):'.ljust(spacing), container['azimuthTimeInterval'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Total_azimuth_band_width (Hz):'.ljust(spacing), container['azimuthBW'][0],'\n')]))
outStream.write(' '.join(['%s%s%s' %('Weighting_azimuth:'.ljust(spacing), str.upper(container['azimuthWind'][0]),'\n')]))
outStream.write(' '.join(['%s%.15f%s' %('Range_time_to_first_pixel (2way) (ms):'.ljust(spacing), (float(container['rangeTimePix'][0])*1000),'\n')]))
outStream.write(' '.join(['%s%.9f%s' %('Range_sampling_rate (computed, MHz):'.ljust(spacing), (float(container['rangeRSR'][0])/1000000),'\n')]))
outStream.write(' '.join(['%s%.9f%s' %('Total_range_band_width (MHz):'.ljust(spacing), (float(container['rangeBW'][0])/1000000),'\n')]))
outStream.write(' '.join(['%s%s%s' %('Weighting_range:'.ljust(spacing), str.upper(container['rangeWind'][0]),'\n')]))


for doppler_index in range(len(container['doppler_azimuth_Time'])):
     # Check for zero doppler azimuth time during this burst. If a first time is found, break the loop, so only one time is saved. (There is only one moment of zero doppler every burst.)
     if hms2sec(container['doppler_azimuth_Time'][doppler_index].split('T')[1]) > hms2sec(container['azimuthTimeStart'][image_num-1].split('T')[1]):
         outStream.write(' '.join(['%s%s%s%s' %('DC_reference_azimuth_time:'.ljust(spacing),\
time.strftime("%d-%b-%Y",time.strptime(container['doppler_azimuth_Time'][doppler_index].split('T')[0],"%Y-%m-%d")),\
' '+container['doppler_azimuth_Time'][doppler_index].split('T')[1],'\n')]))
         outStream.write(' '.join(['%s%s%s' %('DC_reference_range_time:'.ljust(spacing), container['doppler_range_Time'][doppler_index],'\n')]))
         Doppler_param_all=container['dopplerCoeff'][doppler_index]
         Doppler_parameter=Doppler_param_all.split()
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_constant (Hz, early edge):'.ljust(spacing), Doppler_parameter[0],'\n')]))
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_linear (Hz/s, early edge):'.ljust(spacing), Doppler_parameter[1],'\n')]))
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_quadratic (Hz/s/s, early edge):'.ljust(spacing), Doppler_parameter[2],'\n')]))
         break
     # If there is no doppler azimuth time that correspond with the azimuth time of the burst, we use zeros. ??? What will happen then?
     if doppler_index== (len(container['doppler_azimuth_Time'])-1):
         outStream.write(' '.join(['%s%s%s%s' %('DC_reference_azimuth_time:'.ljust(spacing),'0.0',' 0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('DC_reference_range_time:'.ljust(spacing), container['doppler_range_Time'][doppler_index],'\n')]))
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_constant (Hz, early edge):'.ljust(spacing), '0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_linear (Hz/s, early edge):'.ljust(spacing), '0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('Xtrack_f_DC_quadratic (Hz/s/s, early edge):'.ljust(spacing), '0.0','\n')]))

for doppler_index in range(len(container['azimuthFmRate_reference_Azimuth_time'])):
     if hms2sec(container['azimuthFmRate_reference_Azimuth_time'][doppler_index].split('T')[1]) > hms2sec(container['azimuthTimeStart'][image_num-1].split('T')[1]):
         outStream.write('FM_reference_azimuth_time:				%s %s\n'\
 % (time.strftime("%d-%b-%Y",time.strptime(container['azimuthFmRate_reference_Azimuth_time'][doppler_index].\
split('T')[0],"%Y-%m-%d")),container['azimuthFmRate_reference_Azimuth_time'][doppler_index].split('T')[1]))
         outStream.write(' '.join(['%s%s%s' %('FM_reference_range_time:'.ljust(spacing), container['azimuthFmRate_reference_Range_time'][doppler_index],'\n')]))
         if container['azimuthFmRate_c0']:
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_constant_coeff (Hz, early edge):'.ljust(spacing), container['azimuthFmRate_c0'][doppler_index],'\n')]))
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_linear_coeff (Hz/s, early edge):'.ljust(spacing), container['azimuthFmRate_c1'][doppler_index],'\n')]))
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_quadratic_coeff (Hz/s/s, early edge):'.ljust(spacing), container['azimuthFmRate_c2'][doppler_index],'\n')]))
         else:
             FmRate_parameter_all = container['azimuthFmRatePolynomial'][doppler_index]
             FmRate_parameter     = FmRate_parameter_all.split()
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_constant_coeff (Hz, early edge):'.ljust(spacing), FmRate_parameter[0],'\n')]))
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_linear_coeff (Hz/s, early edge):'.ljust(spacing), FmRate_parameter[1],'\n')]))
             outStream.write(' '.join(['%s%s%s' %('FM_polynomial_quadratic_coeff (Hz/s/s, early edge):'.ljust(spacing), FmRate_parameter[2],'\n')]))
         break
     if doppler_index==(len(container['azimuthFmRate_reference_Azimuth_time'])-1):
         outStream.write(' '.join(['%s%s%s%s' %('FM_reference_azimuth_time:'.ljust(spacing),'0.0',' 0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('FM_polynomial_constant_coeff (Hz, early edge):'.ljust(spacing), '0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('FM_polynomial_linear_coeff (Hz/s, early edge):'.ljust(spacing), '0.0','\n')]))
         outStream.write(' '.join(['%s%s%s' %('FM_polynomial_quadratic_coeff (Hz/s/s, early edge):'.ljust(spacing), '0.0','\n')]))


outStream.write('')
outStream.write('*******************************************************************\n')
outStream.write('Datafile:						' +  outputFile_Name_res[:-3] + '.tiff \n' )
outStream.write('Dataformat:						%s\n' % 'tiff')
outStream.write('Number_of_lines_original:				%s\n' % container['imageLines'][0])
outStream.write('Number_of_pixels_original:				%s\n' % container['imagePixels'][0])
outStream.write('Original burst number: %s\n' % image_num)
outStream.write('New burst number: %s\n' % stack_num)
outStream.write('*******************************************************************\n')
outStream.write('* End_readfiles:_NORMAL\n')
outStream.write('*******************************************************************\n')
outStream.write('')
outStream.write('')
outStream.write('*******************************************************************\n')
outStream.write('*_Start_leader_datapoints\n')
outStream.write('*******************************************************************\n')
outStream.write(' t(s)		 X(m)		 Y(m)		 Z(m)\n')
temp=0
for i in range(len(container['orbitTime'])):
   if (hms2sec(container['orbitTime'][i].split('T')[1]) > hms2sec(container['Swath_startTime'][0].split('T')[1])-100)\
    and (hms2sec(container['orbitTime'][i].split('T')[1])< hms2sec(container['Swath_stopTime'][0].split('T')[1])+100):
       temp=temp+1

outStream.write('NUMBER_OF_DATAPOINTS: 			   %d\n' % temp)
outStream.write('')
for i in range(len(container['orbitTime'])):
   if (hms2sec(container['orbitTime'][i].split('T')[1]) > hms2sec(container['Swath_startTime'][0].split('T')[1])-100)\
    and (hms2sec(container['orbitTime'][i].split('T')[1])< hms2sec(container['Swath_stopTime'][0].split('T')[1])+100):
      outStream.write(' %s %s %s %s\n' % (hms2sec(container['orbitTime'][i].split('T')[1]),\
                                  container['orbitX'][i],\
                                  container['orbitY'][i],\
                                  container['orbitZ'][i]))

outStream.write('')
outStream.write('*******************************************************************\n')
outStream.write('* End_leader_datapoints:_NORMAL\n')
outStream.write('*******************************************************************\n')
outStream.close()


##############################################################################
#########*********************Output the Swath information ***************####
#***********Coregistration can be performed at beam level**************#######
#################################################################################


Swath_res_FileName= name_ID + '_iw_' + swath_num + '.res'
if os.path.exists(os.path.join(write_path,Swath_res_FileName)):
    exist = True
    outStream      = open(os.path.join(write_path,Swath_res_FileName),'a')
else:
    exist = False
    outStream      = open(os.path.join(write_path,Swath_res_FileName),'w')

if not exist:
    outStream.write('\n')
    outStream.write('\n for coregistration at beam level in the near future\n')
    outStream.write('Start_process_control\n')
    outStream.write('readfiles:		1\n')
    outStream.write('precise_orbits:         0\n')
    outStream.write('crop:			0\n')
    outStream.write('sim_amplitude:		0\n')
    outStream.write('master_timing:		0\n')
    outStream.write('oversample:		0\n')
    outStream.write('resample:		0\n')
    outStream.write('filt_azi:		0\n')
    outStream.write('filt_range:		0\n')
    outStream.write('NOT_USED:		0\n')
    outStream.write('End_process_control\n')
    outStream.write('\n')
    # ---------------------------------------------------------------------------------------------------------
    outStream.write('\nSentinel_dump_header2doris.py v1.0, doris software,QQ:460249274 (China)\n')
    outStream.write('*******************************************************************\n')
    outStream.write('*_Start_readfiles:\n')
    outStream.write('*******************************************************************\n')
    outStream.write('Scene identification: 				Orbit: %s\n' % (container['orbitABS'][0]))
    outStream.write('Sensor platform mission identifer:         	%s\n' % container['mission'][0])


FirstValidSample_num=" ".join(container['firstValidSample']).split(' ')
FirstValidSample_num_new=list(set(FirstValidSample_num))
FirstValidSample_num_new.sort()
FirstValidSample_num_new.remove('-1')

outStream.write('\n')
outStream.write('Original_swath_datafile_%d:                    %s\n' % (stack_num,xml_file[:-3] + 'tiff'))
outStream.write('Original_burst_number_burst_%d:                %s\n' % (stack_num,image_num))
outStream.write('Number_of_lines_swath_burst_%d: 			    %s\n' % (stack_num,container['numberOfLines_Swath'][0]))
outStream.write('Number_of_pixels_swath_burst_%d: 	            %s\n' % (stack_num,container['numberOfSamples_Swath'][0]))
outStream.write('Number_of_lines_burst_%d: 			            %s\n' % (stack_num,container['imageLines'][0]))
outStream.write('Number_of_pixel_burst_%d: 	                    %s\n' % (stack_num,container['imagePixels'][0]))
outStream.write('Number_of_burst: 	                            %d\n' % stack_num)
outStream.write('First_pixel_burst_%d (w.r.t. original_image): 	%s\n' %  (stack_num, 1))    #start form 1
outStream.write('Last_pixel_burst_%d (w.r.t. original_image): 	%d\n' %  (stack_num, int(container['imagePixels'][0]))) #start form 1
outStream.write('Burst_res_FileName:                            %s\n' % outputFile_Name_res)
outStream.write('First_line_burst_%d (w.r.t. original_image):   %d\n' % (stack_num,container['lastValidSample'][image_num-1].split(' ')[0:int(container['imageLines'][0])/2].count('-1')))
outStream.write('Last_line_burst_%d (w.r.t. original_image):    %d\n' % (stack_num,int(container['imageLines'][0])-container['lastValidSample'][image_num-1].split(' ')[:int(container['imageLines'][0])/2:-1].count('-1')))
