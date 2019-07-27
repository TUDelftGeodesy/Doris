#!/usr/bin/python

#-----------------------------------------------------------------#
# A python code for parsing RS2 XML file into python data structures
# and from there into DORIS res file structure
#
# Author: TUDelft - 2010
# Maintainer: Mahmut Arikan
#
# Developed based on tsx_dump_header2doris.py.
#
# this is rather fair implementation and should be used as a structure
# for the future implementation of XML/GeoTIFF data readers

# NO XPATH version
# 2012.Nov  Fix for ascending mode time parameters Samie Esfahany and Mahmut
# 2013.Sep  Fix for UF and MF acquisition modes and cleanup, Piers van der TOrren
#-----------------------------------------------------------------#

from lxml import etree
import sys
from datetime import datetime

#try:
#    from lxml import etree
#except ImportError:
#    import xml.etree.ElementTree as etree


codeRevision=1.2   # this code revision number

def usage():
    print('INFO    : @(#)Doris InSAR software, $Revision: %s $, $Author: TUDelft $' % codeRevision)
    print()
    print('Usage   : python rs2_dump_header2doris.py rs2_XML_product > outputfile')
    print('                           where rs2_XML_product is the input filename')
    print()
    print('          This software is part of Doris InSAR software package.\n')
    print('(c) 1999-2010 Delft University of Technology, the Netherlands.\n')

try:
    inputFileName  = sys.argv[1]
#    outputFileName = sys.argv[2]
#    outStream      = open(outputFileName,'w')
except:
    print('\nError   : Unrecognized input or missing arguments\n\n')
    usage()
    sys.exit(1)


# Helper functions
def nsmap_none(path, ns='None:'):
    """ Add a namespace to each tag in the given path which doesn't have one.
    """
    def add_ns_if_none(tag):
        if tag in ('', '.', '*') or ':' in tag:
            return tag
        else:
            return ''.join((ns, tag))
    return '/'.join(add_ns_if_none(tag) for tag in path.split('/'))

def hms2sec(hmsString):
    """ input hmsString syntax: XX:XX:XX.xxxxxx
    """
    return int(hmsString[0:2])*3600 + \
        int(hmsString[3:5])*60 + \
        float(hmsString[6:])

# constants
SOL = 299792458.0    # speed of light
dateformat = '%Y-%m-%dT%H:%M:%S.%fZ'

# default namespace
nsmap = {None:"http://www.rsi.ca/rs2/prod/xml/schemas"}
ns = '{' + nsmap[None] + '}'


inTree = etree.parse(inputFileName)

# query syntax for every field
queryList = {
             # mission info
             'mission'             : 'sourceAttributes/satellite',
             # imageFile file
             # fullResolutionImageData pole="HH" # inTree.findall('imageAttributes/fullResolutionImageData')[0].text, and more [1] .. [3]
             'imageFile'           : 'imageAttributes/fullResolutionImageData',
             #'imageLines'          : 'imageAttributes/rasterAttributes/numberOfLines',
             'imageLines'          : 'imageAttributes//numberOfLines',
             'imagePixels'         : 'imageAttributes//numberOfSamplesPerLine',
             'imageLineSpacing'    : 'imageAttributes//sampledLineSpacing',
             'imagePixelSpacing'   : 'imageAttributes//sampledPixelSpacing',
             # volume info
             #'volFile' : 'productComponents/annotation/file/location/filename', # HARDCODED!!! for radarsat-2
             # following smt like Level 1B Product, check manual
             'volID'               : 'productId',
             'volRef'              : 'documentIdentifier',
             # product info
             #productSpec'          : 'generalHeader/referenceDocument',  # TSX
             'productSpec'         : 'documentIdentifier',
             'productVolDate'      : 'imageGenerationParameters//processingTime',
             'productSoftVer'      : 'imageGenerationParameters//softwareVersion',
             'productDate'         : 'sourceAttributes/rawDataStartTime',
             'productFacility'     : 'imageGenerationParameters//processingFacility',
             # scene info
             #'scenePol'           : 'sourceAttributes/radarParameters/acquisitionType',    # Fine Quad Polarization
             'scenePol'            : 'sourceAttributes//polarizations',
             'sceneBeam'           : 'sourceAttributes//beams',
             'sceneBeamMode'       : 'sourceAttributes/beamModeMnemonic',
             'list_sceneLat'       : 'imageAttributes/geographicInformation/geolocationGrid/imageTiePoint/geodeticCoordinate/latitude',
             'list_sceneLon'       : 'imageAttributes/geographicInformation/geolocationGrid/imageTiePoint/geodeticCoordinate/longitude',
             'sceneRecords'        : 'imageGenerationParameters/sarProcessingInformation/numberOfLinesProcessed',
             'antennaLookDir'      : 'sourceAttributes//antennaPointing',
             'missinglines'        : 'sourceAttributes//numberOfMissingLines',
             # orbit info
             'orbitABS'            : 'sourceAttributes/orbitAndAttitude//orbitDataFile',
             'orbitDir'            : 'sourceAttributes//passDirection',
             'list_orbitTime'      : 'sourceAttributes//stateVector/timeStamp',
             'list_orbitX'         : 'sourceAttributes//stateVector/xPosition',
             'list_orbitY'         : 'sourceAttributes//stateVector/yPosition',
             'list_orbitZ'         : 'sourceAttributes//stateVector/zPosition',
             'list_orbitXV'        : 'sourceAttributes//stateVector/xVelocity',
             'list_orbitYV'        : 'sourceAttributes//stateVector/yVelocity',
             'list_orbitZV'        : 'sourceAttributes//stateVector/zVelocity',
             # range
             'list_rangeRSR'       : 'sourceAttributes//adcSamplingRate', # for UF mode there are two subpulses which have to be added together
             'rangeBW'             : 'imageGenerationParameters//rangeLookBandwidth',
             'rangeWind'           : 'imageGenerationParameters//rangeWindow/windowName',
             'rangeWindCoeff'      : 'imageGenerationParameters//rangeWindow/windowCoefficient',
             'rangeTimePix'        : 'imageGenerationParameters//slantRangeTimeToFirstRangeSample',
             # azimuth
             'azimuthPRF'          : 'sourceAttributes//pulseRepetitionFrequency', # for some modes (MF, UF) this value is changed in processing, calculate from other values
             'azimuthBW'           : 'imageGenerationParameters//azimuthLookBandwidth',
             'azimuthWind'         : 'imageGenerationParameters//azimuthWindow/windowName',
             'azimuthWindCoeff'    : 'imageGenerationParameters//azimuthWindow/windowCoefficient',
             'azimuthTimeFirstLine': 'imageGenerationParameters//zeroDopplerTimeFirstLine',
             'azimuthTimeLastLine' : 'imageGenerationParameters//zeroDopplerTimeLastLine',
             # doppler
             'dopplerTime'         : 'imageGenerationParameters//timeOfDopplerCentroidEstimate',
             'dopplerCoeff'        : 'imageGenerationParameters//dopplerCentroidCoefficients',
             # for wavelength computation
             'radarfreq'           : 'sourceAttributes//radarCenterFrequency',
             #  wavelength_computed = (0.000000001*SOL/atof(c8freq)) seems more reliable, BK 03/04
             }


# get variables and parameters from xml
container = {}
for key, value in queryList.items():
    if key.startswith('list_'):
        container[key] = [tag.text for tag in inTree.findall(nsmap_none(value, ns))]
    else:
        container[key] = inTree.findtext(nsmap_none(value, ns))
        if container[key] == None:
            raise Exception('Path {0} not found in XML'.format(value))

container['dopplerCoeff'] = container['dopplerCoeff'].split()

def mean(l):
    return sum(l)/len(l)
container['sceneCenLat'] = mean([float(val) for val in container['list_sceneLat']])
container['sceneCenLon'] = mean([float(val) for val in container['list_sceneLon']])

# sum subpulses for UF like modes
RSR1 = sum(float(val) for val in container['list_rangeRSR'])
# alternative: calculate RSR from given pixel spacing, former way doesn't work correctly for reduced resolution XF images. Difference is a factor 1.0000001
container['rangeRSR'] = SOL/float(container['imagePixelSpacing'])/2
# for backwards compatibility use first method when values are very close to each other
if 0.9999 < RSR1/container['rangeRSR'] < 1.0001:
    container['rangeRSR'] = RSR1

# Calculate PRF
azimuthTimeFirstLine = datetime.strptime(container['azimuthTimeFirstLine'], dateformat)
azimuthTimeLastLine = datetime.strptime(container['azimuthTimeLastLine'], dateformat)
obs_time = (azimuthTimeLastLine - azimuthTimeFirstLine).total_seconds()
# set start time to the first observed line (for ascending the image is flipped)
if obs_time > 0:
    azimuthTimeStart = azimuthTimeFirstLine
else:
    azimuthTimeStart = azimuthTimeLastLine
    obs_time = -obs_time

if container['sceneBeam'] != 'S3': # Hacky fix for S3 merged images
    container['azimuthPRF'] = (float(container['imageLines']) - 1)/obs_time

# ---------------------------------------------------------------------------------------------------------

#print(container['mission'])
#exit()

dummyVar = 'DUMMY'

print('\nrs2_dump_header2doris.py v%s, doris software, 2013\n' % codeRevision)
print('*******************************************************************')
print('*_Start_readfiles:')
print('*******************************************************************')
print('Volume file: 					%s' % 'product.xml') # container['volFile']) # HARDCODED!!! for Radarsat-2
print('Volume_ID: 					%s' % container['volID'])
print('Volume_identifier: 				%s' % container['volRef'])
print('Volume_set_identifier: 				%s' % dummyVar)
print('(Check)Number of records in ref. file: 		%s' % container['sceneRecords'])
print('SAR_PROCESSOR:                                  %s %s' % (str.split(container['productSpec'])[0][:2],container['productSoftVer']))
print('SWATH:                                          %s' % container['sceneBeam'])
print('PASS:                                           %s' % container['orbitDir'])
print('IMAGING_MODE:                                   %s %s' % (container['sceneBeamMode'],container['scenePol']))
print('RADAR_FREQUENCY (Hz):                           %s' % container['radarfreq'])
print('')
print('Product type specifier: 	                %s' % container['mission'])
print('Logical volume generating facility: 		%s' % container['productFacility'])
print('Logical volume creation date: 			%s' % container['productVolDate'])
print('Location and date/time of product creation: 	%s' % container['productDate'])
#print('Scene identification: 				Orbit: %s %s Mode: %s' % (container['orbitABS'].split('_')[0],container['orbitDir'],container['sceneBeamMode']))
print('Scene identification: 				Orbit: %s  %s' % (container['orbitABS'].split('_')[0], azimuthTimeStart.strftime(dateformat)))
print('Scene location: 		                lat: %.4f lon: %.4f' % (float(container['sceneCenLat']),float(container['sceneCenLon'])))
print('')
print('Leader file:                                 	%s' % 'product.xml') # container['volFile']) # HARDCODED!!! for Radarsat-2
print('Sensor platform mission identifer:         	%s' % container['mission'])
print('Scene_centre_latitude:                     	%s' % container['sceneCenLat'])
print('Scene_centre_longitude:                    	%s' % container['sceneCenLon'])
print('Scene_centre_heading:                            %s' % 'Null') # needs to be computed from geoinfo
print('Radar_wavelength (m):                      	%s' % str(SOL/float(container['radarfreq'])))
print('First_pixel_azimuth_time (UTC):			%s' % azimuthTimeStart.strftime('%d-%b-%Y %H:%M:%S.%f'))
print('Pulse_Repetition_Frequency (computed, Hz): 	%s' % container['azimuthPRF'])
print('Total_azimuth_band_width (Hz):             	%s' % float(container['azimuthBW']))
print('Weighting_azimuth:                         	%s %f' % (str.upper(container['azimuthWind']), float(container['azimuthWindCoeff'])))
print('Xtrack_f_DC_constant (Hz, early edge):     	%s' % container['dopplerCoeff'][0])
print('Xtrack_f_DC_linear (Hz/s, early edge):     	%s' % container['dopplerCoeff'][1])
print('Xtrack_f_DC_quadratic (Hz/s/s, early edge): 	%s' % container['dopplerCoeff'][2])
print('Range_time_to_first_pixel (2way) (ms):     	%0.15f' % (float(container['rangeTimePix'])*1000))
print('Range_sampling_rate (computed, MHz):       	%0.6f' % (float(container['rangeRSR'])/1000000))
print('Total_range_band_width (MHz):               	%s' % (float(container['rangeBW'])/1000000))
print('Weighting_range:                            	%s %f' % (str.upper(container['rangeWind']), float(container['rangeWindCoeff'])))
print('')
print('*******************************************************************')
print('Datafile: 					%s' % container['imageFile'])
print('Dataformat: 				%s' % 'GeoTIFF')  # hardcoded!!!
print('Number_of_lines_original: 			%s' % container['imageLines'])
print('Number_of_pixels_original: 	                %s' % container['imagePixels'])
print('*******************************************************************')
print('* End_readfiles:_NORMAL')
print('*******************************************************************')
print('')
print('')
print('*******************************************************************')
print('*_Start_leader_datapoints:  %s ' % container['orbitABS'].split('_')[1])
print('*******************************************************************')
print(' t(s)		X(m)		Y(m)		Z(m)      X_V(m/s)      Y_V(m/s)      Z_V(m/s)')
print('NUMBER_OF_DATAPOINTS: %s' % len(container['list_orbitTime']))

# MA : positions and velocities
for i in range(len(container['list_orbitTime'])):
    print(' %.6f %s %s %s %s %s %s' % (hms2sec(container['list_orbitTime'][i].split('T')[1].strip('Z')),
                                  container['list_orbitX'][i],
                                  container['list_orbitY'][i],
                                  container['list_orbitZ'][i],
                                  container['list_orbitXV'][i],
                                  container['list_orbitYV'][i],
                                  container['list_orbitZV'][i]))

print('')
print('*******************************************************************')
print('* End_leader_datapoints:_NORMAL')
print('*******************************************************************')
# print('\n')
# print('    Current time: %s\n' % time.asctime())
# print('\n')

#EOF
