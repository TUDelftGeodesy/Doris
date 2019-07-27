#!/usr/bin/env python

#-----------------------------------------------------------------#
# A python code for parsing CSK HDF5 file into python data structures
# and from there into DORIS res file structure
#
# Author: TUDelft - 2010
# Maintainer: Prabu Dheenathayalan
#
#-----------------------------------------------------------------#

import numpy, h5py, sys, math, time, string # numpy and h5py  required for HDF5 python support

codeRevision=1.0   # this code revision number

def usage():
    print('INFO    : @(#)Doris InSAR software, $Revision: %s $, $Author: TUDelft $' % codeRevision)
    print()
    print('Usage   : python csk_dump_header2doris.py csk_HDF5_product > OutputFileName')
    print('                               where csk_HDF5_product is the input filename')
    print()
    print('          This software is part of Doris InSAR software package.\n')
    print('(c) 1999-2010 Delft University of Technology, the Netherlands.\n')

try:
  inputFileName  = sys.argv[1]
except:
    print('\nError   : Unrecognized input or missing arguments\n\n')
    usage()
    sys.exit(1)

# accessing the HDF5 product file
f = h5py.File(inputFileName, 'r')
s01 = f.get('/S01')
sbi = s01.get('SBI')
b001 = s01.get('B001')
qlk = s01.get('QLK')


if f.parent.__contains__('/') == False or f.parent.__contains__('/S01') == False or f.parent.__contains__('/S01/SBI') == False or f.parent.__contains__('/S01/B001') == False or f.parent.__contains__('/S01/QLK') == False :
   print('ERROR: Wrong HDF5 format!')

# reading the attributes  
VolumeFile = f.attrs.__getitem__('Product Filename')
Volume_ID = f.attrs.__getitem__('Programmed Image ID') 
Volume_identifier = f.attrs.__getitem__('Product Specification Document') 
Volume_set_identifier='dummy'
NumberOfRecordsInRefFile='dummy'

SAR_PROCESSOR = f.attrs.__getitem__('L1A Software Version')
ProductTypeSpecifier = f.attrs.__getitem__('Product Type')
LogicalVolumeGeneratingFacility = f.attrs.__getitem__('Processing Centre')
LogicalVolumeCreationDate='dummy'
LocationAndDateTimeOfProductCreation = f.attrs.__getitem__('Product Generation UTC') 

SceneIdentification = f.attrs.__getitem__('Orbit Number')
orbitDir = f.attrs.__getitem__('Orbit Direction')
sceneMode = f.attrs.__getitem__('Acquisition Mode')
SceneLocation = 'dummy'
LeaderFile = f.attrs.__getitem__('Product Filename')
SAT_ID=f.attrs.__getitem__('Satellite ID')
SensorPlatformMissionIdentifier = f.attrs.__getitem__('Mission ID')
SceneCentreGeodeticCoordinates = f.attrs.__getitem__('Scene Centre Geodetic Coordinates')

Scene_centre_latitude=SceneCentreGeodeticCoordinates[0]
Scene_centre_longitude=SceneCentreGeodeticCoordinates[1]
Scene_centre_heading = f.attrs.__getitem__('Scene Orientation')
Radar_wavelength = f.attrs.__getitem__('Radar Wavelength')
Radar_frequency = f.attrs.__getitem__('Radar Frequency')
Swath_no = f.attrs.__getitem__('Subswaths Number')
Polarisation = s01.attrs.__getitem__('Polarisation')
ReferenceUTC = f.attrs.__getitem__('Reference UTC')
relZDAFTime = sbi.attrs.__getitem__('Zero Doppler Azimuth First Time') # precision of this variable need to be checked w.r.t tsx

if string.atoi(ReferenceUTC[20:])*10 != 0:
   ZDAFTime_msecsOfDay = round((string.atoi(ReferenceUTC[11:13])*60*60+string.atoi(ReferenceUTC[14:16])*60+string.atoi(ReferenceUTC[17:19])+pow(string.atoi(ReferenceUTC[20:])*10,(-(len(ReferenceUTC)-20)))+relZDAFTime)*1000)
else:
  ZDAFTime_msecsOfDay = round((string.atoi(ReferenceUTC[11:13])*60*60+string.atoi(ReferenceUTC[14:16])*60+string.atoi(ReferenceUTC[17:19])+relZDAFTime)*1000)

ZDAFTime_HH = math.floor(ZDAFTime_msecsOfDay/(60*60*1000))
ZDAFTime_MM = math.floor(ZDAFTime_msecsOfDay % (60*60*1000)/(60*1000))
ZDAFTime_SS = math.floor(ZDAFTime_msecsOfDay%(60*1000)/1000)
ZDAFTime_TTT = round(ZDAFTime_msecsOfDay%1000)

First_pixel_azimuth_time = str('%s %s:%s:%s.%s' % (time.strftime("%d-%b-%Y",time.strptime(ReferenceUTC.split()[0],"%Y-%m-%d")),str(int(ZDAFTime_HH)),str(int(ZDAFTime_MM)),str(int(ZDAFTime_SS)),str(int(ZDAFTime_TTT))))
PRF = s01.attrs.__getitem__('PRF')
Total_azimuth_band_width = s01.attrs.__getitem__('Azimuth Focusing Bandwidth')
Weighting_azimuth = f.attrs.__getitem__('Azimuth Focusing Weighting Function')

Xtrack_f_DC = f.attrs.__getitem__('Centroid vs Range Time Polynomial')
Xtrack_f_DC_constant = Xtrack_f_DC[0]
Xtrack_f_DC_linear = Xtrack_f_DC[1]
Xtrack_f_DC_quadratic = Xtrack_f_DC[2]
Range_time_to_first_pixel = sbi.attrs.__getitem__('Zero Doppler Range First Time')*1000

RSR = s01.attrs.__getitem__('Sampling Rate')/pow(10, 6)
Total_range_band_width = s01.attrs.__getitem__('Range Focusing Bandwidth')/pow(10, 6)
Weighting_range = f.attrs.__getitem__('Range Focusing Weighting Function')
NUMBER_OF_DATAPOINTS = f.attrs.__getitem__('Number of State Vectors')
StateVectorsTimes = f.attrs.__getitem__('State Vectors Times')

if string.atoi(ReferenceUTC[20:])*10 != 0:
   Time_Datapoints = (string.atoi(ReferenceUTC[11:13])*60*60+string.atoi(ReferenceUTC[14:16])*60+string.atoi(ReferenceUTC[17:19])+pow(string.atoi(ReferenceUTC[20:])*10,(-(len(ReferenceUTC)-20)))+StateVectorsTimes) % (60*60*24)
else:
  Time_Datapoints = (string.atoi(ReferenceUTC[11:13])*60*60+string.atoi(ReferenceUTC[14:16])*60+string.atoi(ReferenceUTC[17:19])+StateVectorsTimes) % (60*60*24)

ECEFSatellitePosition = f.attrs.__getitem__('ECEF Satellite Position')
ECEFSatelliteVelocity = f.attrs.__getitem__('ECEF Satellite Velocity')
nrows = len(ECEFSatellitePosition[:, 0])
ncols = len(ECEFSatellitePosition[0, :])
if nrows == 3:
   X_Datapoints=ECEFSatellitePosition[0,:]
   Y_Datapoints=ECEFSatellitePosition[1,:]
   Z_Datapoints=ECEFSatellitePosition[2,:]
   X_Velocity=ECEFSatelliteVelocity[0,:]
   Y_Velocity=ECEFSatelliteVelocity[1,:]
   Z_Velocity=ECEFSatelliteVelocity[2,:]
else:
   X_Datapoints=ECEFSatellitePosition[:,0]
   Y_Datapoints=ECEFSatellitePosition[:,1]
   Z_Datapoints=ECEFSatellitePosition[:,2]
   X_Velocity=ECEFSatelliteVelocity[:,0]
   Y_Velocity=ECEFSatelliteVelocity[:,1]
   Z_Velocity=ECEFSatelliteVelocity[:,2]

Datafile = f.attrs.__getitem__('Product Filename')


if sbi.shape[0]== 2:
   Number_of_lines_original = sbi.shape[2]
   Number_of_pixels_original = sbi.shape[1]
else:
   Number_of_lines_original = sbi.shape[0]
   Number_of_pixels_original = sbi.shape[1]


# ---------------------------------------------------------------------------------------------------------
# print the extracted attributes

dummyVar = 'DUMMY'

print('\ncsk_dump_header2doris.py v%s, doris software, 2010\n' % codeRevision)
print('*******************************************************************')
print('*_Start_readfiles:')
print('*******************************************************************')
print('Volume file:                                     %s' % VolumeFile)
print('Volume_ID:                                       %s' % Volume_ID)
print('Volume_identifier:                               %s' % Volume_identifier)
print('Volume_set_identifier:                           %s' % Volume_set_identifier)
print('(Check)Number of records in ref. file:           %s' % NumberOfRecordsInRefFile)
print('SAR_PROCESSOR:                                   %s L1A %s' % (SensorPlatformMissionIdentifier,SAR_PROCESSOR))
print('SWATH:                                           %s' % Swath_no)
print('PASS:                                            %s' % orbitDir)
print('IMAGING_MODE:                                    %s %s' % (sceneMode,Polarisation))
print('RADAR_FREQUENCY (Hz):                            %s' % Radar_frequency)
print('')
#print('Product type specifier:                          %s' % ProductTypeSpecifier) # returns SCS_B
print('Product type specifier:                          %s %s' % (SensorPlatformMissionIdentifier, ProductTypeSpecifier))
print('Logical volume generating facility:              %s' % LogicalVolumeGeneratingFacility)
print('Logical volume creation date:                    %s' % LogicalVolumeCreationDate)
print('Location and date/time of product creation:      %s' % LocationAndDateTimeOfProductCreation)

#print('Scene identification:                            Orbit: %s %s Mode: %s' % (SceneIdentification,orbitDir,sceneMode))
print('Scene identification:                            Orbit: %s' % SceneIdentification)
print('Scene location:                                  lat: %.4f lon: %.4f' % (float(Scene_centre_latitude),float(Scene_centre_longitude)))
print('')
print('Leader file:                                     %s' % LeaderFile)
#print('Sensor platform mission identifier:               %s' % SensorPlatformMissionIdentifier)
print('Sensor platform mission identifier:              %s' % SAT_ID)
print('Scene_centre_latitude:                           %s' % Scene_centre_latitude)
print('Scene_centre_longitude:                          %s' % Scene_centre_longitude)
print('Scene_centre_heading:                            %s' % Scene_centre_heading)
#print('Radar_wavelength (m):                            0.031') #HARDCODED!!!
print('Radar_wavelength (m):                            %s' % Radar_wavelength) 
print('First_pixel_azimuth_time (UTC):                  %s' % First_pixel_azimuth_time)
print('Pulse_Repetition_Frequency (computed, Hz):       %s' % PRF)
print('Total_azimuth_band_width (Hz):                   %s' % Total_azimuth_band_width)
print('Weighting_azimuth:                               %s' % string.upper(Weighting_azimuth))
print('Xtrack_f_DC_constant (Hz, early edge):           %s' % Xtrack_f_DC_constant)
print('Xtrack_f_DC_linear (Hz/s, early edge):           %s' % Xtrack_f_DC_linear)
print('Xtrack_f_DC_quadratic (Hz/s/s, early edge):      %s' % Xtrack_f_DC_quadratic)
print('Range_time_to_first_pixel (2way) (ms):           %0.15f' % (float(Range_time_to_first_pixel)))
print('Range_sampling_rate (computed, MHz):             %0.6f' % (float(RSR)))
print('Total_range_band_width (MHz):                    %s' % (float(Total_range_band_width)))
print('Weighting_range:                                 %s' % string.upper(Weighting_range))
print('')
print('*******************************************************************')
print('Datafile:                                        %s' % Datafile)
print('Dataformat:                                      %s' % 'HDF5')     # hardcoded!!!
print('Number_of_lines_original:                        %s' % str(Number_of_lines_original))
print('Number_of_pixels_original:                       %s' % str(Number_of_pixels_original))
print('*******************************************************************')
print('* End_readfiles:_NORMAL')
print('*******************************************************************')
print('')
print('')
print('*******************************************************************')
print('*_Start_leader_datapoints')
print('*******************************************************************')
print(' t(s)            X(m)            Y(m)            Z(m)            X_V(m/s)            Y_V(m/s)            Z_V(m/s)')
print('NUMBER_OF_DATAPOINTS:                    %s' % str(NUMBER_OF_DATAPOINTS))
print('')

for i in range(NUMBER_OF_DATAPOINTS):
    print(' %s %s %s %s %s %s %s' % (int(Time_Datapoints[i]),\
                                  X_Datapoints[i],\
                                  Y_Datapoints[i],\
                                  Z_Datapoints[i],\
				  X_Velocity[i],\
                                  Y_Velocity[i],\
                                  Z_Velocity[i]))

print('')
print('*******************************************************************')
print('* End_leader_datapoints:_NORMAL')
print('*******************************************************************')

#EOF
