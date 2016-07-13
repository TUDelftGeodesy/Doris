#!/usr/bin/env python
import os,sys,time
import math
import numpy as np
from numpy import *
import shutil
import re
import scipy as Sci
import scipy.linalg
import scipy.io as sio
from   scipy import interpolate
import struct

from scipy.interpolate import RectBivariateSpline
#import matplotlib.pyplot as plt

import gdal,gdalconst
from gdalconst import *
#from osgeo import gdal,gdalconst
# NOT FETOOLS



################################################################################
def get_parameter(First_param,file_name,format_flag=1):
    Read_contine_flag=0
    orbit_info=""
    value=None
   
    for line in open(file_name):
        if format_flag==1:
            if not (line.find(First_param)):
                index=line.find(' ')
                value=(line[(index+1):].strip(' \n\t')).strip(' ')
                value=value.split(' ')[0]
                return value 
        
        if format_flag==2:
            if not (line.find(First_param)):
                index=line.find(' ')
                value=(line[(index+1):].strip(' \n\t')).strip().split()
                image_line_direction=value[0]
                image_pixel_direction=value[1]
                return image_line_direction,image_pixel_direction 
#------------------------------------------------------------------------------------
def parseFile(filename,dem_x,dem_y):
    # read 1,442,401 (1201x1201) high-endian
    # signed 16-bit words into self.z
    total_size=int64(dem_x)*int64(dem_y)
    #print 'total_size=',total_size
    fi=open(filename,"rb")
    contents=fi.read()
    fi.close()
    format='<'+str(total_size)+'f'
    print 'format=',format 
    dem=struct.unpack(format, contents)
    return dem 


#--------------------------------------------------------------------------------------
def writeFile(filename,dem_x,dem_y,contents):
    # read 1,442,401 (1201x1201) high-endian
    # signed 16-bit words into self.z
    total_size=int64(dem_x)*int64(dem_y)
    print 'total_size=',total_size    
    
    format='f'* (total_size)
    #print 'format=',format 
    dem=struct.pack(format, contents)
    fi=open(filename,"wb")
    fi.write(dem)
    fi.close()
    return True                
#-------------------------------------------------------------------------------------

def freadbk_hgt(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('SRTMHGT')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print 'Could not open'
        sys.exit(1)
    print 'Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
          thisBurstData_file.GetDriver().LongName
    print 'Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
          'x',thisBurstData_file.RasterCount
    print 'Projection is ',thisBurstData_file.GetProjection()
    geotransform = thisBurstData_file.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

    cint_srd=thisBurstData_file.GetRasterBand(1)
    #print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    #print 'thisBurstData=',thisBurstData[130:140,130:140]
    return thisBurstData
###############################################################################
def freadbk(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print 'Could not open'+Path_MFF_HDR
        sys.exit(1)
    print 'Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
          thisBurstData_file.GetDriver().LongName
    print 'Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
          'x',thisBurstData_file.RasterCount
    print 'Projection is ',thisBurstData_file.GetProjection()
    geotransform = thisBurstData_file.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

    cint_srd=thisBurstData_file.GetRasterBand(1)
    #print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    #print 'thisBurstData=',thisBurstData[130:140,130:140]
    return thisBurstData
###############################################################################



Input_dorisin='dem.dorisin'
dem_path = get_parameter('CRD_IN_DEM',Input_dorisin,format_flag=1)
print 'dem_path=',dem_path

dem_format= get_parameter('CRD_IN_FORMAT',Input_dorisin,format_flag=1)
print 'dem_format=',dem_format

dem_size_x,dem_size_y= get_parameter('CRD_IN_SIZE',Input_dorisin,format_flag=2)
print 'dem_size_x,dem_size_y',dem_size_x,dem_size_y

dem_DELTA_x,dem_DELTA_y = get_parameter('CRD_IN_DELTA',Input_dorisin,format_flag=2)
print 'dem_DELTA_x,dem_DELTA_y',dem_DELTA_x,dem_DELTA_y

dem_Upper,dem_Left= get_parameter('CRD_IN_UL',Input_dorisin,format_flag=2)
print 'ddem_Upper,dem_Left=',dem_Upper,dem_Left



############################
#    read/operate data   #
############################
dem_path_orig=dem_path+'.orig'

if os.path.exists(dem_path):
    shutil.move(dem_path,dem_path_orig)  #move back to start again

lon=float(dem_Left)  +np.arange(0,int(dem_size_y))*float(dem_DELTA_y)
lat=float(dem_Upper) -np.arange(0,int(dem_size_x))*float(dem_DELTA_x)
lat=lat[::-1]
print 'lon=',lon
print 'lat=',lat

geoidegm_file='/media/sdb2/wuw/DORIS_TEST/netherlands/aero/'+'geoidegm96grid.mat'

geoidegm=scipy.io.loadmat(geoidegm_file)


f = interpolate.RectBivariateSpline(geoidegm['latbp'][0,:],geoidegm['lonbp'][0,:],array(geoidegm['grid']),kx=1, ky=1 )

#print 'lat=',lat.shape
#print 'lon=',lon.shape
egm96 = f( lat, lon)
egm96 = flipud(egm96)

#print 'egm96.shape=',egm96.shape
Dem_Data=parseFile(dem_path_orig,int64(dem_size_x),int64(dem_size_y))

total_size=int64(dem_size_x) * int64(dem_size_y)

egm96=np.reshape(egm96,total_size)

Dem_Data=Dem_Data+egm96
del egm96

Dem_Data.astype('float32').tofile(dem_path)
