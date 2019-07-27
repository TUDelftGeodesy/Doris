#!/usr/bin/env python
import numpy as np
from numpy import *
import gdal
from gdalconst import *
from scipy import ndimage

def usage():
    print('\nUsage: python get_winpos.py dataFile resFile Nwin outFile                        ')
    print('where   dataFile           is the name of burst you want to deramp                 ')
    print('        resFile            is the .res file of burst                               ')
    print('        Nwin               number of windows to be distributed over the total image')
    print('        outFile            output file name                                        ')
    print('  for example                                                                      ')
    print(' python get_winpos.py 20141003_iw_1_burst_1.raw 20141003_iw_1_burst_1.res 2001 winpos_fine.asc')
    print(' matlab: TU Delft                                                                  ')
    print(' Python: Wu Wenhao   Wuhan QQ:460249274                                            ')
try:
    dataFile   = sys.argv[1]
    resFile    = sys.argv[2]
    Nwin       = sys.argv[3]
    outFile    = sys.argv[4]

except:
    print('Unrecognized input')
    usage()
    sys.exit(1)


################################################################################
def get_parameter(First_param,file_name,format_flag=1,Second_param=None,Third_param=None):
    Read_contine_flag=0
    orbit_info=""
    value=None
    for line in open(file_name):
        if format_flag==1:
            if not (line.find(First_param)):
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                return value

        if format_flag==2:
            if not(line.find(Second_param)):
                Read_contine_flag=1
            if (Read_contine_flag==1) and (not (line.find(First_param))):  #Be careful
                index=line.find(':')
                value=(line[(index+1):].strip(' \n\t'))
                continue
            if Read_contine_flag==1  and (not(line.find(Third_param))):  #Be careful
                Read_contine_flag=0
                return value


        if format_flag==3:
            if not (line.find(First_param)):
                index=line.find(':')
                pixel_time=(line[(index+1):].strip(' \n\t')).split(' ')[1].split(':')
                return pixel_time


        if format_flag==4:

            if not (line.find(First_param)):
                index=line.find(':')
                value=int(line[(index+1):].strip(' \n\t'))
                Read_contine_flag=1
                continue
            if (Read_contine_flag>=1):
                orbit_info=orbit_info+line
                Read_contine_flag=Read_contine_flag+1
                if (Read_contine_flag==(value+1)):
                    return orbit_info
###############################################################################
#thisBurstData = freadbk(['burst' num2str(nBurst)   '/cint_srd.raw'],nofLines1,formatData1, line1, nofLines1,1,nofPixels1);
def freadbk(path_file,line_start=1, Pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print('Could not open'+Path_MFF_HDR)
        sys.exit(1)
    #print('Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
    #      thisBurstData_file.GetDriver().LongName)
    #print('Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
    #      'x',thisBurstData_file.RasterCount)
    #print('Projection is ',thisBurstData_file.GetProjection())
    geotransform = thisBurstData_file.GetGeoTransform()
    #if not geotransform is None:
    #    print('Origin = (',geotransform[0], ',',geotransform[3],')')
    #    print('Pixel Size = (',geotransform[1], ',',geotransform[5],')')

    cint_srd=thisBurstData_file.GetRasterBand(1)
    #print('Band Type=',gdal.GetDataTypeName(cint_srd.DataType))

    if cint_srd.GetOverviewCount() > 0:
            print('Band has ', cint_srd.GetOverviewCount(), ' overviews.')
    thisBurstData= cint_srd.ReadAsArray(int(Pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
###############################################################################


NwinGrid = 5

azSpacing = 20
rSpacing  = 5


resData   = resFile

if get_parameter('First_line (w.r.t. ovs_image)',resData,1): #%oversampled data
  #% First_line (w.r.t. ovs_image):
    l0 = int(get_parameter('First_line (w.r.t. ovs_image)',resData,1));
  #% Last_line (w.r.t. ovs_image):
    lN = int(get_parameter('Last_line (w.r.t. ovs_image)',resData,1));
  #% First_pixel (w.r.t. ovs_image):
    p0 = int(get_parameter('First_pixel (w.r.t. ovs_image)',resData,1));
  #% Last_pixel (w.r.t. ovs_image):
    pN = int(get_parameter('Last_pixel (w.r.t. ovs_image)',resData,1));

    dataFormat = 'cpxfloat32'

else:#original data
    l0 = int(get_parameter('First_line (w.r.t. original_image)',resData,1));
    #% Last_line (w.r.t. original_image):
    lN = int(get_parameter('Last_line (w.r.t. original_image)',resData,1))
    #% First_pixel (w.r.t. original_image):
    p0 = int(get_parameter('First_pixel (w.r.t. original_image)',resData,1))
    #% Last_pixel (w.r.t. original_image):
    pN = int(get_parameter('Last_pixel (w.r.t. original_image)',resData,1))
    dataFormat = 'cpxint16'


# Image size
Nlines = lN-l0+1;
Npixels = pN-p0+1;
print("Nlines =",Nlines)
print("Npixels =",Npixels)


Ngrid = float(Nwin)/NwinGrid;
daz = Nlines*azSpacing;
dr = Npixels*rSpacing;

ratio = float(dr)/daz

Ngrid_az = sqrt(Ngrid/ratio);
Ngrid_r = round(Ngrid_az*ratio);
Ngrid_az = round(Ngrid_az)


Nlines_grid = ceil(Nlines/Ngrid_az)
Nlines_grid_orig = Nlines_grid


Npixels_grid = ceil(Npixels/Ngrid_r)
Npixels_grid_orig = Npixels_grid


RAW_CINT_SRD  = dataFile
Path_MFF_HDR  = dataFile.split('.')[0]+'.hdr';


if dataFormat == 'cpxint16':
   Link_CINT_SRD=dataFile.split('.')[0]+'.j00'
else:
   Link_CINT_SRD=dataFile.split('.')[0]+'.x00'


outStream      = open(Path_MFF_HDR,'w')
outStream.write('IMAGE_FILE_FORMAT = MFF\n')
outStream.write('FILE_TYPE = IMAGE\n')
outStream.write('IMAGE_LINES = %d\n' % int(Nlines))
outStream.write('LINE_SAMPLES = %d\n'% int(Npixels))
outStream.write('BYTE_ORDER = LSB\n')
outStream.write('END\n')
outStream.close()

if (os.path.exists(Link_CINT_SRD)):
    os.remove(Link_CINT_SRD)
RAW_CINT_SRD_ABSOLUTE_PATH=os.path.abspath(RAW_CINT_SRD)
print("RAW_CINT_SRD_ABSOLUTE_PATH=", RAW_CINT_SRD_ABSOLUTE_PATH)
os.symlink(RAW_CINT_SRD_ABSOLUTE_PATH,Link_CINT_SRD)


winpos=np.array([],dtype='int32').reshape(0,2)


for v in range(1,int(Ngrid_az)+1):

    if v==Ngrid_az:
       Nlines_grid = (Nlines-1)%Nlines_grid_orig+1
    else:
       Nlines_grid = Nlines_grid_orig
    ampArray = abs(freadbk(Path_MFF_HDR,int((v-1)*Nlines_grid_orig+1),1,int(Nlines_grid),Npixels ))
    for w in range(1,int(Ngrid_r)+1):
        if w==Ngrid_r:
            Npixels_grid = (Npixels-1)%Npixels_grid_orig+1
        else:
            Npixels_grid = Npixels_grid_orig
        amp = ampArray[:,(w-1)*Npixels_grid_orig:(w-1)*Npixels_grid_orig+Npixels_grid]
        locMaxsInd = amp == ndimage.grey_dilation(amp, size=(5*rSpacing, 5*azSpacing))
        locMaxs = amp[locMaxsInd]
        [az,r] = where(locMaxsInd)
        sortIdx =np.argsort(-locMaxs)
        sortIdx = sortIdx[0:NwinGrid]
        add_winpos=np.array([az[sortIdx]+l0-1+(v-1)*Nlines_grid_orig,r[sortIdx]+p0-1+(w-1)*Npixels_grid_orig]).transpose()
        winpos=np.vstack([winpos,add_winpos])

fidRes        = open(outFile,'w')
cols = winpos.shape[1]
rows = winpos.shape[0]
#print("cols = ",cols)
print("rows = ", rows)
for i_temp in range(0,rows):
    fidRes.write( '%d   %d\n' % (winpos[i_temp,0]+1,winpos[i_temp,1]+1))
fidRes.close()

if (os.path.exists(Link_CINT_SRD)):
    os.remove(Link_CINT_SRD)
if (os.path.exists(Path_MFF_HDR)):
    os.remove(Path_MFF_HDR)
