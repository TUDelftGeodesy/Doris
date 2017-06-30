import os
import numpy as np
import gdal
from gdalconst import *
import sys

def get_ramp(res_file, resampled=0, type='chirp'):
    # Read information
    ################################################################################

    #FM
    t0_FM = np.float64(get_parameter('FM_reference_range_time', res_file,1))
    c0_FM = np.float64(get_parameter('FM_polynomial_constant_coeff (Hz, early edge)', res_file,1))
    c1_FM = np.float64(get_parameter('FM_polynomial_linear_coeff (Hz/s, early edge)', res_file,1))
    c2_FM = np.float64(get_parameter('FM_polynomial_quadratic_coeff (Hz/s/s, early edge)', res_file,1))

    #DC
    azimuthTime_DC = get_parameter('DC_reference_azimuth_time', res_file,3)
    azimuthTime_DC = np.float64(azimuthTime_DC[0])*3600+float(azimuthTime_DC[1])*60+float(azimuthTime_DC[2])

    t0_DC = np.float64(get_parameter('DC_reference_range_time',res_file,1))
    c0_DC = np.float64(get_parameter('Xtrack_f_DC_constant (Hz, early edge)', res_file,1))
    c1_DC = np.float64(get_parameter('Xtrack_f_DC_linear (Hz/s, early edge)', res_file,1))
    c2_DC = np.float64(get_parameter('Xtrack_f_DC_quadratic (Hz/s/s, early edge)', res_file,1))

    Ks = np.float64(get_parameter('Azimuth_steering_rate (deg/s)', res_file,1))

    # Image sampling parameters
    Taz_start = get_parameter('First_pixel_azimuth_time (UTC)', res_file,3)
    Taz_start = np.float64(Taz_start[0])*3600+float(Taz_start[1])*60+float(Taz_start[2])

    Trg_start = np.float64(get_parameter('Range_time_to_first_pixel (2way) (ms)', res_file,1))*1e-3
    fsRg = np.float64(get_parameter('Range_sampling_rate (computed, MHz)', res_file,1))

    dt_az = np.float64(get_parameter('Azimuth_time_interval (s)', res_file,1))
    dt_rg = 1/fsRg/1e6

    # Number of lines
    lNum = int(get_parameter('Number_of_lines_original', res_file,1))

    if resampled == 1:
        l0 = int(get_parameter('First_line (w.r.t. original_master)', res_file,2,'*_Start_resample','* End_resample:_NORMAL'))
        lN = int(get_parameter('Last_line (w.r.t. original_master)', res_file,2,'*_Start_resample','* End_resample:_NORMAL'))
        p0 = int(get_parameter('First_pixel (w.r.t. original_master)', res_file,2,'*_Start_resample','* End_resample:_NORMAL'))
        pN = int(get_parameter('Last_pixel (w.r.t. original_master)', res_file,2,'*_Start_resample','* End_resample:_NORMAL'))
    else:
        l0 = int(get_parameter('First_line (w.r.t. original_image)', res_file, 1))
        lN = int(get_parameter('Last_line (w.r.t. original_image)', res_file, 1))
        p0 = int(get_parameter('First_pixel (w.r.t. original_image)', res_file, 1))
        pN = int(get_parameter('Last_pixel (w.r.t. original_image)', res_file, 1))

    # Get resampled Slv size
    Naz_res = lN-l0+1
    Nrg_res = pN-p0+1

    if resampled == 1:
        # Read the resampled image and slave coordinates in master geometry
        #################################################################################

        Path_MFF_HDR   ='rsmp_orig_slave_pixel'+'.hdr'
        Link_DATA      ='rsmp_orig_slave_pixel'+'.r00'  # the default format should be r00
        Link_rsmp_orig_slave_pixel ='rsmp_orig_slave_pixel.raw'

        if (os.path.isfile(Path_MFF_HDR)):
            os.remove(Path_MFF_HDR)
        if (os.path.isfile(Link_DATA)):
            os.remove(Link_DATA)

        RAW_DATA_ABSOLUTE_PATH=os.path.abspath(Link_rsmp_orig_slave_pixel)
        print "RAW_DATA_ABSOLUTE_PATH=", RAW_DATA_ABSOLUTE_PATH
        os.symlink(RAW_DATA_ABSOLUTE_PATH,Link_DATA)

        outStream      = open(Path_MFF_HDR,'w')
        outStream.write('IMAGE_FILE_FORMAT = MFF\n')
        outStream.write('FILE_TYPE = IMAGE\n')
        outStream.write('IMAGE_LINES = %d\n' % int(Naz_res))
        outStream.write('LINE_SAMPLES = %d\n'% int(Nrg_res))
        outStream.write('BYTE_ORDER = LSB\n')
        outStream.write('END\n')
        outStream.close()

        PixRgGrid = freadbk(Path_MFF_HDR,1, 1,int(Naz_res),int(Nrg_res))
        PixRgGrid = PixRgGrid.astype(np.float64)

        if (os.path.isfile(Path_MFF_HDR)):
            os.remove(Path_MFF_HDR)
        if (os.path.isfile(Link_DATA)):
            os.remove(Link_DATA)
        #################################################################################

        Path_MFF_HDR   ='rsmp_orig_slave_line'+'.hdr'
        Link_DATA      ='rsmp_orig_slave_line'+'.r00'
        Link_rsmp_orig_slave_line ='rsmp_orig_slave_line.raw'

        if (os.path.isfile(Path_MFF_HDR)):
            os.remove(Path_MFF_HDR)
        if (os.path.isfile(Link_DATA)):
            os.remove(Link_DATA)


        RAW_DATA_ABSOLUTE_PATH=os.path.abspath(Link_rsmp_orig_slave_line)
        print "RAW_DATA_ABSOLUTE_PATH=", RAW_DATA_ABSOLUTE_PATH
        os.symlink(RAW_DATA_ABSOLUTE_PATH,Link_DATA)

        outStream      = open(Path_MFF_HDR,'w')
        outStream.write('IMAGE_FILE_FORMAT = MFF\n')
        outStream.write('FILE_TYPE = IMAGE\n')
        outStream.write('IMAGE_LINES = %d\n' % int(Naz_res))
        outStream.write('LINE_SAMPLES = %d\n'% int(Nrg_res))
        outStream.write('BYTE_ORDER = LSB\n')
        outStream.write('END\n')
        outStream.close()

        PixAzGrid = freadbk(Path_MFF_HDR,1, 1,int(Naz_res),int(Nrg_res))
        PixAzGrid=PixAzGrid.astype(np.float64)

        if (os.path.isfile(Path_MFF_HDR)):
            os.remove(Path_MFF_HDR)
        if (os.path.isfile(Link_DATA)):
            os.remove(Link_DATA)

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #% Prepare azimuth and range grids %
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TrgGrid = Trg_start + (PixRgGrid-1) * dt_rg
        TazGrid = (PixAzGrid-1) * dt_az - (lNum/2 * dt_az)

        del PixAzGrid, PixRgGrid

    elif resampled == 0:
        Tvect_rg = Trg_start + np.arange(p0-1,pN) * dt_rg
        Tvect_az = np.arange(l0-1,lN) * dt_az - (lNum/2 * dt_az)
        Tvect_az = Tvect_az[:, None]

        TrgGrid = np.tile(Tvect_rg, (Naz_res, 1))
        TazGrid = np.tile(Tvect_az, (1, Nrg_res))
        
    else:
        print 'variable resampled can only be 0 or 1!'
        return

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% From S-1 steering rate and orbit information %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%%%%%%%%%%%%%%%%%%%%%
    #% Orbit information %
    #%%%%%%%%%%%%%%%%%%%%%

    orbit_number     = int(get_parameter('NUMBER_OF_DATAPOINTS',res_file,1))
    orbit_time       = np.zeros(orbit_number,dtype=np.float64)
    orbit_velocity_x = np.zeros(orbit_number,dtype=np.float64)
    orbit_velocity_y = np.zeros(orbit_number,dtype=np.float64)
    orbit_velocity_z = np.zeros(orbit_number,dtype=np.float64)
    orbit_info       = get_parameter('NUMBER_OF_DATAPOINTS',res_file,4)

    orbit_info=orbit_info.split('\n')

    for row in range(orbit_number):
        orbit_time_position=orbit_info[row]
        orbit_time[row]=float(orbit_time_position.strip().split()[0])
        orbit_velocity_x[row]=float(orbit_time_position.strip().split()[1])
        orbit_velocity_y[row]=float(orbit_time_position.strip().split()[2])
        orbit_velocity_z[row]=float(orbit_time_position.strip().split()[3])
    orbit_velocity = np.sqrt(np.diff(orbit_velocity_x)**2+np.diff(orbit_velocity_y)**2+np.diff(orbit_velocity_z)**2)/np.diff(orbit_time)

    # Compute Nominal DC for the whole burst
    # Compute FM rate along range
    Kfm = c0_FM + c1_FM*(TrgGrid-t0_FM) + c2_FM*(TrgGrid-t0_FM)**2
    Kfm_0 = c0_FM + c1_FM*(Trg_start-t0_FM) + c2_FM*(Trg_start-t0_FM)**2

    # Compute DC along range at reference azimuth time (azimuthTime)
    Df_AzCtr = c0_DC + c1_DC*(TrgGrid-t0_DC) + c2_DC*(TrgGrid-t0_DC)**2
    f_DC_ref_0 = c0_DC + c1_DC*(Trg_start-t0_DC) + c2_DC*(Trg_start-t0_DC)**2
    del TrgGrid

    # From S-1 steering rate and orbit information %
    # Computes sensor velocity from orbits
    C_lambda=np.float64(get_parameter('Radar_wavelength (m)',res_file,1))

    # Frequency rate
    Ks_hz = 2*np.mean(orbit_velocity)/C_lambda*Ks/180*np.pi

    # Time ratio
    alpha_nom = 1 - Ks_hz/Kfm

    # DC Azimuth rate [Hz/s]
    DR_est = Ks_hz/alpha_nom
    del Ks_hz, alpha_nom

    # Reference time
    az_DC = -(Df_AzCtr / Kfm) + (f_DC_ref_0 / Kfm_0)
    del Kfm, Kfm_0
    Taz_vec = TazGrid - az_DC
    del az_DC

    #% Generate inverse chirp %
    if type == 'chirp':
        data = np.exp(1j*2*np.pi*(DR_est/2*Taz_vec+Df_AzCtr)*Taz_vec)
    elif type == 'DC':
        data = Df_AzCtr + Taz_vec * DR_est 
    else:
        print 'Choose either chirp or DC for type'
        return

    return data


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
################################################################################


###############################################################################

def freadbk(path_file,line_start=1, pixels_start=1,nofLines1=None,nofPixels1=None):
    #Driver
    driver=gdal.GetDriverByName('MFF')
    driver.Register()
    gdal.AllRegister()
    thisBurstData_file=gdal.Open(path_file,GA_ReadOnly)
    if thisBurstData_file is None:
        print 'Could not open'+Path_MFF_HDR
        sys.exit(1)
    #print 'Driver: ', thisBurstData_file.GetDriver().ShortName,'/', \
    #      thisBurstData_file.GetDriver().LongName
    #print 'Size is ',thisBurstData_file.RasterXSize,'x',thisBurstData_file.RasterYSize, \
    #      'x',thisBurstData_file.RasterCount
    #print 'Projection is ',thisBurstData_file.GetProjection()
    geotransform = thisBurstData_file.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

    cint_srd=thisBurstData_file.GetRasterBand(1)
    #print 'Band Type=',gdal.GetDataTypeName(cint_srd.DataType)

    if cint_srd.GetOverviewCount() > 0:
            print 'Band has ', cint_srd.GetOverviewCount(), ' overviews.'
    thisBurstData= cint_srd.ReadAsArray(int(pixels_start-1),int(line_start-1),nofPixels1,nofLines1)
    return thisBurstData
##################################################################################
