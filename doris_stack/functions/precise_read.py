# This file contains several functions to read data from precise orbit files of sentinel-1 data.
# Scripts were created by  Wu Wenhao, Wuhan university and adapted by Gert Mulder, TU Delft

import time
import calendar
import numpy as np
import os, sys
from scipy.interpolate import interp1d
import scipy.interpolate as inter
import calendar


def orbit_read(input_EOF_FileName):

    try:
        import xml.etree.cElementTree as etree
    except:
        try:
            from lxml import etree
        except:
            #import xml.etree.ElementTree as etree
            print 'Failed to load lxml.etree or xml.etree.cElementTree'
            sys.exit(1)
    #print(input_EOF_FileName)
    inTree = etree.parse(input_EOF_FileName)
    queryList = {
                # orbit inf
               'Mission'  : './/Earth_Explorer_Header/Fixed_Header/Mission',\
               'Validity_Start': './/Earth_Explorer_Header/Fixed_Header/Validity_Period/Validity_Start',\
               'Validity_Stop': './/Earth_Explorer_Header/Fixed_Header/Validity_Period/Validity_Stop',\
               'orbitABS' : './/Data_Block/List_of_OSVs/OSV/Absolute_Orbit',\
               'orbitTime': './/Data_Block/List_of_OSVs/OSV/UTC',\
               'orbitX'   : './/Data_Block/List_of_OSVs/OSV/X',\
               'orbitY'   : './/Data_Block/List_of_OSVs/OSV/Y',\
               'orbitZ'   : './/Data_Block/List_of_OSVs/OSV/Z',\
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
    return container

#--------------------------------------------------------
def interpolate_orbit(input_orbit_dir, date, input_orbit_type, input_interpolation_method, satellite='S1A'):

    orbit_time = calendar.timegm(time.strptime(date,'%Y-%m-%dT%H:%M:%S.%f'))
    date_start = calendar.timegm(time.strptime(date[:10],'%Y-%m-%d'))

    if input_orbit_type == 'POE':
        input_orbit_dir = os.path.join(input_orbit_dir, 'precise')
    elif input_orbit_type == 'RES':
        input_orbit_dir = os.path.join(input_orbit_dir, 'restituted')

    L = os.listdir(input_orbit_dir)
    Orbit_info = []

    if input_orbit_type == 'POE' and satellite == 'S1A':
        orbit_type = 'S1A_OPER_AUX_POEORB_OPOD_'
    elif input_orbit_type == 'RES' and satellite == 'S1A':
        orbit_type = 'S1A_OPER_AUX_RESORB_OPOD_'
    elif input_orbit_type == 'POE' and satellite == 'S1B':
        orbit_type = 'S1B_OPER_AUX_POEORB_OPOD_'
    elif input_orbit_type == 'RES' and satellite == 'S1B':
        orbit_type = 'S1B_OPER_AUX_RESORB_OPOD_'

    for d in L:
        if d.startswith(orbit_type):
            start_time = calendar.timegm(time.strptime(d[42:57], '%Y%m%dT%H%M%S'))
            end_time = calendar.timegm(time.strptime(d[58:73], '%Y%m%dT%H%M%S'))
            if (start_time < orbit_time) and (end_time > orbit_time):

                meta = orbit_read(os.path.join(input_orbit_dir, d))
                print(d)

                for i in range(len(meta['orbitTime'])):

                    point_time = calendar.timegm(time.strptime(meta['orbitTime'][i][4:-7], '%Y-%m-%dT%H:%M:%S'))
                    if (point_time > orbit_time-290) and (point_time < orbit_time+290):

                        Tuple_orbit=(float(hms2sec(meta['orbitTime'][i][4:].split('T')[1])),\
                                     float(meta['orbitX'][i]), float(meta['orbitY'][i]),\
                                     float(meta['orbitZ'][i]))
                        Orbit_info.append(Tuple_orbit)

    set_list=[]
    Orbit_info=sorted(Orbit_info)

    for element in range(len(Orbit_info)-1):
        temp_element     =Orbit_info[element][0]
        temp_element_add =Orbit_info[element+1][0]
        if int(temp_element) != int(temp_element_add):
            set_list.append(Orbit_info[element])

    Orbit_info=set_list

    orbit_Time=[]
    orbit_X   =[]
    orbit_Y   =[]
    orbit_Z   =[]

    for element in Orbit_info:
        orbit_Time.append(element[0])
        orbit_X.append(element[1])
        orbit_Y.append(element[2])
        orbit_Z.append(element[3])

    if len(orbit_X) == 0:
        return [], orbit_X, orbit_Y, orbit_Z

    del Orbit_info
    orbit_Time=np.array(orbit_Time)
    orbit_X   =np.array(orbit_X)
    orbit_Y   =np.array(orbit_Y)
    orbit_Z   =np.array(orbit_Z)
    if input_interpolation_method=='cubic':
        spl_x=interp1d(orbit_Time,orbit_X,kind='cubic')
        spl_y=interp1d(orbit_Time,orbit_Y,kind='cubic')
        spl_z=interp1d(orbit_Time,orbit_Z,kind='cubic')
    elif input_interpolation_method=='spline':
        spl_x = inter.InterpolatedUnivariateSpline (orbit_Time,orbit_X)
        spl_y = inter.InterpolatedUnivariateSpline (orbit_Time,orbit_Y)
        spl_z = inter.InterpolatedUnivariateSpline (orbit_Time,orbit_Z)

    input_time = np.arange((orbit_time - date_start) - 100, (orbit_time - date_start) + 100).astype(dtype='int32')
    out_orbit_X=spl_x(input_time)
    out_orbit_Y=spl_y(input_time)
    out_orbit_Z=spl_z(input_time)

    return input_time, out_orbit_X, out_orbit_Y, out_orbit_Z

#----------------------------------------------------------------
def hms2sec(hmsString,convertFlag='float'):
    # input hmsString syntax: XX:XX:XX.xxxxxx
    secString = int(hmsString[0:2])*3600 + \
        int(hmsString[3:5])*60 + \
        float(hmsString[6:])
    if convertFlag == 'int' :
        return int(secString)
    elif convertFlag == 'float' :
        return float(secString)
    else:
        return int(secString)
