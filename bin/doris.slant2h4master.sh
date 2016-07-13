#!/bin/bash
#
#
# Author: Mahmut Arikan
#
# TUDelft 2007
#
# Trick Doris to use radarcoded SRTM heights to gecode.

iFile=$1
radarcodedDEM=$2


# Functions
goodfunc(){

grep _Start_slant2h $iFile && echo slant2h info already there, pls check && exit 0;

#m_s=`awk '/INTERFEROGRAM RESULTFILE/{print substr($3,1,11)}' *_*.res`
m_s=master

l0=`grep First_line master.res | awk 'END{print $4}'`
LN=`grep Last_line master.res | awk 'END{print $4}'`
p0=`grep First_pixel master.res | awk 'END{print $4}'`
PN=`grep Last_pixel master.res | awk 'END{print $4}'`
#mlA=`grep Multilookfactor_azimuth_direction master.res | awk 'END{print $2}'`
#mlR=`grep Multilookfactor_range_direction   master.res | awk 'END{print $2}'`
mlA=1
mlR=1

cat << END >> $iFile
*******************************************************************
*_Start_slant2h:
*******************************************************************
Method:                         schwabisch
Data_output_file:                       Outdata/${m_s}_S2H.float
Data_output_format:                     real4
First_line (w.r.t. original_master):    $l0
Last_line (w.r.t. original_master):     $LN
First_pixel (w.r.t. original_master):   $p0
Last_pixel (w.r.t. original_master):    $PN
Multilookfactor_azimuth_direction:      $mlA
Multilookfactor_range_direction:        $mlR
Ellipsoid (name,a,b):                   WGS84 6.37814e+06 6.35675e+06
*******************************************************************
* End_slant2h:_NORMAL
*******************************************************************
END

echo $iFile is updated
tail -n 35 $iFile

# cd Outdata && ln -s ../refDemLP.raw ${m_s}_S2H.float && echo "File: Outdata/${m_s}_S2H.float --> refDemLP.raw"
#cd Outdata && ln -s refdem_hei.raw ${m_s}_S2H.float && echo "File: Outdata/${m_s}_S2H.float --> refDemLP.raw"
cd Outdata && ln -sf ${radarcodedDEM} ${m_s}_S2H.float && echo "File: Outdata/${m_s}_S2H.float --> ${radarcodedDEM}"
}

case $# in
         2) goodfunc ;;
         *) echo -e "Usage: ${0##*/} <interferometric.res> <radarcodedDEM.raw> \n"
           ;;
esac


