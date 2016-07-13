#!/bin/csh -f
###################################################################
# gammaReadfiles.csh
# 1) Process the RAW data to SLC with GAMMA. Provide *.slc.par file.
# 2) Convert dumped header of envisat to a doris resultfile
#    section "readfiles".
#    this includes the orbit in this section.  
# This script is based on the envisat_dumpheader2doris.csh by
# Bert Kampes 16-JUN-2003
# Modified for gamma by Batuhan Osmanoglu. 
###################################################################

set PRG    = `basename "$0"`
set VER    = "v1.1, doris software"
set AUT    = "Batuhan Osmanoglu 2009, Bert Kampes, (c)2003"
# moved later in the file echo " "
# moved later in the file echo "$PRG $VER, $AUT"

# Handle wrong input
if ( $#argv != 3 ) then
cat << __EOFHD
  USAGE:$PRG processParameterFile sensorParameterFile slcFile
              where 
              processParameterFile is the gamma SLC processing  par file (pXXX.slc.par, or XXX.pslcpar)
              sensorParameterFile is the gamma sensor par file (ERS1_ESA.par, or System.par)
              slcfile is the complex SAR image file in 2 bytes integer format (CI2).

  EXAMPLE:
    $PRG p20070105.slc.par System.par 20070105.slc

__EOFHD
exit 1
endif

### Handle input
set PARFILE = $1
set SYSFILE = $2
set SLCFILE = $3
set TMPFILE = gammaReadfiles.tmp

### find out if nawk is available, since awk seems not to know atan2?
set AWK = awk;
echo 1 | /usr/xpg4/bin/awk '{print 1}' >& /dev/null
if ( $status == 0 ) set AWK = '/usr/xpg4/bin/awk'
echo 1 | gawk '{print 1}' >& /dev/null
if ( $status == 0 ) set AWK = 'gawk'
echo 1 | nawk '{print 1}' >& /dev/null
if ( $status == 0 ) set AWK = 'nawk'
echo "Using as awk program: $AWK"

### set LOCALE to decimal.point style
#kampes@newton[18:22]: setenv LC_NUMERIC pl
#kampes@newton[18:22]: echo "5300000000.0" | awk '{print 299792458.0/$1}'
#0,0565646
#kampes@newton[18:22]: setenv LC_NUMERIC C
#kampes@newton[18:22]: echo "5300000000.0" | awk '{print 299792458.0/$1}'
#0.0565646
#we could check first with: locale -ck decimal_point
#setenv LC_ALL POSIX
setenv LC_ALL C
setenv LC_NUMERIC C

### Get parameters from dumped header file.
set today	    = `date`
set sceneDate	    = `$AWK '/^date/{printf "%d%02d%02d\n", $2,$3,$4}' $PARFILE`
set sceneDate	    = `date -d $sceneDate +%d-%b-%Y`
doris -ver >& $TMPFILE
set version	    = `cat $TMPFILE | grep "Software version" | cut -f2 -d:`
set dummy           = "dummy"
set product         = `$AWK '/^title/{print $2}' $PARFILE`
set checkNumLines   = `$AWK '/^azimuth_pixels/{print $2}' $PARFILE`
#Product Type (satellite info is missing in gamma files.)
set productType     = `$AWK '/^sensor_name/{print $2}' $SYSFILE`
set sarProcessor    = "GAMMA"
set frequency	    = `$AWK '/^SAR_center_frequency/{printf "%.6f", $2}' $SYSFILE`
set midLat	    = `$AWK '/^scene_center_latitude/{print $2}' $PARFILE`
set midLon	    = `$AWK '/^scene_center_longitude/{print $2}' $PARFILE`
set pass	    = `$AWK '/^map_coordinate_1:/{lat1=$2};/^map_coordinate_3:/{if ($2>lat1) print "ASCENDING"; else print "DESCENDING"}' $PARFILE` #check first and last latitude.
set wavelength	    = `$AWK '/^SAR_center_frequency/{printf "%.6f", 2.997e+8/$2}' $SYSFILE`
set firstLineTime   = `$AWK '/^raw_data_start_time/{printf "%02d:%02d:%2.6f", $2, $3, $4}' $PARFILE`
set firstLineTimeSec = `hhmmss2sec.py $firstLineTime`
# firstLineTimeSec needs to get corrected for the offset.
set azimuth_offset  = `$AWK '/^azimuth_offset/{print $2}' $PARFILE`
set firstLineTimeSec = `echo ${firstLineTimeSec} ${azimuth_offset} | $AWK '{printf "%.8f", $1+$2};'`
set firstLineTime   = `sec2hhmmss.py $firstLineTimeSec`
#set prf		    = `$AWK '/^start_time/{flt=$2};/^end_time/{llt=$2};/^azimuth_lines/{print $2/(llt-flt)}' $PARFILE`
set prf		    = `$AWK '/^pulse_repetition_frequency/{printf "%.6f", $2}' $PARFILE`
set lastLineTimeSec = `echo $firstLineTimeSec $checkNumLines $prf | $AWK '{printf "%.8f", $1+$2/$3};'`
set lastLineTime    = `sec2hhmmss.py $lastLineTimeSec`
set abw		    = `$AWK '/^pulse_repetition_frequency/{prf=$2};/^azimuth_bandwidth_fraction/{printf "%.6f", prf*$2}' $PARFILE`
set FDC0	    = `$AWK '/^doppler_polynomial/{printf "%.6f", $2}' $PARFILE`
set FDC1	    = `$AWK '/^doppler_polynomial/{printf "%.6f", $3}' $PARFILE`
set FDC2	    = `$AWK '/^doppler_polynomial/{printf "%.6f", $4}' $PARFILE`
# Two Way TravelTime in ms (hence *1e+3)
set TWT		    = `$AWK '/^echo_time_delay/{printf "%.6f", $2*1e+3}' $PARFILE` 
# correct TWT for range extension (hence *1e+3)
#correct for chirp extension (near range extension time changes first pixel time.)
set near_range_extension = `$AWK '/^near_range_extension/{print $2}' $PARFILE`
set RSR		    = `$AWK '/^ADC_sampling_frequency/{printf "%.6f", $2/1e+6}' $SYSFILE`
set TWT 	    = `echo $TWT ${near_range_extension} ${RSR}| $AWK '{printf "%.8f", $1-$2/($3*1e+3)};'`
set RBW		    = `$AWK '/^chirp_bandwidth/{printf "%.6f", $2/1e+6}' $SYSFILE`
set NUMSTATEVECTORS = `$AWK '/^number_of_state_vectors/{print $2}' $PARFILE`
set numlines	    = `$AWK '/^azimuth_pixels/{print $2}' $PARFILE`
set numpixels	    = `$AWK '/^range_pixels/{print $2}' $PARFILE`
### modify some of the values
if ( "$productType" =~ "IS*" ) then
	set productType = "ASAR"
else if ( "$productType" =~ "ERS*" ) then
	set productType = "ERS"
else if ( "$productType" == "PALSAR" ) then
	set productType = "ALOS"
endif	

### create a result section to a tmp file.
cat << __EOFHD
*******************************************************************
*_Start_readfiles:
*******************************************************************
Volume file: $product
Volume_ID:                                      $dummy
Volume_identifier:                              $dummy
Volume_set_identifier:                          $dummy
(Check)Number of records in ref. file:          $checkNumLines
Product type specifier:                         $productType
SAR_PROCESSOR:                                  $sarProcessor
SWATH:						$dummy
PASS:						$pass
RADAR_FREQUENCY (HZ):                           $frequency
Logical volume generating facility:             $dummy
Logical volume creation date:                   $dummy
Location and date/time of product creation:     $dummy
Scene identification:                           ORBIT $dummy
Scene location:                                 FRAME $dummy
Leader file: $product
Sensor platform mission identifer:              $productType
Scene_centre_latitude:                          $midLat
Scene_centre_longitude:                         $midLon
Radar_wavelength (m):                           $wavelength
First_pixel_azimuth_time (UTC):                 $sceneDate $firstLineTime
TIME TO LAST LINE: compute prf:                 $sceneDate $lastLineTime
Pulse_Repetition_Frequency (computed, Hz):      $prf
Total_azimuth_band_width (Hz):                  $abw
Weighting_azimuth:                              $dummy
Xtrack_f_DC_constant (Hz, early edge):          $FDC0
Xtrack_f_DC_linear (Hz/s, early edge):          $FDC1
Xtrack_f_DC_quadratic (Hz/s/s, early edge):     $FDC2
Range_time_to_first_pixel (2way) (ms):          $TWT
Range_sampling_rate (computed, MHz):            $RSR
Total_range_band_width (MHz):                   $RBW
Weighting_range:                                $dummy

*******************************************************************
*_Start_leader_datapoints
*******************************************************************
t(s)            X(m)            Y(m)            Z(m)
NUMBER_OF_DATAPOINTS:   $NUMSTATEVECTORS
__EOFHD

# Dump the state vectors in the required format...
$AWK '/^time_of_first_state_vector/{t=$2};/^state_vector_interval/{dt=$2;c=0};/^state_vector_position/{ printf "%.6f\t%.3f\t%.3f\t%.3f\n", t+(c++)*dt, $2, $3, $4 }' $PARFILE

# Dump the closing section...

cat << __EOFHD
*******************************************************************
* End_leader_datapoints:_NORMAL
*******************************************************************
Datafile:    $SLCFILE
Number_of_lines_original:                       $numlines
Number_of_pixels_original:                      $numpixels
*******************************************************************
* End_readfiles:_NORMAL
*******************************************************************

__EOFHD

rm -rf $TMPFILE

###EOF
#sleep 1

