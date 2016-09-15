#!/bin/csh -f
###################################################################
# envisatdumpheader2doris.csh
# 1) Read envisat file, and dump to TMPFILE such as
#   ASA_IMS_1PNDPA20021025_175208_000000162010_00356_03416_0005.N1
# 2) Convert dumped header of envisat to a doris resultfile
#    section "readfiles".
#    this includes the orbit in this section.  do not run precise
#    orbits step.
# This script depends on the binary envisat_dump_header that must
# be executable on your system.  The source code of envisat_dump_header
# is included in the Doris distribution.
# BK 16-JUN-2003
#   LOCALE settings of awk gave a problem
#   Added more output.  Still report is that range_time_to_first pixel 
#   may be wrong.  (crop GEO seems to crop wrong for Envisat?)
#   added tr -d '"' to remove quotes from strings.  This may not work everywhere?
#%// Bert Kampes, 13-May-2004
#   After bug report, made string with Datafile a bit smaller.
#   I don;t really like this, but OK.
#%// Bert Kampes, 03-Jun-2004
# added Platform heading + Polarization info
# added Orbit velocities
#%// Mahmut Arikan, 29-Jul-2010
# added support for ERS in Envisat Format
#%// Mahmut Arikan, 19-Oct-2010
###################################################################


set PRG    = `basename "$0"`
set VER    = "v1.2, Doris software"
set AUT    = "TUDelft, (c) 2003-2010"
echo " "
echo "$PRG $VER, $AUT"

# Handle wrong input
if ( $#argv != 1 ) then
cat << __EOFHD
  USAGE: $PRG inputfile
              where inputfile is the envisat SLC file.

  EXAMPLE:
    $PRG ASA_IMS_1PNDPA20021025_175208_000000162010_00356_03416_0005.N1

__EOFHD
exit 1
endif


### Handle input
set ASARFILE = $1
set DUMPHEADER = envisat_dump_header
set TMPFILE = envisat_dump_header.log
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



### Run dumpheader program
$DUMPHEADER $ASARFILE > $TMPFILE

if ( $status != 0 ) then
  echo "exit status of $DUMPHEADER not 0"
  echo "continuing, but this may go wrong."
  echo "We experienced this with data processed by ESRIN"
endif



### Get parameters from dumped header file.
set DUMMY           = "dummy"
set PRODUCT         = `$AWK '/^PRODUCT /{print $3}' $TMPFILE`
set PRODUCT_TYPE    = `$AWK '/^PRODUCT /{print substr($3,2,3)}' $TMPFILE`  # MA  ASA --> Envisat ASAR or SAR --> ERS SAR
set PROC_FAC        = `$AWK '/^PROC_CENTER/{print $3}' $TMPFILE`
set PROC_TIME       = `$AWK '/^PROC_TIME/{print $3" "$4}' $TMPFILE`
set TRACK           = `$AWK '/^REL_ORBIT/{print $3}' $TMPFILE`
set ORBIT           = `$AWK '/^ABS_ORBIT/{print $3}' $TMPFILE`
set FRAME           = `$AWK '/^CYCLE/{print $3}' $TMPFILE`
set FIRST_LINE_TIME = `$AWK '/^FIRST_LINE_TIME/{print $3" "$4}' $TMPFILE | cut -f2 -d'"' | cut -f1 -d'"'`
# not used:
set LAST_LINE_TIME  = `$AWK '/^LAST_LINE_TIME/{print $3" "$4}' $TMPFILE`
set FIRST_MID_LAT   = `$AWK '/^FIRST_MID_LAT/{print $3}' $TMPFILE`
set FIRST_MID_LONG  = `$AWK '/^FIRST_MID_LONG/{print $3}' $TMPFILE`
set LAST_MID_LAT    = `$AWK '/^LAST_MID_LAT/{print $3}' $TMPFILE`
set LAST_MID_LONG   = `$AWK '/^LAST_MID_LONG/{print $3}' $TMPFILE`
### center coordinates
set MIDLAT = `echo "$FIRST_MID_LAT $LAST_MID_LAT" | $AWK '{print ($1+$2)/2e6}'`
set MIDLON = `echo "$FIRST_MID_LONG $LAST_MID_LONG" | $AWK '{print ($1+$2)/2e6}'`
set MIDHEA = `$AWK '/^sub_sat_track/{print $3}' $TMPFILE`

set VOLSETID        = `date --date="$FIRST_LINE_TIME"  +%Y%m%d%H%M%S`  # MA

# Decide on satellite
if ( $PRODUCT_TYPE == "SAR" ) then
  set PRODUCT_TYPE = "ERS_SAR"
else
  set PRODUCT_TYPE = "ASAR"
endif
set numpixels       = `$AWK '/^num_samples_per_line/{print $3}' $TMPFILE`
set numlines        = `$AWK '/^num_output_lines/{print $3}' $TMPFILE`
set checknumlines   = `$AWK '/^num_output_lines/{print $3+1}' $TMPFILE`
set wgt_azimuth     = `$AWK '/^filter_az/{print $3}' $TMPFILE`
set wgt_range       = `$AWK '/^filter_window/{print $3}' $TMPFILE`
### remove double quotes:
set wgt_azimuth     = `echo $wgt_azimuth | tr -d '"'`
set wgt_range       = `echo $wgt_range | tr -d '"'`
set RSR             = `$AWK '/^range_samp_rate/{printf "%f\n", $3/1e6}' $TMPFILE`
### After first awk, RBW="{1550000.000000," thus do a cut and to MHz
set RBW             = `$AWK '/^bandwidth.tot_bw_range/{print $3}' $TMPFILE | cut -f2 -d'{' | cut -f1 -d',' | $AWK '{print $1/1e6}'`
### not ok, this is time used in doppler est.
### not ok: set SRT = `$AWK '/slant_range_time/{print $3/1e6}' $TMPFILE`
set SAMPLE          = `$AWK '/^first_line_tie_points.samp_numbers/{print $3}' $TMPFILE | cut -f2 -d'{' | cut -f1 -d','`       # 1
### fast time in [ns], RSR in [MHz], 2-way time.
set T_TO_SAMPLE     = `$AWK '/^first_line_tie_points.slant_range_times/{print $3}' $TMPFILE | cut -f2 -d'{' | cut -f1 -d','`  # 5523270.000000
set SRT = `echo "$T_TO_SAMPLE $SAMPLE $RSR" | $AWK '{printf "%f\n", $1/1e6+($2-1)/($3*1e3)}'` # 5523270.000000 1 19.207680

### also obtain time to last pixel and compute RSR, same for PRF
#RSR = (PIX_N-PIX_1)/DT
#kampes@capone[18:43]: reken 5174/(5810316.50-5540945.5)
#.01920770981286033017


### After first awk, PRF="{14.1223121," thus do a cut 
set PRF             = `$AWK '/^image_parameters.prf_value/{print $3}' $TMPFILE | cut -f2 -d'{' | cut -f1 -d','`
set ABW             = `$AWK '/^to_bw_az/{print $3}' $TMPFILE`
### Doppler parameters
set FDC0            = `$AWK '/^dop_coef/{print $3 $4 $5}' $TMPFILE | cut -f2 -d'{' | cut -f1 -d','`
set FDC1            = `$AWK '/^dop_coef/{print $3 $4 $5}' $TMPFILE | cut -f2 -d'{' | cut -f2 -d','`
set FDC2            = `$AWK '/^dop_coef/{print $3 $4 $5}' $TMPFILE | cut -f2 -d'{' | cut -f3 -d','`
set WAVELENGTH      = `$AWK '/^radar_freq/{print 299792458.0/$3}' $TMPFILE`

### Assume we have 5 state vectors, but check anyway...
set state_x         = `$AWK '/^orbit_state_vectors...x_pos/{printf "%f\n", $3/1e2}' $TMPFILE`
set state_y         = `$AWK '/^orbit_state_vectors...y_pos/{printf "%f\n", $3/1e2}' $TMPFILE`
set state_z         = `$AWK '/^orbit_state_vectors...z_pos/{printf "%f\n", $3/1e2}' $TMPFILE`
set state_xdot      = `$AWK '/^orbit_state_vectors...x_vel/{printf "%f\n", $3/1e2}' $TMPFILE`
set state_ydot      = `$AWK '/^orbit_state_vectors...y_vel/{printf "%f\n", $3/1e2}' $TMPFILE`
set state_zdot      = `$AWK '/^orbit_state_vectors...z_vel/{printf "%f\n", $3/1e2}' $TMPFILE`
### Time seems to be daynumber seconds milliseconds
set state_t_sec     = `$AWK '/^orbit_state_vectors...state_vect_time/{print $4}' $TMPFILE | cut -f2 -d'=' | cut -f1 -d','`
set state_t_dsec    = `$AWK '/^orbit_state_vectors...state_vect_time/{print $5}' $TMPFILE | cut -f2 -d'=' | cut -f1 -d'}'`
set NUMSTATEVECTORS = $#state_x

### SARPROCESSOR: is used by Doris
set SOFTWARE_VER    = `$AWK '/^SOFTWARE_VER/{print $3}'   $TMPFILE`
set SOFTWARE_VER    = "ASAR $SOFTWARE_VER"
### some extra for general information, not used by Doris.
set SWATH           = `$AWK '/^SWATH/{print $3}'          $TMPFILE`
set PASS            = `$AWK '/^PASS/{print $3}'           $TMPFILE`
set SPH_DESCRIPTOR  = `$AWK '/^SPH_DESCRIPTOR/{print $3, $4, $5, $6, $7, $8}' $TMPFILE`
set FREQUENCY       = `$AWK '/^radar_freq/{print $3}'     $TMPFILE`
set POLAR_1         = `$AWK '/^MDS1_TX_RX_POLAR/{print $3}'          $TMPFILE`
set POLAR_2         = `$AWK '/^MDS2_TX_RX_POLAR/{print $3}'          $TMPFILE`
set POLAR           = "$POLAR_1 $POLAR_2"

### create a result section to a tmp file.
cat << __EOFHD


*******************************************************************
*_Start_readfiles:
*******************************************************************
Volume file: $PRODUCT
Volume_ID:                                      $DUMMY
Volume_identifier:                              $DUMMY 
Volume_set_identifier:                          $VOLSETID
(Check)Number of records in ref. file:          $checknumlines
Product type specifier:                         $PRODUCT_TYPE
SAR_PROCESSOR:                                  $SOFTWARE_VER
SWATH:                                          $SWATH
PASS:                                           $PASS
SPH_DESCRIPTOR:                                 $SPH_DESCRIPTOR
IMAGING_MODE:                                   $POLAR
RADAR_FREQUENCY (HZ):                           $FREQUENCY

Logical volume generating facility:             $PROC_FAC
Logical volume creation date:                   $DUMMY
Location and date/time of product creation:     $PROC_TIME
Scene identification:                           ORBIT $ORBIT  TRACK $TRACK
Scene location:                                 FRAME $FRAME
Leader file: $PRODUCT
Sensor platform mission identifer:              ENVISAT-ASAR-SLC
Scene_centre_latitude:                          $MIDLAT
Scene_centre_longitude:                         $MIDLON
Scene_centre_heading:                           $MIDHEA
Radar_wavelength (m):                           $WAVELENGTH
First_pixel_azimuth_time (UTC):                 $FIRST_LINE_TIME
TIME TO LAST LINE: compute prf:                 $LAST_LINE_TIME
Pulse_Repetition_Frequency (computed, Hz):      $PRF
Total_azimuth_band_width (Hz):                  $ABW
Weighting_azimuth:                              $wgt_azimuth
Xtrack_f_DC_constant (Hz, early edge):          $FDC0
Xtrack_f_DC_linear (Hz/s, early edge):          $FDC1
Xtrack_f_DC_quadratic (Hz/s/s, early edge):     $FDC2
Range_time_to_first_pixel (2way) (ms):          $SRT
Range_sampling_rate (computed, MHz):            $RSR
Total_range_band_width (MHz):                   $RBW
Weighting_range:                                $wgt_range

*******************************************************************
*_Start_leader_datapoints
*******************************************************************
t(s)            X(m)            Y(m)            Z(m)
NUMBER_OF_DATAPOINTS:   $NUMSTATEVECTORS
__EOFHD


### t_dsec is given in ms. (place 0's before it, do not use "$t.$dt")
@ i = 1
while ( $i <= $NUMSTATEVECTORS )
  set t = `echo "$state_t_sec[$i] $state_t_dsec[$i]" | awk '{printf "%.8f", $1+$2/1e6}'`
  #echo "$t    $state_x[$i]     $state_y[$i]    $state_z[$i]"           # dump state vectors only
  echo "$t    $state_x[$i]     $state_y[$i]     $state_z[$i]    $state_xdot[$i]     $state_ydot[$i]    $state_zdot[$i]"  # dump velocities too
  @ i++
end



###############################
#%// Bert Kampes, 03-Jun-2004
# report that line is too long:
#Datafile:                                       $ASARFILE
# therefore, suggested to check with
#if ( `echo "$ASARFILE" | awk '{print length($1)}'` > 70 ) then
##or: if ( `echo "$ASARFILE" | wc -c` > 70 ) then
#  echo "Datafile: $ASARFILE"
#else
#  echo "Datafile:                                       $ASARFILE"
#endif
###############################
#does not work
#set DF = "Datafile:                                        $ASARFILE"
#if ( `echo "$ASARFILE" | wc -c` > 70 ) then
  #set DF = "Datafile: $ASARFILE"
#endif
#again check length within 127 total, else create a soft link and warn...

cat << __EOFHD

*******************************************************************
* End_leader_datapoints:_NORMAL
*******************************************************************
#Report that else next line could easily be too long:
Datafile:    $ASARFILE
Dataformat:  ENVISAT 
Number_of_lines_original:                       $numlines
Number_of_pixels_original:                      $numpixels
*******************************************************************
* End_readfiles:_NORMAL
*******************************************************************

__EOFHD





###EOF
