#!/bin/bash -f
##
## baseline.doris.sh
## 
## Made by Petar Marinkovic
## Login   <pmar@cake.lr.tudelft.nl>
## 
## Started on  Tue Jul 18 17:11:53 2006 Petar Marinkovic
## Last update Tue Jul 18 17:24:33 2006 Petar Marinkovic
##
## DESCRIPTION:
##
## CHANGE.LOG:
##
## [BK 27-Oct-2000, 00:00]: initial version in csh ~ baseline.doris
##
## [BK 14-Apr-2003, 06:35]: ??
##
## [BK 18-Jul-2006, 16:00]: script rewrite in bash, doris_v318
## supported, averaging the values over grid
##
## TODO:
##
## [MA 20071022, 16:30]: some more attributes is introduce for comparision with Descw 
##		       : compatiblity with older doris versions improved
##\

PRG=`basename "$0"`
VER="v1.2, DEOS software"
AUT="Bert Kampes, Petar Marinkovic and Mahmut Arikan (c)2007"
# echo "$PRG $VER, $AUT"\\n

### Handle input
if [[ $# < 2 ]]; then
    cat << EOF
    PROGRAM:  ${PRG} -- Obtain baseline parameterization from orbits.
	
    SYNOPSIS:
	${PRG}  master.res  slave.res

    OPTIONS:
	master.res   master result file from doris, precise orbits.
	slave.res    slave result file from doris, precise orbits.

    Basically a dummy run is performed, only to get DUMPBASELINE.
    Then from the standard output the appropriate sections are written.
    The DUMPBASELINE card is used, computing the baseline at a grid 15x10,
    then fitting a polynomial through these values (Bperp, theta and range).
    
    EXAMPLE:
	${PRG} master.res slave.res
EOF
      exit 1
fi

### Correct input.
MASTERFILE=${1}
SLAVEFILE=${2}

if [[ -n ${USER} ]]; then
    USER=$( whoami )
fi

### Create dummy input file(s).
#TMPDIR="/tmp"
TMPDIR=/tmp
#TMPDIR=tmp # MA
DORIS=${DORIS:-$(which doris)} # MA: if defined use externally define DORIS. Also possible to write ${DORIS:=$(which doris)} alone.
DORISIN=${TMPDIR}/${USER}.$$.in
DORISOUT=${TMPDIR}/${USER}.$$.out
LOGFILE=${TMPDIR}/${USER}.$$.log
PRODFILE=${TMPDIR}/${USER}.$$.prod
TMPFILE=${TMPDIR}/${USER}.$$.tmp # [PM]: do not like new temp file

###
cat << EOF > ${DORISIN}
c ***********************************************************
c File created by: $PRG $VER
c Author: $AUT
c  This is an input file for Doris to perform baseline estimation
c  from prompt.
c ***********************************************************
 c
 comment  ___general options___
 c
SCREEN     	info
MEMORY		5
OVERWRITE	OFF
BATCH		ON
LISTINPUT	OFF
DUMPBASELINE    15 10
PROCESS   	coarseorb
LOGFILE     	$LOGFILE
M_RESFILE   	$MASTERFILE
S_RESFILE   	$SLAVEFILE
I_RESFILE   	$PRODFILE
 c
 c
 c
STOP
EOF

# Run Doris, module XXX, use DUMPBASELINE card.
$DORIS $DORISIN > $DORISOUT 2> /dev/null

awk '/Bpar/ {print $5}' ${DORISOUT} > ${TMPFILE}.bpar
awk '/Bperp/ {print $NF}' ${DORISOUT} > ${TMPFILE}.bperp
awk '/Height ambiguity/ {print $NF}' ${DORISOUT} > ${TMPFILE}.hamb
awk '/Look angle|theta \(deg\)/ {print $NF}' ${DORISOUT} > ${TMPFILE}.look

# average values computed on the grid specified by DUMPBASELINE card
BPAR=$( awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' ${TMPFILE}.bpar )
BPERP=$( awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' ${TMPFILE}.bperp )
HAMB=$( awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' ${TMPFILE}.hamb )
LOOK=$( awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' ${TMPFILE}.look )
BTEMP=$( awk '/Btemp/{print $3}' ${PRODFILE} ) ;  # Mahmut additional attributes, mostly based on slave
MISSION=$( awk '{FS=":"}/Product type specifier/{ split($3,a,"."); gsub(/ /,"",$2); if(a[1]== "ERS-1"){print "E1" }else if(a[1]== "ERS-2"){print "E2"}else if($2=="ASAR"){print "N1"}else{print "NA"}}' \
		${SLAVEFILE} )
# MA, I picked slave since you can also do master master to get 0 basline info
#DATE=$( [[ "$MISSION" != "N1" ]] && awk '/Volume_set_identifier/{ print substr($2,1,8)}' ${SLAVEFILE} || awk -F "_" '/Volume file:/{ print substr($3,7)}'  ${SLAVEFILE}  )  # change for ERS due buggy VDF files
DATE=$( awk '/First_pixel_azimuth_time/{ print substr($3,1,11)}' ${SLAVEFILE} ) 
ORBIT=$( awk '/Scene identification:/{printf "%05d", $4}' ${SLAVEFILE} )
mfDC=$( awk '{FS=":"}/Xtrack_f_DC_constant/{gsub(/ /,"",$2); printf "%s\n", $2}' ${MASTERFILE} )
sfDC=$( awk '{FS=":"}/Xtrack_f_DC_constant/{gsub(/ /,"",$2); printf "%s\n", $2}' ${SLAVEFILE} )
dfDC=$( awk -v m="$mfDC" -v s="$sfDC" 'BEGIN{print s-m }' )  # w.r.t master
s_cLON=$( awk '/Scene_centre_longitude/ {printf "%5.3f\n", $2}' ${SLAVEFILE})
s_cLAT=$( awk '/Scene_centre_latitude/  {printf "%5.3f\n", $2}' ${SLAVEFILE})
#INFO=$( awk '/^precise_orbits:/{ print ($2 == "1") ? "preciseorb" : "annotated"}' ${SLAVEFILE} )
INFO=$( awk '/^precise_orbits:/{ yes=$2 }; ( yes == "1" && $2 ~ "Orbit_dir:"){ po= $3  }END{ if( yes == "1" ){ print ( po != "" ) ? po : "preciseorb" }else{ print "annotated"} }' ${SLAVEFILE} ) ; # Orbit_dir tag printed by orbits.precise.sh to .res

echo -e "Mission: \t\t ${MISSION} "
echo -e "Date: \t\t\t ${DATE}"
echo -e "Orbit: \t\t\t ${ORBIT}"
echo -e "Bpar [m]: \t\t ${BPAR}"
echo -e "Bperp [m]: \t\t ${BPERP}"
echo -e "Btemp [days]: \t\t ${BTEMP}"; # MA added
echo -e "fDC [Hz]: \t\t $( awk -v o="${sfDC}" 'BEGIN{printf "%.3f", o}' )"; # MA added
echo -e "dfDC [Hz]: \t\t $( awk -v o="${dfDC}" 'BEGIN{printf "%.3f", o}' )"; # MA added
echo -e "Center Lon [deg]: \t ${s_cLON}"
echo -e "Center Lat [deg]: \t ${s_cLAT}"
echo -e "Look.angl [deg]: \t ${LOOK}"
echo -e "Height.amb [m]: \t ${HAMB}"
echo -e "INFO: \t\t\t $INFO"

### Tidy up.
#rm -rf $TMPDIR
 rm -f $DORISIN
 rm -f $DORISOUT
 rm -f $LOGFILE
 rm -f $PRODFILE
 rm -f $TMPFILE.bpar
 rm -f $TMPFILE.bperp
 rm -f $TMPFILE.hamb
 rm -f $TMPFILE.look

exit 0

### EOF.
