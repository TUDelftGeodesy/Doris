#!/bin/bash
##
## construct_dem.sh
## v1.2, 18-12-2009
## v1.1, 17-12-2008
## v1.0, 22-08-2006
## 
## Made by Freek van Leijen and Zbigniew Perski
## Delft Institute of Earth Observation and Space Systems
## Delft University of Technology
## 
## Started on  Wed May 17 10:00:00 2006 Freek van Leijen
## Previous update Mon May 22 17:36:09 2006 Petar Marinkovic
## Last update Wed June 28 22:28:19 Mahmut Arikan
## added httpSupport (v1.2) on Fri Dec 18 16:39:39 Batuhan Osmanoglu
## 
##
## DESCRIPTION: Downloads, merges and fills voids of SRTM
## data based on coordinates of the area of interest.
## Only basic Linux/Unix commands, wget and GMT are used.
##
## NOTE: Scripts will have a problem with Antarctica and
## when your area crosses the -180/180 meridian.
##
## INPUT: [$1] : project_name
##        [$2] : W
##        [$3] : E
##        [$4] : S
##        [$5] : N
##        [$6] : srtm_code [SRTM1/3] (optional)
##        [$7] : link_1 (optional, default is http://dds.cr.usgs.gov/srtm) 
##        [$8] : link_2 (optional, default is http://dds.cr.usgs.gov/srtm)
##        [$9] : ftp_user   (optional)
##       [$10] : ftp_pass   (optional)
##
## EXAMPLE: ./construct_dem.sh netherlands 3.3 7.3 50.7 53.7 SRTM3
##
## CHANGE.LOG:
##
## MA (some modifications in bash part, keeping some more output files, doris missing prm added)
## FvL fixed bug in input check, explicitely stated ftp_user and ftp_pass
## FvL changed ftp server for SRTM30, update for Doris v4.01, 2008
## MA define variable for awk calls
## Batu added http support, changed folder name to /version2_1/ from version2
##
## TODO Skiping download of never existing tiles, may be a file list for excluded tiles over oceans
##
Revision="v1.2"

AWK=`AWK=$(which nawk 2> /dev/null); [ "$AWK" == "" ] && echo $(which awk) || echo $AWK` # MA awk variable: 1. look for nawk (on old systems) else use awk (gnu awk)

# -- FUNCTIONS declarations [modularity] -----

# ftpcall or wget
downloadFile(){
serverType=`echo ${3}| awk 'BEGIN{ FS=":"};{print $1}'`
case "${serverType}" in 
  ftp )
  ftp_call $1 $2 $3 $4 $5  
  ;;
  http )
  #generate index
  i=`echo ${3}| awk 'BEGIN{ FS="/"};{print $3}'`
  if [ ! -e ${i} ] 
  then 
    generateIndexFile ${3} ${i}
  fi
  tileUrl=`grep "$1/${2}" ${i}| head -n1`
  [ ! -z ${tileUrl} ] && wget ${tileUrl}
  ;;
esac
}

# generates the file list of http site
generateIndexFile(){
servername=$1
filename=$2
#generate temporary folder
str0=$$ #get PID
POS=2  # Starting from position 2 in the string.
LEN=8  # Extract eight characters.
str1=$( echo "$str0" | md5sum | md5sum )
# Doubly scramble:     ^^^^^^   ^^^^^^
tempFolder="${str1:$POS:$LEN}"
#projectFolder=$PWD
mkdir ${tempFolder}
wget -q -k -nd -r -P"./${tempFolder}" -A*.html* "${servername}" #q:quiet, k:convert links, nd:no directories, r:recursive, P:output directory prefix, A:accept pattern
cat ./${tempFolder}/index.html.* | egrep -o "http:.*version2.*SRTM1.*hgt*" | cut -f1 --delimiter=\" >   ${filename}
cat ./${tempFolder}/index.html.* | egrep -o "http:.*version2.*SRTM3.*hgt*" | cut -f1 --delimiter=\" >>  ${filename}
cat ./${tempFolder}/index.html.* | egrep -o "http:.*version2.*SRTM30.*dem*" | cut -f1 --delimiter=\" >> ${filename}
#cat index.html.* | egrep -o "http:.*version1.*SRTM3.*dem*" | cut -f1 --delimiter=\" >> ${projectFolder}/${filename}
rm -rf ./${tempFolder}
[ $? -ne 0 ] && exit 1
}

# downloads DEM patch from ftp site
ftp_call()
{
remotedir=$1
remotefile=$2
ftp_link=$3
ftp_user=$4
ftp_pass=$5
ftp -inv $ftp_link << EOF >> ftp_log 
user $ftp_user $ftp_pass
binary
cd $remotedir
pwd
get $remotefile
bye
EOF
}

#MA more functions
#printing parameters N E S W
pPRM(){ echo -e "Parameters:\n\t\t\tNorth:$north_border\n\tWest:$west\t\t\t\tEast:$east_border
\n\t\t\tSouth:$south  "; }

#MA need by below case statement
ckSRTM(){ 
	   srtm=`echo $srtm | tr a-z A-Z` # change to upper case
	   [[ $srtm != "SRTM1" && $srtm != "SRTM3" ]] && \
	   echo -e "\nNot a valid SRTM version $srtm !!! See usage.\n" && exit 127;
	   #echo -e "\nUsing SRTM Version \033[01;33m$srtm \033[00m";
	   echo -e "\nUsing SRTM Version $srtm";
}

pFTP(){
	   echo -e "\nSRTM   url: $ftp_link_1"
	   echo -e "SRTM30 url: $ftp_link_2"
}


# arguments substitution and checks

# Create/check project directory
project=$1
[ ! -d ${project} ] && mkdir $project

echo $0 $@ > ${project}/construct_dem.command

# FLOOR FIRST(!) integers and then check
west=`echo ${2%%.*}  | ${AWK} '{ if ($1 < 0) $1 = $1 - 1; print $1 }'`
east=`echo ${3%%.*}  | ${AWK} '{ if ($1 < 0) $1 = $1 - 1; print $1 }'`
south=`echo ${4%%.*} | ${AWK} '{ if ($1 < 0) $1 = $1 - 1; print $1 }'`
north=`echo ${5%%.*} | ${AWK} '{ if ($1 < 0) $1 = $1 - 1; print $1 }'`

#FvL
east_border=$[$east+1] 
north_border=$[$north+1]

# MA checking parameters
# input arguments : check
# check on arguments $7 and $8, if specified parsed, otherwise default ftp 
# links used

case $# in
	5) pPRM;  
	   srtm="SRTM3"; # default
	   echo -e "\nUsing default SRTM Version $srtm";
	   #ftp_link_1="e0srp01u.ecs.nasa.gov"; ftp_link_2="e0srp01u.ecs.nasa.gov";
	   ftp_link_1="http://dds.cr.usgs.gov/srtm/";ftp_link_2="http://dds.cr.usgs.gov/srtm/" 
	   ftp_user="anonymous"; ftp_pass="anonymous";
	   pFTP;
	   #exit 127 #debug
	   ;;
	6) pPRM; srtm=$6
	   ckSRTM;
	   #ftp_link_1="e0srp01u.ecs.nasa.gov"; ftp_link_2="e0srp01u.ecs.nasa.gov"; 
	   ftp_link_1="http://dds.cr.usgs.gov/srtm/";ftp_link_2="http://dds.cr.usgs.gov/srtm/" 
	   ftp_user="anonymous"; ftp_pass="anonymous";
	   pFTP;
	   ;;
	7) pPRM; srtm=$6
           ckSRTM;
	   #ftp_link_1=$7; ftp_link_2="e0srp01u.ecs.nasa.gov";
	   ftp_link_1=$7; ftp_link_2="http://dds.cr.usgs.gov/srtm/"
	   ftp_user="anonymous"; ftp_pass="anonymous";
	   pFTP ;;
	8) pPRM; srtm=$6
	   ckSRTM; 
	   ftp_link_1=$7; ftp_link_2=$8;
	   pFTP ;;
	9) pPRM; srtm=$6
	   ckSRTM; 
	   ftp_link_1=$7; ftp_link_2=$8;
           ftp_user=$9; ftp_pass=;
	   pFTP ;;
	10) pPRM; srtm=$6
	   ckSRTM; 
	   ftp_link_1=$7; ftp_link_2=$8;
           ftp_user=$9; ftp_pass=${10}; # careful with ${1} vs ${10}
	   pFTP ;;
	*) echo -e "Doris software,  Revision: $Revision,  Author: TUDelft  ";
	   echo -e "\n[Usage]    : `basename $0` project W E S N SRTM[1|3] <ftp1> <ftp2>" ;
	   echo -e "\n[Example]  : `basename $0` netherlands 3.3 7.3 50.7 53.7 SRTM3 "
	   echo -e "\n[ftp sites]: SRTM: e0srp01u.ecs.nasa.gov & SRTM30: topex.ucsd.edu (defaults)" 
           echo -e "\n[awk ver.] : ${AWK} - $(${AWK} -W version | head -n 1) \n"
             which GMT &> /dev/null || \
	     echo -e "[remark]   : this script requires that GMT package is installed on your system.\n";  # MA
	   echo -e "\nPlease check parameters: $* \n"; 
           exit;;
esac

 
#MA check if west > east, if so exit. BTW  no check for NORTHING in the code
[ ${west%%.*} -gt ${east%%.*} ] && echo -e "\n E:$east can't be less than W:$west! \n" && exit 127
#MA check if west > east, if so exit. BTW  no check for NORTHING in the code
[ ${south%%.*} -gt ${north%%.*} ] && echo -e "\n N:$north can't be less than S:$south! \n" && exit 127

#MA keep some of intermediate files, see end of this script:
echo -e "\nDo you want to keep the following intermediate files:  "
echo -e "${project}/srtm_${project}.grd \n${project}/srtm30_${project}_merged.grd \
\n${project}/final_${project}.grd \t (Y/N) or Ctrl-C to exit? \c"
read keep

#MA SRTM remote dir definitions
case $srtm in
		SRTM1) dir_list="Region_01 Region_02 Region_03 Region_04 Region_05 Region_06 Region_07"
  			 gmt_format=1c ;;
		SRTM3) dir_list="Africa Australia Eurasia Islands North_America South_America"
			 gmt_format=3c ;;
		*) echo "ERROR: You did not specify the srtm version (correctly)"; exit 127 ;;
esac
echo -e " Ftp dirlist:\n \b$dir_list"


# define output files
outfile1=${project}/srtm_${project}.grd
outfile2=${project}/srtm30_${project}_merged.grd
outfile3=${project}/srtm30_${project}.grd
outfile4=${project}/final_${project}.grd
outfile5=${project}/final_${project}.dem
 
# download srtm and merge the tiles
echo ""
echo "--------------------------------------------------"
echo "Downloading srtm and merging the tiles ..."
echo "--------------------------------------------------"
echo ""
countb=1
for ((long=$west; long <= $east; long++))
do
  counta=1
  for ((lat=$south; lat <= $north; lat++))
  do
    long1=$long
    lat1=$lat

    if [ $long1 -lt 0 ] && [ $lat1 -lt 0 ]
    then 
      let "long1 = (( 0 - $long1 ))"
      let "lat1 = (( 0 - $lat1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=S${lat1}W${long1}.hgt

    elif [ $long1 -lt 0 ] && [ $lat1 -ge 0 ]
    then 
      let "long1 = (( 0 - $long1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=N${lat1}W${long1}.hgt

    elif [ $long1 -ge 0 ] && [ $lat1 -lt 0 ]
    then 
      let "lat1 = (( 0 - $lat1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=S${lat1}E${long1}.hgt

    elif [ $long1 -ge 0 ] && [ $lat1 -ge 0 ]
    then
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=N${lat1}E${long1}.hgt
    fi

    echo "Downloading" $file "..."
    if [ ! -e "$file" ]
    then
      # download the tile
      for dir in $dir_list
        do
        remotedir=srtm/version2_1/${srtm}/${dir}
        echo "Checking" ${remotedir} "..."
        downloadFile $remotedir ${file}.zip $ftp_link_1 $ftp_user $ftp_pass

        if [ -e "${file}.zip" ]
        then
          echo "File found"
          unzip -o ${file}.zip
          rm ${file}.zip
          break
        fi
      done
    else
      echo "File already existed."
    fi
    
    # merge the tiles in latitude direction using GMT
    GMT xyz2grd $file -G${project}/tile.grd -I$gmt_format -R$long/$[$long+1]/$lat/$[$lat+1] -ZTLhw
    if [ "$counta" = "1" ]; then
	mv ${project}/tile.grd ${project}/srtm_$long.grd
    else
	GMT grdpaste ${project}/tile.grd ${project}/srtm_$long.grd -G${project}/srtm_$long.grd
    fi
    counta=$[$counta+1]

  done
  
  # merge the tiles in longitude direction using GMT
  if [ "$countb" = "1" ]; then
      mv ${project}/srtm_$long.grd $outfile1
  else
      GMT grdpaste ${project}/srtm_$long.grd $outfile1 -G$outfile1
  fi
  countb=$[$countb+1]

done


# determine strm30 tile(s) (to fill voids)
# based on the 4 corners of the area
# 
# one extra degree is added (or subtracted) to avoid problems
# on the borders of srtm30 (caused by the different sampling
# rate and grid start of srtm1/3 and srtm30)
cnr_long[0]=$[$west-1]
cnr_long[1]=$[$west-1]
cnr_long[2]=$[$east+2]
cnr_long[3]=$[$east+2]
cnr_lat[0]=$[$south-1]
cnr_lat[1]=$[$north+2]
cnr_lat[2]=$[$south-1]
cnr_lat[3]=$[$north+2]

for ((v=0; v<=3; v++))
do
  for ((Xlat=-60; Xlat<=90; Xlat=Xlat+50))
  do
    let "temp = $Xlat-${cnr_lat[${v}]}"
    if [ $temp -ge 0 ] && [ $temp -lt 50 ]
    then
      Xlat1[${v}]=$Xlat
      if [ Xlat1[${v}] != -60 ]
      then
        for ((Xlong=-180; Xlong<=140; Xlong=Xlong+40))
        do
          let "temp = $Xlong-${cnr_long[${v}]}"
          if [ $temp -ge -40 ] && [ $temp -lt 0 ]
          then
            Xlong1[${v}]=$Xlong
          fi
        done
      else
        for ((Xlong=-180; Xlong<=120; Xlong=Xlong+60))
        do
          let "temp = $Xlong-${cnr_long[${v}]}"
          if [ $temp -ge -40 ] && [ $temp -lt 0 ]
          then
            Xlong1[${v}]=$Xlong
          fi
        done
      fi
    fi
  done
done

# determine the unique tile(s)
Xlat=`echo ${Xlat1[*]} | ${AWK} '{for (v = 1; v<=4; v++) print $v}' | sort | uniq | sort -n` 
Xlong=`echo ${Xlong1[*]} | ${AWK} '{for (v = 1; v<=4; v++) print $v}' | sort | uniq | sort -n` 

# download srtm30 (and merge the tiles)
echo ""
echo "--------------------------------------------------"
echo "Downloading srtm30 and merging the tiles ..."
echo "--------------------------------------------------"
echo ""
countb=1
for long in $Xlong
do
  counta=1
  for lat in $Xlat
  do
    long1=$long
    lat1=$lat

    if [ $long1 -lt 0 ] && [ $lat1 -lt 0 ]
    then 
      let "long1 = (( 0 - $long1 ))"
      let "lat1 = (( 0 - $lat1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=W${long1}S${lat1}.DEM
      file2=w${long1}s${lat1}

    elif [ $long1 -lt 0 ] && [ $lat1 -ge 0 ]
    then 
      let "long1 = (( 0 - $long1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=W${long1}N${lat1}.DEM
      file2=w${long1}n${lat1}

    elif [ $long1 -ge 0 ] && [ $lat1 -lt 0 ]
    then 
      let "lat1 = (( 0 - $lat1 ))"
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=E${long1}S${lat1}.DEM
      file2=e${long1}s${lat1}

    elif [ $long1 -ge 0 ] && [ $lat1 -ge 0 ]
    then
      long1=`echo $long1 | ${AWK} '{printf ("%03d",$1)}'`
      lat1=`echo $lat1 | ${AWK} '{printf ("%02d",$1)}'`
      file=E${long1}N${lat1}.DEM
      file2=e${long1}n${lat1}
    fi
    
    echo "Downloading" $file "..."
    if [ ! -e "$file" ]
    then
      # download the tile
      remotedir=srtm/version2_1/SRTM30/$file2
      downloadFile $remotedir ${file2}.dem.zip $ftp_link_2 $ftp_user $ftp_pass

      if [ -e "${file2}.dem.zip" ]
      then
        echo "File found"
        unzip -o ${file2}.dem.zip
        rm ${file2}.dem.zip
      fi
    else
      echo "File already existed."
    fi
  
    # merge the tiles in latitude direction using GMT
    west1=$(echo "$long+0.004166667" | bc)
    east1=$(echo "$long+40-0.004166667" | bc)
    south1=$(echo "$lat-50+0.00416667" | bc)
    north1=$(echo "$lat-0.004166667" | bc)
    GMT xyz2grd $file -G${project}/tile.grd -I30c -R$west1/$east1/$south1/$north1 -ZTLhw
    if [ "$counta" = "1" ]; then
	mv ${project}/tile.grd ${project}/srtm30_$long.grd
        south_first=$south1
    else
        south2=$(echo "$lat-50-0.004166667"| bc)
        GMT grdmath -R$west1/$east1/$south2/$south1 -I30c 0 0 NAN = ${project}/dummy.grd
        GMT grdpaste ${project}/dummy.grd ${project}/tile.grd -G${project}/tile.grd
	GMT grdpaste ${project}/tile.grd ${project}/srtm30_$long.grd -G${project}/srtm30_$long.grd
    fi
    counta=$[$counta+1]

  done
  
  # merge the tiles in longitude direction using GMT
  if [ "$countb" = "1" ]; then
      mv ${project}/srtm30_$long.grd $outfile2
  else
      west2=$(echo "$long-0.004166666" | bc)
      GMT grdmath -R$west2/$west1/$south_first/$north1 -I30c 0 0 NAN = ${project}/dummy.grd
      GMT grdpaste ${project}/dummy.grd ${project}/srtm30_$long.grd -G${project}/srtm30_$long.grd
      GMT grdpaste ${project}/srtm30_$long.grd $outfile2 -G$outfile2
  fi
  countb=$[$countb+1]

done

echo ""
echo "--------------------------------------------------"
echo "Filling voids ..."
echo "--------------------------------------------------"
echo ""

# resample to the same resolution as SRTM data (-I0.000833333) and trimmed to the same size
GMT grdsample $outfile2 -G$outfile3 -I$gmt_format -R${west}/$[${east}+1]/${south}/$[${north}+1]

# define NaN in SRTM (NAN=-32768):
GMT grdmath  $outfile1 -32768 NAN = ${project}/nan.grd

#void fill with AND command (AND 2 NaN if A and B == NaN, B if A == NaN, else A)

GMT grdmath ${project}/nan.grd $outfile3 AND = $outfile4


# write result to binary file (which Doris can read)
GMT grd2xyz $outfile4 -Zf > $outfile5

echo ""
echo "--------------------------------------------------"
echo "Creating output ..."
echo "--------------------------------------------------"
echo ""

# write Doris input lines to file (input.doris_comprefdem)
input_doris=${project}/input.doris_${project}
xmin=`GMT grdinfo ${outfile4} | grep x_min | sed 's/.*x_min: //g' | sed 's/x_max.*//g'`
ymax=`GMT grdinfo ${outfile4} | grep y_max | sed 's/.*y_max: //g' | sed 's/y_inc.*//g'`
xinc=`GMT grdinfo ${outfile4} | grep x_inc | sed 's/.*x_inc: //g' | sed 's/name.*//g'| sed 's/units.*//g'`
yinc=`GMT grdinfo ${outfile4} | grep y_inc | sed 's/.*y_inc: //g' | sed 's/name.*//g'| sed 's/units.*//g'`
Nx=`GMT grdinfo ${outfile4} | grep nx | sed 's/.*nx: //g'`
Ny=`GMT grdinfo ${outfile4} | grep ny | sed 's/.*ny: //g'`
dempath=`pwd`

echo -e "# The processing cards generated by $(basename $0) script." > $input_doris
echo -e "# Using parameters: $@" >> $input_doris
echo -e '# Copy the section(s) that is/are necessary to your processing setup.\n' >> $input_doris
echo "c         ___             ___" >> $input_doris
echo "comment   ___SIM AMPLITUDE___" >> $input_doris
echo "c                            " >> $input_doris
echo "SAM_IN_DEM     $dempath/$outfile5" >> $input_doris
echo -e "SAM_IN_FORMAT   r4 \t\t\t // default is short integer"  >> $input_doris
echo "SAM_IN_SIZE    $Ny $Nx" >> $input_doris
echo "SAM_IN_DELTA   $yinc $xinc" >> $input_doris
echo "SAM_IN_UL      $ymax $xmin"  >> $input_doris
echo "SAM_IN_NODATA  -32768" >> $input_doris
echo -e "SAM_OUT_FILE   master.sam \t // master simulated amplitude" >> $input_doris
echo -e "# SAM_OUT_DEM_LP   master_demhei_lp.raw \t // radarcoded dem to master extend" >> $input_doris
echo -e "# SAM_OUT_THETA_LP  master_theta_lp.raw \t // radarcoded dem to master extend" >> $input_doris
echo " " >> $input_doris
echo " " >> $input_doris
echo "c         ___          ___" >> $input_doris
echo "comment   ___DEM ASSIST___" >> $input_doris
echo "c                            " >> $input_doris
echo "DAC_IN_DEM     $dempath/$outfile5" >> $input_doris
echo -e "DAC_IN_FORMAT   r4 \t\t\t // default is short integer"  >> $input_doris
echo "DAC_IN_SIZE    $Ny $Nx" >> $input_doris
echo "DAC_IN_DELTA   $yinc $xinc" >> $input_doris
echo "DAC_IN_UL      $ymax $xmin"  >> $input_doris
echo "DAC_IN_NODATA  -32768" >> $input_doris
echo " " >> $input_doris
echo " " >> $input_doris
echo "c         ___             ___" >> $input_doris
echo "comment   ___REFERENCE DEM___" >> $input_doris
echo "c                            " >> $input_doris
echo "## CRD_METHOD   DEMINTRPMETHOD" >> $input_doris
echo "CRD_IN_DEM     $dempath/$outfile5" >> $input_doris
echo -e "CRD_IN_FORMAT   r4 \t\t\t // default is short integer"  >> $input_doris
echo "CRD_IN_SIZE    $Ny $Nx" >> $input_doris
echo "CRD_IN_DELTA   $yinc $xinc" >> $input_doris
echo "CRD_IN_UL      $ymax $xmin"  >> $input_doris
echo "CRD_IN_NODATA  -32768" >> $input_doris
# echo -e "CRD_OUT_FILE   master_slave.crd \t // reference dem phase" >> $input_doris


dem_doris=${project}/dem.dorisin
echo "c DEM generated with the following command:" > $dem_doris
echo c $0 $@ >> $dem_doris
echo "c
comment  ___general options___
c
SCREEN          info                          // level of output to standard out
MEMORY          150                             // MB
BEEP            error                            // level of beeping
OVERWRITE                                       // overwrite existing files
BATCH                                           // non-interactive
PROCESS          COMPREFDEM
c                                              //
comment  ___the general io files___            //
c                                              //
LOGFILE         log.out                         // log file
M_RESFILE       master.res                      // parameter file
S_RESFILE       slavedem.res                      // parameter file
I_RESFILE       dem.res                         // parameter file
c
c ___ step comprefdem ___
CRD_OUT_FILE    /dev/null           // do not store, only radarcoded DEM is needed
CRD_OUT_DEM_LP  dem_radar.raw
" >> $dem_doris
echo "CRD_IN_DEM     $dempath/$outfile5" >> $dem_doris
echo -e "CRD_IN_FORMAT   r4 \t\t\t // default is short integer"  >> $dem_doris
echo "CRD_IN_SIZE    $Ny $Nx" >> $dem_doris
echo "CRD_IN_DELTA   $yinc $xinc" >> $dem_doris
echo "CRD_IN_UL      $ymax $xmin"  >> $dem_doris
echo "CRD_IN_NODATA  -32768" >> $dem_doris
echo "c"  >> $dem_doris
echo "STOP"  >> $dem_doris
# ---------------------



# ---------------------
# visualization of SRTM
# ---------------------


#global settings:
GMT gmtdefaults -D >.gmtdefaults
GMT gmtset DOTS_PR_INCH 1200 
GMT gmtset PAPER_MEDIA A2+
GMT gmtset PAGE_ORIENTATION portrait

# parameters
area=-R${west}/$[${east}+1]/${south}/$[${north}+1]
let "central_long = (( ($west+$[${east}+1])/2 ))"
projection=-Jt${central_long}/1:3000000
out=${project}/srtm_${project}.ps


GMT grdgradient $outfile4 -Ne0.6 -A45/315 -G${project}/gradient.grd

echo "-10000 150 10000 150" > ${project}/gray.cpt
#makecpt -T-50/2000/1 -Cz_per.cpt > pl_gtopo30.cpt

GMT grdimage $outfile4 -I${project}/gradient.grd $area $projection -C${project}/gray.cpt -K > $out

GMT pscoast $area $projection -O -Df  -W -I1/1p/0/0/255 -I2/0.5p/0/0/255 -I3/0.25p/0/0/255 -I4/0.25p/0/0/255 -Ic/0.25p/0/0/255 -Ba1g0/0.5g0 -N1/0.35tap >> $out

#MA Clean up
#keep srtm_tr.grd srtm30_tr_merged.grd final_tr.grd 
case $keep in
	y|Y) # Keep intermediate files
	     rm -rf $project/srtm_[0-9][0-9].grd  $project/nan.grd \
	     $project/dummy.grd $project/*.cpt $project/gradient.grd \
	     $project/srtm30_tr.grd $project/tile.grd
	     echo "--------------------------------------------------"
	     echo "Kept intermedia files!"
	     echo "--------------------------------------------------";;
	  *) # Full clean up, just .dem .ps input.doris_comprefdem
	     rm -rf ${project}/*.grd; rm -rf ${project}/*.cpt
	     echo -e "Outputs: \n\t$outfile5 \n\t$input_doris \n\t$out" ;;
esac


echo ""
echo "--------------------------------------------------"
echo "Done!"
echo "--------------------------------------------------"
echo "--------------------Start geoid correction--------"
echo $project
cd   $project
pwd

add_egm96.py
cd ..
echo "--------------------End geoid correction--------"
