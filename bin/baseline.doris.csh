#!/bin/csh -f
# baseline.fast
#%// BK 27-Oct-2000
#%// $Revision: 3.5 $  $Date: 2003/04/14 06:35:15 $
######################################################################
set PRG    = `basename "$0"`
set VER    = "v1.1, FMR software"
set AUT    = "Bert Kampes, (c)2000"
echo "$PRG $VER, $AUT"\\n


### Handle input
if ( $#argv < 2 ) then
#  cat << __EOFHD | more -de
  cat << __EOFHD | more

  PROGRAM:  $PRG -- Obtain baseline parameterization from orbits.

  SYNOPSIS:
    $PRG  master.res  slave.res

  OPTIONS:
    master.res   master result file from doris, precise orbits.
    slave.res    slave result file from doris, precise orbits.

  Basically a dummy run is performed, only to get DUMPBASELINE.
  Then from the standard output the appropriate sections are written.
  The DUMPBASELINE card is used, computing the baseline at a grid 15x10,
  then fitting a polynomial through these values (Bperp, theta and range).

  EXAMPLE:
    $PRG master.res slave.res

__EOFHD
  exit 1
endif

### Correct input.
set MASTERFILE = "$1"
set SLAVEFILE  = "$2"


### Create dummy input file(s).
set TMPDIR     = "/tmp"
set DORISIN    = $TMPDIR/$user.$$.in
set DORISOUT   = $TMPDIR/$user.$$.out
set LOGFILE    = $TMPDIR/$user.$$.log
set PRODFILE   = $TMPDIR/$user.$$.prod


###
cat << __EOFHD >! $DORISIN
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
__EOFHD


### Run Doris, module XXX, use DUMPBASELINE card.
doris $DORISIN > $DORISOUT

### Grep the results from redirected stdout.
### first obtain line number
set L = `grep -n compbaseline $DORISOUT | cut -f1 -d':' | tail -1`
### number of lines of interest
@ NL  = 22
### lastline
@ LL  = $L + $NL
@ NL++

### Display results.
clear
head -n $LL $DORISOUT | tail -n $NL


### Tidy up.
#rm -rf $TMPDIR
rm -f $DORISIN
rm -f $DORISOUT
rm -f $LOGFILE
rm -f $PRODFILE

### EOF.


