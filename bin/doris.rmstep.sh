#!/bin/bash
##
##
## AUTHOR:   Mahmut Arikan
## EMAIL:    M.Arikan@TUDelft.nl
## VERSION:  v1.5
## DATE:     20060612
## UPDATE:   20090519
##
## TUDelft, DEOS Radar Group  - 2006
##
##
## reset doris steps in the result files,
## up to the processing step you want
##

#sed -in -e '/coarse_correl/,/End_process_control/s/1/0/gp' -e '/Start_coarse_correl/,$d' 08956_09958.res
#$d indicates until last line.

case $# in
        2) #    process_name                                              res files
           lineno=$(awk '/_Start_'$1'/{print NR-1}' $2);  # get me the line number of the process line - 1
           [[  "$lineno" == "" ]] && echo -e "No such entry: $1 in $2 .\n"  && exit 127

           if [ "`uname`" == "Darwin" ]; then     # MAC_OSX
             sed -i .bak2 -e '/^'$1'/,/End_process_control/s/1/0/g' -e ''${lineno}',$ d' $2
           else                             # other platforms
             sed -i -e '/^'$1'/,/End_process_control/s/1/0/g' -e ''${lineno}',$ d' $2
           fi

           #sed -i -e '/^'$1'/,/End_process_control/s/1/0/g' -e '/Start_'$1'/,$d' $2 ;;
           ;;
        *) echo -e "\nUSAGE: ${0##*/} <doris_process_name> <.res> \n";
           echo -e "   EX: ${0##*/} coarse_correl   master_slave.res \n ";
           echo -e "       ${0##*/}      resample   slave.res \n ";
           echo -e "Doris process names by order:
 1.  precise_orbits      14. resample
 2.  crop                15. interfero
 3.  sim_amplitude       16. comp_refphase
 4.  master_timing       17. subtr_refphase
 5.  filt_azi            18. comp_refdem
 6.  filt_range          19. subtr_refdem
 7   oversample          20. coherence
 8.  coarse_orbits       21. filtphase
 9.  coarse_correl       22. unwrap
10.  fine_coreg          23. slant2h
11.  timing_error        24. geocoding
12.  dem_assist          25. dinsar
13.  comp_coregpm        26. <extra>     
                          
                                               " ;
echo -e "\n Thank you for using Doris!\n  TU Delft - DEOS Radar Group 2011 \n";
exit 0;;
esac

