#!/usr/bin/env python
import os,sys,time

def usage():
    print('\nUsage: python hhmmss2sec.py time')
    print('  where time in the form HH:MM:SS.sss .')
    print(' ')
    print(' Example ')
    print(' ./hhmmss2sec.py 16:36:40.393 ')
    print('    59800.393000 ')
    print(' ')
    print(' See Also')
    print(' sec2hhmmss.py')
    
try:
    timeOfDay  = sys.argv[1] 
except:
    print('Unrecognized input')
    usage()
    sys.exit(1)

timeOfDay= sys.argv[1].split(':')
hh  = float(timeOfDay[0]);
mm  = float(timeOfDay[1]);
ss  = float(timeOfDay[2]);

secOfDay=hh*3600+mm*60+ss
print( "%f" %(secOfDay))
