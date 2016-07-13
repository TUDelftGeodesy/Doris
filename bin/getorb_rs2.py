#!/usr/bin/env python

#from IPython.Debugger import Tracer; debug_here = Tracer()
#
#-----------------------------------------------------------------#
# A python code for parsing RS2 XML file into python data structures
# and from there into DORIS res file structure
#
# Author: TUDelft - 2011
# Maintainer: Piers Van der Torren
# License: GPL
#
#-----------------------------------------------------------------#

import os, re
from datetime import datetime, timedelta
#import argparse
import optparse

codeRevision=1.0   # this is the code revision number

def orbit_MDA(filename):
    """read MDA orbitfile and return generator with points
    """
    with open(filename) as orbfile:
        line = orbfile.readline().strip()
        while line != '':
            if len(line) == 24 and not line.startswith(';'):
                try:
                    pnt = [ datetime.strptime(line, '%Y-%j-%H:%M:%S.%f') ]
                    pnt.append(orbfile.readline().strip())
                    pnt.append(orbfile.readline().strip())
                    yield pnt
                except ValueError:
                    pass
            line = orbfile.readline().strip()

def orbit_filter(orbit,fromdate,todate):
    """filter orbitpoints and convert to DORIS orbit format
    """
    for pnt in orbit:
        if fromdate < pnt[0] < todate:
            yield('{0} {1} {2}'.format(pnt[0].hour*3600+pnt[0].minute*60+pnt[0].second+pnt[0].microsecond/1000000.0, pnt[1], pnt[2]))

def load_resfile(filename):
    """Load DORIS result file as nested dictionary with sections and key-value pairs
    """
    with open(filename) as resfile:
        resdata = resfile.read()

    # split res string in structures of the form:
    # .name: name of a brocessing block
    # .args: everything after _Start_<name>
    # .data: text between _Start_<name> and End_<name>
    # .status: everything after End_<name>
    regex = re.compile( r'\n[*]*[ \t]*\n[*]_Start_(?P<name>\w+):?[ \t]*' + \
                    r'(?P<args>[^\n]*)\n[*]*[ \t]*\n' + \
                    r'(?P<data>.*?)\n' + \
                    r'[*]*[ \t]*\n[*] End_\1:(?P<status>\w+)\n[*]*[ \t]*', re.DOTALL )
    
    # regex without star lines
    #regex = re.compile( r'\n[*]_Start_(?P<name>\w+):?[ \t]*' + \
    #                r'(?P<args>[^\n]*)\n' + \
    #                r'(?P<data>.*?)\n' + \
    #                r'[*] End_\1:(?P<status>\w+)\n', re.DOTALL )
        
    blocks = {}
    for block in regex.finditer( resdata ):
        blocks[block.group('name')] = block.group('data')
        for block2 in regex.finditer( block.group('data') ):
            blocks[block2.group('name')] = block.group('data')

    # split block.data in two fields:
    # #1=key, excluding the terminating colon
    # #2=value, leading whitespace stripped,
    # plus adjacent non key:value lines.
    regex = re.compile( r'\n(?P<key>[^\n:]+):' + \
                r'[ \t]*(?P<value>[^\n]*[^:]*)' + \
                r'(?=\n[^\n:]+:|\n[*])', re.DOTALL )

    res = {}
    for key in blocks:
        res[key] = {}
        for token in regex.finditer( blocks[key] ):
            pairkey = re.sub('[^a-zA-Z0-9]+', '', token.group('key'))
            res[key][pairkey] = token.group('value').strip()

    return res

def add_to_resfile(filename, section, data, args='', status='_NORMAL'):
    """Add new section to a DORIS result file.
    Write section header, data, and footer, and adjust process control flag.
    """
    with open(filename) as resfile:
        resdata = resfile.read()
    process_control = re.search(r'\n{0}:[ \t]*[01]'.format(section), resdata).span()[1]
    if resdata[process_control-1] != '0':
        raise Exception('already has precise orbits according to process control flags')

    header = '*'*67 + '\n*_Start_' + section + ':' + args + '\n' + '*'*67 + '\n'
    footer = '*'*67 + '\n* End_' + section + ':' + status + '\n' + '*'*67 + '\n'

    resdata = resdata[:process_control-1] + '1' + resdata[process_control:] + \
                        '\n' + header + data + '\n' + footer

    with open(filename, 'w') as resfile:
        resfile.write(resdata)

def clear_resfile_section(filename, section):
    """Clear a section in a DORIS result file.
    Remove the section, and set process control flag to 0 if it exists.
    """
    with open(filename) as resfile:
        resdata = resfile.read()

    process_control = re.search(r'\n{0}:[ \t]*[01]'.format(section), resdata)
    if process_control:
        process_control = process_control.span()[1]
        if resdata[process_control-1] == '1':
            resdata = resdata[:process_control-1] + '0' + resdata[process_control:]

    regex = re.compile( r'\n[*]*[ \t]*\n[*]_Start_{0}:?[ \t]*'.format(section) + \
                    r'(?P<args>[^\n]*)\n[*]*[ \t]*\n' + \
                    r'(?P<data>.*?)\n' + \
                    r'[*]*[ \t]*\n[*] End_{0}:(?P<status>\w+)\n[*]*[ \t]*\n'.format(section), re.DOTALL )

    sectionspan = regex.search(resdata).span()
    resdata = resdata[:sectionspan[0]] + resdata[sectionspan[1]:]

    with open(filename, 'w') as resfile:
        resfile.write(resdata)
    

def indate(str):
    """convert string to date, trying several formats
    """
    for format in ['%y%m%d%H%M%S', '%Y%m%d%H%M%S']:
        try:
            return datetime.strptime(str, format)
        except ValueError:
            pass
    raise ValueError('Date not understood: {0}'.format(str))

#def main():

# parse commandline arguments
parser = optparse.OptionParser(description='getorb-like tool extraction of radarsat orbits.')
parser.add_option('--dir',
                                     help='a directory containing RS2 orbits')
parser.add_option('--extratime', type=int, default=16,
                                     help='extra time before and after acquisition [s]')
parser.add_option('--resfile', 
                                     help='resfile to which to add precise orbits')
args = parser.parse_args()[0]

print('Getting precise orbits for Radarsat2\norbit dir: {dir}\nresfile: {resfile}\nextratime: {extratime}'\
    .format(dir=args.dir,extratime=args.extratime,resfile=args.resfile))
# read res file for neccesary metadata
res = load_resfile( args.resfile )

# get orbit number from line
# Scene identification:        Orbit: 08084    2009-07-02T05:52:34.496664Z 
orbit = int(res['readfiles']['Sceneidentification'][7:12])
# first pixel azimuth time
date = datetime.strptime(res['readfiles']['FirstpixelazimuthtimeUTC'], '%d-%b-%Y %H:%M:%S.%f')

orbitfile = os.path.join( args.dir, '{0:05d}_def.orb'.format(orbit) )
extratime = timedelta(seconds=args.extratime)
fromdate = date - extratime
todate = date + extratime + timedelta(seconds=4) # add 4 seconds as approximate acquisition time

# get orbit section
orbit_points = list(orbit_filter(orbit_MDA(orbitfile), fromdate, todate))

# save the orbit if more than 4 orbit points are found
if len(orbit_points) > 4:
    add_to_resfile( args.resfile, 'precise_orbits',
        't(s) X(m) Y(m) Z(m) X_V(m/s) Y_V(m/s) Z_V(m/s)\n' + \
        'NUMBER_OF_DATAPOINTS: {0}\n'.format(len(orbit_points)) + \
        '\n'.join(orbit_points) ) 

    # remove the leader orbit from the resfile
    clear_resfile_section( args.resfile, 'leader_datapoints')
else:
    print('No matching orbits found, tried orbit file {0}'.format(orbitfile))
