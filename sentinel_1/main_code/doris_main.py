import argparse
from datetime import datetime
from doris_sentinel_1 import DorisSentinel1
"""Doris processing
arguments:  --parameterfilepath, -p
            --startdate -s
            --enddate, -e
            --masterdate, -m
"""

def valid_date(s):
    try:
        return datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

# parse arguments here

parser = argparse.ArgumentParser(description='Doris processing.')
parser.add_argument('--parameterfilepath', '-p', default='./',
                    help='Path to dorisParameter.py file, this file contains case specific parameters')
parser.add_argument('--startdate', '-s', type=valid_date,
                    help='start date of stack to be processed')
parser.add_argument('--enddate', '-e', type=valid_date,
                    help='end date of stack to be processed')
parser.add_argument('--masterdate', '-m', type=valid_date,
                    help='date of master of stack to be processed')

args = parser.parse_args()

#start doris sentinel1 run
doris_sentinel_1 = DorisSentinel1()
doris_sentinel_1.run(args.parameterfilepath, args.startdate, args.enddate, args.masterdate)


