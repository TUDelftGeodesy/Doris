import argparse
from .prepare_datastack import PrepareDatastack

"""Doris prepare datastack
arguments:  --doris_input_file, -i
"""

# parse arguments here

parser = argparse.ArgumentParser(description='Doris prepare datastack.')
parser.add_argument('--doris_input_file', '-i', default='',
                    help='Path to doris input file, this file contains case specific parameters')

args = parser.parse_args()

#start doris sentinel1 run
prepare_data_stack = PrepareDatastack()
prepare_data_stack.prepare(args.doris_input_file)


