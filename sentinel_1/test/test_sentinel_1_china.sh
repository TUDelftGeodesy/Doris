#!/bin/bash

python_source_path='/home/dlevelt/src/Doris_s1_git/sentinel_1'

source activate python27
export PYTHONPATH=$python_source_path/main_code/:$python_source_path/functions/:$PYTHONPATH
python ../main_code/doris_main.py -p //home/dlevelt/Projects//china_s1/ -s 2016-01-13 -e 2016-02-06 -m 2016-01-13

