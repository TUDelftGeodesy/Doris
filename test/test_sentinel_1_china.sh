#!/bin/bash

source_path='/data/src/Doris_s1_git'
python_source_path=$source_path'/sentinel_1'
test_path=$source_path'/test/china_s1'

source activate python27
export PYTHONPATH=$python_source_path/main_code/:$python_source_path/functions/:$PYTHONPATH
python $python_source_path/main_code/doris_main.py -p $test_path -s 2016-01-13 -e 2016-02-06 -m 2016-01-13

