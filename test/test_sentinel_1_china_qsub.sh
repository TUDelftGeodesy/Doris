#!/bin/bash

#PBS -l nodes=1:ppn=4

source_path='/home/dlevelt/src/Doris_s1_git'
python_source_path=$source_path'/sentinel_1'
test_path=$source_path'/test/china_s1'

export PATH=$PATH:/home/everybody/bin/snaphu:/home/dlevelt/src/Doris_s1_git/sar_tools
export PYTHONPATH=/home/fjvanleijen/software/s1_doris/tops-toolbox/tops-reader-python:/home/everybody/python/py_modules/lib64/python2.7/site-packages::/home/everybody/bin/deinsar/deinsar_v0.1.4:$python_source_path/main_code/:$python_source_path/functions/:$PYTHONPATH
python $python_source_path/main_code/doris_main.py -p $test_path -s 2016-01-13 -e 2016-02-06 -m 2016-01-13 2>&1

