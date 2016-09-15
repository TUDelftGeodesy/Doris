import os
import sys

class DorisParameters_Path(object):

    def set(self, doris_parameters_path):
            if(os.path.exists(doris_parameters_path)):
                sys.path.append(doris_parameters_path)
                print 'dorisparameter path: ' + doris_parameters_path
            else:
                print 'dorisparameter path: ' + doris_parameters_path + ' not a valid path'
