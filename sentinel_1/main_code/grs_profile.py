import sys
import time

class GRS_Profile():

    def __init__(self, logfile, verbose):
        self.verbose = verbose
        if(verbose):
            self.start_time = time.localtime()
            self.logfile = logfile + "." + time.strftime("%a, %d %b %Y %H:%M:%S +0000",self.start_time).replace(" ", "_").replace(",", "").replace("+", "")

    def log_time_stamp(self, logtxt):
        if(self.verbose):
            fileHandle = open(self.logfile, 'a')
            message = time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime()) + ' : ' + logtxt + '\n'
            fileHandle.write(message)
            fileHandle.close()
