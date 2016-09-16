'''
	GrsConfig defines paths that are local to the source tree.
	They are copied into DorisParameters for use in Doris python scripts
'''

class GrsConfig(object):

    def __init__(self):
        self.source_path = '/home/dlevelt/src/Doris_s1_git'
        self.doris_path = self.source_path + '/bin/doris'
        self.cpxfiddle_path = self.source_path + '/sar_tools/'
        self.job_handler_script = self.source_path + "/sentinel_1/main_code/jobHandlerScript"
