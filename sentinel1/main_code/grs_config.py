

class GrsConfig(object):

    def __init__(self):
        self.source_path = '/data/src/Doris_s1_git'
        self.doris_path = self.source_path + '/bin/doris'
        self.cpxfiddle_path = self.source_path + '/SARTools/'
        self.job_handler_script = self.source_path + "/sentinel1/main_code/jobHandlerScript"