import os
import time

class Jobs(object):
    """The Jobs class runs a list of jobs in parallel.
       It starts the maximum number of jobs from the list in parallel and monitors for job completion.
       When jobs are finished, new jobs are started until the maximum is reached again.
       This is repeated until all jobs from the list are processed.

       This class executes a jobHandlerScript to execute a job.
       The job handler script sets job start and job finished flags in a flag directory on the system.
       In verbose mode the job handler script prints status info to stdout

       This class generates directories on the system that contain start and finish flags for each job that is run
       Old flag directories are moved to timestamped directories

       Methods:
           Run: executes list of jobs
           """
    def __init__(self, max_jobs, dorisParameters):
        """max_jobs: maximum number of jobs to run simultaniously
           verbose: print status to stdout during execution of job list"""
        self.max_jobs = max_jobs
        self.pid = str(os.getpid())
        self.doris_parameters = dorisParameters
        self.verbose = self.doris_parameters.verbose
        self.flag_dir_root = self.doris_parameters.doris_parallel_flag_dir
        self.between_sleep_time = self.doris_parameters.between_sleep_time
        self.end_sleep_time = self.doris_parameters.end_sleep_time
        self.python_path = os.path.dirname(self.doris_parameters.source_path)
        self.jobs_todo = []
        self.jobs_active = []
        self.jobs_finished = []
        self.flag_dir = ''

    def _new_id(self):
        # static class variable
        Jobs.id = Jobs.id + 1
        return Jobs.id

    def _set_id(self, job_list):
        for job_dict in job_list:
            job_dict['id'] = self._new_id()

    def _create_flag_dir(self):
        self.flag_dir = self.flag_dir_root + "." + time.asctime(time.localtime(time.time())).replace(" ", "_")
        os.mkdir(self.flag_dir)


    def _cleanup_flag_dir(self):
        #
        # cleans stuff from current run, if not verbose
        #
        if(not (self.verbose)):
            os.system("rm -rf " + self.flag_dir)

    def _get_job_id(self, job):
        #
        # returns machine level unique job Id
        return job['path'].replace("/","_") + "." + job['command'].replace("/","_").replace(" ","-") + self.pid


    def _start_jobs(self):
        #
        # starts a number of jobs
        # returns list of unstarted and list of started jobs
        #
        jobs_to_start_count = min((self.max_jobs - len(self.jobs_active)), len(self.jobs_todo))
        for index in range(0, jobs_to_start_count):
            job = self.jobs_todo.pop(0)
            os.chdir(job['path'])
            os.system(self.doris_parameters.job_handler_script + " "
                      + self.flag_dir + " "
                      + str(job['id']) + " "
                      + self._get_job_id(job) + " "
                      + str(self.verbose) + " "
                      + self.python_path + " "
                      + job['command'] + " &")
            self.jobs_active.append(job)
        return

    def _check_active_jobs(self):
        #
        # returns from the list of jobs, the jobs that are started, but not finished
        #
        jobs_active = []
        for job in self.jobs_active:  # find active jobs
            this_job_started = False
            this_job_ready = False
            for file in os.listdir(self.flag_dir):
                if str(job['id']) + ".finished" == file:
                    this_job_ready = True
#                if str(job['id']) + ".started" == file:
#                    this_job_started = True
#            this_job_active = this_job_started & (not (this_job_ready))
#            if (this_job_active):
            if (not (this_job_ready)):
                jobs_active.append(job)
        self.jobs_active = jobs_active
        return

    def run(self, job_list):
        """executes joblist in parallel
        jobList: list of jobs, containing execution path and command to be executed
        """
        self._create_flag_dir()
        self.jobs_todo = job_list
        self._set_id(self.jobs_todo)
        self._start_jobs()
        while len(self.jobs_active):
            if(self.verbose):
                print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + "jobs busy")
            time.sleep(self.between_sleep_time)
            self._check_active_jobs()
            self._start_jobs()
        if (self.verbose):
            print(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + "jobs finished")
        time.sleep(self.end_sleep_time)
        self._cleanup_flag_dir()
