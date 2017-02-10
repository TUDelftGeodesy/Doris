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
        self.flag_dir = self.doris_parameters.doris_parallel_flag_dir
        self.between_sleep_time = self.doris_parameters.between_sleep_time
        self.end_sleep_time = self.doris_parameters.end_sleep_time
        self._cleanup()
        os.mkdir(self.flag_dir)

    def _cleanup(self):
        #
        # cleans stuff from last run
        # moves directory with start/finish flags to timestamped directory
        #
        if(self.verbose):
            os.system("mv "
                      + self.flag_dir
                      + " "
                      + self.flag_dir
                      + "."
                      + time.asctime(time.localtime(time.time())).replace(" ", "_"))
        else:
            os.system("rm -rf "
                      + self.flag_dir)

    def _get_job_id(self, job):
        #
        # returns machine level unique job Id
        #
        id = job[0].replace("/", "") + "." + job[1].replace("/", "").replace(" ", "") + self.pid
        if len(id.split('/')[-1]) > 255:  # If the string becomes to long, we reduce the long paths to max. 40 characters.
            dat = [a[- min(40, len(a)):] + ' ' for a in job[1].split(' ')]
            dat = ' '.join(dat)
            id = job[0].replace("/", "") + "." + dat.replace("/", "").replace(" ", "") + self.pid

        return id

    def _start_jobs(self, jobs, active_jobs, job_count):
        #
        # starts a number of jobs
        # returns list of unstarted and list of started jobs
        #
        jobs_to_start_count = min(job_count, len(jobs))
        for index in range(0, jobs_to_start_count):
            job = jobs.pop(0)
            os.chdir(job[0])
            os.system(self.doris_parameters.job_handler_script + " "
                      + self.flag_dir + " "
                      + self._get_job_id(job) + " "
                      + str(self.verbose) + " "
                      + job[1] + " &")
            active_jobs.append(job)
        return [jobs, active_jobs]

    def _active_jobs(self, jobs):
        #
        # returns from the list of jobs, the jobs that are started, but not finished
        #
        jobs_active = []
        for job in jobs:  # find active jobs
            this_job_started = False
            this_job_ready = False
            for file in os.listdir(self.flag_dir):
                if self._get_job_id(job) + ".finished" == file:
                    this_job_ready = True
                if self._get_job_id(job) + ".started" == file:
                    this_job_started = True
            this_job_active = this_job_started & (not (this_job_ready))
            if (this_job_active):
                jobs_active.append(job)
        return jobs_active

    def run(self, job_list):
        """executes joblist in parallel
        jobList: list of jobs, containing execution path and command to be executed
        """
        [job_list, active_job_list] = self._start_jobs(job_list, [], self.max_jobs)
        while len(active_job_list):
            if(self.verbose):
                print time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + "jobs busy"
            time.sleep(self.between_sleep_time)
            active_job_list = self._active_jobs(active_job_list)
            [job_list, active_job_list] = self._start_jobs(job_list, active_job_list, self.max_jobs - len(active_job_list))
        if (self.verbose):
            print time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + "jobs finished"
        time.sleep(self.end_sleep_time)
