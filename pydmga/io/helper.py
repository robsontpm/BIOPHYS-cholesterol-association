# coding: utf-8
"""
This module gives some helper files to create analysis scripts
of the similar structure to each other.
"""
        

import importlib
import datetime
import os

def load_settings(pathname, base="."):
    settings = importlib.import_module(pathname, ".")    
    return settings


def assure_clean_wd(output_dirpath):
    if os.path.exists(output_dirpath):
        outf_backup = output_dirpath
        counter = 0;
        while os.path.exists(outf_backup):    
            outf_backup = "{0}_{1}.bak".format(output_dirpath, counter)
            counter += 1
        os.rename(output_dirpath, outf_backup)

class Job:
    '''
    Common usage:

    from pydmga.io import helper
    if len(sys.argv) > 1:
        job = helper.Job(sys.argv)
    ...
    job.log("Some info")
    ...
    job.finish()
    '''    
    def __init__(self, setup, confbasename=".", app=False):
        self.end_time = None
        self.start_time = datetime.datetime.now()        
        self.app = app or setup[0]
        self.pathname = setup[1]
        self.jobname = self.pathname
        self.settings = importlib.import_module(self.pathname, confbasename)
        self.workdir = self.settings.OUTPUT_DIRPATH                
        do_backup = (len(setup) < 3) or ("no-backup" not in setup)
        if do_backup:
            assure_clean_wd(self.workdir)
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir) 
        self.log_filepath = os.path.join(self.workdir, "log-{}.txt".format(self.app))
        self.log_file = file(self.log_filepath, "w")             
        self.info = "# Generated with:\n#  \n#     {}\n#\n# Started job '{}' on {}\n".format(
            " ".join(setup), 
            self.jobname, str(self.start_time),
        )
        self.log_file.write(self.info)
        self.user_files = []  

    def is_finished(self):
        return bool(self.end_time)

    def log(self, what, sep='\n'):
        if self.is_finished():
            return
        self.log_file.write("{}{}".format(what, sep))

    def filepath(self, filename, extra_format_dict = {}):
        frmt = self.as_dict()
        frmt.update(extra_format_dict)        
        return os.path.join(self.workdir, filename.format(**frmt))
        
    def finish(self, info=None):
        end_time = datetime.datetime.now()
        if info:
        	self.log(info)	
        self.log("# Finished job '{}' on {}".format(self.jobname, str(end_time)))
        self.log_file.close()
        self.end_time = end_time

    def __del__(self):
        if self.is_finished():
            return
        self.finish()    

    def user_file(self, filename, extra_format_dict = {}, mode="w"):
    	if self.is_finished():
            return None
        self.user_files.append(file(self.filepath(filename, extra_format_dict), mode))
        return self.user_files[-1]

    def as_dict(self):
        return {
            "end_time": self.end_time,
            "start_time": self.start_time,
        	"app": self.app,
		    "pathname": self.pathname,
		    "jobname": self.jobname,
		    "settings": self.settings,
		    "workdir": self.workdir,
        }


# def jobstart(setup):
#     '''
#     Common usage:

#     if len(sys.argv) > 1:
#         jobname, settings, jobstart, workdir, log_file = helper.jobstart(sys.argv)
#     '''
#     pathname = setup[0];
#     do_backup = (len(setup) < 2) or (setup[1] != "no-backup")
#     settings = load_settings(pathname)
#     wd = settings.OUTPUT_DIRPATH
#     if do_backup:
#         assure_clean_wd(wd)
#     os.mkdir(wd)   
#     log_filepath = os.path.join(wd, "log-{}.txt".format(setup[0]))
#     log_file = file(log_filepath, "w")     
#     starttime = datetime.datetime.now()
#     log_file.write("Started job '{}' on {}\n".format(pathname, str(starttime))
#     return (pathname, settings, starttime, wd, log_file)


    
