import abc
import datetime
import os
import shutil
import time

from grocsvs import log


class StepChunk(object):
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def get_steps(options):
        """ """
        raise Exception("this abstract staticmethod needs to be instantiated by subclasses")

    @abc.abstractmethod
    def __init__(self, options, **kwdargs):
        """
        must take as input the options instance and a dict of arguments
        """
        pass

    @abc.abstractmethod
    def outpaths(self, final):
        """
        return a dictionary of names -> output paths; final indicates if the
        paths should be for the temporary, working version (final=False) or the
        completed final version (final=True)
        """
        return

    @abc.abstractmethod
    def run(self):
        """
        actually runs the pipeline step on the arguments, given the options, 
        defined in __init__
        """
        pass

    def __str__(self):
        raise Exception("not implemented")

    @property
    def results_dir(self):
        return self._results_dir(self.options)#return os.path.join(self.options.results_dir, self.__class__.__name__)
    
    @classmethod
    def _results_dir(cls, options):
        # private method to allow classmethod-level access to results_dir
        return os.path.join(options.results_dir, cls.__name__)    

    @property
    def working_dir(self):
        return self._working_dir(self.options)#os.path.join(self.options.working_dir, self.__class__.__name__)

    @classmethod
    def _working_dir(cls, options):
        return os.path.join(options.working_dir, cls.__name__)    


    @property
    def log_path(self):
        return os.path.join(self.options.log_dir, self.__class__.__name__, str(self))

    def start_logging(self):
        self.logger = log.Logger(self.log_path)
        self.logger.log("-- starting logging {} --".format(str(self)))
        self._start_time = time.time()

    def stop_logging(self):
        elapsed = time.time() - self._start_time
        self.logger.log("-> finished running step; time elapsed: {}".format(datetime.timedelta(seconds=elapsed)))
        self.logger.log("-- stopping logging --")
        

    @classmethod
    def clean_all_steps(cls, options):
        import logging
        results_dir = cls._results_dir(options)
        if os.path.exists(results_dir):
            logging.info("Removing {}".format(results_dir))
            shutil.rmtree(results_dir)

        working_dir = cls._working_dir(options)
        if os.path.exists(working_dir):
            logging.info("Removing {}".format(working_dir))
            shutil.rmtree(working_dir)

        log_dir = os.path.join(options.log_dir, cls.__name__)
        if os.path.exists(log_dir):
            logging.info("Removing {}".format(log_dir))
            shutil.rmtree(log_dir)


    def needs_to_run(self):
        """ checks if any of the output files are missing """
        paths = self.outpaths(final=True)
        for name, path in paths.items():
            if not os.path.exists(path):
                return True

        return False

    def finalize(self):
        temp_paths = self.outpaths(final=False)
        final_paths = self.outpaths(final=True)

        for key, working_path in temp_paths.items():
            if working_path != final_paths[key]:
                if not os.path.exists(working_path):
                    raise Exception("{} step failed to produce output: {}".format(
                        self.__class__.__name__, working_path))
                shutil.move(working_path, final_paths[key])
                time.sleep(0.1)

        assert not self.needs_to_run()
