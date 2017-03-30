import logging
import os


from grocsvs import jobmanagers
from grocsvs import utilities

def make_jobmanager(jobmanager_settings, processes, batch_dir):
    jobmanager_classes = {"IPCluster":jobmanagers.IPCluster,
                          "local":    jobmanagers.LocalCluster,
                          "multiprocessing": jobmanagers.MultiprocessingCluster,
                          "admiral": jobmanagers.AdmiralCluster}

    cls = jobmanager_classes[jobmanager_settings.cluster_type]
    return cls(processes, jobmanager_settings, batch_dir)



class Runner(object):
    def __init__(self, options):
        self.options = options
        self.jobmanagers = None


    def run_stage(self, stage, stage_name):
        to_run = []
        all_steps = stage.get_steps(self.options)
        for step_chunk in all_steps:
            if step_chunk.needs_to_run():
                to_run.append(step_chunk)
                
        # should probably have an option to re-use jobmanagers started
        # earlier

        logging.info("{} {} {}".format("="*30, stage_name, "="*30))

        if len(to_run) > 0:
            logging.info("{} chunks to run. Starting...".format(len(to_run)))

            utilities.ensure_dir(to_run[0].results_dir)
            utilities.ensure_dir(to_run[0].working_dir)
            utilities.ensure_dir(os.path.join(self.options.log_dir, to_run[0].__class__.__name__))

            processes = min(self.options.cluster_settings.processes, len(to_run))
            jobmanagers = make_jobmanager(
                self.options.cluster_settings, processes, self.options.working_dir)
            jobmanagers.map(_run_chunk, to_run)

            error_not_done = []
            for step_chunk in all_steps:
                if step_chunk.needs_to_run():
                    error_not_done.append(step_chunk)
            if len(error_not_done) > 0:
                logging.error("Failed to complete the following steps:\n{}".format(
                    "\n".join(step_chunk.log_path for step_chunk in error_not_done)))
                raise Exception(("Stage {} failed to complete. Check the log files " +
                                 "above for more information.").format(stage_name))
                
            logging.info("--> {} completed.\n".format(stage_name))
        else:
            logging.info("--> 0 chunks need to be run. Skipping...\n")


def _run_chunk(chunk):
    try:
        chunk.start_logging()
    except Exception, e:
        print "== Error starting logging =="
        raise

    
    utilities.check_memory(chunk.logger)

    try:
        chunk.run()
        chunk.finalize()
        chunk.stop_logging()
    except Exception, e:
        chunk.logger.exception(e)
        raise
        
    if chunk.needs_to_run():
        chunk.logger.error("Step {} failed to produce output.".format(chunk))
