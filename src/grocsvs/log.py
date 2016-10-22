import datetime
import os
import subprocess
import sys
import traceback

from grocsvs import utilities


class Logger(object):
    """
    This replaces the built-in logging module, which doesn't work reliably 
    across ipyparallel jobs
    """

    def __init__(self, log_path):
        self.log_path = log_path
        self.log_file = None
        self.open_log()

    def open_log(self):
        if self.log_file is None:
            self.log_file = open(self.log_path, "a")

    def log(self, message):
        date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        message_str = "{} - {}\n".format(date_str, message)

        self.log_file.write(message_str)
        self.log_file.flush()

        sys.stderr.write(message_str)
        sys.stderr.flush()

    def exception(self, exception):
        trace_string = traceback.format_exc()
        self.log("="*10 + " Exception " + "="*10)
        for line in trace_string.split("\n"):
            self.log(line)

        self.log(exception)

    def error(self, error_message):
        self.log("ERROR: {}".format(error_message))




def log_command(options, argv):
    main_dir = os.path.dirname(argv[0])

    try:
        with utilities.cd(main_dir):
            with open(os.devnull, "w") as devnull:
                git_hash = subprocess.check_output("git rev-parse HEAD", shell=True, stderr=devnull).strip()
                if subprocess.check_output("git diff", shell=True, stderr=devnull):
                    git_hash += "+"
    except subprocess.CalledProcessError:
        git_hash = ""

    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    debug_str = " DEBUG" if options.debug else ""

    run_log_file = open(os.path.join(options.output_dir, "runlog.txt"), "a")
    run_log_file.write("{} {}{}: {}\n".format(date_str, git_hash, debug_str, " ".join(argv)))
