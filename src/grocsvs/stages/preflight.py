import os
import pkg_resources
import re
import subprocess
import textwrap

from grocsvs import step
from grocsvs import utilities

def indent(text):
    lines = text.splitlines(True)
    return "".join(lines[:1] + ["   "+line for line in lines[1:]])

class PreflightStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield PreflightStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "preflight": os.path.join(directory, "preflight.txt")
        }

        return paths


    def run(self):
        outpath = self.outpaths(final=False)["preflight"]

        results = []
        success = True

        checkers = [self.check_memory,
                    self.check_tabix,
                    self.check_idba]

        for checker in checkers:
            try:
                cur_success, cur_result = checker()
                success &= cur_success
                results.append(indent(cur_result))
                if not cur_success:
                    self.logger.error(cur_result)
            except EOFError:
                success = False
                results.append("Failed to check '{}'".format(checker.__name__))
            
        with open(outpath, "w") as outf:
            for result in results:
                outf.write(result+"\n")

        if not success:
            raise Exception("Preflight failed; please check the log file at {}".format(outpath))


    def check_tabix(self):
        try:
            tabix = self.options.binary("tabix")
        except utilities.BinaryNotFoundError, e:
            return False, "ERROR: " + str(e)

        version = None
        tabix_results = subprocess.Popen(tabix, shell=True, stderr=subprocess.PIPE).stderr.readlines()
        for line in tabix_results:
            line = line.decode("utf-8")
            m = re.match("version: (\S+).*\n", line, flags=re.I)
            if m:
                version = pkg_resources.parse_version(m.group(1))

        if version is None:
            return False, "ERROR: Unable to ascertain tabix version; please make sure version 1.0 or greater is installed"
        elif version < pkg_resources.parse_version("1.0.0"):
            return False, "ERROR: Tabix version {} installed is too old (version 1.0 or greater required)".format(version)

        return True, "Tabix version {} found.".format(version)


    def check_idba(self):
        try:
            idba_ud = self.options.binary("idba_ud")
        except utilities.BinaryNotFoundError, e:
            return False, str(e)

        reads_file_path = os.path.join(self.working_dir, "reads.fa")
        with open(reads_file_path, "w") as reads_file:
            read1 = "ACGTGGCGAG"*25
            read2 = "GGACGACCCA"*25
            reads_file.write(">read1\n{}\n>read2\n{}\n".format(read1, read2))

        command = "{} -r {} -o {}".format(idba_ud, reads_file_path, self.working_dir)
        result = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

        stderr = result.stderr.read()
        if "sequence is too long" in stderr.lower():
            return False, "ERROR: idba_ud failed to assemble long reads; please make sure you use the idba_ud fork " + \
                          "provided at https://github.com/grocsvs/idba/releases/tag/1.1.3g1"

        outpath = os.path.join(self.working_dir, "contig.fa")
        if not os.path.exists(outpath):
            return False, "ERROR: idba_ud failed to produce final assembly; please make sure it is installed and runs\n" +\
                          "correctly; here is the output of the command:\n{}".format(stderr)
        return True, "idba_ud runs and is able to compile reads of length 250."

    def check_memory(self):
        import psutil
        physical_mem_gb = psutil.virtual_memory().total / (1000.**3)
        if physical_mem_gb < 8:
            return False, "ERROR: grocsvs requires at least 8GB of memory and preferably 16GB of memory in order to run\n" +\
                          "but I detect only {:.1f}GB present".format(physical_mem_gb)
        if physical_mem_gb < 16:
            return True, "WARNING: it is suggested to run grocsvs with at least 16GB of memory; I detect only {:.1f}GB present".format(
                physical_mem_gb)
        return True, "Detected {:.1f}GB available. This should be sufficient unless genome coverage is very high.".format(physical_mem_gb)