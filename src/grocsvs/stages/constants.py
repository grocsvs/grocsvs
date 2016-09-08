import json
import os

from grocsvs import step


class ConstantsStep(step.StepChunk):
    """
    store some constants for use later in the pipeline; by storing them in a file,
    we can edit them and run the rest of the pipeline, without mucking around 
    with the stage source files
    """

    @staticmethod
    def get_steps(options):
        steps = []

        step = ConstantsStep(options)
        steps.append(step)

        return steps
        

    def __init__(self, options):
        self.options = options


    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "constants": os.path.join(directory, "constants.json")
        }

        return paths


    def run(self):
        outpath = self.outpaths(final=False)["constants"]

        constants = {
            "readcloud_clustering_dist": 10000, # TODO: this should be higher for gemcode
            "window_size": 10000,
            "min_supporting_fragment_fraction": 0.10,
            "min_supporting_fragments": 40
        }

        json.dump(constants, open(outpath, "w"))


def get_constants(options):
    constants_path = os.path.join(options.results_dir,
                                  "ConstantsStep",
                                  "constants.json")

    if not os.path.exists(constants_path):
        raise Exception("GetConstantsStep needs to be run first! '{}'".format(
            constants_path))

    return json.load(open(constants_path))
