import gzip
import numpy
import os
import pandas
import random

from grocsvs import step
from grocsvs import utilities
from grocsvs import structuralvariants

from grocsvs.stages import cluster_svs

MAX_BARCODES = 200


class BarcodesFromGraphsStep(step.StepChunk):
    """
    Collect barcodes supporting each candidate structural event so that
    we can get all the reads that may support an event and perform assembly
    """

    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            step = BarcodesFromGraphsStep(
                options, sample, dataset)
            yield step


    def __init__(self, options, sample, dataset):
        self.options = options
        self.sample = sample
        self.dataset = dataset


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         str(self.dataset.id)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "sv_barcodes.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id
            )

        paths = {
            "sv_barcodes": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        dist1 = -500
        dist2 = 5000
        outpath = self.outpaths(final=False)["sv_barcodes"]

        self.logger.log("loading...")
        events = self.load_events()

        barcodes_map = {}

        for i, cluster in events.groupby("cluster"):
            barcodes = set()

            for j, event in cluster.iterrows():
                _, _, merged_frags = \
                    structuralvariants.get_supporting_fragments_new(
                        self.options, self.sample, self.dataset,
                        event["chromx"], int(event["x"]),
                        event["chromy"], int(event["y"]),
                        event["orientation"], dist1, dist2,
                        min_reads_per_frag=0)

                cur_bcs = set(merged_frags["bc"])
                if len(cur_bcs) > MAX_BARCODES:
                    print "TOO MANY BARCODES: sampling {} of {} for cluster {}".format(MAX_BARCODES, len(cur_bcs), i)
                    cur_bcs = random.sample(cur_bcs, MAX_BARCODES)
                    
                if len(cur_bcs) == 0:
                    self.logger.log("WARNING: no barcodes found for event {} {}:{}::{}:{}{} dist1={} dist2={}".format(
                        self.sample.name, event.chromx, event.x, event.chromy, event.y, event.orientation, dist1, dist2))
                barcodes.update(cur_bcs)

            barcodes_map[i] = barcodes

        utilities.pickle.dump(barcodes_map, open(outpath, "w"), protocol=-1)

    def load_events(self):
        # path = os.path.join(self.results_dir, "edges.tsv")
        edges_path = cluster_svs.ClusterSVsStep(self.options).outpaths(final=True)["edges"]

        graphs_table = pandas.read_table(edges_path)

        return graphs_table
