import os
import pandas

from grocsvs import datasets as svdatasets
from grocsvs import step
from grocsvs import utilities
from grocsvs import structuralvariants
#from grocsvs.stages import filter_fragments
from grocsvs.stages import call_readclouds
from grocsvs.stages import postassembly_merge


class SupportingBarcodesStep(step.StepChunk):
    """
    """
    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            step = SupportingBarcodesStep(options, sample, dataset)
            yield step


    def __init__(self, options, sample, dataset):
        self.options = options
        self.sample = sample
        self.dataset = dataset

        assert isinstance(self.dataset, svdatasets.TenXDataset)


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        event_barcodes = "event_barcodes.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id)

        events = "events.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id)

        paths = {
            "event_barcodes": os.path.join(directory, event_barcodes),
            "events": os.path.join(directory, events)
        }

        return paths


    def run(self):
        self.logger.log("loading data...")
        #filter_fragments_step = filter_fragments.FilterFragmentsStep(
        #    self.options, self.sample, self.dataset)
        #barcode_map = utilities.pickle.load(
        #    open(filter_fragments_step.outpaths(final=True)["barcode_map"]))

        call_readclouds_step = call_readclouds.CombineReadcloudsStep(
            self.options, self.sample, self.dataset)
        barcode_map = utilities.pickle.load(
            open(call_readclouds_step.outpaths(final=True)["barcode_map"]))
        
        events = self.load_events()
        events_to_barcodes = {}

        self.logger.log("collecting barcodes for {} events...".format(len(events)))        
        for event in events[["chrom", "pos", "orientation"]].itertuples(index=False):
            cur_barcodes = self.get_event_barcodes(event)
            events_to_barcodes[tuple(event)] = set(barcode_map[bc] for bc in cur_barcodes)

        outpaths = self.outpaths(final=False)
        with open(outpaths["event_barcodes"], "w") as f:
            utilities.pickle.dump(events_to_barcodes, f, protocol=-1)

        # we need to save this separately so that we have them ordered properly
        with open(outpaths["events"], "w") as f:
            utilities.pickle.dump(list(events_to_barcodes.keys()), f, protocol=-1)

    def load_events(self):
        path = postassembly_merge.PostAssemblyMergeStep(self.options) \
                                 .outpaths(final=True)["merged_candidates"]

        events = pandas.read_table(path)
        return events
        
    def get_event_barcodes(self, event):
        # TODO: gah dist; should be lower, though, than for discovery
        dist1 = -1000
        dist2 = 5000

        # chrom, position, orientation = event
        frags = structuralvariants.get_fragments_one_side(
            self.options, self.sample, self.dataset,
            event.chrom, event.pos, event.orientation,
            dist1, dist2, strict=True, min_reads_per_frag=0)

        return set(frags["bc"])

