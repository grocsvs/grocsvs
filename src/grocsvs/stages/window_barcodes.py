import numpy
import os

from grocsvs import datasets as svdatasets
from grocsvs import step
from grocsvs import utilities
from grocsvs.stages import call_readclouds

CHUNKSIZE = 5e7

def chunks_for_chrom(options, chrom):
    return int(numpy.ceil(options.reference.chrom_lengths[chrom]/CHUNKSIZE))


class WindowBarcodesStep(step.StepChunk):
    """
    Build a list of all the fragment barcodes overlapping each
    genomic window

    Output files:
        bcwindows.sample.dataset.chrom.pickle - dictionary:
        - barcode_windows - a list of sets of barcode IDs
        - barcod_map - a dict of barcode->barcode ID
        - window_size - the size of the genomic window used; eg 10,000 would 
            mean that window starts were range(0, chrom_length, 10000)
    """

    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            for chrom in options.reference.chroms:
                for chunk in range(chunks_for_chrom(options, chrom)):
                    yield WindowBarcodesStep(
                        options, sample, dataset, chrom, chunk)

    def __init__(self, options, sample, dataset, chrom, chunk):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chrom = chrom
        self.chunk = chunk

        assert isinstance(self.dataset, svdatasets.TenXDataset)


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id, 
                         self.chrom,
                         str(self.chunk)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "bcwindows.{}.{}.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id,
            self.chrom,
            self.chunk)

        paths = {
            "bcwindows": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        import logging
        logging.info("running!")

        window_size = self.options.constants["window_size"]
        outpath = self.outpaths(final=False)["bcwindows"]


        self.logger.log("Loading barcode map...")
        # call_readclouds_step = call_readclouds.FilterFragmentsStep(
        input_step = call_readclouds.CombineReadcloudsStep(
            self.options, self.sample, self.dataset)
        barcode_map = utilities.pickle.load(
            open(input_step.outpaths(final=True)["barcode_map"]))

        chrom_length = self.options.reference.chrom_lengths[self.chrom]
        start = int(self.chunk*CHUNKSIZE)
        end = int(min((self.chunk+1)*CHUNKSIZE, chrom_length))


        self.logger.log("Running chunk: {}:{:,}-{:,}".format(self.chrom, start, end))

        fragments = call_readclouds.load_fragments(
            self.options, self.sample, self.dataset,
            self.chrom, start, end, min_reads_per_frag=0)

        barcode_windows = get_barcode_windows(
            fragments, barcode_map, window_size, chrom_length, start, end)


        self.logger.log("Saving results...")
        result = {
            "barcode_windows": barcode_windows,
            # "barcode_map":     barcode_map,
            "nbcs": len(barcode_map),
            "window_size":     window_size
        }
        utilities.pickle.dump(result, open(outpath, "w"), protocol=-1)



def get_barcode_windows(fragments, barcode_map, window_size, chrom_length, start, end):
    window_starts = range(start, end, window_size)

    barcode_windows = []

    for start in window_starts:
        end = start + window_size
        overlap = utilities.frags_overlap_same_chrom(fragments, start, end)
        barcodes = set(barcode_map[bc] for bc in overlap["bc"])
        barcode_windows.append(barcodes)

    return barcode_windows
