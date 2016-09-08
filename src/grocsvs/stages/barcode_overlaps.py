import itertools
# import gzip
import h5py
import logging
import numpy
import os
import scipy.stats


from grocsvs import datasets as svdatasets
from grocsvs import step
from grocsvs import utilities

from grocsvs.stages import window_barcodes



class BarcodeOverlapsStep(step.StepChunk):
    """
    Compare two sets of genomic windows, counting the number of overlapping
    barcodes between each pair of windows

    Output files:
        bcoverlaps.hist.sample.dataset.chromx.chromy.npy.gz - numpy array:
        - shape = n x m, where n=len(bcwindowsx) and m=len(bcwindowsy)
        - each value is the number of barcodes that are in both 
        bcwindowsx[i] and bcwindowsy[j]

        bcoverlaps.p.sample.dataset.chromx.chromy.npy.gz - numpy array:
        - similar to the hist above, except contains binomial p-values for 
        every pair of windows
    """

    @staticmethod
    def get_steps(options):
        chroms = options.reference.chroms
        for sample, dataset in options.iter_10xdatasets():
            for chromx, chromy in itertools.product(chroms, chroms):
                if options.reference.compare_chroms(chromx, chromy) < 0: continue

                # if chromx == chromy: continue


                for chunk in BarcodeOverlapsStep.chunks_for_chroms(
                        options, sample, dataset, chromx, chromy):
                    yield chunk

    @staticmethod
    def chunks_for_chroms(options, sample, dataset, chromx, chromy):
        nchunksx = window_barcodes.chunks_for_chrom(options, chromx)
        nchunksy = window_barcodes.chunks_for_chrom(options, chromy)

        for chunkx in range(0, nchunksx):
            for chunky in range(0, nchunksy):
                yield BarcodeOverlapsStep(
                    options, sample, dataset,
                    chromx, chromy, chunkx, chunky)


    def __init__(self, options, sample, dataset, chromx, chromy, chunkx, chunky):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chromx = chromx
        self.chromy = chromy
        self.chunkx = chunkx
        self.chunky = chunky

        assert isinstance(self.dataset, svdatasets.TenXDataset)

    def __str__(self):
        return ".".join([self.__class__.__name__,
            self.sample.name,
            self.dataset.id,
            self.chromx,
            self.chromy,
            str(self.chunkx),
            str(self.chunky)
            ])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "{}.{}.{}.{}.{}.{}.hdf5".format(
            self.sample.name,
            self.dataset.id,
            self.chromx,
            self.chromy,
            self.chunkx,
            self.chunky)

        paths = {
            "bcoverlaps": os.path.join(directory,
                                 "bcoverlaps.{}".format(file_name))
        }

        return paths


    def run(self):
        outpaths = self.outpaths(final=False)
        bcwindowsx, bcwindowsy, nbcs = self.get_bcwindows()

        hist, p = get_overlaps(bcwindowsx, bcwindowsy, nbcs, self.logger)

        # because marginally significant results are sparse,
        # this dramatically reduces space being used after compression
        p[p>0.001] = 1

        f = h5py.File(outpaths["bcoverlaps"], "w")
        f.create_dataset("hist", data=hist, compression="gzip")
        f.create_dataset("p", data=p, compression="gzip")


    def get_bcwindows(self):
        input_pathx = window_barcodes.WindowBarcodesStep(
            self.options,
            self.sample,
            self.dataset,
            self.chromx,
            self.chunkx).outpaths(final=True)

        input_pathy = window_barcodes.WindowBarcodesStep(
            self.options,
            self.sample, 
            self.dataset, 
            self.chromy,
            self.chunky).outpaths(final=True)

        bcwindowsx = utilities.pickle.load(open(
            input_pathx["bcwindows"]))
        bcwindowsy = utilities.pickle.load(open(
            input_pathy["bcwindows"]))

        nbcs = bcwindowsx["nbcs"]
        assert nbcs == bcwindowsy["nbcs"]

        return (bcwindowsx["barcode_windows"],
                bcwindowsy["barcode_windows"],
                nbcs)


def get_overlaps(bcwindowsx, bcwindowsy, nbcs, logger):
    width = len(bcwindowsx)
    height = len(bcwindowsy)

    hist = numpy.zeros((height, width))
    hist[:] = numpy.nan

    ps = hist.copy()
    ps[:,:] = numpy.nan
    
    nbcs = float(nbcs)
    
    for i in range(height):
        # if i % 100 == 0:
            # logger.log("{} of {} ({:.1f}%)".format(i, height, i/float(height)*100.0))

        y = bcwindowsy[i]
        prob = len(y) / nbcs
        
        counts = []
        n = []
        for j in range(width):
            x = bcwindowsx[j]
            curcounts = len(x.intersection(y))
            hist[i,j] = curcounts
            
            counts.append(curcounts)
            n.append(len(x))
        p = scipy.stats.binom.sf(counts, n, prob)
        ps[i, :] = p
    return hist, ps
