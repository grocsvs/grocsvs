# import gzip
import h5py
import itertools
import numpy
import os
import pandas

from grocsvs import step
from grocsvs import utilities
from grocsvs import structuralvariants

from grocsvs.stages import barcode_overlaps



class SVCandidateRegionsStep(step.StepChunk):
    """
    Takes the matrix of barcode overlaps and p-values and finds clusters
    of significant window-pairs

    Output files:
        svs.sample.dataset.chromx.chromy.pickle - list
        - each item is the x,y coordinates of a significant window-pair
    """

    @staticmethod
    def get_steps(options):
        chroms = options.reference.chroms

        for sample, dataset in options.iter_10xdatasets():
            for chromx, chromy in itertools.product(chroms, chroms):
                if options.reference.compare_chroms(chromx, chromy) >= 0:
                    yield SVCandidateRegionsStep(
                            options, sample, dataset, chromx, chromy)

    def __init__(self, options, sample, dataset, chromx, chromy):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chromx = chromx
        self.chromy = chromy


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         self.chromx,
                         self.chromy])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "sv_regions.{}.{}.{}.{}.tsv".format(
            self.sample.name,
            self.dataset.id,
            self.chromx,
            self.chromy
            )

        paths = {
            "sv_regions": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        outpath = self.outpaths(final=False)["sv_regions"]

        self.logger.log("loading...")
        hist, p = self.get_barcode_overlaps()

        self.logger.log("getting SV regions...")
        # TODO: learn these constants
        window_size = self.options.constants["window_size"]
        offdiag_dist = int(self.max_distance_between_events()/window_size)

        self.logger.log("Using offdiagonal distance of {} (with window size {})".format(
            offdiag_dist, window_size))
        sv_regions = get_sv_regions(
            hist, 
            p, 
            self.chromx, 
            self.chromy, 
            offdiag_dist=10, 
            clustering_dist=1,
            overlap_min=25,
            window_size=self.options.constants["window_size"],
            logger=self.logger)

        if len(sv_regions) > 0:
            # utilities.pickle.dump(sv_regions, open(outpath, "w"), protocol=-1)
            sv_regions.to_csv(outpath, sep="\t", index=False)
        else:
            # touch
            open(outpath, "w")

    def max_distance_between_events(self, pct=80):
        sample_info = self.options.sample_info(self.sample.name)
        dataset_info = sample_info[self.dataset.id]

        lengths = numpy.concatenate(
            dataset_info["frag_length_distributions"])

        return int(numpy.percentile(lengths, pct))

    def get_barcode_overlaps(self):
        window_size = self.options.constants["window_size"]

        chromx_length = self.options.reference.chrom_lengths[self.chromx]
        chromy_length = self.options.reference.chrom_lengths[self.chromy]

        hist = numpy.zeros((chromy_length/window_size+1,
                            chromx_length/window_size+1))

        p = numpy.zeros((chromy_length/window_size+1,
                         chromx_length/window_size+1))

        hist[:] = numpy.nan
        p[:] = numpy.nan

        for chunk in barcode_overlaps.BarcodeOverlapsStep.chunks_for_chroms(
            self.options, self.sample, self.dataset, self.chromx, self.chromy):

            inpaths = chunk.outpaths(final=True)

            cur_file = h5py.File(inpaths["bcoverlaps"])
            cur_hist = numpy.array(cur_file["hist"])
            cur_p =    numpy.array(cur_file["p"])

            hist[chunk.chunky*5000:(chunk.chunky+1)*5000,
                 chunk.chunkx*5000:(chunk.chunkx+1)*5000] = cur_hist

            p[chunk.chunky*5000:(chunk.chunky+1)*5000,
              chunk.chunkx*5000:(chunk.chunkx+1)*5000] = cur_p

        print "COUNTs:", numpy.isnan(hist).sum(), (~numpy.isnan(hist)).sum()
        return hist, p


# TODO: refactor the below code
# it should:
# - use cutoffs based on the data for the *_dist and overlap_min params
# - figure out if we need to adjust the p-value cutoff?
# - 


def get_sv_regions(hist, p, chrom1, chrom2, offdiag_dist, clustering_dist,
           overlap_min, window_size, logger):
    mask = (hist > overlap_min) & (p < 1e-7)
    coords = numpy.where(mask)

    print hist.shape, mask.sum()
    svs_raw = structuralvariants.offdiag(coords, offdiag_dist=offdiag_dist)
    if chrom1 == chrom2:
        svs_raw = [coord for coord in svs_raw if coord[0]>coord[1]]

    logger.log("raw: {}".format(len(svs_raw)))

    if not len(svs_raw):
        return []

    svs_clustered = structuralvariants.do_grid_clustering(
        svs_raw, clustering_dist)

    logger.log("clustered: {}".format(len(svs_clustered)))
    sv_regions = []
    for cluster in svs_clustered:
        ys, xs = zip(*cluster)
        startx, endx = min(xs), max(xs)+1
        starty, endy = min(ys), max(ys)+1
        # print "^^^", startx, endx, starty, endy
        # sv_regions.append({"startx":startx*window_size,
        #                    "endx":endx*window_size,
        #                    "starty":starty*window_size,
        #                    "endy":endy*window_size,
        #                    "chromx": chrom1,
        #                    "chromy": chrom2})

        region_positions = itertools.product(range(starty, endy), range(startx, endx))
        region_positions = zip(*structuralvariants.offdiag(zip(*region_positions), offdiag_dist))

        # print "*"*1000, p[region_positions]

        best_position = numpy.where(p[region_positions] == p[region_positions].min())[0][0]

        bestp =     p[region_positions][best_position]
        bestcount = hist[region_positions][best_position]
        # print best_position, bestp, bestcount
        sv_regions.append([chrom1,
                           chrom2,
                           startx*window_size,
                           endx*window_size,
                           starty*window_size,
                           endy*window_size,
                           bestp,
                           bestcount])

    sv_regions = pandas.DataFrame(sv_regions, 
        columns=["chromx", "chromy", "startx", "endx", "starty", "endy", "p", "count"])
    logger.log("final regions: {}".format(len(sv_regions)))

    # from biorpy import r

    # r.png("clusters.{}.{}.png".format(chrom1, chrom2), width=5000, height=5000)

    # r.plot([0], xlim=[0,hist.shape[1]], ylim=[0,hist.shape[0]])

    # colors = ["red","blue","orange","green", "gray", "black", "purple"]
    # for i, cluster in enumerate(svs_clustered):
    #     if i % 100 == 0:
    #         print i
    #     xs, ys = zip(*cluster)
    #     startx, endx = min(xs), max(xs)+1
    #     starty, endy = min(ys), max(ys)+1

    #     r.rect(startx, starty, endx, endy, col=colors[i%len(colors)])

    #     # r.points(*zip(*cluster))
    # r.devoff()


    return sv_regions

