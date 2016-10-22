import collections
import numpy
import os
import pysam

from grocsvs import datasets as svdatasets
from grocsvs import step
from grocsvs import utilities
from grocsvs.stages import call_readclouds

# def get_sample_info_steps(options):
#     steps = []

#     for sample_name, sample in options.samples.items():
#         step = SampleInfoStep(
#             options, sample)
#         steps.append(step)

#     return steps


class SampleInfoStep(step.StepChunk):
    """

    Output files:
    """

    @staticmethod
    def get_steps(options):
        steps = []

        for sample_name, sample in options.samples.items():
            step = SampleInfoStep(
                options, sample)
            steps.append(step)

        return steps

    def __init__(self, options, sample):
        self.options = options
        self.sample = sample


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "sample_info.{}.pickle".format(
            self.sample.name)

        paths = {
            "sample_info": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        outpath = self.outpaths(final=False)["sample_info"]

        results_by_dataset = {}
        for dataset in self.sample.datasets:
            self.logger.log("Getting sample info for {}".format(dataset.id))
            results_by_dataset[dataset.id] = self.get_dataset_info(dataset)

        utilities.pickle.dump(results_by_dataset,
            open(outpath, "w"), protocol=-1)

    def get_dataset_info(self, dataset):
        if isinstance(dataset, svdatasets.TenXDataset):
            return self.get_10xdataset_info(dataset)

        if (isinstance(dataset, svdatasets.ShortFragDataset) or
            isinstance(dataset, svdatasets.MatePairDataset)):
            return self.get_short_frag_info(dataset)

    def get_10xdataset_info(self, dataset):
        self.logger.log(" loading fragments...")
        fragments = call_readclouds.load_fragments(
            self.options, self.sample, dataset)

        info = {
            "dataset_type": "10X",
            # "frag_length_quantiles": fragments["obs_len"].quantile(
                # numpy.arange(0, 1, 0.01)),
            "bam_coverage": get_bam_coverages(
                self.options, self.sample, dataset),
            "good_bc_count": len(fragments["bc"].unique()),
            "good_frags_count": fragments.shape[0],
            "frag_length_distributions": get_bg_fragments_distributions(
                fragments, self.options.reference, self.logger),
            "barcode_read_totals": get_barcode_read_totals(fragments)
        }

        info["frag_length_info"] = get_frag_lengths(fragments)
        info["physical_depths"] = [len(d) for d in info["frag_length_distributions"]]
        info["coverage_of_fragments"] = get_read_coverage_of_fragments(
            self.options, self.sample, dataset, fragments)

        # this assumes the median coverage is diploid, assumes most events are 
        # going to be haploid, and that not all fragments spanning the 
        # junctions will be recovered
        info["min_supporting_fragments"] = \
            numpy.median(info["physical_depths"]) * 0.20



        return info

    def get_short_frag_info(self, dataset):
        info = {}
        info["insert_sizes"], info["orientation"] = sample_insert_sizes(
            self.options, self.sample, dataset)
        return info

    
def get_bg_fragments_distributions(frags, reference, logger):
    distributions = []
    skip_length = 1e6
    for chrom in reference.chroms:
        logger.log("BG Frag Length Distribution: {}".format(chrom))
        chrom_length = reference.chrom_lengths[chrom]
        curfrags = frags.loc[frags["chrom"]==chrom]

        if chrom_length < 2*skip_length: continue
        for pos in numpy.arange(skip_length, chrom_length - skip_length, skip_length):
            e = utilities.frags_overlap_same_chrom(curfrags, pos, pos)
            if e.shape[0] < 5:
                continue
            distributions.append(e["obs_len"].values)

    return distributions

def n50(values, fraction=0.5):
    values = values.copy()
    values.sort()

    cumsum = numpy.cumsum(values)
    sum_ = cumsum[-1]
    
    i = numpy.where(cumsum>=fraction*sum_)[0][0]
    
    return values[i]

def get_frag_lengths(frags):
    lengths = frags.obs_len.values

    print lengths
    # lengths = numpy.concatenate(lengths)

    percentiles = numpy.arange(0,101,5)
    length_quantiles = dict(zip(percentiles/100.0, numpy.percentile(lengths, percentiles)))

    lengths_sampled = numpy.random.choice(lengths, 10000)

    results = {
        "quantiles": length_quantiles,
        "sampled": lengths_sampled,
        "mean": numpy.mean(lengths),
        "std": numpy.std(lengths),
        "N25": n50(lengths, 0.25),
        "N50": n50(lengths, 0.50),
        "N75": n50(lengths, 0.75)
    }

    return results


def get_search_regions(reference, window_skip, window_size, skip_at_ends):
    # get the chromosomes and move X, Y, M/MT to the end

    for chrom in reference.chroms:
        length = reference.chrom_lengths[chrom]
        if length >= (skip_at_ends*2 + window_size):
            for start in range(skip_at_ends, length-skip_at_ends, window_skip):
                yield (chrom, start, start+window_size)

def sample_insert_sizes(options, sample, dataset):
    inserts = []

    orientations = collections.Counter()

    to_sample = int(1e7)

    window_skip = int(1e6)
    window_size = int(1e5)
    skip_at_ends = int(1e6)

    min_mapq = 55

    max_insert_size = 20000

    bam = pysam.AlignmentFile(dataset.bam)

    for chrom, start, end in get_search_regions(options.reference, 
                                                window_skip, window_size, skip_at_ends):
        if chrom not in bam.references:
            chrom = chrom.replace("chr", "")
            if chrom not in bam.references:
                raise Exception("Unknown chrom: {}; make sure all your samples (10X + frag) all use the same reference")
            
        for read in bam.fetch(chrom, start, end):
            if not read.is_paired:
                raise Exception("Unpaired read found!")
                # orientations["unpaired"] += 1
                # readLengths.append(len(read.seq))
                # continue
                
            if not read.is_read1:
                continue
            
            if not read.is_proper_pair:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.tid != read.rnext:
                continue
            if read.mapq < min_mapq:
                continue
            
            if abs(read.isize) < max_insert_size:
                inserts.append(abs(read.isize))

            cur_orientation = (read.is_reverse, read.mate_is_reverse)
            if read.reference_start > read.next_reference_start:
                cur_orientation = not cur_orientation[0], not cur_orientation[1]
            orientations[cur_orientation] += 1

            if to_sample % 100000 == 0:
                print to_sample
            to_sample -= 1
            if to_sample <= 0:
                break
        if to_sample <= 0:
            break

    orientation = max(orientations, key=lambda x: orientations[x])
    switch = {False: "+", True:"-"}
    orientation = switch[orientation[0]] + switch[orientation[1]]

    print len(inserts), max(inserts), inserts.count(max(inserts))
    inserts = removeOutliers(inserts)
    print len(inserts), max(inserts), (inserts==max(inserts)).sum()
    inserts = collections.Counter(inserts)

    return inserts, orientation


def get_bam_coverages(options, sample, dataset):
    window_skip = int(1e6)
    window_size = 1000
    skip_at_ends = int(1e6)

    bam = pysam.AlignmentFile(dataset.bam)

    print "getting coverages"
    coverages = []
    count = 0
    for chrom, start, end in get_search_regions(options.reference, 
                                                window_skip, window_size, skip_at_ends):
        coverages.append(bam.count(chrom, start, end))

        if count % 1000 == 0:
            print count
        count += 1
    return coverages


def get_read_coverage_of_fragments(options, sample, dataset, fragments):
    bam = pysam.AlignmentFile(dataset.bam)

    coverages = []
    lengths = []

    sampled_frags = fragments.sample(1000)

    print "sampling long fragments, getting coverage by short reads..."
    for i, (_, frag) in enumerate(sampled_frags.iterrows()):
        if i % 10 == 0:
            print i
        positions = set()
        
        for read in bam.fetch(frag["chrom"], frag["start_pos"], frag["end_pos"]):
            if read.has_tag("BX") and read.get_tag("BX") == frag["bc"]:
                positions.update(read.get_reference_positions())
        
        positions = [x for x in positions if frag["start_pos"]<=x<=frag["end_pos"]]
        coverages.append(len(positions))
        lengths.append(frag["end_pos"]-frag["start_pos"])

    result = {}
    result["coverages"] = numpy.array(coverages)
    result["lengths"] = numpy.array(lengths)
    result["mean_cr"] = numpy.mean(result["coverages"] / result["lengths"].astype(float))

    return result




def removeOutliers(data, m = 10.):
    """ a method of trimming outliers from a list/array using 
    outlier-safe methods of calculating the center and variance;
    only removes the upper tail, not the lower tail """
    if len(data) < 2:
        return data
        
    data = numpy.array(data)
    d_abs = numpy.abs(data - numpy.median(data))
    d = data - numpy.median(data)
    mdev = numpy.median(d_abs)
    s = d/mdev if mdev else 0.
    return data[s<m]

def get_barcode_read_totals(frags):
    return frags.groupby("bc")["num_reads"].sum()
