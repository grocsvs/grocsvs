import collections
import gzip
import numpy
import os
import pandas
import pysam
import cStringIO as StringIO
import subprocess

from grocsvs import step
from grocsvs import utilities

Readcloud = collections.namedtuple("Readcloud", "chrom start_pos end_pos bc num_reads obs_len hap")


class CombineReadcloudsStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            yield CombineReadcloudsStep(options, sample, dataset)

    def __init__(self, options, sample, dataset):
        self.options = options
        self.sample = sample
        self.dataset = dataset

    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        readclouds = "readclouds.{}.{}.tsv.gz".format(
            self.sample.name,
            self.dataset.id
            )

        barcode_map_file = "barcode_map.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id
            )

        paths = {
            "barcode_map": os.path.join(directory, barcode_map_file),
            "readclouds": os.path.join(directory, readclouds),
            "index": os.path.join(directory, readclouds+".tbi")

        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)

        self.logger.log("Loading read clouds...")

        readclouds = []
        for i, inpath in enumerate(self.get_input_paths()):
            try:
                readclouds.append(pandas.read_table(inpath, compression="gzip"))
            except pandas.io.common.EmptyDataError:
                self.logger.log("No read clouds found in {}; skipping".format(inpath))
        readclouds = pandas.concat(readclouds)
        goodbcs = get_good_barcodes(readclouds)
        readclouds = readclouds.loc[readclouds["bc"].isin(goodbcs)]

        good_barcodes = readclouds["bc"].unique()
        barcode_map = get_barcode_map(good_barcodes)


        self.logger.log("Writing barcode map to file...")
        with open(outpaths["barcode_map"], "w") as outf:
            utilities.pickle.dump(barcode_map, outf, protocol=-1)


        self.logger.log("Writing readclouds to file...")

        tmp_readclouds_path = outpaths["readclouds"][:-3]
        readclouds.to_csv(tmp_readclouds_path, sep="\t", index=False)

        bgzip_cmd = "bgzip {}".format(tmp_readclouds_path)
        bgzip_proc = subprocess.Popen(bgzip_cmd, shell=True)
        bgzip_proc.wait()


        self.logger.log("Indexing readclouds file...")
        # define the chrom, start and end columns; and indicate that
        # the first (header) row should be skipped
        tabix_cmd = "tabix -s 1 -b 2 -e 3 -S 1 {}".format(outpaths["readclouds"])
        subprocess.check_call(tabix_cmd, shell=True)


    def get_input_paths(self):
        paths = []
        for chrom in self.options.reference.chroms:
            input_step = CallReadcloudsStep(self.options, self.sample, self.dataset, chrom)
            paths.append(input_step.outpaths(final=True)["readclouds"])

        return paths


class EstimateReadCloudParamsStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            yield EstimateReadCloudParamsStep(options, sample, dataset)

    def __init__(self, options, sample, dataset):
        self.options = options
        self.sample = sample
        self.dataset = dataset

    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        inter_read_distances = "inter_read_distances.{}.{}.pickle".format(
            self.sample.name,
            self.dataset.id
            )

        paths = {
            "inter_read_distances": os.path.join(directory, inter_read_distances),
        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)

        inter_read_distances = sample_inter_read_distances(self.dataset.bam)
        result = {}
        result["sampled_inter_read_distances"] = numpy.random.choice(
            inter_read_distances, int(0.5e6), replace=False)
        result["read_cloud_clustering_distance"] = \
            max(10000, int(numpy.ceil(numpy.percentile(inter_read_distances, 99) / 5000)) * 5000)

        self.logger.log("{} {} {}".format(*numpy.percentile(inter_read_distances, [50, 99, 99.9])))
        self.logger.log(result["read_cloud_clustering_distance"])
                        
        with open(outpaths["inter_read_distances"], "w") as outf:
            utilities.pickle.dump(result, outf, protocol=-1)
            
def sample_inter_read_distances(bam_path, window_size=0.5e6, skip_size=5e7):
    bam = pysam.AlignmentFile(bam_path)

    window_size = int(window_size)
    skip_size = int(skip_size)
    if skip_size < window_size:
        skip_size = window_size
        
    distances = []
    
    for chrom, length in zip(bam.references, bam.lengths):
        if length < 2*skip_size+2*window_size: continue
            
        print chrom
        for start in range(skip_size, length-skip_size, skip_size):
            bc_last_pos = {}

            for read in bam.fetch(chrom, start, start+window_size):
                if read.is_secondary or read.is_supplementary or read.is_unmapped or read.is_read2:
                    continue
                if not read.has_tag("BX"):
                    continue
                bc = read.get_tag("BX")
                if bc in bc_last_pos:
                    d = read.pos - bc_last_pos[bc]
                    if d < 60000:
                        distances.append(d)
                        if len(distances) > 10e6: # that's a plenty big sample!
                            return distances
                bc_last_pos[bc] = read.pos

    if len(distances) < 25 and skip_size > window_size:
        new_skip_size = skip_size / 100
        return sample_inter_read_distances(bam_path, window_size, new_skip_size)
    
    return distances



class CallReadcloudsStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        for sample, dataset in options.iter_10xdatasets():
            for chrom in options.reference.chroms:
                yield CallReadcloudsStep(
                    options, sample, dataset, chrom)

        
    def __init__(self, options, sample, dataset, chrom):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chrom = chrom


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         self.chrom])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        readclouds = "readclouds.{}.{}.{}.tsv.gz".format(
            self.sample.name,
            self.dataset.id,
            self.chrom
            )

        paths = {
            "readclouds": os.path.join(directory, readclouds)
        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)

        input_step = EstimateReadCloudParamsStep(self.options, self.sample, self.dataset)
        input_path = input_step.outpaths(final=True)["inter_read_distances"]
        info = utilities.pickle.load(open(input_path))
        
        max_dist = info["read_cloud_clustering_distance"]
        
        self.logger.log("Using {} for distance between readclouds".format(max_dist))

        bam_path = self.dataset.bam
        readclouds = call_readclouds(bam_path, self.chrom, max_dist)

        readclouds.to_csv(outpaths["readclouds"], sep="\t", index=False, compression="gzip")


def call_readclouds(bam_path, chrom, max_dist):
    """
    we'll load everything into memory so we can easily sort; this isn't strictly
    necessary if we were to start running out of memory
    """

    detector = ReadcloudDetector(bam_path, max_dist)

    readcloud_iter = detector.get_read_clouds(chrom)
    dataframe = pandas.DataFrame(readcloud_iter)

    if len(dataframe) > 1:
        dataframe = dataframe.sort_values(["start_pos", "end_pos"])

    return dataframe




def is_good_read(read):
    if not read.has_tag("BX"):
        return False
    # if read.mapq < 30:
        # return False
    if read.is_duplicate:
        return False
    if read.is_unmapped:
        return False
    return True


class ReadcloudDetector(object):
    def __init__(self, bam_path, max_dist, min_end_mapq=30, require_one_mapq=60):
        self.bam = pysam.AlignmentFile(bam_path)
        self.max_dist = max_dist
        self.min_end_mapq = min_end_mapq
        self.require_one_mapq = require_one_mapq

        self._cur_chrom = None
        self._barcodes_to_reads = {}
        self._barcodes_to_haplotypes = {}
        self._recently_observed_bcs = set()

        self.cur_start = 0

        
    def get_read_clouds(self, chrom):
        assert self._cur_chrom is None, "can only run on one chrom at once"
        self._cur_chrom = chrom

        self.cur_start = 0
        _progress = int(1e6)
        for i, read in enumerate(self.bam.fetch(chrom)):
            if i % _progress == 0:
                print "{:>15} {:15,.0f}".format(i, read.pos)

            # if i > 2e6:
            #     break

            if (read.pos - self.cur_start) > self.max_dist:
                self.cur_start = read.pos
                for read_cloud in self._flush():
                    yield read_cloud
            read_cloud = self._add_read(read)
            if read_cloud is not None:
                yield read_cloud

        for read_cloud in self._flush():
            yield read_cloud

        self._cur_chrom = None



    def _make_read_cloud(self, bc):
        frag_reads = self._barcodes_to_reads.pop(bc)

        if max(read.mapq for read in frag_reads) < self.require_one_mapq:
            return None

        while frag_reads[-1].mapq < self.min_end_mapq:
            frag_reads.pop(-1)

        count = len(frag_reads)

        # This is how 10X counts; not really sure why, since by this approach,
        # we could have a 100kb read cloud with only a single read counted
        # sum(1 for read in frag_reads if read.mapq >= self.require_one_mapq)
        # if count == 0:
        #     return None

        start = min(read.pos for read in frag_reads)
        end = max(read.aend for read in frag_reads)
        obs_len = end - start
        haplotype = self._barcodes_to_haplotypes[bc]
        if haplotype is not None:
            haplotype = haplotype[1]

        read_cloud = Readcloud(self._cur_chrom, start, end, bc, count, obs_len, haplotype)
        return read_cloud

    def _flush(self):
        to_flush = set(self._barcodes_to_reads) - self._recently_observed_bcs

        for bc in to_flush:
            read_cloud = self._make_read_cloud(bc)

            # if bc == "AGTGAAAAGCTTGGTT-1": print read_cloud
            
            if read_cloud:
                yield read_cloud

        self._recently_observed_bcs = set()

    def _add_read(self, read):
        if not is_good_read(read):
            return None

        bc = read.get_tag("BX")

        # if bc == "AGTGAAAAGCTTGGTT-1":
        #     print read.pos, read.mapq, len(self._barcodes_to_reads.get(bc, []))
            
        if bc in self._barcodes_to_reads:
            previous_read = self._barcodes_to_reads[bc][-1]
            if (read.pos - previous_read.aend) > self.max_dist:
                read_cloud = self._make_read_cloud(bc)
                self._barcodes_to_reads[bc] = [read]
                self._barcodes_to_haplotypes[bc] = get_read_haplotype(read)
                self._recently_observed_bcs.add(bc)
                return read_cloud
            else:
                self._barcodes_to_reads[bc].append(read)
                self._recently_observed_bcs.add(bc)
        elif read.mapq > self.min_end_mapq:
            self._barcodes_to_reads[bc] = [read]
            self._recently_observed_bcs.add(bc)
            self._barcodes_to_haplotypes[bc] = get_read_haplotype(read)
                
        return None

def get_read_haplotype(read):
    try:
        hp = read.get_tag("HP")
        ps = read.get_tag("PS")
        return (ps, hp)
    except KeyError:
        return None


def get_good_barcodes(fragments, proportion=0.90):
    """
    return the top barcodes which together comprise 90% of reads
    """
    read_counts = fragments.groupby("bc").sum()["num_reads"].copy()

    read_counts.sort_values(inplace=True, ascending=False)
    cutoff = proportion * read_counts.sum()
    cutoff = numpy.where(read_counts.cumsum() >= cutoff)[0][0]

    return sorted(read_counts.index[:cutoff])



def get_barcode_map(barcodes):
    barcodes = sorted(barcodes)
    return dict(zip(barcodes, range(len(barcodes))))


def load_fragments(options, sample, dataset, chrom=None, start=None, end=None, usecols=None, 
                   min_reads_per_frag=1):
    if start is not None:
        if start < 0:
            raise Exception("start coord is negative: {}:{}-{}".format(chrom, start, end))
    if end is not None:
        if start >= end:
            raise Exception("end coord is before start: {}:{}-{}".format(chrom, start, end))

    readclouds_path = os.path.join(
        options.results_dir,
        "CombineReadcloudsStep",
        "readclouds.{}.{}.tsv.gz".format(sample.name, dataset.id))

    tabix = pysam.TabixFile(readclouds_path)
    
    if chrom is not None and chrom not in tabix.contigs:
        print("MISSING:", chrom)
        return pandas.DataFrame(columns="chrom start_pos end_pos bc num_reads obs_len hap".split())
    
    if usecols is not None and "num_reads" not in usecols:
        usecols.append("num_reads")
        
    s = StringIO.StringIO("\n".join(tabix.fetch(chrom, start, end)))
    readclouds = pandas.read_table(s, header=None, names=Readcloud._fields, usecols=usecols)
    readclouds["chrom"] = readclouds["chrom"].astype("string")
    
    if min_reads_per_frag > 0:
        readclouds = readclouds.loc[readclouds["num_reads"]>min_reads_per_frag]

    return readclouds

    
