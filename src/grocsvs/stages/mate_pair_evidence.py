import numpy
import os
import pandas
import pysam

from grocsvs import datasets as svdatasets
from grocsvs import step

from grocsvs.stages import pair_evidence



class MatePairEvidenceStep(step.StepChunk):
    """
    """
    @staticmethod
    def get_steps(options):
        for sample_name, sample in options.samples.items():
            tenx_dataset = sample.get_10x_dataset()
            for short_frag_dataset in sample.datasets:
                if isinstance(short_frag_dataset, svdatasets.MatePairDataset):
                    evidence_steps = pair_evidence.PairEvidenceStep.get_steps_for_dataset(
                        options, sample, tenx_dataset)

                    for evidence_step in evidence_steps:
                        yield MatePairEvidenceStep(
                            options, sample, short_frag_dataset, evidence_step)


    def __init__(self, options, sample, dataset, evidence_step):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.evidence_step = evidence_step

        assert isinstance(self.dataset, svdatasets.MatePairDataset)


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         str(self.evidence_step.chunk)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "mate_pair_evidence.{}.{}.{}.tsv".format(
            self.sample.name,
            self.dataset.id,
            self.evidence_step.chunk)

        paths = {
            "mate_pair_evidence": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        self.logger.log("running!")

        path = self.evidence_step.outpaths(final=True)["evidence"]
        evidence = pandas.read_table(path)
        evidence["chromx"] = evidence["chromx"].astype("string")
        evidence["chromy"] = evidence["chromy"].astype("string")
        
        sample_info = self.options.sample_info(self.sample.name)
        dataset_info = sample_info[self.dataset.id]

        insert_size_dict = dataset_info["insert_sizes"]
        bam = pysam.AlignmentFile(self.dataset.bam)
        short_frag_evidence = []

        print "@@@@", get_max_insert_size(insert_size_dict)
        # for i, j in numpy.transpose(numpy.where(binom_ps<1e-4)):
        for i, row in evidence.iterrows():
            # if common_counts[i,j] < 20: continue
            if row["shared"] < 20: continue
            chromx, x, orientationx = row["chromx"], int(row["new_x"]), row["orientation"][0]
            chromy, y, orientationy = row["chromy"], int(row["new_y"]), row["orientation"][1]

            cur_evidence = get_evidence(chromx, x, chromy, y, orientationx, orientationy, 
                bam, get_max_insert_size(insert_size_dict), 5000)

            short_frag_evidence.append(
                [chromx, row["original_x"], orientationx,
                 chromy, row["original_y"], orientationy]
                 + list(cur_evidence))

        columns = ["chromx", "x", "orientationx", "chromy", "y", "orientationy",
                   "common", "total"]
        df = pandas.DataFrame(short_frag_evidence, columns=columns)

        outpath = self.outpaths(final=False)["mate_pair_evidence"]
        df.to_csv(outpath, sep="\t", index=False)

def get_max_insert_size(insert_size_dict):
    isd = numpy.array([insert_size_dict.get(x, 0) for x in range(max(insert_size_dict)+1)])

    cum_isize_counts = numpy.cumsum(isd)
    max_ins = numpy.where(cum_isize_counts >= 0.9*cum_isize_counts[-1])[0][0]

    return max_ins


def _coords(position, orientation, dist1, dist2):
    if orientation == "+":
        start = position-dist2
        end = position-dist1
    elif orientation == "-":
        start = position+dist1
        end = position+dist2

    start = max(0, start)
    end = max(0, end)

    return start, end

def get_evidence(chromx, x, chromy, y, orientationx, orientationy, bam, max_ins, imprecision, backwards=False, min_mapq=55):
    extend = max_ins * 2 + imprecision
    startx, endx = _coords(x, orientationx, -imprecision, extend)
    starty, endy = _coords(y, orientationy, -imprecision, extend)
    
    flip = {"+":"-", "-":"+"}
    if backwards:
        orientationx = flip[orientationx]
        orientationy = flip[orientationy]
        
    return _get_evidence(chromx, startx, endx, chromy, starty, endy, orientationx, orientationy, bam, min_mapq)

def _get_evidence(chromx, startx, endx, chromy, starty, endy, orientationx, orientationy, bam, min_mapq=55):
    readsx = get_reads(bam, chromx, startx, endx, orientationx, min_mapq)
    readsy = get_reads(bam, chromy, starty, endy, orientationy, min_mapq)

    common = len(readsx.intersection(readsy))
    union = len(readsx.union(readsy))

    return common, union


def get_reads(bam, chrom, start, end, orientation, min_mapq=55):
    reads = set()

    for read in bam.fetch(chrom, start, end):
        if read.is_reverse == (orientation=="+"):
            if (read.mapq >= min_mapq) and not (read.is_secondary or read.is_supplementary):
                reads.add(read.query_name)
    return reads
