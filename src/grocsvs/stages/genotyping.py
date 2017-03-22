import collections
import numpy
import os
import pandas
import pybedtools
import pysam

from grocsvs import step
from grocsvs import structuralvariants
from grocsvs.stages import final_clustering
from grocsvs.stages import refine_grid_search_breakpoints

pandas.options.display.max_rows = 500

CHUNKSIZE = 10
def nchunks(events):
    return numpy.ceil(len(events)/float(CHUNKSIZE)).astype(int)


class MergeGenotypesStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield MergeGenotypesStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "genotypes": os.path.join(directory, "genotypes.tsv")
        }

        return paths

    def run(self):
        events = load_events(self.options)
        genotypes = []

        for chunk in range(nchunks(events)):
            inpath = GenotypingStep(self.options, chunk).outpaths(final=True)["genotypes"]
            genotypes.append(pandas.read_table(inpath))

        genotypes = pandas.concat(genotypes, ignore_index=True)
        genotypes["chromx"] = genotypes["chromx"].astype("string")
        genotypes["chromy"] = genotypes["chromy"].astype("string")
        
        counts = genotypes.groupby("cluster").count()["total"]

        for row in genotypes.itertuples():
            print row.Index, row.cluster, counts[row.cluster]
            genotypes.loc[row.Index, "cluster_size"] = int(counts[row.cluster])

        genotypes["dist"] = (genotypes["x"] - genotypes["y"]).abs()
        genotypes.loc[genotypes["chromx"]!=genotypes["chromy"], "dist"] = numpy.nan

        genotypes.to_csv(self.outpaths(final=False)["genotypes"], sep="\t", index=False)



class GenotypingStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        events = load_events(options)

        for chunk in range(nchunks(events)):
            yield GenotypingStep(options, chunk)

    def __init__(self, options, chunk):
        self.options = options
        self.chunk = chunk

    def __str__(self):
        return ".".join([self.__class__.__name__, str(self.chunk)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "genotypes.{}.tsv".format(self.chunk)

        paths = {
            "genotypes": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        clustered_events = load_events(self.options, self.chunk)
        breakpoints = clustered_events.loc[clustered_events["kind"]=="breakpoint"]

        self.logger.log("Genotyping...")
        genotypes = genotype_breakpoints(breakpoints, self.options)

        columns = [c for c in self.get_columns() if c in genotypes]
        genotypes = genotypes[columns]

        self.logger.log("Detecting segmental duplications...")
        try:
            frag_length_filter(genotypes, self.options)
        except:
            self.logger.log("warning: could not accurately apply fragment length filter;"
                            " this can happen when there is not genome-wide data availabe")

        segdup_detector(genotypes, self.options)        
        compare_to_blacklist(genotypes, self.options)
        apply_filters(genotypes, self.options)

        outpath = self.outpaths(final=False)["genotypes"]
        genotypes.to_csv(outpath, sep="\t", index=False)


    def get_columns(self):
        columns = ["chromx", "x", "chromy", "y", "orientation", "cluster", 
                   "kind", "assembled", "p", "shared", "total"]

        for sample, dataset in self.options.iter_10xdatasets():
            columns.extend(map(lambda x: "{}_{}".format(sample.name, x), ["p_resampling", "shared", "total", "x_hap0", "x_hap1", "y_hap0", "y_hap1"]))

        return columns


def load_events(options, chunk=None):
    final_clustering_step = final_clustering.FinalClusterSVsStep(options)
    clustered_events_path = final_clustering_step.outpaths(final=True)["edges"]
    clustered_events = pandas.read_table(clustered_events_path)
    clustered_events["chromx"] = clustered_events["chromx"].astype("string")
    clustered_events["chromy"] = clustered_events["chromy"].astype("string")

    if chunk is not None:
        clustered_events = clustered_events.iloc[chunk*CHUNKSIZE:((chunk+1)*CHUNKSIZE)]
    # if self.options.debug:
    #     self.logger.log("LOADING ONLY 100 EVENTS....")
    #     clustered_events = clustered_events.iloc[:100]

    return clustered_events

# def get_good_bc_counts_by_dataset(options):
#     good_bc_counts_by_dataset = {}

#     for sample, dataset in options.iter_10xdatasets():
#         sample_info = options.sample_info(sample.name)
#         dataset_info = sample_info[dataset.id]
#         good_bc_counts_by_dataset[dataset.id] = dataset_info["good_bc_count"]

#     return good_bc_counts_by_dataset
        

def genotype_breakpoints(breakpoints, options):
    # TODO: gah
    dist1 = -500
    dist2 = 5000

    genotypes = []
    good_bc_counts_by_dataset, barcode_frequencies_by_dataset = \
        refine_grid_search_breakpoints.get_barcode_info(options)

    for i, breakpoint in enumerate(breakpoints.itertuples()):
        if i % 10 == 0:
            print i, "of", len(breakpoints)
        genotype = refine_grid_search_breakpoints.quantify_breakpoint(
            breakpoint.chromx, breakpoint.x, 
            breakpoint.chromy, breakpoint.y,
            breakpoint.orientation,
            options, good_bc_counts_by_dataset, 
            barcode_frequencies_by_dataset,
            dist1, dist2,
            with_phasing=True)

        genotype = genotype.rename({"new_x":"x", "new_y":"y"})
        genotype["kind"] = "breakpoint"
        genotype["cluster"] = breakpoint.cluster
        genotype["assembled"] = breakpoint.assembled

        genotypes.append(genotype)

    genotypes = pandas.DataFrame(genotypes)

    return genotypes
    

# def genotype_facing(clustered_events, options):
#     facing = clustered_events.loc[clustered_events["kind"]=="facing"]
#     import graphing

#     n1y = graphing.Node(facing["chromx"], facing["x"], facing["orientation"][0])
#     n2x = graphing.Node(facing["chromy"], facing["y"], facing["orientation"][1])

def compare_to_blacklist(genotypes, options, distance=10000):
    bad = collections.defaultdict(set)

    for blacklist_path in options.blacklists:
        blacklist_bed = ensure_file_is_bed(blacklist_path)

        for which in "xy":
            svs_bed = _convert_svs_to_bed(genotypes, which)

            c = svs_bed.closest(blacklist_bed, d=True) \
                .to_dataframe(names=["chrom","start","end","event_id",
                                     "black_chrom","black_start","black_end","black_type",
                                     "distance"])

            close = c.loc[c["distance"].abs() < distance]
            for row in close.itertuples():
                bad[row.event_id].add(row.black_type)

    genotypes["blacklist"] = ""
    genotypes["quality"] = "PASS"
    for event_id, reasons in bad.items():
        genotypes.loc[event_id, "blacklist"] = ",".join(sorted(set(reasons)))
        genotypes.loc[event_id, "quality"] = "FAIL"

def ensure_file_is_bed(path):
    if path.endswith(".bedpe"):
        table = pandas.read_table(path)
        columns = ["chromx", "startx", "endx", "chromy", "starty", "endy", "name"]
        i = 0
        while len(columns) < len(table.columns):
            columns.append("extra_{}".format(i))

        table.columns = columns
        bedx = pandas.DataFrame()
        bedx["chrom"] = table["chromx"]
        bedx["start"] = table["startx"]
        bedx["end"] =   table["endx"]
        bedx["name"] =  table["name"]

        bedy = pandas.DataFrame()
        bedy["chrom"] = table["chromy"]
        bedy["start"] = table["starty"]
        bedy["end"] =   table["endy"]
        bedy["name"] =  table["name"]

        bed = pandas.concat([bedx, bedy], ignore_index=True)

        return pybedtools.BedTool.from_dataframe(bed).sort()
    return pybedtools.BedTool(path).sort()

def _convert_svs_to_bed(table, col):
    svs = pandas.DataFrame()
    svs["chrom"] = table["chrom{}".format(col)]
    svs["start"] = table[col]
    svs["end"] = table[col] + 1
    svs["cluster"] = table.index

    svs_bed = pybedtools.BedTool.from_dataframe(svs).sort()

    return svs_bed


## Segmental duplication detection
def segdup_detector(genotypes, options):
    nearby_snv_counts = []

    for i, event in genotypes.iterrows():
        cur_nearby_snvs_x = count_nearby_snvs(options, event["chromx"], event["x"])
        cur_nearby_snvs_y = count_nearby_snvs(options, event["chromy"], event["y"])

        nearby_snv_counts.append(cur_nearby_snvs_x+cur_nearby_snvs_y)

    genotypes["nearby_snvs"] = nearby_snv_counts


def count_nearby_snvs(options, chrom, pos):
    nearby_snvs = 0
    start = pos - 500
    end = pos + 500

    for sample, dataset in options.iter_10xdatasets():
        ref_counts = []
        non_ref_counts = []

        cur_ref_counts, cur_non_ref_counts = count_ref_reads(
            dataset.bam, options.reference, chrom, start, end)

        # ref_counts.append(cur_ref_counts)
        # non_ref_counts.append(cur_non_ref_counts)

        # ref_counts = numpy.sum(ref_counts, axis=0)
        # non_ref_counts = numpy.sum(non_ref_counts, axis=0)
    
        cur_nearby_snvs = (cur_non_ref_counts/(cur_non_ref_counts+cur_ref_counts).astype(float) > 0.20).sum()
        nearby_snvs = max(nearby_snvs, cur_nearby_snvs)
        
    return nearby_snvs


def count_ref_reads(bampath, reference, chrom, start, end):
    ref_counts = numpy.zeros(end-start)
    non_ref_counts = numpy.zeros(end-start)

    bam = pysam.AlignmentFile(bampath)

    # stepper = "all" skips dupe, unmapped, secondary, and qcfail reads
    start = max(0, start)
    
    for col in bam.pileup(chrom, start, end, truncate=True, stepper="all"):
        refnuc = reference.fasta[chrom][col.reference_pos].upper()
        nuc_counts = collections.Counter()
        
        for read in col.pileups:
            if read.query_position is None:
                # nuc_counts["indel"] += 1
                pass
            else:
                nuc_counts[read.alignment.query_sequence[read.query_position]] += 1
    
        ref_counts[col.reference_pos-start] = nuc_counts[refnuc]
        non_ref_counts[col.reference_pos-start] = sum(nuc_counts.values()) - nuc_counts[refnuc]

    return ref_counts, non_ref_counts


def count_nearby_Ns(options, genotypes):
    nearby_Ns = []

    for i, event in genotypes.iterrows():
        chromx, startx, endx = event["chromx"], event["x"]-5000, event["x"]+5000
        chromy, starty, endy = event["chromy"], event["y"]-5000, event["y"]+5000

        startx = max(0, startx)
        starty = max(0, starty)
        
        seqx = options.reference.fasta[chromx][startx:endx].upper()
        seqy = options.reference.fasta[chromy][starty:endy].upper()

        nearby_Ns.append(max(seqx.count("N"), seqy.count("N")))

    return numpy.array(nearby_Ns)


def apply_filters(genotypes, options):
    # TODO: should probably more intelligently combine these filters to produce 
    # some sort of quality score
    # TODO: should learn this from the data (otherwise may not work with 
    # non-human data)
    genotypes.loc[genotypes["nearby_snvs"] >= 15, "quality"] = "FAIL"

    nearby_Ns = count_nearby_Ns(options, genotypes)
    bad_nearby_N = (nearby_Ns > 50)
    genotypes.loc[bad_nearby_N, "quality"] = "FAIL"

    bad_nearby_N_vector = (",N=" + pandas.Series(nearby_Ns[bad_nearby_N].astype(str))).values
    genotypes.loc[bad_nearby_N, "blacklist"] += bad_nearby_N_vector

    print genotypes.loc[bad_nearby_N, "blacklist"]


def frag_length_filter(genotypes, options):
    dist1 = -500
    dist2 = 5000

    flgs = {}

    for sample, dataset in options.iter_10xdatasets():
        sample_info = options.sample_info(sample.name)
        dataset_info = sample_info[dataset.id]
        flgs[sample.name] = FragLengthGenotyper(dataset_info["frag_length_distributions"])

    # for n, cluster in genotypes.groupby("cluster"):
    for event in genotypes.itertuples():
        frag_length_filter = False

        for sample, dataset in options.iter_10xdatasets():
            fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
                options, sample, dataset,
                event.chromx, event.x, event.chromy, event.y, event.orientation,
                dist1, dist2)

            if len(merged) < 5:
                continue
            
            lengths = calc_frag_lengths(event.x, event.y, event.orientation, merged)
            frag_length_filter |= flgs[sample.name].genotype(lengths)

            genotypes.loc[event.Index, "{}_lengths_50".format(sample.name)] = numpy.percentile(lengths, 50)
            genotypes.loc[event.Index, "{}_lengths_90".format(sample.name)] = numpy.percentile(lengths, 90)

        genotypes.loc[event.Index, "frag_length_passes"] = frag_length_filter
    print genotypes




def calc_frag_lengths(x, y, orientation, merged_frags):
    """
    given a merged fragments table, calculate the fragment lengths assuming the 
    fragments span the provided x and y breakpoints; adds frag lengths inplace
    """
    if orientation[0] == "+":
        partlen_x = x - merged_frags["start_pos_x"]
    else:
        partlen_x = merged_frags["end_pos_x"] - x
        
    if orientation[1] == "+":
        partlen_y = y - merged_frags["start_pos_y"]
    else:
        partlen_y = merged_frags["end_pos_y"] - y

    merged_len = partlen_x + partlen_y
    return merged_len




class FragLengthGenotyper(object):
    def __init__(self, distributions):
        self._set_distributions(distributions)

        self.quantiles = {}
        for quantile in [0.2, 0.5, 0.8]:
            self.quantiles[quantile] = [numpy.percentile(d, quantile*100) 
                                        for d in self.distributions]
        
        self.cutoffs = {}
        for quantile in sorted(self.quantiles):
            self.cutoffs[quantile] = numpy.percentile(
                self.quantiles[quantile], [2.5, 97.5])
        
        self.median_median = numpy.median(self.quantiles[0.5])
        
    def _set_distributions(self, distributions):
        """
        do some filtering on the fragment length distributions to remove 
        outliers
        """

        medians = [numpy.median(d) for d in distributions]
        sizes = [len(d) for d in distributions]

        size_cutoffs = numpy.percentile(sizes, [2.5, 97.5])
        distributions = [d for d in distributions 
                         if (size_cutoffs[0]<len(d)<size_cutoffs[1])]

        median_cutoffs = numpy.percentile(medians, [2.5, 97.5])
        distributions = [d for d in distributions 
                         if (median_cutoffs[0]<numpy.median(d)<median_cutoffs[1])]

        self.distributions = distributions

    def genotype(self, frag_lengths, lenient=False):
        geno = True
        # which = {}
        # exceptions = []
        
        for quantile in sorted(self.quantiles):
            curquantile = numpy.percentile(frag_lengths, quantile*100)

            curcutoffs = self.cutoffs[quantile].copy()
            if quantile != 0.5:
                curcutoffs *= [0.95, 1.05]
            if lenient:
                curcutoffs *= [0.95, 1.05]
            print quantile, curquantile, curcutoffs
            if not (curcutoffs[0]<curquantile<curcutoffs[1]):
                geno = False
                # which[quantile] = "F"
                # exceptions.append((curquantile, quantile))
            # else:
                # which[quantile] = "T"
        return geno#, which, exceptions

    def genotype_after_shifting(self, frag_lengths, lenient=False):
        offset = self.get_offset_at_medians(frag_lengths)
        shifted = frag_lengths + offset
        return self.genotype(shifted, lenient=lenient)
    
    def get_offset_at_medians(self, frag_lengths):
        median = numpy.median(frag_lengths)
        offset = self.median_median - median
        return offset
