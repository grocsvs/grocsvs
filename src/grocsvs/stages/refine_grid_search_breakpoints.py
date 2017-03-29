# FOR GRID SEARCH CANDIDATES

import itertools
import numpy
import os
import pandas
import scipy.stats

from grocsvs import step
from grocsvs import structuralvariants

from grocsvs.stages import sv_candidates


class CombineRefinedBreakpointsStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield CombineRefinedBreakpointsStep(options)
        
    def __init__(self, options):
        self.options = options


    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "refined_pairs": os.path.join(directory, "refined_pairs.tsv")
        }

        return paths


    def run(self):
        inputs = []
        
        chroms = self.options.reference.chroms

        for chromx, chromy in itertools.product(chroms, chroms):
            if self.options.reference.compare_chroms(chromx, chromy) < 0: continue
            input_step = RefineGridSearchBreakpointsStep(self.options, chromx, chromy)
            inpath = input_step.outpaths(final=True)["refined_pairs"]

            if os.stat(inpath).st_size > 0:
                inputs.append(pandas.read_table(inpath))
        if len(inputs) == 0:
            raise Exception("No candidate SVs discovered.")

        combined = pandas.concat(inputs)
        combined["chromx"] = combined["chromx"].astype("string")
        combined["chromy"] = combined["chromy"].astype("string")

        combined.to_csv(self.outpaths(final=False)["refined_pairs"], sep="\t", index=False)
        

class RefineGridSearchBreakpointsStep(step.StepChunk):
    """
    Takes the rough grid search (aka barcode overlaps) candidates, then performs
    breakpoint refinement to find the best potential breakpoint in the expected
    orientation
    """

    # TODO: refactor so that this and the final refine breakpoints steps share 
    # most of the refinement code
    @staticmethod
    def get_steps(options):
        chroms = options.reference.chroms

        for chromx, chromy in itertools.product(chroms, chroms):
            if options.reference.compare_chroms(chromx, chromy) < 0: continue
            yield RefineGridSearchBreakpointsStep(options, chromx, chromy)
        
    def __init__(self, options, chromx, chromy):
        self.options = options
        self.chromx = chromx
        self.chromy = chromy


    def __str__(self):
        return ".".join([self.__class__.__name__,
                     self.chromx,
                     self.chromy])

        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "refined_pairs": os.path.join(directory, "refined_pairs.{}.{}.tsv".format(self.chromx, self.chromy))
        }

        return paths


    def run(self):
        outpath = self.outpaths(final=False)["refined_pairs"]

        events = self.load_events()
        if len(events) > 0:
            refined = refine_events(events, self.options, self.logger)

            refined.to_csv(outpath, sep="\t", index=False)
        else:
            open(outpath, "w")

    def load_events(self):
        significant_events = []
        cur_events = []
        
        for sample, dataset in self.options.iter_10xdatasets():
            if self.options.reference.compare_chroms(self.chromx, self.chromy) < 0: continue

            input_step = sv_candidates.SVCandidatesStep(self.options, sample, dataset, self.chromx, self.chromy)
            inpath = input_step.outpaths(final=True)["svs"]
            if os.stat(inpath).st_size > 0:
                cur_events.append(pandas.read_table(inpath))

        if len(cur_events) > 0:
            significant_events = combine_nearby_events(pandas.concat(cur_events))
            significant_events = significant_events[["chromx", "x", "chromy", "y", "orientation"]]
            significant_events["chromx"] = significant_events["chromx"].astype("string")
            significant_events["chromy"] = significant_events["chromy"].astype("string")
            return significant_events
        else:
            return []
 

def combine_nearby_events(table, max_distance=5000):
    """
    2d-clustering of breakpoints (ie pairs of breakENDs)
    """
    if len(table) == 0:
        return table
    
    combined_tables = []

    table = table.reset_index(drop=True)
    for orientation, cur_table in table.groupby("orientation"):
        # it's already a copy, but this will suppress a warning
        cur_table = cur_table.copy()
        points = [(row.x, row.y, row.Index) for row in cur_table.itertuples()]
        clusters = structuralvariants.do_free_clustering(points, max_dist=5000)

        cur_table["cluster"] = 0
        for i, cluster in enumerate(clusters):
            for point in cluster:
                cur_table.loc[point[2], "cluster"] = i

        cur_combined = cur_table.groupby("cluster").aggregate(
            {"chromx": lambda x:x.iloc[0],
             "chromy": lambda x:x.iloc[0],
             "x": numpy.mean,
             "y": numpy.mean,
             "orientation": lambda x:x.iloc[0],             
             })
        combined_tables.append(cur_combined)

    combined_table = pandas.concat(combined_tables, ignore_index=True)
    
    combined_table["x"] = combined_table["x"].astype(int)
    combined_table["y"] = combined_table["y"].astype(int)

    return combined_table


def refine_events(events, options, logger):
    # TODO: gah
    refinement_dist1 = -20000
    refinement_dist2 = 20000
    refinement_extend = 20000

    quantification_dist1 = -500
    quantification_dist2 = 5000

    good_bc_counts_by_dataset, barcode_frequencies_by_dataset = get_barcode_info(options)

    results = []
    count = 0
    for i, event in events.iterrows():
        print ">>>", i, event.dtypes
        logger.log("{}:{}::{}:{}{}".format(event["chromx"], event["x"],
                                           event["chromy"], event["y"],
                                           event["orientation"]))
        if count % 10 == 0:
            logger.log("{} of {} ({:.0%})".format(count, len(events), count/float(len(events))))
        count += 1

        # First get better breakpoints
        refined = refine_breakpoint(
            event["chromx"], event["x"], 
            event["chromy"], event["y"], 
            event["orientation"], options,
            refinement_dist1, refinement_dist2, refinement_extend)

        if refined is None:
            continue
        newx, newy = refined

        # Next quantify the event based on the better breakpoint loci
        quantification = quantify_breakpoint(
            event["chromx"], newx, 
            event["chromy"], newy, 
            event["orientation"],
            options, good_bc_counts_by_dataset,
            barcode_frequencies_by_dataset,
            quantification_dist1, quantification_dist2)

        quantification["original_x"] = event["x"]
        quantification["original_y"] = event["y"]
        results.append(quantification)

    return pandas.DataFrame(results)



def get_shared_frags(options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2):
    fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
        options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2)

    bcx = set(fragsx["bc"])
    bcy = set(fragsy["bc"])

    common_barcodes = bcx.intersection(bcy)

    shared_fragsx = fragsx.loc[fragsx["bc"].isin(common_barcodes)]
    shared_fragsy = fragsy.loc[fragsy["bc"].isin(common_barcodes)]

    return shared_fragsx, shared_fragsy

def refine_breakpoint(chromx, x, chromy, y, orientation, 
                      options, dist1, dist2, extend):
    shared_fragsx = []
    shared_fragsy = []

    # because all we're concerned with for refinement is the fragments
    # with common barcodes across the breakpoint, we'll do refinement
    # with all datasets without worrying if a given dataset supports
    # the event
    for sample, dataset in options.iter_10xdatasets():
        cur_fragsx, cur_fragsy = get_shared_frags(
            options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2)
        shared_fragsx.append(cur_fragsx)
        shared_fragsy.append(cur_fragsy)

    shared_fragsx = pandas.concat(shared_fragsx)
    shared_fragsy = pandas.concat(shared_fragsy)

    if len(shared_fragsx) < 1:
        return None

    breakpointx = get_breakpoint(shared_fragsx, x, orientation[0], extend)
    breakpointy = get_breakpoint(shared_fragsy, y, orientation[1], extend)

    return breakpointx, breakpointy


def get_breakpoint(frags, pos, orientation, extend=20000):
    density = numpy.zeros(extend*2)
    for i, frag in frags.iterrows():
        curstart = max(frag["start_pos"]-(pos-extend), 0)
        curend = min(frag["end_pos"]-(pos-extend), len(density))

        density[curstart:curend] += 1

    peaks = numpy.where(density>(0.9*density.max()))[0]
    if orientation == "+":
        peak = peaks[0]
    elif orientation == "-":
        peak = peaks[-1]
    else:
        raise Exception("unknown orientation: {}".format(orientation))
        
    diff = density[peak] - density
    dist = numpy.sqrt(numpy.abs(numpy.arange(len(density))-peak))
    score = numpy.ma.masked_array(diff / dist.astype(float), mask=False)
    score.mask[numpy.isnan(score)] = True
    
    if orientation == "+":
        score.mask[numpy.arange(0, peak)] = True
    elif orientation == "-":
        score.mask[numpy.arange(peak, len(score))] = True
    else:
        raise Exception("unknown orientation: {}".format(orientation))

    breakpoint = numpy.where(score==score.max())[0][0]
    breakpoint += pos - extend
    return breakpoint


def get_barcode_info(options):
    good_bc_counts_by_dataset = {}
    barcode_frequencies_by_dataset = {}

    for sample, dataset in options.iter_10xdatasets():
        sample_info = options.sample_info(sample.name)
        dataset_info = sample_info[dataset.id]
        good_bc_counts_by_dataset[dataset.id] = dataset_info["good_bc_count"]

        sample_info = options.sample_info(sample.name)
        dataset_info = sample_info[dataset.id]
        barcode_frequencies = dataset_info["barcode_read_totals"]
        barcode_frequencies /= barcode_frequencies.sum().astype(float)
        barcode_frequencies = barcode_frequencies.values
        barcode_frequencies_by_dataset[dataset.id] = barcode_frequencies

    return good_bc_counts_by_dataset, barcode_frequencies_by_dataset


def quantify_breakpoint(chromx, x, chromy, y, orientation, 
                        options, good_bc_counts_by_dataset, barcode_frequencies_by_dataset,
                        dist1, dist2, with_phasing=False):
    
    cur_result = {}
    cur_result["chromx"] = chromx
    cur_result["new_x"] = x
    cur_result["chromy"] = chromy
    cur_result["new_y"] = y
    cur_result["orientation"] = orientation

    cur_result["shared"] = 0
    cur_result["total"] = 0

    for sample, dataset in options.iter_10xdatasets():
        barcode_frequencies = barcode_frequencies_by_dataset[dataset.id]

        fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
            options, sample, dataset,
            chromx, x, chromy, y, orientation, dist1, dist2,
            with_phasing=with_phasing)

        bcx = set(fragsx["bc"])
        bcy = set(fragsy["bc"])


        common_barcodes = bcx.intersection(bcy)
        total_barcodes = bcx.union(bcy)
        cur_result["{}_total".format(sample.name)] = len(total_barcodes)

        if len(common_barcodes) < 1:
            continue

        good_bc_count = good_bc_counts_by_dataset[dataset.id]
        contingency_table = numpy.array([[len(common_barcodes), len(bcx-bcy)],
                                         [len(bcy-bcx), good_bc_count-len(total_barcodes)]])
        p_fisher = scipy.stats.fisher_exact(contingency_table, alternative="greater")[1]
        p_resampling = structuralvariants.score_event(
            len(bcx), len(bcy), len(common_barcodes), barcode_frequencies, resamples=100)

        cur_result["{}_shared".format(sample.name)] = len(common_barcodes)

        cur_result["{}_p_fisher".format(sample.name)] = p_fisher
        cur_result["{}_p_resampling".format(sample.name)] = p_resampling


        if with_phasing:
            cur_result["{}_x_hap0".format(sample.name)] = (merged["hap_x"] == 0).sum()
            cur_result["{}_x_hap1".format(sample.name)] = (merged["hap_x"] == 1).sum()
            cur_result["{}_y_hap0".format(sample.name)] = (merged["hap_y"] == 0).sum()
            cur_result["{}_y_hap1".format(sample.name)] = (merged["hap_y"] == 1).sum()

        # TODO: constants should be constant across steps!
        if (p_resampling < 1e-4) and (len(common_barcodes)/float(len(total_barcodes)) > 0.10):
            cur_result["shared"] += len(common_barcodes)
            cur_result["total"] += len(total_barcodes)

    cur_result["p_resampling"] = min(cur_result.get("{}_p_resampling".format(sample_name), 1.0)
                          for sample_name in options.samples)

    return pandas.Series(cur_result)
