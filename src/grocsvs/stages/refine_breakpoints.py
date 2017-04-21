# FOR FRAG ENDS

import numpy
import os
import pandas
# import scipy.stats

from grocsvs import step
from grocsvs import structuralvariants

from grocsvs.stages import pair_evidence
from grocsvs.stages import refine_grid_search_breakpoints
from grocsvs.stages import walk_assemblies


class RefineBreakpointsWithAssembliesStep(step.StepChunk):
    """
    """

    @staticmethod
    def get_steps(options):
        input_steps = pair_evidence.PairEvidenceStep.get_steps(options)

        chunks = set()
        for input_step in input_steps:
            chunks.add(input_step.chunk)

        for chunk in sorted(chunks):
            yield RefineBreakpointsWithAssembliesStep(options, chunk)
        
    def __init__(self, options, chunk):
        self.options = options
        self.chunk = chunk


    def __str__(self):
        return ".".join(map(str, [self.__class__.__name__, self.chunk]))
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "refined_pairs": os.path.join(directory, "refined_pairs.{}.tsv".format(self.chunk))
        }

        return paths


    def run(self):
        events, assembly_evidence = self.load_events()
        refined = refine_events(events, assembly_evidence, self.options, self.logger)
        self.save_refined_events(refined)


    def load_events(self):
        input_steps = pair_evidence.PairEvidenceStep.get_steps(self.options)

        significant_events = []
        for input_step in input_steps:
            if input_step.chunk == self.chunk:
                inpath = input_step.outpaths(final=True)["evidence"]
                significant_events.append(pandas.read_table(inpath))

        print significant_events
        significant_events = pandas.concat(significant_events)
        
        significant_events = significant_events[["chromx", "x", "chromy", "y", "orientation", "facing"]]
        significant_events["x"] = significant_events["x"].astype(int)
        significant_events["y"] = significant_events["y"].astype(int)
        significant_events["chromx"] = significant_events["chromx"].astype("string")
        significant_events["chromy"] = significant_events["chromy"].astype("string")
        significant_events = significant_events.drop_duplicates()

        assemblies_path = walk_assemblies.WalkAssembliesStep(self.options) \
                              .outpaths(final=True)["walk_assemblies"]
        assembly_evidence = pandas.read_table(assemblies_path)
        assembly_evidence["chromx"] = assembly_evidence["chromx"].astype("string")
        assembly_evidence["chromy"] = assembly_evidence["chromy"].astype("string")

        return significant_events, assembly_evidence

 
    def save_refined_events(self, refined):
        outpath = self.outpaths(final=False)["refined_pairs"]

        print "::", len(refined)
        if len(refined) > 0:
            refined = refined.loc[(refined["new_x"]>=0) & (refined["new_y"]>=0)]
            print "::", len(refined)

        refined.to_csv(outpath, sep="\t", index=False)

        if len(refined) == 0:
            with open(outpath, "w") as outf:
                outf.write("<empty>\n")


def is_assembled(event, assembly_evidence):
    # TODO: would probably be simpler if we made a copy of event and did this using a 
    # table merge
    if ((assembly_evidence["chromx"]==event["chromx"]) &
        (assembly_evidence["chromy"]==event["chromy"]) &
        (assembly_evidence["x"]==event["x"]) &
        (assembly_evidence["y"]==event["y"]) &
        (assembly_evidence["orientationx"]==event["orientation"][0]) &
        (assembly_evidence["orientationy"]==event["orientation"][1])).any():
        return True
    if ((assembly_evidence["chromy"]==event["chromx"]) &
        (assembly_evidence["chromx"]==event["chromy"]) &
        (assembly_evidence["y"]==event["x"]) &
        (assembly_evidence["x"]==event["y"]) &
        (assembly_evidence["orientationy"]==event["orientation"][0]) &
        (assembly_evidence["orientationx"]==event["orientation"][1])).any():
        return True
    return False

def refine_events(events, assembly_evidence, options, logger):
    # TODO: gah
    refinement_dist1 = -2000
    refinement_dist2 = 5000
    refinement_extend = 20000

    quantification_dist1 = -500
    quantification_dist2 = 5000

    good_bc_counts_by_dataset, barcode_frequencies_by_dataset = \
        refine_grid_search_breakpoints.get_barcode_info(options)

    results = []
    count = 0
    ass_count = 0
    for i, event in events.iterrows():
        #if count % 10 == 0:
        logger.log("{} of {} ({:.0%})".format(count, len(events), count/float(len(events))))
        count += 1

        if event["facing"] != "breakpoint":
            continue

        refined = None

        assembled = is_assembled(event, assembly_evidence)
        if assembled:
            ass_count += 1
            refined = event["x"], event["y"]
        else:
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
        quantification = refine_grid_search_breakpoints.quantify_breakpoint(
            event["chromx"], newx, 
            event["chromy"], newy, 
            event["orientation"],
            options, good_bc_counts_by_dataset,
            barcode_frequencies_by_dataset,
            quantification_dist1, quantification_dist2)

        # TODO: standardize or at least document the x, new_x, original_x and clustered_x
        # columns as they appear in different files
        quantification["original_x"] = event["x"]
        quantification["original_y"] = event["y"]
        quantification["assembled"] = assembled

        results.append(quantification)

    print "*"*100, ass_count, len(results)
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


# def quantify_breakpoint(chromx, x, chromy, y, orientation, 
#                         options, good_bc_counts_by_dataset, barcode_frequencies_by_dataset,
#                         dist1, dist2):
    
#     cur_result = {}
#     cur_result["chromx"] = chromx
#     cur_result["new_x"] = x
#     cur_result["chromy"] = chromy
#     cur_result["new_y"] = y
#     cur_result["orientation"] = orientation

#     cur_result["shared"] = 0
#     cur_result["total"] = 0

#     for sample, dataset in options.iter_10xdatasets():
#         barcode_frequencies = barcode_frequencies_by_dataset[dataset.id]

#         fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
#             options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2)

#         bcx = set(fragsx["bc"])
#         bcy = set(fragsy["bc"])

#         common_barcodes = bcx.intersection(bcy)
#         if len(common_barcodes) < 1:
#             continue

#         total_barcodes = bcx.union(bcy)

#         good_bc_count = good_bc_counts_by_dataset[dataset.id]
#         contingency_table = numpy.array([[len(common_barcodes), len(bcx-bcy)],
#                                          [len(bcy-bcx), good_bc_count-len(total_barcodes)]])
#         p_fisher = scipy.stats.fisher_exact(contingency_table, alternative="greater")[1]

#         p_resampling = structuralvariants.score_event(
#             len(bcx), len(bcy), len(common_barcodes), barcode_frequencies, resamples=100)

#         cur_result["{}_shared".format(sample.name)] = len(common_barcodes)
#         cur_result["{}_total".format(sample.name)] = len(total_barcodes)
#         cur_result["{}_p_fisher".format(sample.name)] = p_fisher
#         cur_result["{}_p_resampling".format(sample.name)] = p_resampling

#         # TODO: constants should be constant across steps!
#         if (p_resampling < 1e-4) and (len(common_barcodes)/float(len(total_barcodes)) > 0.05):
#             cur_result["shared"] += len(common_barcodes)
#             cur_result["total"] += len(total_barcodes)

#     cur_result["p_resampling"] = min(cur_result.get("{}_p_resampling".format(sample_name), 1.0)
#                           for sample_name in options.samples)

#     return pandas.Series(cur_result)

def get_breakpoint(frags, pos, orientation, extend=20000):
    density = numpy.zeros(extend*2)
    for i, frag in frags.iterrows():
        curstart = max(frag["start_pos"]-(pos-extend), 0)
        curend = min(frag["end_pos"]-(pos-extend), len(density))

        density[int(curstart):int(curend)] += 1

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
