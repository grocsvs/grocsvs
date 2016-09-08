import glob
import json
import logging
import numpy
import os
import pandas
import scipy.stats

from grocsvs import datasets as svdatasets
from grocsvs import step
from grocsvs import structuralvariants
from grocsvs import utilities

from grocsvs.stages import supporting_barcodes
from grocsvs.stages import sv_candidates


def get_events(options, sample, dataset, with_barcodes=False):
    which = "event_barcodes" if with_barcodes else "events"

    path = supporting_barcodes.SupportingBarcodesStep(options, sample, dataset)\
                              .outpaths(final=True)[which]

    return utilities.pickle.load(open(path))


class PairEvidenceStep(step.StepChunk):
    """
    """
    @staticmethod
    def get_steps(options, sample=None, dataset=None):
        for sample, dataset in options.iter_10xdatasets():
            for step in PairEvidenceStep.get_steps_for_dataset(
                options, sample, dataset):
                yield step

    @staticmethod
    def get_steps_for_dataset(options, sample, dataset):
        events = get_events(options, sample, dataset)

        nchunks = min(100, len(events)/5)
        #nchunks = 1
        for i in range(nchunks):
            step = PairEvidenceStep(options, sample, dataset, i, nchunks)
            yield step

    def __init__(self, options, sample, dataset, chunk, nchunks):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chunk = chunk
        self.nchunks = nchunks

        assert isinstance(self.dataset, svdatasets.TenXDataset)


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         str(self.chunk)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        evidence = "pair_evidence.{}.{}.{}.tsv".format(
            self.sample.name, self.dataset.id, self.chunk)

        common = "common.{}.{}.{}.tsv".format(
            self.sample.name, self.dataset.id, self.chunk)
        total = "total.{}.{}.{}.tsv".format(
            self.sample.name, self.dataset.id, self.chunk)

        paths = {
            "evidence": os.path.join(directory, evidence),
            "common" : os.path.join(directory, common),
            "total" : os.path.join(directory, total),
        }

        return paths


    def run(self):
        self.logger.log("running!")

        events_to_barcodes = get_events(
            self.options, self.sample, self.dataset, with_barcodes=True)
        events = get_events(self.options, self.sample, self.dataset)
        
        good_bc_count = utilities.get_good_bc_count(self)
        evidence = get_pair_evidence(
            events, events_to_barcodes, self.chunk, self.nchunks,
            good_bc_count)

        self.logger.log("1")
        common_df = pandas.DataFrame(evidence[0], index=events, columns=events)
        common_df.to_csv(self.outpaths(final=False)["common"])
        self.logger.log("2")
        total_df = pandas.DataFrame(evidence[1], index=events, columns=events)
        total_df.to_csv(self.outpaths(final=False)["total"])
        self.logger.log("3")
        significant_events = get_significant_events(evidence, events)
        self.logger.log("4")
        # evidence = refine_events(events, evidence,
        #                          self.options, self.sample, self.dataset, good_bc_count)
        # TODO: combine nearby events, probably after combining tables from all chunks

        outpath = self.outpaths(final=False)["evidence"]
        significant_events.to_csv(outpath, sep="\t", index=False)
        self.logger.log("5")
        # with open(outpath, "w") as f:
        #     utilities.pickle.dump(evidence, f, protocol=-1)



def get_mat_coords(events, chunk, nchunks):
    length = len(events)
    chunk_size = (length**2 / nchunks) + 1
    from_ = chunk_size*chunk
    to_ = chunk_size*(chunk+1)

    return chunk_size, length, from_, to_

def get_pair_evidence(events, events_to_barcodes, chunk, nchunks, good_bc_count):
    chunk_size, length, from_, to_ = get_mat_coords(events, chunk, nchunks)

    common_counts = numpy.empty((length, length))
    common_counts[:] = numpy.nan

    total_counts = common_counts.copy()
    binom_ps = common_counts.copy()

    # this gobbedlygook is to get it to parallelize evenly even
    # though we're only filling the upper diagonal of the matrix
    for i in range(length):
        barcodesx = events_to_barcodes[events[i]]
        row_common = []
        row_totals = []
        row_start = None

        for j in range(length):
            if j <= i: continue
            n = i * length + j
            if n < from_: continue
            if n >= to_: break

            if row_start is None:
                row_start = j

            barcodesy = events_to_barcodes[events[j]]

            common = len(barcodesx.intersection(barcodesy))
            total = len(barcodesx.union(barcodesy))
            common_counts[i,j] = common
            total_counts[i,j] = total

            row_common.append(common)
            row_totals.append(total)

        if row_start is not None:
            row_expected = len(barcodesx)/float(good_bc_count)
            curps = scipy.stats.binom.sf(row_common, row_totals, row_expected)
            binom_ps[i, row_start:row_start+len(curps)] = curps

    return common_counts, total_counts, binom_ps

def _facing(chromx, x, orientationx, chromy, y, orientationy):
    # TODO: fix constant
    if chromx == chromy and abs(x-y)<100000:
        if (((x < y) and orientationx=="-" and orientationy=="+") or
            ((x > y) and orientationx=="+" and orientationy=="-")):
            if abs(x-y)<1000:
                return "bad"
            return "facing"

    return "breakpoint"


def get_significant_events(evidence, events):
    common_counts, total_counts, binom_ps = evidence

    significant_events = []
    #significant_coords = numpy.transpose(numpy.where( (binom_ps<1e-4) | (common_counts>100) ))
    significant_coords = numpy.transpose(numpy.where( ((common_counts/(total_counts.astype(float)) > 0.01) & (common_counts >= 15)) ))
    
    for i, j in significant_coords:
        # if common_counts[i,j] < 15: continue

        chromx, x, orientationx = events[i]
        chromy, y, orientationy = events[j]

        facing = _facing(chromx, x, orientationx, chromy, y, orientationy)    
        if facing == "bad":
            continue

        significant_events.append(
            [chromx, x, chromy, y, orientationx+orientationy,
             common_counts[i,j], total_counts[i,j], binom_ps[i,j], facing])

    columns = ["chromx", "x", "chromy", "y", "orientation",
               "shared", "total", "p", "facing"]
    significant_events = pandas.DataFrame(
        significant_events, columns=columns)

    return significant_events
# def refine_events(events, evidence, options, sample, dataset, good_bc_count):
#     # TODO: gah
#     refinement_dist1 = -5000
#     refinement_dist2 = 15000
#     refinement_extend = 20000

#     quantification_dist1 = -500
#     quantification_dist2 = 5000

#     common_counts, total_counts, binom_ps = evidence

#     evidence = []
#     xx = numpy.transpose(numpy.where(binom_ps<1e-4))
#     for i, j in xx:
#         if common_counts[i,j] < 15: continue

#         chromx, x, orientationx = events[i]
#         chromy, y, orientationy = events[j]


#         facing = _facing(chromx, x, orientationx, chromy, y, orientationy)
#         if facing in ["facing", "bad"]:
#             newx, newy = x, y
#             cur_common_counts = common_counts[i,j]
#             cur_total_counts = total_counts[i,j]
#             p = binom_ps[i,j]
#         else:
#             # First get better breakpoints
#             refined = refine_breakpoint(
#                 chromx, x, 
#                 chromy, y, 
#                 orientationx+orientationy, 
#                 options, sample, dataset, 
#                 refinement_dist1, refinement_dist2, refinement_extend)

#             if refined is None:
#                 continue
#             newx, newy = refined

#             # Next quantify the event based on the better breakpoint loci
#             quantification = quantify_breakpoint(
#                 chromx, newx, 
#                 chromy, newy, 
#                 orientationx+orientationy, 
#                 options, sample, dataset, 
#                 good_bc_count, quantification_dist1, quantification_dist2)
#             # if abs(200903945-x) < 1e6 or abs(200903945-y) < 1e6:
#             #     print "HERE"*6
#             #     print "::", refined
#             #     print "..", quantification

#             if quantification is None:
#                 continue
#             cur_common_counts, cur_total_counts, p = quantification

#         evidence.append(
#             [chromx, x, newx, chromy, y, newy, orientationx+orientationy,
#              cur_common_counts, cur_total_counts, p, facing])

#     columns = ["chromx", "original_x", "new_x",
#                "chromy", "original_y", "new_y", "orientation",
#                "shared", "total", "p", "facing"]

#     return pandas.DataFrame(evidence, columns=columns)



# def refine_breakpoint(chromx, x, chromy, y, orientation, 
#                       options, sample, dataset, 
#                       dist1, dist2, extend):

#     fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
#         options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2)

#     bcx = set(fragsx["bc"])
#     bcy = set(fragsy["bc"])

#     common_barcodes = bcx.intersection(bcy)
#     if len(common_barcodes) < 1:
#         return None

#     shared_fragsx = fragsx.loc[fragsx["bc"].isin(common_barcodes)]
#     shared_fragsy = fragsy.loc[fragsy["bc"].isin(common_barcodes)]

#     breakpointx = sv_candidates.get_breakpoint(shared_fragsx, x, orientation[0], extend)
#     breakpointy = sv_candidates.get_breakpoint(shared_fragsy, y, orientation[1], extend)

#     return breakpointx, breakpointy


# def quantify_breakpoint(chromx, x, chromy, y, orientation, 
#                         options, sample, dataset, 
#                         good_bc_count, dist1, dist2):

#     fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
#         options, sample, dataset, chromx, x, chromy, y, orientation, dist1, dist2)

#     bcx = set(fragsx["bc"])
#     bcy = set(fragsy["bc"])

#     common_barcodes = bcx.intersection(bcy)
#     if len(common_barcodes) < 1:
#         return None

#     total_barcodes = bcx.union(bcy)

#     contingency_table = numpy.array([[len(common_barcodes), len(bcx-bcy)],
#                                      [len(bcy-bcx), good_bc_count-len(total_barcodes)]])
#     p = scipy.stats.fisher_exact(contingency_table, alternative="greater")[1]

#     return len(common_barcodes), len(total_barcodes), p
