import os
import pandas

from grocsvs import step
from grocsvs.stages import refine_grid_search_breakpoints
from grocsvs.stages import walk_assemblies


class PostAssemblyMergeStep(step.StepChunk):
    """
    """

    @staticmethod
    def get_steps(options):
        yield PostAssemblyMergeStep(options)
        
    def __init__(self, options):
        self.options = options


    def __str__(self):
        return self.__class__.__name__
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "merged_candidates": os.path.join(directory, "merged_candidates.tsv")
        }

        return paths

    def run(self):
        outpath = self.outpaths(final=False)["merged_candidates"]

        events = self.load_events()

        events = combine_nearby_breakends(events)
        events.to_csv(outpath, sep="\t", index=False)
 

    def load_events(self):
        walk_assemblies_step = walk_assemblies.WalkAssembliesStep(self.options)
        assembled_events_path = walk_assemblies_step.outpaths(final=True)["walk_assemblies"]
        assembled_events = pandas.read_table(assembled_events_path)
        assembled_events["assembled"] = True

        refine_breakpoints_step = refine_grid_search_breakpoints.CombineRefinedBreakpointsStep(
            self.options)
        refine_breakpoints_path = refine_breakpoints_step.outpaths(final=True)["refined_pairs"]
        refine_breakpoints_events = pandas.read_table(refine_breakpoints_path)
        refine_breakpoints_events["assembled"] = False
        refine_breakpoints_events["x"] = refine_breakpoints_events["new_x"]
        refine_breakpoints_events["y"] = refine_breakpoints_events["new_y"]

        refine_breakpoints_events["orientationx"] = refine_breakpoints_events["orientation"].str[0]
        refine_breakpoints_events["orientationy"] = refine_breakpoints_events["orientation"].str[1]

        events = pandas.concat([assembled_events, refine_breakpoints_events])

        return events


def combine_nearby_breakends(events, distance=5000):
    """
    1d clustering, prioritizing assembled breakpoint coords
    """

    breakends = []

    positions = get_positions(events)

    for (chrom, orientation), cur_events in positions.groupby(["chrom", "orientation"]):
        cur_events = cur_events.sort_values("pos")
        groups = ((cur_events["pos"]-cur_events["pos"].shift()) > distance).cumsum()

        for i, cur_group in cur_events.groupby(groups):
            if cur_group["assembled"].any():
                cur_combined = cur_group.loc[cur_group["assembled"]].copy()
                cur_combined["assembled"] = True
            else:
                cur_orientations = cur_group["orientation"].unique()
                cur_combined = pandas.DataFrame({"orientation":cur_orientations})
                cur_combined["chrom"] = chrom
                cur_combined["pos"] = int(cur_group["pos"].mean())
                cur_combined["assembled"] = False

            breakends.append(cur_combined)

    return pandas.concat(breakends, ignore_index=True)



def get_positions(evidence):
    """
    gets a list of all (chrom,pos) pairs from x and y coords
    """

    x = evidence[["chromx", "x", "orientationx", "assembled"]].copy()
    x.columns = ["chrom", "pos", "orientation", "assembled"]

    y = evidence[["chromy", "y", "orientationy", "assembled"]].copy()
    y.columns = ["chrom", "pos", "orientation", "assembled"]


    positions = pandas.concat([x,y], ignore_index=True)\
                      .sort_values(["chrom","orientation", "pos"])\
                      .drop_duplicates()\
                      .reset_index(drop=True)
    return positions
