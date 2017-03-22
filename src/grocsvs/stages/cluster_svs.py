import numpy
import os
import pandas

from grocsvs import step
from grocsvs import graphing
from grocsvs.stages import refine_grid_search_breakpoints

pandas.options.display.width = 150



class ClusterSVsStep(step.StepChunk):
    """
    Takes the matrix of barcode overlaps and p-values and finds clusters
    of significant window-pairs

    Output files:
        svs.sample.dataset.chromx.chromy.pickle - list
        - each item is the x,y coordinates of a significant window-pair
    """

    @staticmethod
    def get_steps(options):
        yield ClusterSVsStep(options)


    def __init__(self, options):
        self.options = options

    def __str__(self):
        return self.__class__.__name__
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "breakpoints": os.path.join(directory, "breakpoints.tsv"),
            "edges": os.path.join(directory, "edges.tsv"),
            "graphs": os.path.join(directory, "graphs")
        }

        return paths


    def run(self):
        outpaths = self.outpaths(final=False)

        self.logger.log("Loading data...")
        self.evidence, short_frag_support, mate_pair_support = load_evidence(
            self.options)

        self.evidence["p"] = self.evidence["p_resampling"]
        
        self.evidence = precluster(self.evidence)
        # self.evidence = add_facing(self.evidence)

        print self.evidence

        self.logger.log("Building graph, pruning...")
        graph = graphing.build_graph(self.evidence)#common_counts, total_counts, binom_ps, events)
        graph = self.prune_graph(graph)

        graphs = graphing.get_subgraphs(graph)
        print "graphs:", len(graphs)
        total_breakpoints = 0
        for graph in graphs:
            total_breakpoints += sum(1 for n1,n2,data in graph.edges(data=True) if data["kind"]=="breakpoint")
        print total_breakpoints

        graphing.visualize_graphs(outpaths["graphs"], graphs, self.evidence)

        table = graphing.graphs_to_table(graphs)
        table.to_csv(outpaths["edges"], sep="\t", index=False)

        table.loc[table["kind"]=="breakpoint"].to_csv(
            outpaths["breakpoints"], sep="\t", index=False)


    def prune_graph(self, graph):
        pruned = graphing.pick_best_edges(graph)
        fixed = graphing.fix_breakpoints(self.options, pruned)

        cleaned = graphing.cleanup_fixed_graph(fixed)

        return cleaned


def load_evidence(options):
    short_frag_support = []
    mate_pair_support = []
    path = refine_grid_search_breakpoints.CombineRefinedBreakpointsStep(options)\
                             .outpaths(final=True)["refined_pairs"]
    evidence = pandas.read_table(path)

    evidence["chromx"] = evidence["chromx"].astype("string")
    evidence["chromy"] = evidence["chromy"].astype("string")

    evidence["p"] = evidence["p_resampling"]
    
    return evidence, short_frag_support, mate_pair_support

def get_positions(evidence, which="new", with_orientation=False):
    """
    gets a list of all (chrom,pos) pairs from x and y coords
    """
    x = evidence[["chromx", "{}_x".format(which)]].copy()
    x.columns = ["chrom", "pos"]

    y = evidence[["chromy", "{}_y".format(which)]].copy()
    y.columns = ["chrom", "pos"]

    if with_orientation:
        x["orientation"] = evidence["orientation"].str[0]
        y["orientation"] = evidence["orientation"].str[1]

    positions = pandas.concat([x,y], ignore_index=True)\
                      .sort_values(["chrom","pos"])\
                      .drop_duplicates()\
                      .reset_index(drop=True)
    return positions

def mean_(vals):
    return int(numpy.mean(vals))
def precluster(evidence, distance=3000):
    """
    we need to cluster putative breakpoints in 1d space so that the graphs
    can be built properly
    """

    positions = get_positions(evidence)
    transforms = []

    for chrom, p in positions.groupby("chrom"):
        p = p.copy()
        p["group"] = ((p["pos"]-p["pos"].shift()) > distance).cumsum()
        p["clustered"] = p.groupby("group")["pos"].transform(mean_)
        transforms.append(p[["chrom","pos","clustered"]])
    transform = pandas.concat(transforms)

    cur_evidence = evidence.copy()
    cur_evidence = pandas.merge(cur_evidence, transform[["chrom", "pos", "clustered"]], 
                                left_on=["chromx", "new_x"], right_on=["chrom","pos"], how="left")\
                         .drop(["chrom", "pos"], axis=1)

    cur_evidence = cur_evidence.rename(columns={"clustered": "clustered_x"})

    cur_evidence = pandas.merge(cur_evidence, transform[["chrom", "pos", "clustered"]], 
                                left_on=["chromy", "new_y"], right_on=["chrom","pos"], how="left")\
                         .drop(["chrom", "pos"], axis=1)
    cur_evidence = cur_evidence.rename(columns={"clustered": "clustered_y"})
    cur_evidence

    return cur_evidence
