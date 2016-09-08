import os
import pandas

from grocsvs import step
from grocsvs import graphing

from grocsvs.stages import walk_assemblies
from grocsvs.stages import refine_breakpoints


pandas.options.display.width = 150



class FinalClusterSVsStep(step.StepChunk):
    """
    """

    @staticmethod
    def get_steps(options):
        yield FinalClusterSVsStep(options)

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
            "graphs": os.path.join(directory, "graphs")#,
            # "raw": os.path.join(directory, "raw")
        }

        return paths


    def run(self):
        outpaths = self.outpaths(final=False)

        self.logger.log("Loading data...")
        evidence = load_evidence(self.options)
        assembly_evidence = load_assembly_evidence(self.options)

        self.logger.log("Building graph, pruning...")
        graph = graphing.build_graph(evidence)
        graph = self.prune_graph(graph, assembly_evidence)

        graphs = graphing.get_subgraphs(graph)
        print "graphs:", len(graphs)
        total_breakpoints = 0
        for graph in graphs:
            print graph.nodes()
            total_breakpoints += sum(1 for n1,n2,data in graph.edges(data=True) 
                                     if data["kind"]=="breakpoint")
        print total_breakpoints

        table = graphing.graphs_to_table(graphs)
        table.to_csv(outpaths["edges"], sep="\t", index=False)

        table.loc[table["kind"]=="breakpoint"].to_csv(
            outpaths["breakpoints"], sep="\t", index=False)

        graphing.visualize_graphs(outpaths["graphs"], graphs, evidence)
        # graphing.visualize_frags(outpaths["raw"], graphs, self.options)


    def prune_graph(self, graph, assembly_evidence):
        pruned = pick_assembled_breakpoints(graph, assembly_evidence)
        pruned = graphing.pick_best_edges(pruned)
        fixed = graphing.fix_breakpoints(self.options, pruned)

        cleaned = graphing.cleanup_fixed_graph(fixed)

        return cleaned


def load_evidence(options):
    evidences = []

    input_steps = refine_breakpoints.RefineBreakpointsWithAssembliesStep.get_steps(options)

    for input_step in input_steps:
        path = input_step.outpaths(final=True)["refined_pairs"]
        try:
            evidence = pandas.read_table(path)
            evidences.append(evidence)
        except ValueError:
            print "empty file", path
            print open(path).read()

    evidence = pandas.concat(evidences).dropna(axis=1, how="all")
    evidence = evidence.rename(columns={"original_x":"clustered_x",
                                        "original_y":"clustered_y",
                                        "p_resampling":"p"})

    print evidence
    return evidence

def load_assembly_evidence(options):
    path = walk_assemblies.WalkAssembliesStep(options) \
                          .outpaths(final=True)["walk_assemblies"]
    evidence = pandas.read_table(path)

    return evidence


def pick_assembled_breakpoints(graph, assembly_evidence):
    pruned = graph.copy()

    for event in assembly_evidence.itertuples():
        nodex = graphing.Node(event.chromx, event.x, event.orientationx)
        nodey = graphing.Node(event.chromy, event.y, event.orientationy)

        if not pruned.has_edge(nodex, nodey):
            continue
        data = pruned[nodex][nodey]
        if data["ratio"] < 0.05 or data["shared"] < 10:
            continue

        data["assembled"] = True

        to_delete = set()
        for n1,n2,data in pruned.edges_iter((nodex, nodey), data=True):
            if (data["kind"]=="breakpoint") and set([n1,n2]) != set([nodex, nodey]):
                to_delete.add((n1,n2))
        pruned.remove_edges_from(to_delete)
    
    return pruned

