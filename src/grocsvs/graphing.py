import collections
import os
import networkx
import numpy
import pandas
import pysam
import scipy.stats

from grocsvs import structuralvariants
from grocsvs import datasets as svdatasets
from grocsvs import utilities
from grocsvs.stages import call_readclouds
 
Node = collections.namedtuple("Node", "chrom position orientation")


def build_graph(evidence):
    graph = networkx.Graph()

    # TODO: constants: p, min_shared

    for i, row in evidence.iterrows():

        if row["p"] > 1e-4: continue
        if row["total"] < 20: continue
        if row["shared"]/float(row["total"]) < 0.01: continue
        if row["shared"] < 20: continue

        chromx, x, orientationx = row["chromx"], int(row["clustered_x"]), row["orientation"][0]
        chromy, y, orientationy = row["chromy"], int(row["clustered_y"]), row["orientation"][1]


        from_ = Node(chromx, x, orientationx)
        to_ = Node(chromy, y, orientationy)

        kind = _facing(chromx, x, orientationx, chromy, y, orientationy)
        if kind == "bad": continue
        if kind == "facing": continue

        graph.add_edge(from_, to_,
            kind = kind,
            shared = row["shared"],
            total = row["total"],
            ratio = row["shared"]/float(row["total"]),
            p = row["p"],
            chrom_x = chromx,
            chrom_y = chromy,
            orientation = orientationx+orientationy,
            new_x = int(row["new_x"]),
            new_y = int(row["new_y"]),
            id = i)#,

    return graph

def _facing(chromx, x, orientationx, chromy, y, orientationy):
    # TODO: fix constant
    if chromx == chromy and abs(x-y)<100000:
        if (((x < y) and orientationx=="-" and orientationy=="+") or
            ((x > y) and orientationx=="+" and orientationy=="-")):
            if abs(x-y)<1000:
                return "bad"
            return "facing"

    return "breakpoint"



def pick_best_edges(graph):
    pruned = graph.copy()
    # prune_reasons = {}
    remaining_nodes = pruned.nodes()

    # prune by taking best edge amongst nodes with trans-degree > 1
    # and removing all other trans-edges

    while len(remaining_nodes) > 0:
        best_edge = None
        best_ratio = None
        to_delete = set()

        for n1,n2,data in pruned.edges(remaining_nodes, data=True):
            if data["kind"] != "breakpoint": #continue
                to_delete.add((n1,n2))

            if ((best_edge is None) or (data["ratio"]>best_ratio)):
                best_edge = n1,n2
                best_ratio = data["ratio"]
        
        if best_edge is None:
            break
        for n1,n2,data in pruned.edges(best_edge, data=True):
            if (data["kind"]=="breakpoint") and set([n1,n2]) != set(best_edge):
                to_delete.add((n1,n2))
        pruned.remove_edges_from(to_delete)
        remaining_nodes = [node for node in remaining_nodes if node not in best_edge]
    
    return pruned


def cleanup_fixed_graph(pruned):
    # remove cycles
    while True:
        cycles = networkx.cycle_basis(pruned)
        if len(cycles) == 0:
            break

        to_delete = None
        worst_val = None
        for cycle in cycles:
            cur_worst_val = None
            cur_worst = None
            for a,b,data in pruned.edges(cycle, data=True):
                cur_val = data["p"]
                if cur_worst_val is None or cur_val < cur_worst_val:
                    cur_worst_val = cur_val
                    cur_worst = (a,b)
            if worst_val is None or cur_worst_val < worst_val:
                worst_val = cur_worst_val
                to_delete = cur_worst

        pruned.remove_edge(*to_delete)

    # remove all cis-edges at the ends of subgraphs
    degrees = pruned.degree()
    to_delete = []
    for node, degree in dict(degrees).items():
        if degree == 1:
            edge = list(pruned.edges([node], data=True))[0]
            if edge[2]["kind"]=="facing":
                to_delete.append(node)
    pruned.remove_nodes_from(to_delete)

    # remove unconnected nodes
    pruned.remove_nodes_from([node for (node, degree) in dict(pruned.degree()).items() if degree==0])

    return pruned

def get_subgraphs(graph):
    subgraphs = []
    for subgraph in networkx.connected_components(graph):
        if len(subgraph) > 0:
            subgraphs.append(graph.subgraph(subgraph))
    return subgraphs


def get_node_label(node):
    if node.orientation == "+":
         return "]{}:{:,}".format(node.chrom, int(node.position))
    else:
        return "{}:{:,}[".format(node.chrom, int(node.position))

def graphs_to_table(graphs):
    tables = []

    for i, graph in enumerate(graphs):
        column_names = ["chromx", "x", "chromy", "y", "orientation", 
                        "support_shared", "support_total", "p",
                        "kind", "cluster", "assembled"]
        rows = []

        for nodex, nodey, data in graph.edges(data=True):
            chromx, x, orientationx = nodex
            chromy, y, orientationy = nodey
            orientation = orientationx+orientationy
            support_shared = data["shared"]
            support_total = data["total"]
            p = data["p"]
            kind = data["kind"]
            assembled = True if "assembled" in data else False
            # new_x = data["new_x"]
            # new_y = data["new_y"]

            rows.append([chromx, x, chromy, y, orientation,  
                         support_shared, support_total, p,
                         kind, i, assembled])

        cur_table = pandas.DataFrame(rows, columns=column_names)
        tables.append(cur_table)

    table = pandas.concat(tables).reset_index(drop=True)

    return table



def fix_breakpoints(options, graph):
    fixed = networkx.Graph()

    # good_bc_count = pair_evidence.get_good_bc_count(self)

    good_bc_counts_by_dataset = {}

    for sample, dataset in options.iter_10xdatasets():
        sample_info = options.sample_info(sample.name)
        dataset_info = sample_info[dataset.id]
        good_bc_counts_by_dataset[dataset.id] = dataset_info["good_bc_count"]


    for n1, n2, data in graph.edges(data=True):
        if data["kind"] != "breakpoint": continue

        from_ = Node(data["chrom_x"], data["new_x"], data["orientation"][0])
        to_ = Node(data["chrom_y"], data["new_y"], data["orientation"][1])

        fixed.add_edge(from_, to_, **data)

    to_delete = []
    for node, degree in dict(fixed.degree()).items():
        if degree > 1:
            edges = list(fixed.edges([node], data=True))
            edges.sort(key=lambda x: x[2]["ratio"], reverse=True)
            for edge in edges[1:]:
                to_delete.append(edge)
    fixed.remove_edges_from(to_delete)

    nodes = list(fixed.nodes())
    print "::", len(nodes)
    dist1 = -500
    dist2 = 5000

    def edge_for_node(_node):
        _edges = list(fixed.edges(_node))
        if len(_edges) < 1:
            return None
        for n1,n2 in _edges:
            if _edges[0][0] == _node:
                return _edges[0][1]
            elif _edges[0][1] == _node:
                return _edges[0][0]

    for i, n1x in enumerate(nodes[:-1]):
        for n2x in nodes[i+1:]:
            if _facing(n1x.chrom, n1x.position, n1x.orientation,
                       n2x.chrom, n2x.position, n2x.orientation) == "facing":

                if abs(n1x.position - n2x.position) < dist2:
                    print "too close!", n1x.chrom, n1x.position, n2x.position
                    continue

                n1y = edge_for_node(n1x)
                n2y = edge_for_node(n2x)

                if n1y is None or n2y is None:
                    continue

                quantification = compare_breakpoint_pair_pairs(
                    options, n1x, n1y, n2x, n2y, good_bc_counts_by_dataset, dist1, dist2)

                if quantification is None:
                    continue

                common_counts, total_counts, p = quantification
                ratio = common_counts/float(total_counts)
                if ratio < 0.01 or common_counts < 20 or p > 1e-4:
                    continue
                fixed.add_edge(n1x, n2x,
                    kind="facing",
                    shared=common_counts,
                    total=total_counts,
                    ratio=ratio,
                    p=p)

    return fixed


def compare_breakpoint_pair_pairs(options, n1x, n1y, n2x, n2y, 
                                  good_bc_counts_by_dataset, dist1, dist2):
    comparisons = []
    for sample, dataset in options.iter_10xdatasets():
        cur_result = _compare_for_dataset(
            options, sample, dataset, n1x, n1y, n2x, n2y, 
            good_bc_counts_by_dataset[dataset.id], dist1, dist2)

        if cur_result is not None:
            comparisons.append(cur_result)

    if len(comparisons) == 0:
        return None

    best = min(comparisons, key=lambda x: x[2])
    return best

        
def _compare_for_dataset(options, sample, dataset, n1x, n1y, n2x, n2y, good_bc_count, dist1, dist2):
    barcodes1 = _barcodes_for_breakpoint(options, sample, dataset, n1x, n1y, dist1, dist2)
    barcodes2 = _barcodes_for_breakpoint(options, sample, dataset, n2x, n2y, dist1, dist2)

    if barcodes1 is None or barcodes2 is None:
        return None
    total_barcodes = barcodes1.union(barcodes2)
    common_barcodes = barcodes1.intersection(barcodes2)

    contingency_table = numpy.array([[len(common_barcodes), len(barcodes1-barcodes2)],
                                     [len(barcodes2-barcodes1), good_bc_count-len(total_barcodes)]])
    p = scipy.stats.fisher_exact(contingency_table, alternative="greater")[1]

    return len(common_barcodes), len(total_barcodes), p

def _barcodes_for_breakpoint(options, sample, dataset, nodex, nodey, dist1, dist2):
    fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
        options, sample, dataset, 
        nodex.chrom, nodex.position, nodey.chrom, nodey.position, 
        nodex.orientation+nodey.orientation, dist1, dist2)

    bcx = set(fragsx["bc"])
    bcy = set(fragsy["bc"])

    common_barcodes = bcx.intersection(bcy)
    if len(common_barcodes) < 1:
        return None
    
    return common_barcodes


def visualize_graphs(outdir, graphs, evidence, file_label=""):
    try:
        utilities.ensure_dir(outdir)

        supported = 0
        missing = 0
        breakpoints = 0

        for i, graph in enumerate(graphs):
            print "visualize", i, graph
            graph = graph.copy()

            for n1,n2,data in graph.edges(data=True):
                data["label"] = "{}/{}={:.2g};{:.2g}".format(
                    int(data["shared"]),
                    int(data["total"]),
                    data["shared"]/float(data["total"]),
                    data["p"])
                if data["kind"]=="facing":
                    data["label"] = "[{}]".format(data["label"])
                    data["style"] = "dashed"
                elif data["kind"] == "breakpoint":
                    breakpoints += 1
                elif data["kind"] == "weak":
                    data["fontsize"] = 11
                    data["color"] = "gray"

                if "assembled" in data:
                    data["color"] = "orange"

            for node in graph.nodes():
                graph.node[node]["label"] = get_node_label(node)

            dot = networkx.nx_agraph.to_agraph(graph)
            if len(dot.edges())> 1000:
                print("  skipping")
                continue
            dot.draw("{}/temp{}{}.pdf".format(outdir, file_label, i), prog="dot")

        print("Supported:", supported, "Missing:", missing, "Total breakpoints:", breakpoints)
    except:
        pass


def breakend_from_label(node):
    if node.startswith("]"):
        return (node[1:].split(":")[0], int(node.split(":")[1].replace(",","")), "+")
    elif node.endswith("["):
        return (node.split(":")[0], int(node[:-1].split(":")[1].replace(",","")), "-")
    
def visualize_frags(outdir, graphs, options):
    from rpy2.robjects import r

    utilities.ensure_dir(outdir)

    for i, graph in enumerate(graphs):
        r.pdf(os.path.join(outdir, "fragments.cluster_{}.pdf".format(i)))

        for component in networkx.connected_components(graph):
            subgraph = graph.subgraph(component)
            
            ends = [node for node,degree in subgraph.degree_iter() if degree==1]
            breakends = [node for node in list(networkx.shortest_simple_paths(subgraph, ends[0], ends[1]))[0]]
            # breakends = [breakend_from_label(node) for node in breakends]
            breakends = breakends[:-1:2] + breakends[-1:]
            # print ")"*100, breakends

            for sample, dataset in sorted(options.iter_10xdatasets()):
                plot_frags(breakends, options, sample, dataset)
        # plot_frags(breakpoints, options, sample, dataset)
        r["dev.off"]()

def visualize_frag_cluster(breakpoints, options):
    """
    breakpoints are (chromx, x, chromy, y, orientation)
    """

    graph = networkx.Graph()
    for breakpoint in breakpoints:
        chromx, x, chromy, y, orientation = breakpoint
        graph.add_edge((chromx, x, orientation[0]), (chromy, y, orientation[1]))

    ends = [node for node,degree in graph.degree_iter() if degree==1]
    breakends = [node for node in list(networkx.shortest_simple_paths(graph, ends[0], ends[1]))[0]]
    # breakends = [breakend_from_label(node) for node in breakends]
    breakends = breakends[:-1:2] + breakends[-1:]
    # print ")"*100, breakends

    for sample, dataset in sorted(options.iter_10xdatasets()):
        plot_frags(breakends, options, sample, dataset)






scale = 1e6
extend = 100000

def get_frags_for_breakends(breakends, islast, options, sample, dataset):
    if len(breakends) == 2:
        assert breakends[0][0] == breakends[1][0]
        assert breakends[0][2] != breakends[1][2]

    dist1 = -500
    dist2 = 5000
    
    chrom = breakends[0][0]
    start = min(b[1] for b in breakends)
    end = max(b[1] for b in breakends)
    
    window_start = max(1, start - extend)
    window_end = end + extend
    
    #cur_frags = filter_fragments.get_frags_with_phasing(
    cur_frags = call_readclouds.load_fragments(
                    options, sample, dataset, 
                    chrom, window_start, window_end,
                    min_reads_per_frag=0)
    cur_frags["supporting"] = False
    cur_frags["event"] = str(breakends)
    
#     for breakend in breakends:
    if len(breakends) == 1:
        position = breakends[0][1]
        orientation = breakends[0][2]
        if orientation == "+":
            condition = (position-cur_frags["end_pos"]).between(dist1, dist2)
        elif orientation == "-":
            condition = (cur_frags["start_pos"]-position).between(dist1, dist2)
        cur_frags["supporting"] |= condition
    else:
        between = (cur_frags["start_pos"] > (start + dist1)) & \
                  (cur_frags["end_pos"] < (end - dist1))
        
        cur_frags["supporting"] = between
    
    # print "  ", len(cur_frags), cur_frags["supporting"].sum()

#     if breakends[0][2] == "+":
    if is_reversed(breakends, islast):
        cur_frags = cur_frags.sort_values("end_pos", ascending=False)
    else:
        cur_frags = cur_frags.sort_values("start_pos")

    return cur_frags

def get_frags(event, options, sample, dataset):#, reverse_order=False):
    breakpoints_to_frags = collections.OrderedDict()
    bc_counts = collections.defaultdict(set)
    
    # if len(event) == 2:
    #     ends = sorted(event.ends)
    # if reverse_order:
    #     ends = ends[::-1]
    ends = event.ends
    breakend_sets = [(ends[0],)]
    breakend_sets.extend(event.facing)
    breakend_sets.append((ends[1],))
    
    for i, breakends in enumerate(breakend_sets):
        islast = (i == len(breakend_sets)-1)

        cur_frags = get_frags_for_breakends(breakends, islast, options, sample, dataset)
        breakpoints_to_frags[breakends] = cur_frags
        for row in cur_frags.itertuples():
            if row.supporting:
                bc_counts[row.bc].add(row.event)
    
    good_bcs = set(bc for bc in bc_counts if len(bc_counts[bc])>1)
    
    bcs_to_rows = {}
    cur_row = 0
            
    for breakends in breakend_sets:
        cur_frags = breakpoints_to_frags[breakends]
        cur_frags = cur_frags.loc[cur_frags.bc.isin(good_bcs)]
        breakpoints_to_frags[breakends] = cur_frags.copy()
        
        for row in cur_frags.itertuples():
            if row.supporting and (row.bc not in bcs_to_rows):
                bcs_to_rows[row.bc] = cur_row
                cur_row += 1
            
    return breakpoints_to_frags, bcs_to_rows

def is_reversed(breakpoint, islast):
    isreversed = False
    if len(breakpoint)==2:
        if breakpoint[0][1] > breakpoint[1][1]:
            isreversed = True
    elif breakpoint[0][2]=="-":
        isreversed = True

    if islast:
        isreversed = not isreversed
    return isreversed

def do_plotting(breakpoints_to_frags, bcs_to_rows, options, sample, dataset,
                coverage_normalizations=None, normalize_coverage_to=None):

    # TODO: biorpy required
    # from biorpy import r

    from rpy2 import robjects as ro
    oldpar = ro.r.par(mfrow=ro.IntVector([2, len(breakpoints_to_frags)]),
                      mar=ro.FloatVector([5,4,4,0]))
    
    windows = []

    for i, (breakpoint, frags) in enumerate(breakpoints_to_frags.items()):
        if i > 0:
            ro.r.par(mar=ro.FloatVector([5,2,4,0]))
        chrom, position, orientation = breakpoint[0]
        islast = (i == len(breakpoints_to_frags)-1)
        isreversed = is_reversed(breakpoint, islast)
        
        if len(breakpoint) > 1:
            position = int((breakpoint[0][1] + breakpoint[1][1]) / 2)
        
        breakpoint_lines = numpy.array([b[1] for b in breakpoint])
        windows.append((chrom, position, isreversed, breakpoint_lines))

        xlim = numpy.array([position-extend, position+extend])/scale
        if isreversed:
            xlim = xlim[::-1]
        
        ylim = [0,0]
        if len(bcs_to_rows) > 0:
            ylim = [0,max(bcs_to_rows.values())]

        ro.r.plot(ro.FloatVector([0]),
                  type="n", bty="n",
                  xlim=ro.FloatVector(xlim),
                  ylim=ro.FloatVector(ylim),
                  xlab="{} (mb)".format(chrom), ylab="Barcode")
        
        ro.r.abline(v=ro.FloatVector(breakpoint_lines/scale), lty=2, col="gray")

        if len(bcs_to_rows) == 0:
            continue

        ypos = numpy.array([bcs_to_rows[bc] for bc in frags["bc"]])
        # colors = numpy.where(cur_frags["supporting"], "black", "gray")
        support = []
        for row in frags.itertuples():
            if not row.supporting:
                support.append("background")
            elif numpy.isnan(row.hap):
                support.append("nohap")
            elif row.hap == 1:
                support.append("hap0")
            elif row.hap == 2 or row.hap == 0:
                support.append("hap1")
            else:
                raise Exception("can't interpret hap: {}".format(row))

        frags = frags.copy()
        frags["support"] = support
        frags["ypos"] = ypos
        
        support_types = ["background", "nohap", "hap0", "hap1"]
        colors = ["lightgray", "black"]

        major_color = "orange" #"#F2694C"
        minor_color = "cyan" #"#6AADD5"
        hap_counts = frags["hap"].value_counts()
        
        if hap_counts.get(0, 0) > hap_counts.get(1, 0):
            colors.extend([major_color, minor_color])
        else:
            colors.extend([minor_color, major_color])
        
        for support_type, color in zip(support_types, colors):
            if len(frags) == 0:
                continue
            cur_frags = frags.loc[frags["support"]==support_type]

            ro.r.segments(
                ro.FloatVector(cur_frags["start_pos"].values/scale), ro.FloatVector(cur_frags["ypos"]),
                ro.FloatVector(cur_frags["end_pos"].values/scale), ro.FloatVector(cur_frags["ypos"]),
                col=color)

    plot_coverage(windows, options, sample, dataset, coverage_normalizations, normalize_coverage_to)

    ro.r.par(oldpar)
    
def plot_frags(cluster, options, sample, dataset, coverage_normalizations=None, normalize_coverage_to=None, reverse_order=False):
    if len(cluster) > 1:
        assert cluster.kind.str.contains("facing").any()

    complex_event = ComplexEvent(cluster, reverse_order)

    breakpoints_to_frags, bcs_to_rows = get_frags(
        complex_event, options, sample, dataset)#, reverse_order=reverse_order)

    for bnd, frags in breakpoints_to_frags.items():
        frags = frags.loc[frags.bc.isin(bcs_to_rows)]
        print bnd, len(frags[["chrom", "start_pos", "end_pos", "supporting"]])
    print len(bcs_to_rows)
    
    do_plotting(breakpoints_to_frags, bcs_to_rows, options, sample, dataset,
                coverage_normalizations, normalize_coverage_to)

class ComplexEvent(object):
    def __init__(self, cluster, reverse=False):
        self.graph = networkx.Graph()
        self.reverse = reverse
        for row in cluster.itertuples():
            nodex = (row.chromx, row.x, row.orientation[0])
            nodey = (row.chromy, row.y, row.orientation[1])

            self.graph.add_edge(nodex, nodey, **{"kind":row.kind, "id":row.Index})
            
    @property
    def ends(self):
        ends = sorted([node for node,degree in dict(self.graph.degree()).items() if degree==1])
        if len(ends) > 2:
            print "*"*100, ends
        if self.reverse:
            ends = ends[::-1]
        return ends
    @property
    def nodes_in_order(self):
        ends = self.ends
        return list(networkx.shortest_simple_paths(self.graph, ends[0], ends[1]))[0]
    @property
    def edges(self):
        nodes_in_order = self.nodes_in_order
        return [nodes_in_order[i:i+2] for i in range(len(nodes_in_order)-1)]
    @property
    def breakpoints(self):
        return [tuple(e) for e in self.edges if self.graph[e[0]][e[1]]["kind"]=="breakpoint"]
    @property
    def facing(self):
        return [tuple(e) for e in self.edges if self.graph[e[0]][e[1]]["kind"]=="facing"]

def plot_coverage(windows, options, sample, dataset, coverage_normalizations, normalize_coverage_to):
    """
    coverage_normalizations - dictionary of constants to use to normalize coverage profiles (optional)
    normalize_coverage_to - path of normal bam (optional)
    """
    # from biorpy import r
    from rpy2 import robjects as ro

    if not isinstance(dataset, svdatasets.ShortFragDataset):
        for d in sample.datasets:
            if isinstance(d, svdatasets.ShortFragDataset):
                dataset = d
    cur_bam = pysam.AlignmentFile(dataset.bam)
    normal_bam = None
    if normalize_coverage_to is not None:
        normal_bam = pysam.AlignmentFile(normalize_coverage_to)

    step = 1000
    scale = 1e6

    event_coverages = []
    xs = []
    ymax = 0

    for window in windows:
        chrom, position, isreversed, breakpoint_lines = window
        start, end = position-extend, position+extend
        start = max(0, start)

        cur_coverages = []
        cur_x = []

        for i in range(start, end, step):
            cur_count = cur_bam.count(chrom, i, i+step)
            if coverage_normalizations is not None:
                cur_count /= coverage_normalizations[sample.name]

            if normal_bam is not None:
                normal_count = normal_bam.count(chrom, i, i+step)
                if coverage_normalizations is not None:
                    normal_count /= coverage_normalizations["normal"]

                cur_count /= float(normal_count)

                # make in terms of copies (diploid, etc)
                cur_count *= 2.0

            cur_coverages.append(cur_count)
            cur_x.append(i)

        cur_coverages = pandas.Series(cur_coverages)
        cur_coverages = cur_coverages.rolling(5, center=True).median()
        ymax = max(ymax, cur_coverages.max())

        event_coverages.append(cur_coverages)
        xs.append(cur_x)

    ylim = [0, ymax]

    for i, window in enumerate(windows):
        if i == 0:
            # TODO: fix 
            ro.r.par(mar=ro.FloatVector([5,4,4,0]))
        else:
            ro.r.par(mar=ro.FloatVector([5,2,4,0]))
        chrom, position, isreversed, breakpoint_lines = window
        start, end = position-extend, position+extend


        xlim = [start, end]
        if isreversed:
            xlim = xlim[::-1]

        ro.r.plot(
            ro.FloatVector(numpy.array(xs[i])/scale),
            ro.FloatVector(event_coverages[i].values),
            xlim=ro.FloatVector(numpy.array(xlim)/scale),
            ylim=ro.FloatVector(numpy.array(ylim)),
            type="n",
            bty="n", xlab="{} (mb)".format(chrom), ylab="Copy number", main="")

        ro.r.abline(v=ro.FloatVector(breakpoint_lines/scale), lty=2, col="gray")
        ro.r.abline(h=ro.FloatVector(numpy.arange(0,10)), lty=2, col="gray")

        ro.r.lines(ro.FloatVector(numpy.array(xs[i])/scale), ro.FloatVector(event_coverages[i].values), lwd=2)
