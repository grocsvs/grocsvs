import collections
import networkx
import numpy
import os
import pandas
import pysam
import re
import scipy.stats

from grocsvs import step
from grocsvs import structuralvariants
from grocsvs import utilities

from grocsvs.stages import assembly
from grocsvs.stages import cluster_svs
from grocsvs.stages import call_readclouds

BAM_CMATCH     = 0
BAM_CINS       = 1
BAM_CDEL       = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5




class WalkAssembliesStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        yield WalkAssembliesStep(options)


    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        walk_assemblies = "walk_assemblies.tsv"

        graphs = "walk_assemblies.graphs"

        paths = {
            "walk_assemblies": os.path.join(directory, walk_assemblies),
            "graphs": os.path.join(directory, graphs)
        }

        return paths

    def run(self):
        edges_path = cluster_svs.ClusterSVsStep(self.options).outpaths(final=True)["edges"]

        clusters = pandas.read_table(edges_path)

        assembled = []

        utilities.ensure_dir(self.outpaths(final=False)["graphs"])
        for cluster_number, cluster in clusters.groupby("cluster"):
            self.logger.log(cluster_number)
            try:
                cur_assembled = self.walk(cluster_number, cluster)
                assembled.append(cur_assembled)
            except IOError:
                print "not found", cluster_number

        # TODO: deal with empty list
        # TODO: normalize coordinates according to reference.compare_chroms()
        assembled = pandas.concat(assembled, ignore_index=True)
        assembled["x"] = assembled["x"].astype(int)
        assembled["y"] = assembled["y"].astype(int)

        print self.options.reference.chroms
        print assembled["chromx"].unique()
        print assembled["chromy"].unique()
        
        outpath = self.outpaths(final=False)["walk_assemblies"]
        assembled.to_csv(outpath, sep="\t", index=False)


    def walk(self, event_name, cluster):        
        assembly_directory = assembly.AssemblyStep(self.options, event_name)\
                                     .outpaths(final=True)["assembly_dir"]

        bam, contigs = self.get_contigs(assembly_directory)
    
        cluster = cluster.loc[cluster["kind"]=="breakpoint"]
        cluster = cluster.sort_values("p")

        chains = set()

        for i, event in cluster.iterrows():
            if event["x"] < 0 or event["y"] < 0:
                continue


            # print "starting from:\n", event
            chain = get_chain(bam, event["chromx"], event["x"], contigs)
            # print "chain:", chain
            chains.update(chain)

            # print "starting from:\n", event
            chain = get_chain(bam, event["chromy"], event["y"], contigs)
            # print "chain:", chain

            chains.update(chain)

        # self.analyze_chains(chains, event_name)

        cur_edges = pandas.DataFrame(
            list(chains), columns=["chromx", "x", "orientationx", "chromy", "y", "orientationy", "contig"])

        cur_edges = cur_edges.loc[(cur_edges["chromx"].isin(self.options.reference.chroms)) & 
                                  (cur_edges["chromy"].isin(self.options.reference.chroms))]
        
        return cur_edges

    def _barcodes_for_breakpoint(self, chromx, x, orientationx, 
                                       chromy, y, orientationy, dist1, dist2):
        # TODO: refactor to re-use same version as cluster_svs
        fragsx, fragsy, merged = structuralvariants.get_supporting_fragments_new(
            self.options, self.sample, self.dataset, 
            chromx, x, chromy, y,
            orientationx+orientationy, dist1, dist2)

        bcx = set(fragsx["bc"])
        bcy = set(fragsy["bc"])

        common_barcodes = bcx.intersection(bcy)
        if len(common_barcodes) < 1:
            return None
        
        return common_barcodes

    def _compare_breakpoint_pair_pairs(self, breakpoint1, breakpoint2, good_bc_count, dist1, dist2):
        chrom1x, pos1x, orientation1x, chrom1y, pos1y, orientation1y, _ = breakpoint1
        chrom2x, pos2x, orientation2x, chrom2y, pos2y, orientation2y, _ = breakpoint2

        # TODO: refactor to re-use same version as cluster_svs
        barcodes1 = self._barcodes_for_breakpoint(
            chrom1x, pos1x, orientation1x, chrom1y, pos1y, orientation1y, dist1, dist2)
        barcodes2 = self._barcodes_for_breakpoint(
            chrom2x, pos2x, orientation2x, chrom2y, pos2y, orientation2y, dist1, dist2)

        if barcodes1 is None or barcodes2 is None:
            return None
        total_barcodes = barcodes1.union(barcodes2)
        common_barcodes = barcodes1.intersection(barcodes2)

        contingency_table = numpy.array([[len(common_barcodes), len(barcodes1-barcodes2)],
                                         [len(barcodes2-barcodes1), good_bc_count-len(total_barcodes)]])
        p = scipy.stats.fisher_exact(contingency_table, alternative="greater")[1]

        return len(common_barcodes), len(total_barcodes), p



    def get_contigs(self, assembly_directory):
        bam_path = os.path.join(assembly_directory, "contigs.sorted.bam")

        names_to_reads = collections.defaultdict(list)

        self.logger.log(bam_path)
        bam = pysam.AlignmentFile(bam_path)
        for read in bam.fetch():
            if read.reference_length > 45:
                names_to_reads[read.query_name].append(read)

        return bam, names_to_reads


    # def analyze_chains(self, chains, event_name):
    #     dist1 = -500
    #     dist2 =  5000

    #     good_bc_count = utilities.get_good_bc_count(self)
        
    #     full_graph = networkx.Graph()
    #     barcode_supported_graph = networkx.Graph()


    #     for chain in chains:
    #         for breakpoint in chain:
    #             chromx, x, orientationx, chromy, y, orientationy, _ = breakpoint

    #             nodex = get_node_label(chromx, x, orientationx)
    #             nodey = get_node_label(chromy, y, orientationy)

    #             if not full_graph.has_edge(nodex, nodey):
    #                 quant = quantify_breakpoint(
    #                     chromx, x,
    #                     chromy, y,
    #                     orientationx+orientationy,
    #                     self.options, self.sample, self.dataset,
    #                     good_bc_count, dist1, dist2)
    #                 if quant is None: continue
    #                 cur_common_counts, cur_total_counts, p = quant

    #                 print breakpoint, cur_common_counts, cur_total_counts, p
    #                 ratio = cur_common_counts/float(cur_total_counts)
    #                 label = "{}/{}={:.2g};{:.2g}".format(
    #                     int(cur_common_counts),
    #                     int(cur_total_counts),
    #                     ratio,
    #                     p)
                    
    #                 if ratio > 0.08 and p < 1e-4 and cur_total_counts > 10:
    #                     barcode_supported_graph.add_edge(nodex, nodey, label=label)
    #                 full_graph.add_edge(nodex, nodey, label=label)
    #         for j, breakpoint1 in enumerate(chain[:-1]):
    #             breakpoint2 = chain[j+1]

    #             quant = self._compare_breakpoint_pair_pairs(
    #                 breakpoint1, breakpoint2, good_bc_count, dist1, dist2)
    #             if quant is None: continue
    #             common_counts, total_counts, p = quant

    #             ratio = common_counts/float(total_counts)

    #             node1y = get_node_label(*breakpoint1[3:6])
    #             node2x = get_node_label(*breakpoint2[:3])

    #             label = "[{}/{}={:.2g};{:.2g}]".format(
    #                 int(common_counts),
    #                 int(total_counts),
    #                 ratio,
    #                 p)

    #             if ratio > 0.08 and p < 1e-4 and total_counts > 10:
    #                 barcode_supported_graph.add_edge(node1y, node2x, label=label, style="dashed")
    #             full_graph.add_edge(node1y, node2x, label=label, style="dashed")


    #     print ":: ALL:", full_graph.edges(data=True)
    #     print ":: SUPPORTED:", barcode_supported_graph.edges(data=True)

    #     outdir = self.outpaths(final=False)["graphs"]
    #     barcode_supported_dot = networkx.nx_agraph.to_agraph(barcode_supported_graph)
    #     barcode_supported_dot.draw("{}/barcode_supported.{}.pdf".format(outdir, event_name), prog="dot")

    #     full_dot = networkx.nx_agraph.to_agraph(full_graph)
    #     full_dot.draw("{}/full_graph.{}.pdf".format(outdir, event_name), prog="dot")

    #     def breakend_from_label(node):
    #         if node.startswith("]"):
    #             return (node[1:].split(":")[0], int(node.split(":")[1].replace(",","")), "+")
    #         elif node.endswith("["):
    #             return (node.split(":")[0], int(node[:-1].split(":")[1].replace(",","")), "-")
            
    #     from rpy2.robjects import r

    #     r.pdf("{}/raw{}.pdf".format(outdir, event_name))
    #     for component in networkx.connected_components(barcode_supported_graph):
    #         subgraph = barcode_supported_graph.subgraph(component)
            
    #         ends = [node for node,degree in subgraph.degree_iter() if degree==1]
    #         breakends = [node for node in list(networkx.shortest_simple_paths(subgraph, ends[0], ends[1]))[0]]
    #         breakends = [breakend_from_label(node) for node in breakends]
    #         breakends = breakends[:-1:2] + breakends[-1:]
    #         # print ")"*100, breakends
    #         plot_frags(breakends, self.options, self.sample, self.dataset)

    #     r["dev.off"]()



def sort_by_ref_pos(q1, r1, q2, r2, l=None):
    if r1 < r2:
        return q1, r1, q2, r2
    else:
        return q2, r2, q1, r1 

def get_query_positions(read):
    cigar = read.cigar
    
    qstart = 0
    rstart = read.pos
    if cigar[0][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
        qstart += cigar[0][1]
        
    qend = 0
    rend = read.pos
    
    for op, length in cigar:
        if op in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CINS]:
            qend += length
        elif op == BAM_CMATCH:
            qend += length
            rend += length
        elif op in [BAM_CDEL, BAM_CREF_SKIP]:
            rend += length
    
    total_length = qend
    
    if cigar[-1][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
        qend -= cigar[-1][1]
    if read.is_reverse:
        qstart, rstart, qend, rend = total_length-qend, rend, total_length-qstart, rstart
    return qstart, rstart, qend, rend, total_length


def is_clipped(read, min_clip_length=40):
    q1, r1, q2, r2, l = get_query_positions(read)
    x = ""
    
    # print read.query_name, q1, q2, l
    if q1 > min_clip_length:
        if r2 > r1:
            x += ">"
        else:
            x += "<"
    if q2 < l-min_clip_length:
        if r2 > r1:
            x += "<"
        else:
            x += ">"
            
    return x, min(r1,r2), max(r1,r2)

def covered_positions(start, end, reads):
    length = end-start
    positions = [False]*length
    for read in reads:
        for pos in read.get_reference_positions():
            pos = pos-start
            if 0 <= pos < length:
                positions[pos] = True
    return positions

def walk_local(bam, chrom, pos, direction, offset=4000, window_size=10000):
    reads = set()
    
    if direction == "-":
        offset = -offset
        window_size = -window_size
        
    cur_start = pos
    cur_end = pos + window_size
#     cur_start, cur_end = sorted([cur_start, cur_end])
    
    while True:
        if min(cur_start, cur_end) < abs(offset)*2:
            return []

        cur_reads = list(bam.fetch(chrom, min(cur_start, cur_end), max(cur_start, cur_end)))
        if len(cur_reads) == 0 or max(read.reference_length for read in cur_reads) < 40:
            break
            
        reads.update(cur_reads)
        cur_start += offset
        cur_end += offset
    
    if len(reads) == 0:
        return []

    reads = sorted(reads, key=lambda x: x.pos)
    if direction == "+":
        start_index = 0
        end_index = len(reads)-1
        increment = 1
    elif direction == "-":
        start_index = len(reads)-1
        end_index = 0
        increment = -1
    
    positions = set()
    for read in reads:
        positions.update(read.get_reference_positions())
    positions = numpy.array(sorted(positions))
    # print len(numpy.diff(positions))
    # print numpy.where(numpy.diff(positions)>abs(offset))[0]

    # print numpy.split(positions, numpy.where(numpy.diff(positions)>abs(offset))[0]+1)
    split = numpy.split(positions, numpy.where(numpy.diff(positions)>abs(offset))[0]+1)
    # for s in split:
        # print "  ", direction, s[0], s[-1]

    if direction == "+":
        positions = split[0]
    else:
        positions = split[-1]
    # print "POSITIONS:", positions[0], positions[-1]

    filtered_reads = [read for read in reads if (read.reference_end>=positions[0] and
                                                 read.reference_start<=positions[-1])]

    if direction == "-":
        filtered_reads = sorted(filtered_reads, key=lambda x: x.reference_start)
    elif direction == "+":
        filtered_reads = sorted(filtered_reads, key=lambda x: x.reference_end)

    return filtered_reads


def are_matched(read1, orientation1, read2, orientation2):
    if reads_overlap(read1, read2):
        # TODO: should figure out if we can improve the assemblies
        # so that we can detect small events
        return False, None, None

    qstart1, rstart1, qend1, rend1 = sort_by_ref_pos(*get_query_positions(read1))
    qstart2, rstart2, qend2, rend2 = sort_by_ref_pos(*get_query_positions(read2))
    
    # print ":: {} {:,} {} {:,}||  {} {:,} {} {:,}".format(qstart1, rstart1, qend1, rend1, qstart2, rstart2, qend2, rend2)
    switched = False
    if (qstart1+qend1) > (qstart2+qend2):
        # print "switch"
        switched = True
        qstart1, rstart1, qend1, rend1, qstart2, rstart2, qend2, rend2 = qstart2, rstart2, qend2, rend2, qstart1, rstart1, qend1, rend1
        orientation1, orientation2 = orientation2, orientation1
    
    if orientation1 == "+":
        qadj1 = qend1
        radj1 = rend1
    elif orientation1 == "-":
        qadj1 = qstart1
        radj1 = rstart1
        
    if orientation2 == "+":
        qadj2 = qend2
        radj2 = rend2
    elif orientation2 == "-":
        qadj2 = qstart2
        radj2 = rstart2
        
    # print qadj1, qadj2
    
    if switched:
        radj1, radj2 = radj2, radj1
    if abs(qadj1 - qadj2) < 20:
        return True, radj1, radj2
    return False, None, None


def reads_overlap(read1, read2):
    if read1.reference_id != read2.reference_id:
        return False
    if read1.reference_end < read2.reference_start:
        return False
    if read2.reference_end < read1.reference_start:
        return False
    return True

def _get_chain(bam, chrom, pos, direction, reads_by_name, previous=None):
    # direction is "-" if we're walking from the right to the left, ie looking
    # for a "-" orientation breakpoint; and "+" if we're walking from left to right
    # print "\n** _get_chain", chrom, pos, direction
    reads = walk_local(bam, chrom, pos, direction)
    
    if len(reads) == 0:
        return []
    
    if direction == "+":
        read = reads[-1]
        clipping = is_clipped(read)
        if "<" in clipping[0]:
            end = clipping[2]
        else:
            return []
    elif direction == "-":
        read = reads[0]
        clipping = is_clipped(read)
        if ">" in clipping[0]:
            end = clipping[1]
        else:
            return []
    
    # print "last read:", read.reference_start, read.reference_end
    for other_read in reads_by_name[read.query_name]:
        if other_read == read: continue
        is_match = False
        
        for other_orientation in "+-":
            match = are_matched(read, direction, other_read, other_orientation)
            if match[0]:
                is_match = True
                other_end = match[2]
                break
        
        if is_match:
            next_chrom = other_read.reference_name
            next_start = (other_read.pos+other_read.reference_end)/2
            next_orientation = {"+":"-", "-":"+"}[other_orientation]

            # if reads_overlap(read, other_read):
            #     return []
 
            bad_chrom = re.search(r"gl|hap|un|_", next_chrom, re.IGNORECASE)
            # print ">>>>>>>>>>>", chrom, end, direction, next_chrom, next_start, next_orientation
            if bad_chrom is not None:
                return []

            if previous is None:
                previous = set()

            if (next_chrom, next_start, next_orientation) in previous:
                return []
            previous.add((next_chrom, next_start, next_orientation))

            next_chain = _get_chain(bam, next_chrom, next_start, next_orientation, 
                reads_by_name, previous)

            chain = [(chrom, end, direction, next_chrom, other_end, other_orientation, read.query_name)]

            if next_chain is not None:
                chain.extend(next_chain)

            return chain
        
    return []

def get_chain(bam, chrom, pos, reads_by_name):
    chain = []
    # print "orientation:", "-"*30

    left_chain = _get_chain(bam, chrom, pos+5000, "-", reads_by_name)
    for link in left_chain[::-1]:
        chain.append(tuple(list(link[3:6])+list(link[:3])+[link[6]]))
    
    # print "orientation:", "+"*30
    right_chain = _get_chain(bam, chrom, max(0, pos-5000), "+", reads_by_name)
    if right_chain is not None:
        chain.extend(right_chain)
    return chain


# def visualize_chain(chain):
#     graph = networkx.Graph()
#     prev = None

#     for link in chain:
#         chromx, posx, orientationx = link[:3]
#         chromy, posy, orientationy = link[3:6]

#         graph.add_edge(
#             get_node_label(chromx, posx, orientationx),
#             get_node_label(chromy, posy, orientationy)
#             )


def get_node_label(chrom, position, orientation):
    if orientation == "+":
         return "]{}:{:,}".format(chrom, int(position))
    else:
        return "{}:{:,}[".format(chrom, int(position))

