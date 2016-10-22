import collections
import numpy
import pandas
from scipy import stats

from grocsvs.stages import call_readclouds


class StructuralVariant(object):
    """
    Represent an arbitrary structural variant as a pair of breakpoints
    """
    def __init__(self, chromx, chromy, x, y, orientation):
        self.chromx = chromx
        self.chromy = chromy
        self.x = x
        self.y = y
        self.orientation = orientation

        self.support = {}

    def __repr__(self):
        return "{}:{}::{}:{}({})".format(
            self.chromx, self.x,
            self.chromy, self.y,
            self.support)

    def get_sv_info(self, fragsx, fragsy, ext_dist=100000, winsize=500):
        # TODO: learn this parameter
        startx, endx = self.x-ext_dist, self.x+ext_dist
        starty, endy = self.y-ext_dist, self.y+ext_dist

        fragsx, fragsy, merged_frags = merge_fragments(
            startx, endx,
            starty, endy, 
            fragsx, fragsy)


        merged_frags = merged_frags.loc[(merged_frags["chrom_x"]!=merged_frags["chrom_y"]) | 
                                        (merged_frags["start_pos_x"]!=merged_frags["start_pos_y"])]

        mat = overlap_matrix(
            merged_frags, 
            startx, endx, 
            starty, endy, 
            winsize=winsize)

        breakpoint = find_breakpoint(
            mat, startx, starty, endx, endy)

        return fragsx, fragsy, merged_frags, mat, breakpoint

    def calculate_characteristics(self, merged_frags, filtered_merged_frags, frag_length_genotyper):
        supporting_frags = filtered_merged_frags.shape[0]
        
        query_x = """~((start_pos_x > {}) | ({} > end_pos_x))"""
        query_x = query_x.format(self.x+12000, self.x-12000)
        query_y = """~((start_pos_y > {}) | ({} > end_pos_y))"""
        query_y = query_y.format(self.y+12000, self.y-12000)

        overlapping_frags = merged_frags.query("({}) & ({})".format(query_x, query_y))
        # print sv
        # print merged_frags
        # print overlapping_frags

        total_frags = overlapping_frags.shape[0]
        if total_frags > 0:
            supporting_fraction = supporting_frags/float(total_frags)
        else:
            supporting_fraction = 0
        merged_lengths = filtered_merged_frags["merged_len"]
        if len(merged_lengths) > 0:
            frag_length_dist_does_match = frag_length_genotyper.genotype(
                merged_lengths, lenient=False)
        else:
            frag_length_dist_does_match = False

        characteristics = {
            "supporting": supporting_frags,
            "overlapping": overlapping_frags.shape[0],
            "supporting_fraction": supporting_fraction,
            "frag_length_distributions_match": frag_length_dist_does_match
        }
        
        return characteristics

            
        
def merge_svs(svs):
    if len(svs) == 1:
        return svs[0]
    raise Exception("need to merge the breakpoints here (pick the best? mean?)")

    sv = svs[0]

    for othersv in svs[1:]:
        for key in othersv.support:
            if not key in sv.support:
                sv.support[key] = {}
            sv.support[key].update(othersv.support[key])

    return sv


def offdiag(coords, offdiag_dist):
    filtered_coords = []
    for x, y in zip(*coords):
        if abs(x-y) > offdiag_dist:
            filtered_coords.append((x,y))
    return filtered_coords

def closest(a, b):
    min_dist = None
    aidx = None
    bidx = None
    # equivalent to zip()'ping and taking the first two items
    bb = numpy.transpose(b)[:2].astype(int)
    for i, aa in enumerate(a):
        dist = numpy.sqrt((bb[0]-aa[0])**2 + (bb[1]-aa[1])**2)
        cur_min_dist = dist.min()
        if cur_min_dist < min_dist or min_dist is None:
            min_dist = cur_min_dist
            aidx = i
            bidx = numpy.where(dist==cur_min_dist)[0][0]
                       
    return aidx, bidx, min_dist
                    
def do_free_clustering(points, max_dist):
    points = points[:]
    clusters = []
    
    while len(points) > 1:
        p = points.pop()
        curcluster = [p]
        
        while True:
            i, j, d = closest(curcluster, points)
            if d is not None and d < max_dist:
                curcluster.append(points.pop(j))
                if len(points) == 0:
                    break
            else:
                break
        clusters.append(curcluster)

    if len(points) > 0:
        clusters.append([points.pop()])
     
    return clusters

def do_grid_clustering(points, max_dist):
    clusters = numpy.zeros((max(p[0] for p in points)+1, max(p[1] for p in points)+1))
    assignments = collections.defaultdict(list)
    n = 0
    collisions = 0
    
    for i,j in points:
        neighbors = numpy.unique(clusters[max(0, i-max_dist):i+max_dist+1,
                                          max(0, j-max_dist):j+max_dist+1])
        
        neighbors = sorted(neighbors[neighbors>0])
        
        if len(neighbors) == 0:
            n += 1
            cur_cluster = n
        elif len(neighbors) == 1:
            cur_cluster = neighbors[0]
        else:
            cur_cluster = neighbors[0]
            for neighbor in neighbors[1:]:
                for idx in assignments[neighbor]:
                    clusters[idx] = cur_cluster
                assignments[cur_cluster].extend(assignments.pop(neighbor))
            collisions += 1
            
        assignments[cur_cluster].append((i,j))
        clusters[i,j] = cur_cluster

    # y, x = zip(*points)
    # cluster = clusters[y,x]

    # return zip(x, y, cluster)
    return assignments.values()
        

def get_fragments_one_side(options, sample, dataset, chrom, position,
                             orientation, dist1, dist2, strict=False,
                             min_reads_per_frag=1, with_phasing=False):
    """
    strict : if strict is True, only include those frags starting/ending
        within the search region; otherwise, include all frags overlapping
        the search region
    """
    if orientation == "+":
        start = position-dist2
        end = position-dist1
    elif orientation == "-":
        start = position+dist1
        end = position+dist2

    start = max(0, start)
    
    if with_phasing:
        frags = call_readclouds.get_frags_with_phasing(
            options, sample, dataset, chrom, start, end,
            min_reads_per_frag=min_reads_per_frag)
    else:
        frags = call_readclouds.load_fragments(
            options, sample, dataset, chrom, start, end,
            min_reads_per_frag=min_reads_per_frag)


    # TODO: figure out whether we want to be using strict or not...
    if strict:
        if orientation == "+":
            condition = (position-frags["end_pos"]).between(dist1, dist2)
        elif orientation == "-":
            condition = (frags["start_pos"]-position).between(dist1, dist2)
        frags = frags.loc[condition]

    return frags

def get_supporting_fragments_new(options, sample, dataset, chromx, x, chromy, y, 
                             orientation, dist1, dist2, min_reads_per_frag=1, with_phasing=False):
    
    fragsx = get_fragments_one_side(options, sample, dataset, 
                                    chromx, x, orientation[0], 
                                    dist1, dist2, strict=True, 
                                    min_reads_per_frag=min_reads_per_frag,
                                    with_phasing=with_phasing)
    fragsy = get_fragments_one_side(options, sample, dataset, 
                                    chromy, y, orientation[1], 
                                    dist1, dist2, strict=True,
                                    min_reads_per_frag=min_reads_per_frag,
                                    with_phasing=with_phasing)


    merged = pandas.merge(fragsx, fragsy, on="bc", suffixes=["_x", "_y"], how="outer")

    # filter out fragments that span from region x to region y
    merged = merged.loc[(merged["chrom_x"]!=merged["chrom_y"]) | 
                        (merged["start_pos_x"]!=merged["start_pos_y"])]

    fragsx = fragsx.loc[fragsx["bc"].isin(merged["bc"])]
    fragsy = fragsy.loc[fragsy["bc"].isin(merged["bc"])]

    # converts to inner join
    # print "=" * 100
    # print "BEFORE"
    # print merged    
    merged = merged.dropna(subset=["chrom_x", "chrom_y"])
    # print "AFTER"
    # print merged
    # print "=" * 100

    return fragsx, fragsy, merged




def get_supporting_fragments(options, sample, dataset, chromx, x, chromy, y, 
                             orientation, dist1, dist2):
    """
    returns frags overlapping and sharing barcodes across x and y, as well as
    those starting/ending near the breakpoints (using dist1/dist2 for both)

    dist1 should proportional to the variance in observed fragment ends 
        (relative to the actual ends)
    dist2 should be proportional to the uncertainty in breakpoint position
    """

    fragsx = get_fragments_one_side(
        options, sample, dataset, chromx, x, orientation[0], dist1, dist2)
    fragsy = get_fragments_one_side(
        options, sample, dataset, chromy, y, orientation[1], dist1, dist2)

    merged = pandas.merge(fragsx, fragsy, on="bc", suffixes=["_x", "_y"])

    # TODO: this filtering step should be unnecessary if the strict option
    # is used above to get_fragments_one_side()
    filtered_frags = filter_merged(merged, x, y, orientation, dist1, dist2)

    return fragsx, fragsy, merged, filtered_frags


def merge_fragments(startx, endx, starty, endy, frags_chromx, frags_chromy):
    """
    take the fragments from each region and join them on the barcode
    """
    query = """~((start_pos > {}) | ({} > end_pos))"""
    qx = query.format(endx, startx)
    qy = query.format(endy, starty)

    fragsx = frags_chromx.query(qx)
    fragsy = frags_chromy.query(qy)

    merged = pandas.merge(fragsx, fragsy, on="bc", suffixes=["_x", "_y"])

    return fragsx, fragsy, merged

def filter_merged(merged, x, y, orientation, dist1, dist2):
    """
    filter the merged fragments table for only those fragments with endpoints 
    close to both breakpoints; assumes that the fragments should extend from 
    each breakpoint only in a single direction
    """

    if orientation == "++":
        condition = (((x-merged["end_pos_x"]).between(dist1, dist2)) &
                     ((y-merged["end_pos_y"]).between(dist1, dist2)))
    elif orientation == "+-":
            condition = (((x-merged["end_pos_x"]).between(dist1, dist2)) &
                         ((merged["start_pos_y"]-y).between(dist1, dist2)))
    elif orientation == "-+":
        condition = (((merged["start_pos_x"]-x).between(dist1, dist2)) &
                     ((y-merged["end_pos_y"]).between(dist1, dist2)))
    elif orientation == "--":
            condition = (((merged["start_pos_x"]-x).between(dist1, dist2)) &
                         ((merged["start_pos_y"]-y).between(dist1, dist2)))
    else:
        raise Exception("unknown orientation: '{}'".format(orientation))

    filtered = merged.loc[condition].copy()

    return filtered

def overlap_matrix(merged, startx, endx, starty, endy, winsize=250):
    """
    get a 2d histogram of fragment overlaps as a matrix

    "merged" is the result of merge_fragments()
    """
    # quick fix for the diagonal
    # merged = merged.loc[(merged["start_pos_x"]!=merged["start_pos_y"]) |
    #                     (merged["end_pos_x"]!=merged["end_pos_y"]) |
    #                     (merged["chrom_x"]!=merged["chrom_y"])]

    lenx, leny = (endx-startx)/winsize, (endy-starty)/winsize
    mat = numpy.zeros([leny, lenx])    
    
    for i, row in merged.iterrows():
        xa = numpy.floor((max(startx, row["start_pos_x"]) - startx) / float(winsize))
        xb = numpy.ceil((min(endx, row["end_pos_x"]) - startx) / float(winsize))

        ya = numpy.floor((max(starty, row["start_pos_y"]) - starty) / float(winsize))
        yb = numpy.ceil((min(endy, row["end_pos_y"]) - starty) / float(winsize))
        
        mat[ya:yb, xa:xb] += 1
    
    mat = numpy.ma.masked_array(mat, mask=False)
    # if (merged.shape[0]>0 and (merged["chrom_x"].iloc[0] == merged["chrom_y"].iloc[0]
    #     and (starty<=startx<=endy or starty<=endx<=endy))):
    #     for col in range(startx, endx, winsize):
    #         for row in range(starty, endy, winsize):
    #             if col > row:
    #                 mat.mask[(row-starty)/winsize, (col-startx)/winsize] = True
            
    return mat[::-1,:]

def find_breakpoint(mat, x1, y1, x2, y2):
    margx = mat.max(axis=0)
    breakpointx, orientationx = _find_breakpoints(margx)

    margy = mat.max(axis=1)[::-1]
    breakpointy, orientationy = _find_breakpoints(margy)
    
    if breakpointy is None:
        return None
    if breakpointx is None:
        return None

    curx = int(breakpointx/float(len(margx)) * (x2-x1) + x1)
    cury = int(breakpointy/float(len(margy)) * (y2-y1) + y1)
    
    # print "*"*100
    # print y1, y2, cury
    # print "*"*100

    o = {True:"+", False:"-"}
    orientation = o[orientationx]+o[orientationy]

    return curx, cury, orientation

    
def _find_breakpoints(x):
    if x.max() < 1:
        return None, None

    m = numpy.where((x==x.max()))[0][0]
    diff = x.max() - x
    dist = numpy.sqrt(numpy.abs(numpy.arange(len(x)) - m))

    n = diff / dist.astype(float)
    n[dist==0] = 0

    p = numpy.where(n==n.max())[0][0]
    return p, m<p


# from biorpy import r

def score_event(nx, ny, ncommon, barcode_frequencies, resamples=100):
    """
    perform a resampling test, based on the number of fragments in each region
    and the barcode frequency distribution (the latter because some barcodes
    have more DNA in them than others)
    """

    samples = []
    for _ in range(resamples):
        s1 = numpy.random.choice(len(barcode_frequencies), nx, replace=False, p=barcode_frequencies)
        s2 = numpy.random.choice(len(barcode_frequencies), ny, replace=False, p=barcode_frequencies)
        common = len(set(s1).intersection(set(s2)))
        samples.append(common)

    # make a one-sided one-sample Wilcoxon test
    statistic, pvalue = stats.wilcoxon(numpy.array(samples)-ncommon)
    if numpy.mean(samples) > ncommon:
        pvalue = 1 - pvalue/2.0
    else:
        pvalue = pvalue/2.0

    # result = r["wilcox.test"](samples, mu=ncommon, alternative="less")
    # pvalue2 = result.rx2("p.value")[0]

    # print "::::::", nx, ny, ncommon, (numpy.array(samples)>ncommon).sum(), pvalue, pvalue2
    return pvalue
