import itertools
import numpy
import os
import pandas
from scipy import ndimage

from grocsvs import step
from grocsvs import structuralvariants
from grocsvs import utilities

from grocsvs.stages import call_readclouds
from grocsvs.stages import sv_candidate_regions



class SVCandidatesStep(step.StepChunk):
    """
    this finds SV candidates within candidate regions -- ie, we're honing in on 
    the actual coordinates

    the basic approach is to find pileups of fragment ends that are substantially
    higher than the number we'd expect given the barcode counts on the margins 
    (ie x and y bins); this is repeated for the four possible orientation
    combinations
    """

    @staticmethod
    def get_steps(options):
        chroms = options.reference.chroms
        for sample, dataset in options.iter_10xdatasets():
            for chromx, chromy in itertools.product(chroms, chroms):
                if options.reference.compare_chroms(chromx, chromy) < 0: continue

                yield SVCandidatesStep(
                    options, sample, dataset, chromx, chromy)


    def __init__(self, options, sample, dataset, chromx, chromy):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chromx = chromx
        self.chromy = chromy


    def __str__(self):
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         self.chromx,
                         self.chromy])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        svs_file_name = "svs.{}.{}.{}.{}.tsv".format(
            self.sample.name,
            self.dataset.id,
            self.chromx,
            self.chromy
            )

        report_file_name = "sv_report.{}.{}.{}.{}.tsv".format(
            self.sample.name,
            self.dataset.id,
            self.chromx,
            self.chromy
            )

        paths = {
            "svs": os.path.join(directory, svs_file_name),
            "report": os.path.join(directory, report_file_name)
        }

        return paths


    def run(self):
        outpaths = self.outpaths(final=False)

        self.logger.log("loading...")
        sv_regions = self.load_candidate_regions()
        # print sv_regions

        self.logger.log("getting SV breakpoints...")
        svs = []
        reports = []
        for i, (_,sv_region) in enumerate(sv_regions.iterrows()):
            if i % 100 == 0:
                print i, "of", len(sv_regions)#, sv_region
            cur_svs = self.get_svs(sv_region.to_dict())
            svs.append(cur_svs)
            cur_report = sv_region.copy()
            cur_report["found"] = len(cur_svs)

            reports.append(cur_report)

        if len(svs) > 0:
            svs = pandas.concat(svs, ignore_index=True)
            svs.to_csv(outpaths["svs"], index=False, sep="\t")
        else:
            # touch
            open(outpaths["svs"], "w")

        report = pandas.DataFrame(reports)
        report.to_csv(outpaths["report"], sep="\t", index=False)

    def load_candidate_regions(self):    
        input_step = sv_candidate_regions.SVCandidateRegionsStep(
            self.options, 
            self.sample, 
            self.dataset, 
            self.chromx, 
            self.chromy)

        inpaths = input_step.outpaths(final=True)    
        path = inpaths["sv_regions"]


        if os.stat(path).st_size == 0:

        # if True:
            # print "SKIPPING "*1000
            candidate_regions = pandas.DataFrame()
        else:
            candidate_regions = pandas.read_table(path)

            candidate_regions["chromx"] = candidate_regions["chromx"].astype("string")
            candidate_regions["chromy"] = candidate_regions["chromy"].astype("string")

        if self.chromx == self.chromy:
            diagonal = []
            chrom_length = self.options.reference.chrom_lengths[self.chromx]
            distance = input_step.max_distance_between_events() * 2

            for start in range(0, chrom_length, distance/2):
                cur_region = [self.chromx, self.chromy,
                              start, start+distance,
                              start, start+distance]
                diagonal.append(cur_region)

            diagonal = pandas.DataFrame(diagonal, 
                columns=["chromx", "chromy", "startx", "endx", "starty", "endy"])

            candidate_regions = pandas.concat([candidate_regions, diagonal])

        
        return candidate_regions


    def get_svs(self, sv_region):
        breakpoints = []
        window_size = 10000
        columns = ["chromx", "x", "chromy", "y", "orientation"]

        sv_region = sv_region.copy()
        sv_region["startx"] -= 20000
        sv_region["endx"] += 20000
        sv_region["starty"] -= 20000
        sv_region["endy"] += 20000

        sv_region["startx"] = max(0, sv_region["startx"])
        sv_region["starty"] = max(0, sv_region["starty"])

        fragsx, fragsy, merged, mat, pc = self.get_sv_data(sv_region)
        if len(merged) < 20:
            # print "skipping", sv_region, pc, len(merged)
            return pandas.DataFrame(columns=columns)

        print sv_region
        # print ":::::::::::", pc, len(merged)
        # print fragsx
        # print fragsy
        # print merged

        # print sv_region
        mats = self.get_fragment_ends_matrices(merged, sv_region, window_size)
        bg_mats = get_bg_mats(fragsx, fragsy, sv_region, window_size)
        
        

        rolling = 0
        for o in bg_mats.keys():
            cur_breakpoints = get_svs(mats[o], bg_mats[o], sv_region, window_size, rolling=rolling)
            for bp in cur_breakpoints:
                breakpoints.append([sv_region["chromx"], bp[0], sv_region["chromy"], bp[1], o])

        breakpoints = pandas.DataFrame(breakpoints, columns=columns)




        # if len(breakpoints) > 1:
        # # from biorpy import r, plotting
        #     import rpy2.robjects as ro
        #     from rpy2.robjects import numpy2ri
        #     numpy2ri.activate()

        #     ro.r.pdf("temp.{},{}.pdf".format(sv_region["startx"], sv_region["starty"]))
        #     for orientation in sorted(mats):
        #         print mats[orientation]
        #         print bg_mats[orientation]

        #         x1 = sv_region["startx"]
        #         x2 = sv_region["endx"]
        #         y1 = sv_region["starty"]
        #         y2 = sv_region["endy"]
        #         scale = 1e6

        #         try:
        #             utilities.plot_matrix_as_image(mats[orientation], x1/scale, y1/scale, x2/scale, y2/scale, 
        #                 main="fg {:.1f}".format(mats[orientation].max()))
        #         except:
        #             ro.r.plot(numpy.array([0]), main="fg")
        #             pass

        #         try:
        #             utilities.plot_matrix_as_image(bg_mats[orientation], x1/scale, y1/scale, x2/scale, y2/scale,
        #                 main="bg {:.1f}".format(bg_mats[orientation].max()))
        #         except:
        #             ro.r.plot(numpy.array([0]), main="bg")
        #             pass

        #     ro.r["dev.off"]()

        
        return breakpoints



    def get_sv_data(self, sv_region):
        fragsx = call_readclouds.load_fragments(
            self.options, self.sample, self.dataset, 
            sv_region["chromx"], sv_region["startx"], sv_region["endx"],
            min_reads_per_frag=0)
        fragsy = call_readclouds.load_fragments(
            self.options, self.sample, self.dataset, 
            sv_region["chromy"], sv_region["starty"], sv_region["endy"],
            min_reads_per_frag=0)
        
        fragsx["idx"] = fragsx.index
        fragsy["idx"] = fragsy.index

        _,_,merged = structuralvariants.merge_fragments(
                sv_region["startx"], sv_region["endx"],
                sv_region["starty"], sv_region["endy"], 
                fragsx, fragsy)
        
        precount = len(merged)
        non_spanning = ((merged["start_pos_x"]!=merged["start_pos_y"]) |
                        (merged["end_pos_x"]!=merged["end_pos_y"]) |
                        (merged["chrom_x"]!=merged["chrom_y"]))

        merged = merged.loc[non_spanning]

        mat = structuralvariants.overlap_matrix(
                merged, 
                sv_region["startx"], sv_region["endx"], 
                sv_region["starty"], sv_region["endy"], 10000)
        
        return fragsx, fragsy, merged, mat, precount

    def get_fragment_ends_matrices(self, merged, sv_region, window_size):
        mats = {}
        selectors = {"+":"end_pos_{}",
                     "-":"start_pos_{}"}
        
        binsx = numpy.arange(sv_region["startx"], sv_region["endx"]+window_size*2, window_size)
        binsy = numpy.arange(sv_region["starty"], sv_region["endy"]+window_size*2, window_size)
        
        for orientationx in "+-":
            for orientationy in "+-":
                fx = merged[selectors[orientationx].format("x")].values
                fy = merged[selectors[orientationy].format("y")].values
                hist = numpy.histogram2d(fy, fx, (binsy, binsx))[0]
                mats[orientationx+orientationy] = hist
        return mats


def get_bg_mat(bcsx, bcsy):
    mat = numpy.empty((len(bcsy), len(bcsx)))
    for i in range(len(bcsx)):
        for j in range(len(bcsy)):
            mat[j,i] = len(bcsx[i].union(bcsy[j]))
    return mat

def get_bg_mats(fragsx, fragsy, sv_region, window_size):
    bg_mats = {}
    selectors = {"+":"end_pos", "-":"start_pos"}

    binsx = numpy.arange(sv_region["startx"], sv_region["endx"]+window_size, window_size)
    binsy = numpy.arange(sv_region["starty"], sv_region["endy"]+window_size, window_size)
    
    for orientationx in "+-":
        binx = numpy.digitize(fragsx[selectors[orientationx]], binsx)-1
        gx = fragsx.groupby(binx)
        bcsx = [set(gx.get_group(k)["bc"]) if k in gx.groups else set() for k in range(len(binsx))]
        
        for orientationy in "+-":
            biny = numpy.digitize(fragsy[selectors[orientationy]], binsy)-1
            gy = fragsy.groupby(biny)
            bcsy = [set(gy.get_group(k)["bc"]) if k in gy.groups else set() for k in range(len(binsy))]
    
            bg_mats[orientationx+orientationy] = get_bg_mat(bcsx, bcsy)
            
    return bg_mats

def get_svs(mat, bg_mat, sv_region, window_size, rolling=0):
    if rolling > 0:
        weights = numpy.ones((rolling*2+1, rolling*2+1))    
        mat = ndimage.convolve(mat, weights, mode="constant")
        bg_mat = ndimage.convolve(bg_mat, weights, mode="constant")

    norm = mat/bg_mat
    norm[numpy.isnan(norm)] = 0
    norm = numpy.ma.masked_array(norm, mask=False)

    breakpoints = []

    while not norm.mask.all():
        where = numpy.where(norm==norm.max())
        where = (where[0][0], where[1][0])

        # TODO: constants
        print "*"*100
        print "DEBUG "*10
        print "MAT"
        print mat
        print mat.dtype
        print "NORM"
        print norm
        print norm.dtype
        print "MAX:", norm.max()
        print "CALC", numpy.where(norm==norm.max())
        print "WHERE:", where
        
        is_good = (mat[where] > 25 and norm[where] > 0.05)
        
        if is_good:
            breakpoint = (where[1]*window_size + sv_region["startx"],
                          where[0]*window_size + sv_region["starty"])

            breakpoints.append(breakpoint)

            # TODO: constant for extend; this determines the closest
            # any two breakpoints can be from one another
            norm.mask[get_matrix_rect(norm, where, 10)] = True
        else:
            break

    return breakpoints

def get_matrix_rect(mat, pos, extend):
    top = max(0, pos[0]-extend)
    bottom = min(mat.shape[0], pos[0]+extend+1)
    left = max(0, pos[1]-extend)
    right = min(mat.shape[1], pos[1]+extend+1)

    return numpy.s_[top:bottom, left:right]
