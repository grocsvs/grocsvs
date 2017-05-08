import collections
import numpy
import os
import pandas

from grocsvs import graphing
from grocsvs import step
from grocsvs import structuralvariants
from grocsvs import utilities

from grocsvs.stages import final_clustering
from grocsvs.stages import genotyping
from grocsvs.stages import call_readclouds


try:
    from rpy2 import robjects as ro
    from rpy2.robjects import numpy2ri
    numpy2ri.activate()
except ImportError, e:
    ro = None

class VisualizeStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        if ro is None:
            print(" ** rpy2 not installed correctly; skipping visualization step ** ")
            return

        edges = load_edges(options)

        for cluster in edges["cluster"].unique():
            yield VisualizeStep(options, cluster)

    def __init__(self, options, cluster):
        self.options = options
        self.cluster = cluster

    def __str__(self):
        return ".".join([self.__class__.__name__, str(self.cluster)])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "event_{}.pdf".format(self.cluster)

        paths = {
            "visualization": os.path.join(directory, file_name)
        }

        return paths


    def run(self):
        # open(self.outpaths(final=False)["visualization"], "w")

        edges = load_edges(self.options)
        cluster = edges.loc[edges["cluster"]==self.cluster]
        # breakpoints = get_cluster_breakpoints(self.options, self.cluster)

        from rpy2.robjects import r
        r.pdf(self.outpaths(final=False)["visualization"])

        for sample, dataset in sorted(self.options.iter_10xdatasets()):
            graphing.plot_frags(cluster, self.options, sample, dataset)
            for i, row in cluster.iterrows():
                plot_triangles(row, self.options)
        # print "::", breakpoints
        # graphing.visualize_frag_cluster(breakpoints, self.options)


        r["dev.off"]()


def load_edges(options):
    clustering_step = final_clustering.FinalClusterSVsStep(options)
    edges = pandas.read_table(
        clustering_step.outpaths(final=True)["edges"])

    return edges
    # genotyping_step = genotyping.MergeGenotypesStep(options)
    # genotypes = pandas.read_table(
    #     genotyping_step.outpaths(final=True)["genotypes"])
    # return genotypes


def get_cluster_breakpoints(options, cluster):
    # genotyping_step = genotyping.GenotypingStep(options)
    # genotypes = pandas.read_table(
    #     genotyping_step.outpaths(final=True)["genotypes"])

    clustering_step = final_clustering.FinalClusterSVsStep(options)
    edges = pandas.read_table(
        clustering_step.outpaths(final=True)["edges"])
    clustered_edges = edges.loc[edges["cluster"]==cluster]

    breakpoints = []
    for row in clustered_edges.itertuples():
        breakpoints.append((row.chromx, row.x, row.chromy, row.y, row.orientation))

    return breakpoints

def plot_triangles(event, options):
    samples_to_mats = collections.OrderedDict()

    extend_distance = 100000

    startx = event.x - extend_distance
    endx   = event.x + extend_distance
    starty = event.y - extend_distance
    endy   = event.y + extend_distance

    startx = max(0, startx)
    starty = max(0, starty)

    for sample, dataset in sorted(options.iter_10xdatasets()):
        sv = structuralvariants.StructuralVariant(
            event.chromx, event.chromy, event.x, event.y, event.orientation)

        fragsx = call_readclouds.load_fragments(
                options, sample, dataset, event.chromx, startx, endx)
        fragsy = call_readclouds.load_fragments(
                options, sample, dataset, event.chromy, starty, endy)

        _, _, merged_frags, mat, breakpoint = sv.get_sv_info(
                    fragsx, fragsy, ext_dist=extend_distance, winsize=100)

        samples_to_mats[sample.name] = mat


    scale = 1e6

    cur_max = int(max(mat.max() for mat in samples_to_mats.values()) + 1)

    for sample_name, mat in samples_to_mats.items():
        sample = options.samples[sample_name]
        dataset = sample.get_10x_dataset()

        palette = get_palette(mat, 0, cur_max)
        mat = convert(mat, palette)

        ro.r.layout(numpy.array([[1,2],[0,0]]), widths=[9,1], heights=[10,0])
        oldpar = ro.r.par(mar=numpy.array([5,5,3,0]))

        main = "{} {}:{:,}::{}:{:,}{}".format(sample_name, event.chromx, event.x, event.chromy, event.y, event.orientation)
        utilities.plot_matrix_as_image(mat, startx/scale, starty/scale, endx/scale, endy/scale, main=main,
            xlab="{} (MB)".format(event.chromx), ylab="{} (MB)".format(event.chromy))

        p = [ro.r.rgb(*palette(i)) for i in range(cur_max)]
        ro.r.par(mar=numpy.array([5,2,3,0.5]))
        color_bar(p, cur_max, nticks=5)


def color_bar(lut, max_, nticks=11, ticks=None, title=""):
    min_ = 0

    scale = float((len(lut)-1))/(max_-min_)
    if ticks is None:
        tick_max = max_/10*10
        ticks = numpy.linspace(min_, tick_max, nticks)

    ro.r.plot(numpy.array([0,10]), numpy.array([min_,max_]), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    ro.r.axis(2, ticks, las=1)
    for i in range(1, len(lut)):
        y = (i-1)/scale + min_
        ro.r.rect(0,y,10,(y+1/scale)*1.025, col=lut[i], border=ro.NA_Logical)


def convert(matrix, fn):
    converted = []
    for row in range(matrix.shape[0]):
        currow = []
        for col in range(matrix.shape[1]):
            currow.append(fn(matrix[row,col]))
        converted.append(currow)
    return numpy.array(converted)


def get_palette(mat=None, min_=None, max_=None):
    if max_ is None:
        max_ = mat.max()
    if min_ is None:
        min_ = mat.min()

    print min_, max_
    diff = float(max_ - min_ - 10)

    def color(x):
        if x < 10:
            # return (0.8,0.8,0.8)
            v = (x + 5) / 15.0
            return (v,v,1)
        frac = (x-min_-10)/diff
        return (1, 1-frac, 1-frac)
    return color