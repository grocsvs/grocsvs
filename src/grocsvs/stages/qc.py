import collections
import numpy
import os
import pandas

from grocsvs import step


class QCStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        # for sample, dataset in options.iter_10xdatasets():
        yield QCStep(options)#, sample, dataset)
        # yield step
        
    def __init__(self, options):
        self.options = options


    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "report": os.path.join(directory, "qc_report.tsv"),
            "visual": os.path.join(directory, "report.pdf")
        }

        return paths

    def run(self):
        outpaths = self.outpaths(final=False)

        qc = collections.OrderedDict()
        tenx_sample_infos = collections.OrderedDict()

        for sample_name, sample in self.options.samples.items():
            sample_info = self.options.sample_info(sample_name)
            tenx_dataset = sample.get_10x_dataset()
            tenx_sample_info = sample_info[tenx_dataset.id]

            qc[sample_name] = get_report(tenx_sample_info)
            tenx_sample_infos[sample_name] = tenx_sample_info

        visualize_report(tenx_sample_infos, outpaths["visual"])

        qc = pandas.DataFrame(qc)
        print qc

        # this prints a nicely formatted output table, rather than tab-delimited
        with open(outpaths["report"], "w") as outf:
            outf.write(qc.to_string())



def get_report(tenx_sample_info):
    qc = collections.OrderedDict()

    quantiles = tenx_sample_info["frag_length_info"]["quantiles"]

    qc["Fragment lengths"] = ""
    qc[" - 10%"] = "{:,}".format(int(quantiles[0.1]))
    qc[" - 25%"] = "{:,}".format(int(quantiles[0.25]))
    qc[" - 50%"] = "{:,}".format(int(quantiles[0.5]))
    qc[" - 75%"] = "{:,}".format(int(quantiles[0.75]))
    qc[" - 90%"] = "{:,}".format(int(quantiles[0.9]))
    qc[" - 95%"] = "{:,}".format(int(quantiles[0.95]))

    qc[" - mean"] = "{:,}".format(int(tenx_sample_info["frag_length_info"]["mean"]))
    qc[" - std"] = "{:,}".format(int(tenx_sample_info["frag_length_info"]["std"]))

    qc[" - N50"] = "{:,}".format(int(tenx_sample_info["frag_length_info"]["N50"]))

    depths = tenx_sample_info["physical_depths"]
    qc["Physical depth (Cf)"] = "{:.1f} +/- {:.1f}".format(numpy.mean(depths), numpy.std(depths))

    C_Rs = tenx_sample_info["coverage_of_fragments"]
    C_Rs = C_Rs["coverages"] / C_Rs["lengths"].astype(float)

    qc["Read coverage per fragment position (Cr)"] = "{:.2f} +/- {:.2f} (median={:.2f})".format(
        numpy.mean(C_Rs), numpy.std(C_Rs), numpy.median(C_Rs))

    qc["Good barcodes"] = "{:,}".format(tenx_sample_info["good_bc_count"])


    return pandas.Series(qc)


def visualize_report(tenx_sample_infos, outpath):
    try:
        from rpy2 import robjects as ro
        from grocsvs import plotting
        
        
        ro.r.pdf(outpath)

        frag_lengths = [numpy.log10(sample_info["frag_length_info"]["sampled"]) for sample_info in tenx_sample_infos.values()]
        max_ = max(numpy.percentile(cur_frag_lengths, 99.5) for cur_frag_lengths in frag_lengths)

        plotting.ecdf(frag_lengths, tenx_sample_infos.keys(), xlim=[0, max_], main="Fragment length distribution",
            xlab="Fragment length (log10)", legendWhere="bottomright")

        bc_counts = dict((name, sample_info["good_bc_count"]) for name, sample_info in tenx_sample_infos.items())
        plotting.barPlot(bc_counts, main="Number of high-quality barcodes", ylim=[0, 1.1*max(bc_counts.values())])

        # oldpar = r.par(mfrow=[min(3, len(tenx_sample_infos)), 2])
        oldpar = ro.r.par(mfrow=[2,1])

        for name, tenx_sample_info in tenx_sample_infos.items():
            C_Rs = tenx_sample_info["coverage_of_fragments"]
            C_Rs = C_Rs["coverages"] / C_Rs["lengths"].astype(float)

            ro.r.hist(ro.FloatVector(C_Rs), breaks=50, xlab="Fragment coverage by short-reads (C_R)", main=name)
            ro.r.hist(ro.FloatVector(tenx_sample_info["physical_depths"]), breaks=100, xlab="Coverage by long fragments (C_F)",
                main=name)

        ro.r.par(oldpar)

        ro.r["dev.off"]()
    except:
        with open(outpath, "w") as f:
            f.write("[the visual report requires rpy2 to be correctly installed]")
    
