import collections
import numpy
import os
import pandas
import pysam

from grocsvs import step
from grocsvs import datasets as svdatasets
from grocsvs.stages import genotyping
from grocsvs.stages import mate_pair_evidence
from grocsvs.stages import export_vcf

pandas.options.display.max_rows = 250
pandas.options.display.width = 250


def cwp(count, total):
    """
    count with percentage
    """
    if total == 0:
        return "0/0"
    
    return "{} ({:.0%})".format(count, count/float(total))


class PostprocessingStep(step.StepChunk):
    """
    """

    @staticmethod
    def get_steps(options):
        yield PostprocessingStep(options)

    def __init__(self, options):
        self.options = options

    def __str__(self):
        return ".".join([self.__class__.__name__])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        paths = {
            "summary": os.path.join(directory, "summary.txt"),
            "classes": os.path.join(directory, "classes.txt"),
            "vcf":     os.path.join(directory, "svs.vcf"),
        }

        if len(self.get_mate_pair_datasets()) > 0:
            paths["validation"] = os.path.join(directory, "validation.tsv")

        return paths


    def run(self):
        results = self.load_results()
        genotype_data = results.loc[results["kind"]=="breakpoint"]

        # genotype_data = genotype_data.iloc[:25]
        
        # self.frag_length_filter(genotype_data)

        self.convert_to_vcf(genotype_data)
        short_frag_comparison = self.compare_to_frag_libraries(genotype_data)
        #validation = self.validate(genotype_data)
        # self.logger.log("skipping validation!")
        # validation = None
        self.summarize(genotype_data, short_frag_comparison)

        self.summarize_classes(genotype_data, short_frag_comparison)

                
    def convert_to_vcf(self, genotypes):
        outpath = self.outpaths(final=False)["vcf"]

        with open(outpath, "w") as outfile:
            for line in export_vcf.convert_to_vcf(self.options, genotypes):
                outfile.write(line+"\n")
                print(line)


    def get_mate_pair_datasets(self):
        # figure out which samples have mate-pair datasets
        mp_datasets = {}
        for sample_name, sample in self.options.samples.items():
            for dataset in sample.datasets:
                if isinstance(dataset, svdatasets.MatePairDataset):
                    mp_datasets[sample_name] = dataset
        return mp_datasets

    def get_short_frag_datasets(self):
        # figure out which samples have mate-pair datasets
        sf_datasets = {}
        for sample_name, sample in self.options.samples.items():
            for dataset in sample.datasets:
                if isinstance(dataset, svdatasets.ShortFragDataset):
                    sf_datasets[sample_name] = dataset
        return sf_datasets

    def compare_to_frag_libraries(self, genotype_data):
        mp_datasets = self.get_mate_pair_datasets()
        sf_datasets = self.get_short_frag_datasets()

        if len(mp_datasets) + len(sf_datasets) == 0:
            return None
        
        # check each event for support in each mate-pair dataset        
        results = {}
        for i, (_, event) in enumerate(genotype_data.iterrows()):
            if i % 10 == 0:
                print i
            key = "{}:{}::{}:{}{}".format(*event[["chromx", "x", "chromy", "y", "orientation"]])
            results[key] = {}

            if len(mp_datasets) > 0:
                results[key].update(self.do_comparisons(event, mp_datasets, frag_type="mp"))

            if len(sf_datasets) > 0:
                results[key].update(self.do_comparisons(event, sf_datasets, frag_type="sf"))

            # results[key]["assembled"] = event["assembled"]
            results[key]["chromx"] = event["chromx"]
            results[key]["x"] = event["x"]
            results[key]["chromy"] = event["chromy"]
            results[key]["y"] = event["y"]
            results[key]["orientation"] = event["orientation"]
            
            # if i > 31:
            #     break

        df = pandas.DataFrame(results).T
        df = pandas.merge(genotype_data, df, on=["chromx","x","chromy","y","orientation"], how="outer")
        print df

        outpath = self.outpaths(final=False)["validation"]
        df.to_csv(outpath, sep="\t")
                
        return df

    
    def do_comparisons(self, event, frag_datasets, frag_type):
        samples_to_results = {}

        if frag_type == "mp":
            max_ins = 10000
            backwards = False
            if self.options.flip_matepair_orientation:
                backwards = True
        elif frag_type == "sf":
            max_ins = 1000
            backwards = True

        imprecision = 0

        for name, sample in self.options.samples.items():
            if not name in frag_datasets: continue
            mpdataset = frag_datasets[name]

            bam = pysam.AlignmentFile(mpdataset.bam)

            # only validate the event if it's supposedly present in the 10X data
            # this is only helpful for estimating false positives, not false negatives
            # but it probably wouldn't be the best way to do that anyways

            # TODO: contstant
            # if event["{}_shared".format(sample.name)]/float(event["{}_total".format(sample.name)]) > 0.05:
            try:
                result = mate_pair_evidence.get_evidence(
                    event["chromx"], event["x"], event["chromy"], event["y"],
                    event["orientation"][0], event["orientation"][1],
                    bam, max_ins, imprecision, backwards=backwards)
            except ValueError:
                # TODO: should handle different chrom names (eg "1" vs "chr1") better here
                chromx = str(event["chromx"]).replace("chr", "")
                chromy = str(event["chromy"]).replace("chr", "")
                result = mate_pair_evidence.get_evidence(
                    chromx, event["x"], chromy, event["y"],
                    event["orientation"][0], event["orientation"][1],
                    bam, max_ins, imprecision, backwards=backwards)
                
            samples_to_results["{}_shared_{}".format(sample.name, frag_type)] = result[0]
            samples_to_results["{}_total_{}".format(sample.name, frag_type)] = result[1]

        return samples_to_results



    def summarize(self, breakpoints, validation):
        # TODO: this shouldn't be part of the pipeline, really
        if "SarcomaC1" in breakpoints:
            # SOMATIC
            print "** SOMATIC **"
            breakpoints = breakpoints.loc[breakpoints["SarcomaC1"]!=True]


        summary = collections.OrderedDict()

        summary["total breakpoints"] = len(breakpoints)

        breakpoints_per_event = breakpoints["cluster"].value_counts().value_counts()
        summary["breakpoints"] = ""
        for n, count in breakpoints_per_event.sort_index().iteritems():
            summary["{:>18}".format("n={}".format(n))] = count

        assembled = breakpoints.loc[breakpoints["assembled"]]
        summary["assembled"] = cwp(len(assembled), len(breakpoints))

        intra = breakpoints.loc[breakpoints["chromx"]==breakpoints["chromy"]]
        summary["intrachromosomal"] = cwp(len(intra), len(breakpoints))

        distances = ((intra["x"]-intra["y"])/1e3).abs()
        summary[" - median distance"] = "{:,.0f} kb".format(distances.median())
        summary[" -    min distance"] = "{:,.0f} kb".format(distances.min())
        summary[" -    max distance"] = "{:,.0f} kb".format(distances.max())

        count = breakpoints[self.options.samples.keys()].sum(axis=1)
        private = breakpoints.loc[count==1]
        shared = breakpoints.loc[count>1]

        summary["private"] = cwp(len(private), len(breakpoints))
        private_intra = (private["chromx"]==private["chromy"]).sum()
        summary[" - pintrachromosomal"] = cwp(private_intra, len(private))

        summary["shared"] = cwp(len(shared), len(breakpoints))
        shared_intra = (shared["chromx"]==shared["chromy"]).sum()
        summary[" - sintrachromosomal"] = cwp(shared_intra, len(shared))

        validation_columns = []
        validation_samples = []
        for sample_name in sorted(self.options.samples):
            summary[sample_name] = breakpoints[sample_name].sum()

            if validation is not None:
                validation_column = "{}_shared_mp".format(sample_name)
                if validation_column in validation:
                    validation_columns.append(validation_column)
                    validation_samples.append(sample_name)

                validation["{}_present".format(sample_name)] = (validation["{}_shared".format(sample_name)]/validation["{}_total".format(sample_name)].astype(float)>0.10)
            
        if validation is not None:
            print "VALIDATION COLUMNS:", validation_columns

            to_validate = validation.loc[validation[["{}_present".format(name) for name in validation_samples]].any(axis=1)]
            
            validated = (to_validate[validation_columns]>=10).any(axis=1)
            
            summary["validated"] = cwp(validated.sum(), len(validated))

            to_validate_blacklist = to_validate.loc[to_validate["quality"]=="PASS"]
            validated_blacklist = (to_validate_blacklist[validation_columns]>=10).any(axis=1)
            summary[" - passing blacklist"] = cwp(validated_blacklist.sum(), len(validated_blacklist))

            to_validate_simple = to_validate.loc[(to_validate["frag_length_passes"])]
            validated_simple = (to_validate_simple[validation_columns]>=10).any(axis=1)
            summary[" - simple"] = cwp(validated_simple.sum(), len(validated_simple))

            to_validate_assembled = to_validate.loc[to_validate["assembled"]]
            validated_assembled = (to_validate_assembled[validation_columns]>=10).any(axis=1)
            summary[" - assembled"] = cwp(validated_assembled.sum(), len(validated_assembled))

            to_validate_assembled_or_blacklist = to_validate.loc[to_validate["assembled"] | (to_validate["quality"]=="PASS")]
            validated_assembled_or_blacklist = (to_validate_assembled_or_blacklist[validation_columns]>=10).any(axis=1)
            summary[" - assembled|pass blacklist"] = cwp(validated_assembled_or_blacklist.sum(), len(validated_assembled_or_blacklist))

            to_validate_assembled_and_simple = to_validate.loc[to_validate["assembled"] & (to_validate["frag_length_passes"])]
            validated_assembled_and_simple = (to_validate_assembled_and_simple[validation_columns]>=10).any(axis=1)
            summary[" - assembled&simple"] = cwp(validated_assembled_and_simple.sum(), len(validated_assembled_and_simple))

            to_validate_assembled_or_nsegdups = to_validate.loc[to_validate["assembled"] | (to_validate["nearby_snvs"]<10)]
            to_validate_assembled_or_nsegdups = to_validate_assembled_or_nsegdups.loc[~to_validate_assembled_or_nsegdups["blacklist"].str.contains("N=").astype(bool)]
            validated_assembled_or_nsegdups = (to_validate_assembled_or_nsegdups[validation_columns]>=10).any(axis=1)
            summary[" - assembled|non-segdups"] = cwp(validated_assembled_or_nsegdups.sum(), len(validated_assembled_or_nsegdups))

        summary = pandas.Series(summary)

        outpath = self.outpaths(final=False)["summary"]
        summary.to_csv(outpath, sep="\t")
            
        print summary


        chrom_pair_counts = collections.defaultdict(collections.Counter)
            
        for i, row in breakpoints.iterrows():
            key = "{},{}".format(*sorted((row.chromx, row.chromy)))

            if row[self.options.samples.keys()].sum() == 1:
                chrom_pair_counts["private"][key] += 1
                cur = row[self.options.samples.keys()]
                cur = cur.index[cur.values.astype(bool)][0]
                chrom_pair_counts[cur][key] += 1
            elif row[self.options.samples.keys()].sum() > 1:
                chrom_pair_counts["shared"][key] += 1

        for key, group in breakpoints.groupby("cluster"):
            chrom_pairs = set()
            for i, row in group.iterrows():
                chrom_pairs.add(tuple(sorted((row.chromx, row.chromy))))

            if len(group) > 1:
                for i, row in group.iterrows():
                    key = "{},{}".format(*sorted((row.chromx, row.chromy)))
                    chrom_pair_counts["complex"][key] += 1
            else:
                pair = chrom_pairs.pop()
                chrom_pair_counts["simple"]["{},{}".format(*pair)] += len(group)
            
        chrom_pair_counts = pandas.DataFrame(chrom_pair_counts)
        chrom_pair_counts.index = pandas.MultiIndex.from_tuples([k.split(",") for k in chrom_pair_counts.index])
        
        if not "complex" in chrom_pair_counts:
            chrom_pair_counts["complex"] = 0

        chrom_pair_counts["complex_frac"] = chrom_pair_counts["complex"].astype(float)/(chrom_pair_counts["complex"]+chrom_pair_counts["simple"])
        print chrom_pair_counts.sort_index()


    def summarize_classes(self, genotypes, short_frag_comparison):
        classes = collections.OrderedDict()

        for col in ["chromx", "x", "chromy", "y", "orientation"]:
            classes[col] = genotypes[col]

        for sample_name in self.options.samples:
            cur_afs = genotypes["{}_shared".format(sample_name)] \
                          / genotypes["{}_total".format(sample_name)].astype(float)

            cur_classes = pandas.Series(index=cur_afs.index)
            cur_classes[(cur_afs < 0.02) | cur_afs.isnull()] = 0
            cur_classes[cur_afs >= 0.02] = 1

            classes[sample_name] = cur_classes.fillna(9).astype(int)
            classes["{}_shared".format(sample_name)] = genotypes["{}_shared".format(sample_name)]
            classes["{}_total".format(sample_name)] = genotypes["{}_total".format(sample_name)]
            classes["{}_p_resampling".format(sample_name)] = genotypes["{}_p_resampling".format(sample_name)]

            for short_frag_type in ["mp", "sf"]:
                shared_key = "{}_shared_{}".format(sample_name, short_frag_type)
                total_key = "{}_total_{}".format(sample_name, short_frag_type)
                if short_frag_comparison is not None and shared_key in short_frag_comparison.columns:
                    classes[shared_key] = short_frag_comparison[shared_key]
                    classes[total_key] = short_frag_comparison[total_key]

            
        classes = pandas.DataFrame(classes)
        classes_summary = classes[self.options.samples.keys()].astype(str).apply(lambda x: "".join(x), axis=1)
        #classes.insert(6, "classes", classes_summary)
        classes["classes"] = classes_summary

        classes.to_csv(self.outpaths(final=False)["classes"], sep="\t")

        
    def load_results(self):
        # input_step = final_clustering.FinalClusterSVsStep(self.options)
        input_step = genotyping.MergeGenotypesStep(self.options)
        inpath = input_step.outpaths(final=True)["genotypes"]

        results = pandas.read_table(inpath)

        min_fraction = 0.10
        max_p = 1e-6

        for sample_name in self.options.samples:
            p = results["{}_p_resampling".format(sample_name)]
            ratio = results["{}_shared".format(sample_name)]/results["{}_total".format(sample_name)].astype(float)

            results[sample_name] = ((p < max_p) & (ratio >= min_fraction))

        print results
        return results
