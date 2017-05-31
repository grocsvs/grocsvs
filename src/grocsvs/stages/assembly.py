import collections
import glob
import logging
import os
import pandas
import shutil
import subprocess

from grocsvs import step
from grocsvs import utilities

from grocsvs.stages import collect_reads_for_barcodes
from grocsvs.stages import cluster_svs


class AssemblyStep(step.StepChunk):
    @staticmethod
    def get_steps(options):
        edges_path = cluster_svs.ClusterSVsStep(options).outpaths(final=True)["edges"]
        graphs_table = pandas.read_table(edges_path)
        graphs_table["chromx"] = graphs_table["chromx"].astype("string")
        graphs_table["chromy"] = graphs_table["chromy"].astype("string")

        clusters = sorted(graphs_table["cluster"].unique())

        for cluster in clusters:
            yield AssemblyStep(options, cluster)


    def __init__(self, options, cluster):
        self.options = options
        self.cluster = cluster

    def __str__(self):
        return ".".join([self.__class__.__name__,
                         str(self.cluster)])

    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        assembly_dir = "assembly.{}".format(self.cluster)

        paths = {
            "assembly_dir": os.path.join(directory, assembly_dir)
        }

        return paths

    def run(self):
        self.assembly_dir = self.outpaths(final=False)["assembly_dir"]
        utilities.ensure_dir(self.assembly_dir)

        # TODO: how many reads is reasonable here?
        max_reads = 5e6

        fasta_path = self.combine_fastas(max_reads)
        contigs_path = self.run_assembly(fasta_path)
        self.align_contigs(contigs_path)

    def combine_fastas(self, max_reads):
        self.logger.log("Combining input fastas...")

        input_file_paths = []

        for sample, dataset in self.options.iter_10xdatasets():
            steps = collect_reads_for_barcodes.CollectReadsForBarcodesStep.get_steps_for_dataset(
                self.options, sample, dataset)

            
            for cur_step in steps:
                self.logger.log("step: "+str(cur_step))
                input_directory = cur_step.outpaths(final=True)["event_fastqs"]
                self.logger.log(input_directory)
                input_file_path = os.path.join(input_directory, "{}.fa".format(self.cluster))
                input_file_paths.append(input_file_path)

        outpath = os.path.join(self.assembly_dir, 
                               "cluster_{}.fa".format(self.cluster))

        combine_process = subprocess.Popen("xargs cat | sort -k1,1 -k2,2n -S 3G",
                                           shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        combine_process.stdin.write("\n".join(input_file_paths))
        combine_process.stdin.close()

        #print input_file_paths
        with open(outpath, "w") as outf:
            reads_to_ids = collections.defaultdict(list)
            counts = 0
            
            for line in combine_process.stdout:
                #self.logger.log("LINE: "+line)
                if counts < 10:
                    print line.strip()
                fields = line.split()
                
                #some fastqs split across two lines? 
                if len(fields) < 2: continue
                
                outf.write("\n".join([">"+fields[0], fields[2], ""]))
                reads_to_ids[fields[0]].append(line)
                counts += 1

                if max_reads is not None and counts >= 2*max_reads:
                    self.logger.log("Number of reads ({}) exceeds maximum; stopping early".format(
                        counts))
                    break
                    
        for id_ in reads_to_ids:
            if len(reads_to_ids[id_]) != 2:
                self.logger.log("Pair error: {}".format(reads_to_ids[id_]))
                    
                
        self.logger.log("Found {:,} reads for a total of {:,} read pairs".format(
            counts, counts/2))

        if counts == 0:
            raise Exception("Error: event '{}' not supported by any reads".format(self.cluster))

        return outpath

    def _combine_fastas(self, input_file_names, outpath, max_reads=None):
        total_lines = 0
        with open(outpath, "w") as outf:
            for input_file_name in input_file_names:
                inpath = os.path.join(
                    input_file_name,
                    "{}.fa".format(self.cluster))

                inf = open(inpath)
                for line in inf:
                    outf.write(line)

                    total_lines += 1
                    if max_reads is not None and (total_lines/4) >= max_reads:
                        self.logger.log("Hit maximum number of reads; stopping early")
                        return total_lines
        return total_lines

    def run_assembly(self, fasta_path):
        # we're parallelizing across clusters so specify only a single thread
        # here; it's possible if there are very few clusters that we could
        # speed things up by using more threads, but idba_ud isn't actually
        # that well parallelized so we're probably not losing a lot here
        contigs_path = os.path.join(self.assembly_dir, "contig.fa")

        command = "{idba} -r {fasta} -o {outdir} --num_threads 1"
        command = command.format(
            idba=self.options.binary("idba_ud"),
            fasta=fasta_path,
            outdir=self.assembly_dir)

        self.logger.log("Running assembly command: '{}'".format(command))
        cmd = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
        retcode = cmd.wait()
        
        if not os.path.exists(contigs_path):
            if retcode != 0 and  "invalid insert distance" in cmd.stderr.read():
                self.logger.log("idba_ud internal error; this is probably a bug in idba_ud, so we have to bail here")
                open(contigs_path, "w")
            else:
                error_message = "assembly failed to produce contig.fa"
                self.logger.log(error_message)
                raise Exception()
        elif retcode != 0:
            self.logger.log("something funky happened while running idba_ud (got error code "
                            "{}) but there's not much we can do so continuing".format(retcode))
            
        return contigs_path


    def align_contigs(self, contigs_path):
        sorted_bam_path = os.path.join(self.assembly_dir, "contigs.sorted.bam")

        bwa_command = "bwa mem -t1 {index} {contigs}"
        samtools_command = "samtools sort -T asdf -O bam - > {sorted_bam}"
        bwa_command = bwa_command + " | " + samtools_command
        bwa_command = bwa_command.format(
            index=self.options.bwa_index,
            contigs=contigs_path,
            sorted_bam=sorted_bam_path)

        self.logger.log("Aligning contigs against reference using bwa: "
                        "'{}'".format(bwa_command))
        subprocess.check_call(bwa_command, shell=True)

        error_message = "bwa alignment failed to produce expected bam file"
        assert os.path.exists(sorted_bam_path), error_message

        index_command = "samtools index {sorted_bam}".format(
            sorted_bam=sorted_bam_path)
        subprocess.check_call(index_command, shell=True)

        error_message = "failed to create bam index"
        assert os.path.exists(sorted_bam_path+".bai"), error_message

        return sorted_bam_path
