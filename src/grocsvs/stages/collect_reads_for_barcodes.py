import collections
import errno
import glob
import gzip
import os
import pysam

from grocsvs import step
from grocsvs import utilities
from grocsvs import reference

from grocsvs.stages import barcodes_from_graphs

    


class CollectReadsForBarcodesStep(step.StepChunk):
    """
    """

    @classmethod
    def get_steps(cls, options):
        for sample, dataset in options.iter_10xdatasets():
            for step in cls.get_steps_for_dataset(options, sample, dataset):
                yield step 
                
    @staticmethod
    def get_steps_for_dataset(options, sample, dataset):
        bam = pysam.AlignmentFile(dataset.bam)
        chroms_to_lengths = dict(zip(bam.references, bam.lengths))

        ##limit to chromosomes (>1mb)
        all_chroms_to_lengths = dict(zip(bam.references, bam.lengths))
        chroms_to_lengths=dict()
        #for c in dict(zip(bam.references, bam.lengths)):
        #    if all_chroms_to_lengths[c] > 1e6:
        #        chroms_to_lengths[c] = all_chroms_to_lengths[c]
        large_references = tuple([c for c in chroms_to_lengths])
                

        
        #for (chrom, start, end) in reference.split_genome(bam.references, chroms_to_lengths, 50e6):
        for (chrom, start, end) in reference.split_genome(large_references, chroms_to_lengths, 50e6):
            step = CollectReadsForBarcodesStep(
                options, sample, dataset, chrom, start, end)
            yield step
            
        past_end_step = CollectReadsForBarcodesStep(
                options, sample, dataset, "unmapped", 0, 0)
        yield past_end_step
        
                
    def __init__(self, options, sample, dataset, chrom, start, end):
        self.options = options
        self.sample = sample
        self.dataset = dataset
        self.chrom = chrom
        self.start = start
        self.end = end
        past_end = False
        if self.chrom == "unmapped":
            past_end = True
        self.past_end = past_end

    def __str__(self):
        past_end_str = "+" if self.past_end else ""
        
        return ".".join([self.__class__.__name__,
                         self.sample.name, 
                         self.dataset.id,
                         self.chrom,
                         str(self.start)])
        
    def outpaths(self, final):
        directory = self.results_dir if final \
                    else self.working_dir

        file_name = "event_fastqs.{}.{}.{}.{}".format(
            self.sample.name,
            self.dataset.id,
            self.chrom,
            self.start,
            )

        paths = {
            "event_fastqs": os.path.join(directory, file_name+".fa"),
            "counts": os.path.join(directory, file_name+".counts")
        }

        return paths

    def ensure_dir(self):
        self.events_dir = self.outpaths(final=False)["event_fastqs"]

        try:
            os.makedirs(self.events_dir)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise


    def run(self):
        self.ensure_dir()
        clusters, barcodes_map = load_barcodes_map(self.options, self.sample, self.dataset)
        
#        clusters_to_out_fastqs = self.get_output_fastas(clusters)
        clusters_to_out_fastqs = self.get_output_fasta_paths(clusters)

        count = 0

        self.logger.log("Running on {}:{}-{}".format(self.chrom, self.start, self.end))
        
        for count, (barcode, read) in enumerate(self.fetch()):
            if self.past_end:
                assert read.reference_id == -1
                
            if count % 1000000 == 0:
                self.logger.log("{} {:,} {:.1%}".format(count, read.pos, (read.pos-self.start)/float(self.end-self.start+1)))
                self.logger.log(str(read.reference_id))

            if barcode in barcodes_map:
                for cluster in barcodes_map[barcode]:
#                    out_fasta = clusters_to_out_fastqs[cluster]
                    out_fasta = open(clusters_to_out_fastq_paths[cluster],'w')

                    seq = read.seq if not read.is_reverse else utilities.revcomp(read.seq)
                    seq = seq.replace('\n','')
                    order = "0" if read.is_read1 else "1"
                    out_fasta.write("\t".join([read.query_name, order, seq])+"\n")
                    out_fasta.close()

        with open(self.outpaths(final=False)["counts"], "w") as counts_file:
            counts_file.write("{}\n".format(count))
            

    def fetch(self):
        bam = pysam.AlignmentFile(self.dataset.bam)
        
        if self.past_end:
            chrom = bam.references[-1]
            fetch = bam.fetch(chrom, bam.lengths[-1]-100) # should seek to end
            for i, read in enumerate(fetch):
                if i % 100 == 0: print i
            fetch = bam.fetch(until_eof=True)
        else:
            fetch = bam.fetch(self.chrom, self.start, self.end)

        for read in fetch:
            if read.pos < self.start and not read.pos < 0:
                continue
            if read.pos >= self.end and not self.past_end:
                return
            if read.is_secondary or read.is_supplementary:
                continue
            if self.past_end and read.reference_id != -1:
                continue

            if read.has_tag("BX"):
                barcode = read.get_tag("BX")
            elif read.has_tag("RX"):
                barcode = read.get_tag("RX")
            else:
                continue
                #self.logger.log("ERROR: read has no BX nor RX tags: {}".format(str(read)))
                #sys.exit(1)
                
            barcode = barcode.split(",")[0]#.split("-")[0]

            yield barcode, read
            
    def get_output_fastas(self, clusters):
        clusters_to_out_fastqs = {}
        for cluster in clusters:
            cluster_name = "{}.fa".format(cluster)
            clusters_to_out_fastqs[cluster] = \
                open(os.path.join(self.events_dir, cluster_name), "w")

        return clusters_to_out_fastqs

    def get_output_fasta_paths(self, clusters):
        clusters_to_out_fastqs = {}
        for cluster in clusters:
            cluster_name = "{}.fa".format(cluster)
            clusters_to_out_fastqs[cluster] = os.path.join(self.events_dir, cluster_name)

        return clusters_to_out_fastqs

    

def load_barcodes_map(options, sample, dataset):
    inpaths = []

    input_step = barcodes_from_graphs.BarcodesFromGraphsStep(
        options, sample, dataset)
    inpath = input_step.outpaths(final=True)["sv_barcodes"]
    
    sv_barcodes = utilities.pickle.load(open(inpath))
    barcodes_map = collections.defaultdict(set)
    
    for cluster, barcodes in sv_barcodes.items():
        for barcode in barcodes:
            barcodes_map[barcode].add(cluster)
    
    return sorted(sv_barcodes.keys()), barcodes_map
