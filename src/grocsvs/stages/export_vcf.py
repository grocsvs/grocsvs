import datetime
import numpy

import grocsvs

VCF_HEADER = \
"""##fileformat=VCFv4.3
##fileDate={date}
##source=grocsvs_v{version}
##reference=file://{reference}
{chromosomes}
##INFO=<ID=ASSEMBLED,Number=0,Type=Flag,Description="Event is supported by a read-cloud-based sequence assembly across the breakpoints"> 
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> 
##FORMAT=<ID=BS,Number=1,Type=int,Description="Number of barcodes shared between loci"> 
##FORMAT=<ID=BT,Number=1,Type=int,Description="Total number of barcodes between the loci (union)"> 
{samples}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_columns}
"""

CHROMOSOME = "##contig=<ID={name},length={length}>"
SAMPLE = "##SAMPLE=<ID={name},URL=file://{path}>"
FORMAT = "GT:BS:BT"

def get_header(options):
    chromosomes = []
    for chrom, length in options.reference.chrom_lengths.items():
        chromosomes.append(CHROMOSOME.format(name=chrom, length=length))
    chromosomes = "\n".join(chromosomes)

    samples = []
    sample_columns = []
    for sample_name, sample in options.samples.items():
        samples.append(SAMPLE.format(name=sample_name, path=sample.get_10x_dataset().bam))
        sample_columns.append(sample_name)
    samples = "\n".join(samples)
    sample_columns = "\t".join(sample_columns)

    header = VCF_HEADER.format(
            date=datetime.datetime.now().strftime("%Y%m%d"),
            version=grocsvs.__version__,
            reference=options.reference.path,
            chromosomes=chromosomes,
            samples=samples,
            sample_columns=sample_columns
        )

    return header


def get_alt_breakend(options, chromx, x, chromy, y, orientation):
    # note that vcf uses one-based coords so we add 1 to all positions

    ref = options.reference.fasta[chromx][x].upper()
    if orientation == "+-":
        alt = "{}[{}:{}[".format(ref, chromy, y+1)
    elif orientation == "++":
        alt = "{}]{}:{}]".format(ref, chromy, y+1)
    elif orientation == "-+":
        alt = "]{}:{}]{}".format(chromy, y+1, ref)
    elif orientation == "--":
        alt = "[{}:{}[{}".format(chromy, y+1, ref)

    return ref, alt

def get_sample_info(sample, event):
    sample_info = []
    shared = event["{}_shared".format(sample.name)]
    total = event["{}_total".format(sample.name)]
    if shared > 5:
        sample_info.append("0/1")
    else:
        sample_info.append("0/0")

    sample_info.append(shared)
    sample_info.append(total)

    return ":".join(map(str, sample_info))

def get_filters(event):
    filters = []

    try:
        if "n=" in event["blacklist"]:
            filters.append("NEAR_REF_GAP")
    except TypeError:
        pass

    if len(filters) == 0:
        return "PASS"

    return ";".join(filters)


def get_info(event, mateid=None):
    info = []

    info.append("SVTYPE=BND")

    if mateid is not None:
        info.append("MATEID={}".format(mateid))

    if event["assembled"]:
        info.append("ASSEMBLED")

    if len(info) == 0:
        return "."

    return ";".join(info)

def convert_event(options, event_id, event):
    chromx, x, chromy, y, orientation = event["chromx"], event["x"], event["chromy"], event["y"], event["orientation"]

    for i in [0,1]:
        if i == 1:
            chromx, x, chromy, y = chromy, y, chromx, x
            orientation = orientation[1]+orientation[0]

        this_id = "{}:{}".format(event_id, i)
        other_id = "{}:{}".format(event_id, 1-i)
        
        ref, alt_breakend = get_alt_breakend(options, chromx, x, chromy, y, orientation)
        info = get_info(event, mateid=other_id)
        filters = get_filters(event)

        fields = [chromx, x+1, this_id, ref, alt_breakend, 0, filters, info, FORMAT]

        for sample_name, sample in options.samples.items():
            fields.append(get_sample_info(sample, event))

        yield "\t".join(map(str, fields))

def convert_to_vcf(options, results):
    header = get_header(options)

    for line in header.splitlines():
        yield line

    for event_id, result in results.iterrows():
        for line in convert_event(options, event_id, result):
            yield line
