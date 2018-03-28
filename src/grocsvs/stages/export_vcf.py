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
##INFO=<ID=EVENT,Number=1,Type=String,Description="Cluster ID used to specify when multiple breakpoints belong to a single complex event">
##FILTER=<ID=NOLONGFRAGS,Description="No samples found with long fragments supporting event (likely segdup)">
##FILTER=<ID=NEARBYSNVS,Description="Very high number of possible SNVs detected in vicinity of breakpoint (likely segdup)">
##FILTER=<ID=ASSEMBLYGAP,Description="Nearby gap in the reference assembly found; these gaps can appear to mimic SVs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (currently either 1 when breakpoint is present in sample or 0 if not)">
##FORMAT=<ID=BS,Number=1,Type=Integer,Description="Number of barcodes supporting the event">
##FORMAT=<ID=BT,Number=1,Type=Integer,Description="Total number of barcodes, calculated as the union of barcodes at each site">
##FORMAT=<ID=F50,Number=1,Type=Float,Description="Median fragment length supporting the breakpoint; estimated as the total across all supporting segments">
##FORMAT=<ID=F90,Number=1,Type=Float,Description="90th percentile of fragment lengths supporting the breakpoint; estimated as the total across all supporting segments">
##FORMAT=<ID=PR,Number=1,Type=Float,Description="p-value for the event, calculated by resampling from the number of supporting and total barcodes">
##FORMAT=<ID=H1x,Number=1,Type=String,Description="Number of supporting barcodes assigned to haplotype 1 (left side of breakpoint)">
##FORMAT=<ID=H2x,Number=1,Type=String,Description="Number of supporting barcodes assigned to haplotype 2 (left side of breakpoint)">
##FORMAT=<ID=H1y,Number=1,Type=String,Description="Number of supporting barcodes assigned to haplotype 1 (right side of breakpoint)">
##FORMAT=<ID=H2y,Number=1,Type=String,Description="Number of supporting barcodes assigned to haplotype 2 (right side of breakpoint)">
{samples}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_columns}
"""

CHROMOSOME = "##contig=<ID={name},length={length}>"
SAMPLE = "##SAMPLE=<ID={name},URL=file://{path}>"
FORMAT = ["GT", "BS", "BT", "F50", "F90", "PR", "H1x", "H2x", "H1y", "H2y"]

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
    sample_info = {}
    sample_info["BS"] = event["{}_shared".format(sample.name)]
    sample_info["BT"] = event["{}_total".format(sample.name)]
    if sample_info["BS"] > 5:
        sample_info["GT"] = "1"
    else:
        sample_info["GT"] = "0"

    if "{}_lengths_50".format(sample.name) in event:
        sample_info["F50"] = event["{}_lengths_50".format(sample.name)]
        sample_info["F90"] = event["{}_lengths_90".format(sample.name)]
    else:
        sample_info["F50"] = "."
        sample_info["F90"] = "."

    sample_info["PR"] = event["{}_p_resampling".format(sample.name)]

    sample_info["H1x"] = event["{}_x_hap0".format(sample.name)]
    sample_info["H2x"] = event["{}_x_hap1".format(sample.name)]
    sample_info["H1y"] = event["{}_y_hap0".format(sample.name)]
    sample_info["H2y"] = event["{}_y_hap1".format(sample.name)]

    sample_info = [sample_info[key] for key in FORMAT]

    return ":".join(map(str, sample_info))

def get_filters(event):
    filters = []

    try:
        if "N=" in event["blacklist"]:
            filters.append("ASSEMBLYGAP")
    except TypeError:
        pass

    if event["nearby_snvs"] >= 15:
        filters.append("NEARBYSNVS")

    if "frag_length_passes" in event and not event["frag_length_passes"]:
        filters.append("NOLONGFRAGS")

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

    info.append("EVENT={}".format(event["cluster"]))
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

        fields = [chromx, x+1, this_id, ref, alt_breakend, 0, filters, info, ":".join(FORMAT)]

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
