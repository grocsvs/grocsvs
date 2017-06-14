import argparse
import collections
import errno
import os
import six
import sys
import tqdm

import _bcs
import xopen

def ensure_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

def load_gem_barcodes(bcs_path):
    bcs_file = xopen.xopen(bcs_path)
    bcs = []
    for line in bcs_file:
        bcs.append(line.strip())
    bcs = set(bcs)
    return bcs

def load_sample_barcodes(barcode_set):
    try:
        return getattr(_bcs, barcode_set)
    except AttributeError:
        print("Unable to find barcode_set '{}'".format(barcode_set))
        sys.exit(1)

def fastq_iter(inpath):
    fastq_file = xopen.xopen(inpath)

    while True:
        lines = [fastq_file.readline().strip() for i in range(4)]
        if len(lines[-1]) < 1:
            break
        yield lines

def interleaved_fastq_iter(inpath1, inpath2):
    fastq1_iter = fastq_iter(inpath1)
    fastq2_iter = fastq_iter(inpath2)

    for reads1, reads2 in six.moves.zip_longest(fastq1_iter, fastq2_iter):
        read_name, comment = reads1[0].split(" ")
        sample_bc = comment.split(":")[-1]
        if "N" in sample_bc: continue
        yield sample_bc, reads1, reads2

def get_out_fastqs(sample_names, outdir, un=None):
    samples_to_out_fastqs = {}
    if un is not None:
        sample_names = sample_names.copy()
        sample_names.add(un)
        
    for sample_name in sample_names:
        samples_to_out_fastqs[sample_name] = []
        for i in [1,2]:
            outpath = os.path.join(outdir, "{}_{}.fastq.gz".format(sample_name, i))
            outfile = xopen.xopen(outpath, "w")

            samples_to_out_fastqs[sample_name].append(outfile)

    return samples_to_out_fastqs


def demultiplex(fastq1_path, fastq2_path, sample_barcodes_to_samples, gem_barcodes, outdir, un=None):
    gem_barcode_length = 16
    samples_to_out_fastqs = get_out_fastqs(set(sample_barcodes_to_samples.values()), outdir, un)

    total = 0
    good = 0
    bad_gem = 0
    bad_sample = 0
    bc_counts = collections.Counter()

    t_iter = tqdm.tqdm(interleaved_fastq_iter(fastq1_path, fastq2_path))
    for sample_bc, reads1, reads2 in t_iter:
        total += 1
        if total % int(1e6) == 0:
            t_iter.write("{:>13,} good:{:.1%} bad_gem:{:.1%} bad_sample:{:.1%}".format(total, good/float(total), bad_gem/float(total), bad_sample/float(total)))

        bc_counts[sample_bc] += 1
        if sample_bc not in sample_barcodes_to_samples:
            bad_sample += 1
            if un is None:
                continue

        gem_bc = reads1[1][:gem_barcode_length]
        if gem_bc not in gem_barcodes:
            bad_gem += 1
            continue

        out_reads1 = reads1[:]
        out_reads1[0] = out_reads1[0].split(" ")[0] + " BX:Z:" + gem_bc
        out_reads1[1] = out_reads1[1][gem_barcode_length:]
        out_reads1[3] = out_reads1[3][gem_barcode_length:]

        out_reads2 = reads2[:]
        out_reads2[0] = out_reads1[0]

        sample_name = sample_barcodes_to_samples.get(sample_bc, un)
        samples_to_out_fastqs[sample_name][0].write("\n".join(out_reads1)+"\n")
        samples_to_out_fastqs[sample_name][1].write("\n".join(out_reads2)+"\n")

        good += 1

    print("{:>13,} good:{:.1%} bad_gem:{:.1%} bad_sample:{:.1%}".format(total, good/float(total), bad_gem/float(total), bad_sample/float(total)))

    # Sanity check that the user specified all well-used sample barcodes
    # this is only a warning because the user may want to ignore a sample on purpose
    # eg only analyze 3/4 samples in a run
    for bc, count in bc_counts.items():
        if bc not in sample_barcodes_to_samples and count > (total / float(len(sample_barcodes_to_samples)) * 0.1):
            print("Got many counts of the following unspecified sample barcode: {} (n={})".format(bc, count))


def main():
    parser = argparse.ArgumentParser(description="Demultiplexes pooled 10x Chromium fastq files by sample barcode and moves "
                                                 "droplet (aka GEM) barcodes into fastq comment. "
                                                 "Note that a list of proper barcode names can be found in the "
                                                 "_bcs.py file distributed with this script.")
    
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("fastq1", help="read end 1 fastq file (ie the first read out of the read pairs); "
                                              "optionally, this can be a .gz, .xz or .bz2 compressed file")
    required_args.add_argument("fastq2", help="read end 2 fastq file")

    required_args.add_argument("-s", "--sample", action="append",
                               help="this tells us how to demultiplex reads by sample; must be of the format "
                                    "barcode=sample_name, where barcode is the 10x barcode name and sample_name "
                                    "is a label for this sample, (eg 'SI_GA_A1=strain211g'); this argument should "
                                    "be specified once for each sample")
    required_args.add_argument("-b", "--gem_barcodes", help="a file listing the proper droplet (aka gem) barcodes; "
                                                            "optionally, this can be a compressed file, for example "
                                                            "the '4M-with-alts-february-2016.txt.xz' file distributed "
                                                            "with this script which describes the chromium barcodes")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-o", "--outdir", help="output directory; default is the current working directory")
    optional_args.add_argument("-u", "--unmatched", help="an output file to put all reads without a barcode match")
    

    args = parser.parse_args()

    if args.outdir is None:
        args.outdir = os.getcwd()
    else:
        ensure_dir(args.outdir)
    print(args)

    if args.sample is None and not args.unmatched:
        raise Exception("Must specify at least one barcode=sample_name assignment using the '-s' option AND/OR an --unmatched file")
    if args.sample is None:
        args.sample = []
    if args.gem_barcodes is None:
        raise Exception("Must specify the file listing correct GEM barcodes using the '-b' option")
            
    sample_barcodes_to_samples = {}
    for sample in args.sample:
        barcode_set, sep, sample_name = sample.partition("=")
        assert len(sample_name) > 1, "invalid sample specifier '{}'; should be of format 'barcode_set=sample_name'".format(sample)

        sample_barcodes = load_sample_barcodes(barcode_set)
        for bc in sample_barcodes:
            assert not bc in sample_barcodes_to_samples
            sample_barcodes_to_samples[bc] = sample_name

    print(sample_barcodes_to_samples)

    print("Loading the list of proper droplet barcodes...")
    gem_barcodes = load_gem_barcodes(args.gem_barcodes)

    print("Running demultiplexing; this may take a while!")
    demultiplex(args.fastq1, args.fastq2, sample_barcodes_to_samples, gem_barcodes, args.outdir, args.unmatched)

if __name__ == '__main__':
    main()
