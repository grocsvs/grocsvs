Simple 10x Chromium Demultiplexing and Mapping
----------------------------------------------

The following are a set of scripts and directions for simplified demultiplexing and mapping of 10x Genomics Chromium data. The input to this pipeline is a standard fastq file produced by your sequencing facility -- this replaces the 10x Long Ranger bcl conversion pipeline (aka ``longranger mkfastq``) with a far simpler 2-step pipeline.


Demultiplexing
==============

Use the demux.py script to extract the inline droplet barcodes from the Chromium fastq files, optionally demultiplexing samples from multiple different Chromium preps (eg different biological samples prepared separately then pooled into a single Illumina lane).

There are several required arguments:

- the ``fastq1`` and ``fastq2`` files you've obtained from the sequencer/sequencing core
- the ``-s`` argument is used to specify which sample barcodes correspond to which biological samples; each Chromium barcode is actually made up of four different sequences, so just use the name of the barcode set you received from 10x, for example ``SI_GA_C1``; the list of possible options can be viewed in the accompanying ``_bcs.py`` file; you'll also specify the sample label for each sample here, for example ``SI_GA_C1=sample1``
- the ``-b`` argument points to a file listing the valid Chromium droplet barcodes used (these are the "GEM" or gel emulsion droplet barcodes, which in Chromium are at the very beginning of the first read); you can just specify the accompanying ``4M-with-alts-february-2016.txt.xz`` file unless you have a newer set of barcodes

If the ``-o`` option specifying the output directory is not specified, the demultiplexed fastq files will be placed in the current working directory.

Please note that, in order to run ``demux.py`` you need to either be in the ``grocsvs/simple_demux_map`` directory or to add that directory to your ``$PYTHONPATH`` environment variable.

For example, you might run the following command:

.. code-block:: bash

    cd .../path/to/grocsvs/simple_demux_map
    python demux.py -s SI_GA_A1=SampleA -s SI_GA_B1=SampleB -b 4M-with-alts-february-2016.txt.xz -o results reads_1.fastq.gz reads_2.fastq.gz

This will create ``SampleA_1.fastq.gz``, ``SampleA_2.fastq.gz``, ``SampleB_1.fastq.gz`` and ``SampleB_2.fastq.gz`` in the ``results`` directory.

Mapping
=======

After extracting the droplet barcodes and demultiplexing, you can run the standard `bwa mem <https://github.com/lh3/bwa>`_ software on each pair of fastq files using the ``-C`` option which transfers the fastq comment (which is where we put the droplet barcodes in the first step above) into the bam tag.

If you ran the ``demux.py`` code above, you could use the following command to map the SampleA fastqs, sorting the resulting output and converting to bam format:

.. code-block:: bash

    bwa mem -C /path/to/index.fa results/SampleA_1.fastq.gz SampleA_2.fastq.gz | samtools -bS - | samtools sort -T results/temp_sorting -o results/SampleA.sorted.bam
    samtools index results/SampleA.sorted.bam

The resulting bam file is valid input to the `GROC-SVs <https://github.com/grocsvs/grocsvs>`_ pipeline.
