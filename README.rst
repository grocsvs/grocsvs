GROC-SVs
--------

Genome-wide Reconstruction of Complex Structural Variants, or GROC-SVs, is a software pipeline for identifying structural variants, performing sequence assembly at the breakpoints, and reconstructing the complex structural variants using the long-fragment information from the 10x Genomics platform.


Installation
============

Prerequisites: the following programs must be installed prior to running GROC-SVs:

* `idba_ud <https://github.com/grocsvs/idba/releases/tag/1.1.3g1>`_ -- please use this version, as the version distributed by the original author does not support paired reads longer than 128 bp
* `samtools and htslib <http://www.htslib.org/download/>`_ -- version 1.0 or later of the ``samtools``, ``bgzip``, and ``tabix`` programs must all be in your ``$PATH``
* `bwa-mem <https://github.com/lh3/bwa/releases>`_
* `graphviz <http://www.graphviz.org/Download..php>`_ (optional)

We recommend setting up a `virtualenv <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ prior to installing GROC-SVs (or using `virtualenvwrapper <http://www.simononsoftware.com/virtualenv-tutorial-part-2/>`_):

.. code-block:: bash

    sudo pip install virtualenv
    virtualenv grocsvs_env

The virtualenv can be activated by running the following command:

.. code-block:: bash

    source grocsvs_env/bin/activate

Then, to install grocsvs:

.. code-block:: bash

    cd /path/to/grocsvs
    pip install .

To test that GROC-SVs is installed correctly, you can simply run ``grocsvs`` from the commandline, which should show help text without any error messages.


Running GROC-SVs
================

Overview:

1. extract barcodes and align your 10x sequencing data to the reference genome
2. setup a ``configuration.json`` file describing your samples and your compute (eg cluster) setup
3. run grocsvs


1. Extract barcodes and align reads
"""""""""""""""""""""""""""""""""""

There are two options to align 10x Genomics data to the reference genome for downstream use by GROC-SVs. The simplest option is to use the accompanying `simple_demux_map`_ script followed by read alignment using ``bwa mem``. Note that this will extract the 10x droplet barcodes for use by GROC-SVs (optionally demultiplexing pooled samples) but does not perform `barcode-aware read alignment <http://genome.cshlp.org/content/25/10/1570>`_.

.. _simple_demux_map: simple_demux_map/

The second option is to use the 10x Genomics `longranger align <http://support.10xgenomics.com/genome-exome/software>`_ pipeline, which can optionally perform barcode-aware alignment. While not necessary, the full 10x longranger pipeline may be run, which adds phasing information that GROC-SVs can include in its analysis.


2. Setup a configuration file
"""""""""""""""""""""""""""""

Because GROC-SVs can analyze multiple samples jointly, and each sample can involve multiple input files, GROC-SVs uses a configuration file to specify inputs and settings for each run. See the examples directory for an example ``configuration.json`` file.

The configuration file is in the `JSON <http://www.json.org>`_ format, and contains the following three parts:

1. reference genome paths
2. sample information
3. compute cluster settings

**Reference genome** The following paths must be defined:

* ``ref_fasta``: path to the reference genome FASTA file (hg19.fa, GRCh38.fa, etc)
* ``bwa_index``: path to the genome index used by bwa-mem

In addition, the following optional paths may be specified:

* ``blacklists``: a list of paths of blacklist regions, either `bed <https://genome.ucsc.edu/FAQ/FAQformat.html>`_ or `bedpe <http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format>`_ format
* ``binaries``: a hash containing paths for any of the following binaries: ``idba_ud``, ``bwa``, ``samtools``. If these are in your ``$PATH``, there is no need to specify them in the configuration file.

**Sample information** Samples is a hash with key specifying sample name, and the value is a list of datasets. Each sample must have one 10x dataset specified, and may optionally specify a separate standard Illumina short-frag library or a mate-pair library (these are used for validation and comparison only).

Each dataset is defined as a hash. 10x datasets must define the following items:

* ``bam``: the path of the bam file produced by longranger. From the root longranger output directory, this is typically the file ``$ROOT/outs/possorted_bam.bam``.
* ``id``: this is a name used to identify the dataset; typically, it would be something like "``sample_name_10x``"
* ``type``: this should be the string "``TenXDataset``"

**Tumor/normal, trio and other multi-sample analyses** GROC-SVs performs variant calling jointly on all samples specified in the configuration.json file, and no additional arguments are required to indicate the biological meaning of the sample. See the Output section below for descriptions of the "genotypes.tsv" and "classes.txt" files, which can be filtered in order to obtain events that are somatic (ie private to the tumor, and not present in the germline sample) or de novo in the child (ie private to the child and not present in either parent).

**Compute cluster settings** This has defines the compute environment being used to perform the analysis. A standard cluster setup looks like this:

.. code-block:: json

    "cluster_settings": {
        "cluster_type": "IPCluster",
        "processes": 128,
        "cluster_options": {
            "scheduler": "slurm",
            "queue": "normal",
            "extra_params": {"mem":16}
        }
    }

Where ``processes`` specifies the maximum number of separate jobs (1 processor per job) to allow. ``scheduler`` may be any of the clusters supported by `ipython-cluster-helper <https://github.com/roryk/ipython-cluster-helper>`_. Currently, these are Platform LSF ("lsf"), Sun Grid Engine ("sge"), Torque ("torque"), and SLURM ("slurm").

To run in parallel on a single machine, use ``cluster_type":"multiprocessing"`` and specify the desired number of ``processes``.

To override the cluster options in the configuration.json file, use ``--local`` to specify single-core mode or ``--multiprocessing`` to specify running in parallel using all cores on a single machine.

3. Run GROC-SVs
"""""""""""""""

To run GROC-SVs, use the ``grocsvs /path/to/experiment/configuration.json`` command. If you are using a virtualenv, remember to run ``source grocsvs_env/bin/activate`` to activate the virtualenv prior to running ``grocsvs``. 

The output will be placed in the directory containing configuration (in this case, in ``/path/to/experiment/``), so make sure this filesystem has enough space for the analysis (~40GB per sample). GROC-SVs typically requires about 12-16 GB of memory in order to run, though this depends on your samples. If you have less than 16 GB of memory available on your machine, a warning will be output but the pipeline will continue to run as best as it can.

Note that the ``grocsvs`` command will continue running until all steps have completed. The ``grocsvs`` command itself is lightweight, and so can be run from a head node on your cluster.

Logging output for each step will be put in ``/path/to/experiment/logs``. The final results will be put in ``/path/to/experiment/results``.


Output
""""""

Final results of interest might be:

* ``results/MergeGenotypesStep/genotypes.tsv``: the structural variant calls, including coordinates, information on which samples are positive for each event, which events together form complex events, and some filtering information (eg blacklist annotations provided above, genome gaps, etc) to remove potential false-positives
* ``results/QCStep/qc_report.tsv``: some basic quality control statistics, including fragment lengths and number of barcodes per sample
* ``results/AssemblyStep/assembly.i``: the sequence assemblies for event ``i``; in this directory, ``contigs.sorted.bam`` contain the contigs aligned back to the reference genome (this file may be viewed with `IGV <https://www.broadinstitute.org/igv/>`_)
* ``results/FinalClusterSVsStep/edges.tsv``: full information relating breakpoints in complex structural variants
* ``results/PostprocessingStep/classes.txt``: this file includes a simple presence/absence call for each structural variant for each sample, denoted as a 0 for absence and a 1 for presence. For example, if your tumor sample were the first sample, and the matched normal sample were the second sample, a "10" would indicate a somatic event and a "11" would indicate a germline event. These classes are determined using a simple allele-frequency cutoff which in our experience has been quite robust. More statistically motivated filters can be established by filtering on the p-values for each sample, which are indicated in this file as "sarcoma_p_resampling" if your sample name were "sarcoma" (note that missing p-values should be treated as 1).


Docker (and example dataset)
============================

A docker image is available for grocsvs. If you wish to download and run grocsvs on an example dataset (~1.3GB required), you can run the following commands:

.. code-block:: bash
    
    # use 'curl -O' if you're on a mac without wget
    wget http://mendel.stanford.edu/public/noah/grocsvs_example.tar.gz 
    tar -xzf grocsvs_example.tar.gz

Assuming `docker <https://docs.docker.com/engine/installation/>`_ is installed, the following command can be used to analyze the example data from within docker (make sure you are in the same directory where you downloaded and extracted grocsvs_example.tar.gz):

.. code-block:: bash

    docker run -v `pwd`:/data -w /data/grocsvs_example/ grocsvs/grocsvs-docker grocsvs configuration.json --local

This requires ~16GB of memory to run and will take ~1 hour to complete. If you are running docker for Mac, please make sure that your virtual machine has access to at least 16GB of memory.

The output can be found in ``grocsvs_example/results``.

Comparison to Long Ranger Pipeline
==================================

Briefly, GROC-SVs was designed to detect and characterize complex structural variants such as those frequently found in cancer or in orphan diseases. The Long Ranger software available from 10x Genomics can also perform SV detection using inferred long-fragment sequence information, but is more well-suited to analysis of individual germline genomes. Note that both GROC-SVs and Long Ranger are being actively developed, and so some features may migrate between packages.

GROC-SVs:

* performs sequence assembly of structural variants
* reconstructs large-scale complex structural variants
* is designed for multi-sample analyses (tumor/normal, or trios) - this is important when identifying somatic or de novo germline events, as analyzing multiple samples separately can result in false negative calls in the control or parent samples


Troubleshooting
===============

The ``grocsvs /path/to/experiment/configuration.json`` command may be run multiple times to resume the pipeline.

If you are having trouble installing or running grocsvs, the docker file (see above) may help you diagnose the issue.

If an error arises, the output from ``grocsvs`` or the log files may be informative.

**ShortSequence: Sequence is too long.** If you get this error during assembly, please make sure you are using `the grocsvs fork of idba_ud <https://github.com/grocsvs/idba/releases/tag/1.1.3g1>`_.


Please submit issues on the `github page for grocsvs <https://github.com/grocsvs/grocsvs/issues>`_.

