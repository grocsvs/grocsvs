GROC-SVs
--------

Genome-wide Reconstruction of Complex Structural Variants, or GROC-SVs, is a software pipeline for identifying structural variants, performing sequence assembly at the breakpoints, and reconstructing the complex structural variants using the long-fragment information from the 10x Genomics platform.


Installation
============

Prerequisites: the following programs must be installed prior to running GROC-SVs:

* `idba_ud <https://github.com/loneknightpy/idba>`_
* `samtools <http://www.htslib.org/download/>`_ version 1.0 or later
* `bwa-mem <https://github.com/lh3/bwa/releases>`_

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

1. run `longranger <http://support.10xgenomics.com/genome-exome/software>`_ on your 10x sequencing data
2. setup a ``configuration.json`` file describing your samples and your compute (eg cluster) setup
3. run grocsvs


1. Run longranger align
"""""""""""""""""""""""

GROC-SVs uses the read alignments produced by the 10x software pipeline, so ``longranger align`` must be run prior to ``grocsvs``. Optionally, the full 10x longranger pipeline may be run, which adds phasing information that GROC-SVs can include in its analysis.


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

Note that for SLURM clusters, the preferred batch submission method is using ``admiral``, which is installed by default with grocsvs. See the tumor-normal configuration in the examples directory.


3. Run GROC-SVs
"""""""""""""""

To run GROC-SVs, use the ``grocsvs /path/to/experiment/configuration.json`` command. If you are using a virtualenv, remember to run ``source grocsvs_env/bin/activate`` to activate the virtualenv prior to running ``grocsvs``. 

The output will be placed in the directory containing configuration (in this case, in ``/path/to/experiment/``), so make sure this filesystem has enough space for the analysis (~40GB per sample).

Note that the ``grocsvs`` command will continue running until all steps have completed. The ``grocsvs`` command itself is lightweight, and so can be run from a head node on your cluster.

Logging output for each step will be put in ``/path/to/experiment/logs``. The final results will be put in ``/path/to/experiment/results``.


Output
""""""

Final results of interest might be:

* ``results/MergeGenotypesStep/genotypes.tsv``: the structural variant calls, including coordinates, information on which samples are positive for each event, which events together form complex events, and some filtering information (eg blacklist annotations provided above, genome gaps, etc) to remove potential false-positives
* ``results/QCStep/qc_report.tsv``: some basic quality control statistics, including fragment lengths and number of barcodes per sample
* ``results/AssemblyStep/assembly.i``: the sequence assemblies for event ``i``; in this directory, ``contigs.sorted.bam`` contain the contigs aligned back to the reference genome (this file may be viewed with `IGV <https://www.broadinstitute.org/igv/>`_)
* ``results/FinalClusterSVsStep/edges.tsv``: full information relating breakpoints in complex structural variants




Troubleshooting
===============

The ``grocsvs /path/to/experiment/configuration.json`` command may be run multiple times to resume the pipeline.

If an error arises, the output from ``grocsvs`` or the log files may be informative.

Please submit issues on the `github page for grocsvs <https://github.com/grocsvs/grocsvs/issues>`_.

