import logging
import os
import pandas
import pysam
import cStringIO as StringIO
import sys

from grocsvs import utilities
from grocsvs.utilities import get_key

class Dataset(object):
    def serialize(self):
        d = self.__dict__.copy()
        d["type"] = self.__class__.__name__
        del d["sample"]
        return d

    @staticmethod
    def deserialize(sample, dict_):
        if not isinstance(dict_, dict):
            print "samples must be of type 'dict', not '{}': '{}'".format(type(dict_).__name__, dict_)
            sys.exit(1)
            
        dict_ = dict_.copy()
        dataset_types = [TenXDataset, ShortFragDataset, MatePairDataset]
        dataset_types = dict((x.__name__, x) for x in dataset_types)

        type_ = get_key(dict_, "type", error_msg="sample") # just for type-checking
        if not type_ in dataset_types:
            print "ERROR: Sample type must be one of '{}'; got '{}' instead".format(dataset_types.keys(), type_)
            sys.exit(1)
            
        dataset_class = dataset_types[dict_.pop("type")]

        dict_["sample"] = sample
        #try:
        return dataset_class(**dict_)
        #except TypeError:
        #    print "MISSING FIELD FOR SAMPLE:", sample.name, dataset_class
        #    print "  Fields provided:", dataset_class.__class__.__name__, dict_
        #    sys.exit(1)



class TenXDataset(Dataset):
    def __init__(self, **kwdargs):#sample, bam, fragments, phased_fragments, id, sorted_fastqs=None):
        self.sample = get_key(kwdargs, "sample", None, error_msg="TenXDataset")
        self.bam = os.path.realpath(get_key(kwdargs, "bam", error_msg="TenXDataset"))
        #self.fragments = get_key(kwdargs, "fragments", error_msg="TenXDataset")
        #self.phased_fragments = get_key(kwdargs, "phased_fragments", error_msg="TenXDataset")
        #self.sorted_fastqs = get_key(kwdargs, "sorted_fastqs", default=None, error_msg="TenXDataset")

        self.id = get_key(kwdargs, "id", error_msg="TenXDataset")

        self.validate()
        
        
    def validate(self):
        assert os.path.exists(self.bam), "missing bam file '{}' for sample '{}' and dataset '{}'".format(
                                         self.bam, self.sample.name, self.id)


    # @staticmethod
    # def from_longranger_dir(self, longranger_dir):
    #     fragments = os.path.join(longranger_dir,
    #         "PHASER_SVCALLER_CS/PHASER_SVCALLER/_REPORTER/"
    #         "REPORT_SINGLE_PARTITION/fork0/files/fragments.h5")

    #     bam = os.path.join(longranger_dir,
    #         "PHASER_SVCALLER_CS/PHASER_SVCALLER/ATTACH_PHASING/"
    #         "fork0/files/phased_possorted_bam.bam")

    #     phased_fragments = os.path.join(longranger_dir,
    #         "10XSARCOMAC1/PHASER_SVCALLER_CS/PHASER_SVCALLER/"
    #         "_SNPINDEL_PHASER/PHASE_SNPINDELS/fork0/files/"
    #         "fragment_phasing.tsv.gz")

    #     self.validate()

    #     return TenXDataset(bam, fragments, phased_fragments)

    # def load_phased_fragments(self, chrom=None, start=None, end=None):

    #     columns = ["chrom", "start_pos", "end_pos", "phase_set", "ps_start", 
    #                "ps_end", "bc", "h0", "h1", "hmix", "unkn"]

    #     try:
    #         tabix = pysam.TabixFile(self.phased_fragments)
    #         s = StringIO.StringIO("\n".join(tabix.fetch(chrom, start, end)))
    #         frags = pandas.read_table(s)
    #         frags.columns = columns
    #     except (IOError, ValueError):
    #         frags = pandas.DataFrame(columns=columns)
        
    #     return frags

    # def load_fragments(self, chrom=None, start=None, end=None):
    #     tabix = pysam.TabixFile()

        # try:
        #     fragments = utilities.read_data_frame(self.fragments)
        #     goodbcs = utilities.get_good_barcodes(fragments)
        #     fragments = fragments.loc[fragments["bc"].isin(goodbcs)]
        #     # fragments = fragments.loc[fragments["num_reads"]>5]
        #     if chrom is not None:
        #         fragments = fragments.loc[fragments["chrom"]==chrom]

        #     return fragments
        # except:
        #     logging.exception("Unable to load fragments from fragments file "
        #         "'{}'".format(self.fragments))
        #     raise


class ShortFragDataset(Dataset):
    def __init__(self, sample, bam, id):
        self.sample = sample
        self.bam = os.path.realpath(bam)
        self.id = id


class MatePairDataset(Dataset):
    def __init__(self, sample, bam, id):
        self.sample = sample
        self.bam = os.path.realpath(bam)
        self.id = id
