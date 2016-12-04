import collections
import json
import os
import sys

from grocsvs.reference import Reference
from grocsvs import utilities
from grocsvs.utilities import get_key
from grocsvs.stages import constants
from grocsvs import datasets as svdatasets


class Sample(object):
    def __init__(self, name):
        self.name = name
        self.datasets = []
        if len(self.datasets) != len(set(d.id for d in self.datasets)):
            raise Exception("Dataset IDs must be unique")

    def serialize(self):
        return [dataset.serialize() for dataset in self.datasets]

    def get_10x_dataset(self):
        """
        This assumes that there is a single 10X dataset for each sample
        """
        tenx_datasets = []
        for dataset in self.datasets:
            if isinstance(dataset, svdatasets.TenXDataset):
                tenx_datasets.append(dataset)
        if len(tenx_datasets) == 0:
            raise Exception("No 10X dataset specified for sample {}".format(self.name))
        if len(tenx_datasets) > 1:
            raise Exception("Currently only support one 10X dataset per sample {}".format(self.name))
        return tenx_datasets[0]

    def get_shortfrag_dataset(self):
        """
        This assumes that there is a single short-frag dataset for each sample
        """
        for dataset in self.datasets:
            if isinstance(dataset, svdatasets.ShortFragDataset):
                return dataset
        raise Exception("No short-frag dataset specified for sample {}".format(self.name))


class ClusterSettings(object):
    def __init__(self):
        self.cluster_type = "local"
        self.processes = utilities.cpu_count_physical()
        self.cluster_options = {}

    @staticmethod
    def deserialize(options_dict):
        settings = ClusterSettings()

        if "processes" in options_dict:
            settings.processes = options_dict["processes"]
        if "cluster_type" in options_dict:
            settings.cluster_type = options_dict["cluster_type"]
        if "cluster_options" in options_dict:
            settings.cluster_options = options_dict["cluster_options"]

        return settings

    def serialize(self):
        return {
            "processes": self.processes,
            "cluster_type": self.cluster_type,
            "cluster_options": self.cluster_options
        }



    
class Options(object):
    def __init__(self, options_path, debug=False):
        self.options_path = options_path
        self._output_dir = os.path.dirname(self.options_path)

        self.samples = {}
        self.ref_fasta = None
        self.bwa_index = None
        self.binaries = {}
        self.blacklists = []
        self._reference = None
        self._constants = None

        self.cluster_settings = ClusterSettings()

        self.debug = debug


    def serialize(self, ):
        samples = dict((name, self.samples[name].serialize())
                       for name in self.samples)

        d = {"samples": samples,
             "ref_fasta": self.ref_fasta,
             "blacklists": self.blacklists,
             "cluster_settings": self.cluster_settings.serialize(),
             "bwa_index": self.bwa_index,
             "binaries": self.binaries
        }

        return d

    @staticmethod
    def deserialize(options_dict, options_path):
        samples = collections.OrderedDict()
        sample_names = get_key(options_dict, "samples", dict)
        sorted_names = sorted(sample_names.keys(), key=utilities.natural_sort_key)

        for sample_name in sorted_names:
            sample_info = get_key(options_dict["samples"], sample_name, list)
            cursample = Sample(sample_name)
            for dataset in sample_info:
                cursample.datasets.append(
                    svdatasets.Dataset.deserialize(cursample, dataset))

            cursample.get_10x_dataset() # make sure a 10X dataset was defined
            
            samples[sample_name] = cursample

        options = Options(options_path)
        options.samples = samples
        options.ref_fasta = get_key(options_dict, "ref_fasta")
        options.blacklists = get_key(options_dict, "blacklists", list, default=[])

        options.bwa_index = get_key(options_dict, "bwa_index")
        options.binaries = get_key(options_dict, "binaries", dict, default={})

        options.cluster_settings = ClusterSettings.deserialize(
            options_dict.get("cluster_settings", {}))

        return options

    def iter_10xdatasets(self):
        for sample_name, sample in self.samples.items():
            for dataset in sample.datasets:
                if isinstance(dataset, svdatasets.TenXDataset):
                    yield sample, dataset
    @property
    def output_dir(self):
        return self._output_dir
    
    @property
    def results_dir(self):
        return os.path.join(self.output_dir, "results")

    @property
    def working_dir(self):
        return os.path.join(self.output_dir, "working")

    @property
    def log_dir(self):
        return os.path.join(self.output_dir, "logs")

    @property
    def reference(self):
        if self._reference is None:
            self._reference = Reference(self.ref_fasta, self.debug)
        return self._reference
    
    @property
    def constants(self):
        if self._constants is None:
            self._constants = constants.get_constants(self)
        return self._constants
    
    def sample_info(self, sample_name):
        """
        Gets some pre-calculated values for the given sample (eg frag length 
        distributions)
        """
        info_path = os.path.join(
            self.results_dir, 
            "SampleInfoStep",
            "sample_info.{}.pickle".format(sample_name))

        if not os.path.exists(info_path):
            print info_path
            raise Exception("SampleInfoStep needs to be run first!")

        return utilities.pickle.load(open(info_path))

    def binary(self, name):
        """
        Checks to see if a path has been specified for an external binary,
        otherwise just return the name of the binary to try running it
        if it's in $PATH
        """
        
        bin_path = self.binaries.get(name, name)
        if utilities.which(bin_path) is None:
            raise Exception("Failed to locate binary '{}'; please make sure it is in ".format(name) + 
                            "your $PATH or add it to the configuration.json file")
        return bin_path


    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, mode=True):
        self._reference = None
        self._debug = mode

    def __str__(self):
        d = self.serialize()
        d["debug"] = self.debug
        return json.dumps(d, sort_keys=True, indent=4)

    def __getstate__(self):
        """
        allows pickling of Options instances, necessary for ipyparallel
        """
        state = self.__dict__.copy()
        state["_reference"] = None
        state["_constants"] = None

        return state


def validate_options(options):
    for sample in options.samples.values():
        tenx_count = 0
        for dataset in sample.datasets:
            if isinstance(dataset, svdatasets.TenXDataset):
                tenx_count += 1

        if tenx_count > 1:
            raise Exception(
                "Not yet implemented -- multiple 10X datasets per sample")
