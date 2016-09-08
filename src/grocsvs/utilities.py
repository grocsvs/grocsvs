import collections
import errno
import h5py
import os
import numpy
import pandas as pandas
import pyfaidx
import re
import string 
import sys

try:
    import cPickle as pickle
except ImportError:
    import pickle

assert pickle.HIGHEST_PROTOCOL >= 2, \
    "A relatively recent version of pickle is required"


def ensure_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

class cd: 
  def __init__(self, newPath):
    self.newPath = newPath

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)
    

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]    


def get_key(options_dict, key, type_=basestring, default="error", error_msg="configuration"):
    if default == "error" and key not in options_dict:
        print "CONFIG ERROR: {} key '{}' is missing".format(error_msg, key)
        sys.exit(1)
    value = options_dict.get(key, default)
    if type_ is not None and not isinstance(value, type_):
        print "CONFIG ERROR: {} key '{}' should be type '{}', not '{}'".format(
            error_msg, key, type_.__name__, type(value).__name__)
        sys.exit(1)
    return value


comp = string.maketrans('ATCGatcg','TAGCtagc')
def revcomp(seq):
    return seq[::-1].translate(comp)

    
###############################################################################
###############################################################################

# based on hdf5.py from 10X Genomics

LEVEL_GROUP = "_levels"

def has_levels(ds):
    """ Determine if a data column is leveled """
    if "levels" in ds.attrs.keys():
        return True

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            return True

    return False

def get_levels(ds):
    """ Get the level index for a dataset """
    if "levels" in ds.attrs.keys():
        return ds.attrs["levels"][:]

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            return level_ds[:]

    return None

def get_column_intersection(column_names, columns):
    if len(columns) > 0:
        column_names = sorted(list(set(columns) & set(column_names)))
        if len(column_names) == 0:
            raise Exception("No valid column specifications.")
    return column_names

def read_data_frame(fn, query_cols=[]):
    """ Load a pandas DataFrame from an HDF5 file. If a column list is 
    specified, only load the matching columns """

    with h5py.File(fn, "r") as f:
        column_names = f.attrs.get("column_names")
        column_names = get_column_intersection(column_names, query_cols)

        df = pandas.DataFrame()

        # Add the columns progressively to save memory
        for name in column_names:
            ds = f[name]
            if has_levels(ds):
                indices = ds[:]
                uniques = get_levels(ds)
                # This method of constructing of Categorical avoids copying the 
                # indices array which saves memory for big datasets
                df[name] = pandas.Categorical(indices, categories=uniques, 
                                              ordered=False, fastpath=True)
            else:
                df[name] = pandas.Series(ds[:])

        return df

###############################################################################
###############################################################################


def get_good_barcodes(fragments, proportion=0.90):
    """
    return the top barcodes which together comprise 90% of reads
    """
    read_counts = fragments.groupby("bc").sum()["num_reads"].copy()

    read_counts.sort_values(inplace=True, ascending=False)
    cutoff = proportion * read_counts.sum()
    cutoff = numpy.where(read_counts.cumsum() >= cutoff)[0][0]

    return sorted(read_counts.index[:cutoff])

def get_good_bc_count(step):
    sample_info = step.options.sample_info(step.sample.name)
    dataset_info = sample_info[step.dataset.id]
    good_bc_count = dataset_info["good_bc_count"]
    return good_bc_count
    

def cpu_count_physical():
    """
    tries to get the number of physical (ie not virtual) cores
    """
    try:
        import psutil
        return psutil.cpu_count(logical=False)
    except:
        import multiprocessing
        return multiprocessing.cpu_count()


def frags_overlap_same_chrom(frags, start, end):
    """
    get the fragments overlapping the interval [start, end], assuming
    all fragments in the input table are already on the correct chromosome
    """
    f = frags.loc[((frags["start_pos"] < start) & (frags["end_pos"] > start)) |
                  ((frags["start_pos"] < end) & (frags["end_pos"] > end))]
    return f



###############################################################################
###############################################################################




def plot_matrix_as_image(mat, x1=None, y1=None, x2=None, y2=None, maxval=None, main="", xlab="", ylab=""):
    """ adds the image to the current plot if x1, x2, y1 and y2 are defined;
    otherwise, create a new image with the dimensions of the matrix """
    from rpy2.robjects import r

    r.plot(numpy.array([0]),
           xlim=numpy.array([x1,x2]),
           ylim=numpy.array([y1,y2]),
           type="n", bty="n",
           main=main, xlab=xlab, ylab=ylab)

    if maxval is None:
        maxval = mat.max()
    r.rasterImage(r["as.raster"](mat, max=maxval), x1, y1, x2, y2)


