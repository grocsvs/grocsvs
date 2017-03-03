from setuptools import setup, find_packages

def get_version(path):
    """ Parse the version number variable __version__ from a script. """
    import re
    string = open(path).read()
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(
    name = 'grocsvs',
    version = get_version("src/grocsvs/__init__.py"),

    packages = find_packages('src'),
    package_dir = {"": "src"},

    entry_points = {
        'console_scripts' : ["grocsvs = grocsvs.main:main"]
    },

    install_requires = ["admiral", "h5py", "networkx", "pandas", "pybedtools", 
                        "pyfaidx", "pysam>=0.10.0", "scipy", "ipython-cluster-helper",
                        "pygraphviz"],

)
