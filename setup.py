from setuptools import setup, find_packages
from setuptools.command.install import install as InstallCommand

class Install(InstallCommand):
    """ Customized setuptools install command which uses pip. """

    def run(self, *args, **kwargs):
        import pip
        pip.main(['install', '.'])
        InstallCommand.run(self, *args, **kwargs)

setup(
    name="4D-CHAINS",
    author="Thomas Evangelidis",
    author_email="tevang3@gmail.com",
    maintainer="Thomas Evangelidis",
    maintainer_email="tevang3@gmail.com",
    description="\n*** 4D-CHAINS: complete protein backbone and sidechain chemical shift assignment from two 4D NMR spectra. ***\n",
    long_description="To be added...",
    url="Under construction.",
    license="4D-CHAINS is free ONLY for non-commercial usage.",
    version="1.0",
    platforms="Unix, Mac (to be supported)",
    cmdclass={
        'install': Install,
    },
    packages=find_packages(where='.',  exclude=()),
    #package_dir={'':'dev'},
    install_requires=['scipy', 'numpy', 'argparse', 'cluster', 'scoop', 'tabulate', 'ordereddict', 'ete3', 'multiprocessing']
)

# This is the full list of the required python packages. Those that come along the standard python2.7 installation were excluded from above:
# ['scipy', 'numpy', 'sys', 're', 'os', 'cPickle', 'traceback', 'shutil', 'bz2', 'math', 'argparse', 'cluster', 'itertools', 'operator', 'scoop', 'tabulate', 'ordereddict', 'collections', 'ftplib', 'ete3', 'time', 'subprocess', 'multiprocessing', 'copy', 'csv']
