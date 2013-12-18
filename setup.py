from distutils.core import setup

setup(
    name='MetaDBSCAN',
    version='0.0.1',
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['metadbscan', 'metadbscan.plots'],
    scripts=['bin/metadbscan'],
    url='http://pypi.python.org/pypi/metadbscan/',
    license='GPL3',
    description='Bin contigs into population genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.6.1",
        "scipy >= 0.10.1",
        "matplotlib >= 1.1.0",
        "pysam >= 0.7.4",
        "sqlite3 >= 2.6.0"],
)
