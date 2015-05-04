from distutils.core import setup

setup(
    name='dbb',
    version='1.0.4',
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['dbb', 'dbb.plots'],
    scripts=['bin/dbb'],
    package_data={'': ['data/*.txt']},
    url='http://pypi.python.org/pypi/dbb/',
    license='GPL3',
    description='Bin scaffolds into population genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.6.1",
        "scipy >= 0.10.1",
        "matplotlib >= 1.3.0",
        "pysam >= 0.7.4, <0.8.0"],
)
