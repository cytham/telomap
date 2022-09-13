from setuptools import setup, find_packages
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

exec(open("telomap/version.py").read())

setup(
    name='telomap',
    version=__version__,
    packages=find_packages(),
    package_data={'telomap.ref': ['*.fa']},
    include_package_data=True,
    scripts=['telomap/telomap'],
    url='https://github.com/cytham/telomap',
    download_url='https://github.com/cytham/telomap/releases',
    license='gpl-3.0',
    author='Tham Cheng Yong',
    author_email='chengyong.tham@u.nus.edu',
    description='A tool for analyzing telobait-captured long-read telomere sequencing data',
    keywords=['telomap', 'telomere', 'pacbio', 'ont', 'oxford nanopore', 'long read'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['pandas>=1.4.2', 'numpy>=1.22.4', 'scipy>=1.8.1', 'matplotlib>=3.5.2',
                      'pysam>=0.19.0', 'seaborn>=0.11.2', 'natsort>=8.1.0', 'biopython>=1.79', 'scikit-learn>=1.1.2'],
    python_requires='>=3.6',
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
