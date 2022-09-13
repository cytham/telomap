Work in progress

## Telomap - A tool for analyzing telobait-captured long-read telomere sequencing data

[![Build Status](https://app.travis-ci.com/cytham/telomap.svg?branch=master)](https://app.travis-ci.com/github/cytham/telomap)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/telomap)](https://pypi.org/project/telomap/)
[![PyPI versions](https://img.shields.io/pypi/v/telomap)](https://pypi.org/project/telomap/)
[![Conda](https://img.shields.io/conda/v/bioconda/telomap)](https://anaconda.org/bioconda/telomap)
[![Github release](https://img.shields.io/github/v/release/cytham/telomap?include_prereleases)](../../releases)
[![PyPI license](https://img.shields.io/pypi/l/telomap)](./LICENSE.txt)

Telomap is a tool for downstream analysis on telobait-captured long-read telomere sequencing data.

### Basic information

* Identify long-reads that contain telomeric sequences
* De-multiplex barcoded samples
* Map telomeric reads to chromosome ends based on T2T-CHM13 assembly
* Analyze the distribution of non-canonical telomeric motifs (Telomere variant sequences (TVS))
* Create QC and analytical plots

## Getting Started

### Quick run

```bash
telomap [reads.fa] [capture_oligo.fa] [barcodes.fa] [data_type] [no_cores] [working_directory]
```

| Argument | Comment |
| :--- | :--- |
| reads.fa | Long-read file in fasta/fastq/bam/pacbio-bam formats |
| capture_oligo.fa | A FASTA file containing names and sequences of capture oligos |
| barcodes.fa | A FASTA file containing names and sequences of sample barcode |
| data_type | Type of long-read data (fasta/fastq/bam/pacbio-bam) |
| no_cores | Number of cores to specify to hasten subtelomere clustering |
| working_directory | Working directory |

#### Output

| Output file | Comment |
| :--- | :--- |
| ${sample}.telomap.tsv | Main output table |
| ${sample}.telomap.anchors.tsv | Table containing details of subtelomeric anchor sequences for chromosome-end mapping |

For more information, see [wiki](https://github.com/cytham/telomap/wiki).

### Operating system

* Linux (x86_64 architecture, tested in Ubuntu 16.04)  

### Installation

There are three ways to install Telomap:

#### Option 1: Conda (Recommended)

```bash
conda install -c bioconda telomap
```

#### Option 2: PyPI (See dependencies below)

```bash
pip install telomap
```

#### Option 3: GitHub (See dependencies below)

```bash
git clone https://github.com/cytham/telomap.git 
cd telomap
pip install .
```

## Versioning

See [CHANGELOG](./CHANGELOG.txt)

## Citation

Not available yet

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

This project is licensed under GNU General Public License - see [LICENSE.txt](./LICENSE.txt) for details.

## Limitations

* Current chromosomal end mapping is incomplete
