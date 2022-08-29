
# WORK IN PROGRESS

## Telomap - A tool for analyzing telobait-captured long-read telomere sequencing data'

Telomap is a tool for handling telobait-captured long-read telomere sequencing data.

### Basic information:
* Identify long-reads that contain telomeric sequences
* De-multiplex barcoded samples
* Map telomeric reads to chromosome ends based on T2T-CHM13 assembly
* Analyze the distribution of non-canonical telomeric motifs (Telomere variant sequences (TVS))
* Create QC and analytical plots

## Getting Started

### Quick run

```
telomap [reads.fa] [capture_oligo.fa] [barcodes.fa] [data_type] [no_cores] [working_directory]
```

| Argument | Comment |
| :--- | :--- |
| reads | Long-read file in fasta/fastq/bam/pacbio-bam formats |
| capture_oligo.fa | A FASTA file containing names and sequences of capture oligos |
| barcodes.fa | A FASTA file containing names and sequences of sample barcode |
| data_type | Type of long-read data (fasta/fastq/bam/pacbio-bam) |
| no_cores | Number of cores to specify to hasten subtelomere clustering|
| working_directory | Working directory |

#### Output
| Output file | Comment |
| :--- | :--- |
| ${sample}.telomap.tsv | Main output table |
| ${sample}.telomap.anchors.tsv | Table containing details of subtelomeric anchor sequences for chromosome-end mapping |

For more information, see [wiki](https://github.com/cytham/telomap/wiki).

### Operating system: 
* Linux (x86_64 architecture, tested in Ubuntu 16.04)  

### Installation:
There are three ways to install Telomap:
#### Option 1: Conda (Recommended)
```
conda install -c bioconda telomap
```
#### Option 2: PyPI (See dependencies below)
```
pip install telomap
```
#### Option 3: GitHub (See dependencies below)
```
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
