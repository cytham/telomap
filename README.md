## Telomap - A tool for analyzing telobait-captured long-read telomere sequencing data

[![Build Status](https://app.travis-ci.com/cytham/telomap.svg?branch=main)](https://app.travis-ci.com/github/cytham/telomap)

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
| no_cores | Number of cores for subtelomere clustering |
| working_directory | Working directory |

### Generating capture oligo and barcode FASTA files
The capture oligo sequence refers to the sequence after the barcode on the telobait, while the barcode sequence refers to the start of the telobait until the end of the barcode.

See example below:
```
# Telobait 
      5'-ATCGANNNNNNNNGATGCCAGATGCACGGAGCA...-Biotin-3'
         |||||||||||||||||||||||||||||||||
3'-CCAATCTAGCTNNNNNNNNCTACGGTCTACGTGCCTCGT...-5'

where NNNNNNNN is the sample barcode

# capture_oligo.fa
>capture_oligo
GATGCCAGATGCACGGAGCA

# barcodes.fa
>sample1
ATCGACGGTTCAA
>sample2
ATCGAGCTGGATT
>sample3
ATCGATAACTCGG
```
#### Output

| Output file | Comment |
| :--- | :--- |
| ${sample}.telomap.tsv | Main output table |
| ${sample}.telomap.anchors.tsv | Table containing details of subtelomeric anchor sequences for chromosome-end mapping |

### Operating system

* Linux (x86_64 architecture, tested in Ubuntu 16.04)  

### Installation

```bash
git clone https://github.com/cytham/telomap.git 
cd telomap
pip install .
```

## Versioning

See [CHANGELOG](./CHANGELOG.txt)

## Citation

Tham CY, Poon L, Yan T, Koh JYP, Ramlee MK, Teoh VSI, Zhang S, Cai Y, Hong Z, Lee GS, Liu J, Song HW, Hwang WYK, Teh BT, Tan P, Xu L, Koh AS, Osato M, Li S. High-throughput telomere length measurement at nucleotide resolution using the PacBio high fidelity sequencing platform. Nat Commun. 2023 Jan 17;14(1):281. doi: 10.1038/s41467-023-35823-7. PMID: 36650155; PMCID: PMC9845338. https://www.nature.com/articles/s41467-023-35823-7

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

This project is licensed under GNU General Public License - see [LICENSE.txt](./LICENSE.txt) for details.

## Limitations

* Current chromosomal end mapping is incomplete
