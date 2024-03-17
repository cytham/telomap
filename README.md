## New feature: Telomap version 0.1.0 onwards now supports input of whole-genome sequencing (WGS) long-read data

## Telomap - Telomere capture tool for WGS or telobait long-read sequencing data

[![Build Status](https://app.travis-ci.com/cytham/telomap.svg?branch=main)](https://app.travis-ci.com/github/cytham/telomap)

Telomap is a tool to identify and analyse telomeric reads within WGS or telobait-captured long-read sequencing data.

### Basic information

* Identify long-reads that contain telomeric sequences
* De-multiplex barcoded samples (for telobait-capture data)
* Map telomeric reads to chromosome ends based on T2T-CHM13 assembly
* Analyze the distribution of non-canonical telomeric motifs (Telomere variant sequences (TVS))
* Generate QC and analytic plots

## Getting Started

### Installation

```bash
git clone https://github.com/cytham/telomap.git 
cd telomap
pip install .
```

### Quick run

```bash
telomap [options] [RUN_MODE] [READS] [WORK_DIRECTORY]

# For WGS mode with 24 threads
telomap -t 24 wgs /path/to/reads.fastq.gz /path/to/work_dir

# For telobait mode with 24 threads
telomap -t 24 -c /path/to/capture_oligo.fa -b /path/to/barcodes.fa telobait /path/to/reads.fastq.gz /path/to/work_dir
# See below for instructions to create capture_oligo.fa and barcodes.fa files
```

### Full usage

```
usage: telomap [options] [RUN_MODE] [READS] [WORK_DIRECTORY]

Telomap is a tool for analyzing telomeric reads from WGS or telobait-capture long-read sequencing data

required arguments:
  [RUN_MODE]            run modes:
                        wgs (whole genome sequencing mode)
                        telobait (telobait capture mode)
  [READS]               path to input reads (fasta/fastq/bam/pacbio-bam)
  [WORK_DIRECTORY]      path to work directory

optional arguments:
  -c path, --capoligo path
                        path to capture oligo fasta file
  -b path, --barcodes path
                        path to barcodes fasta file
  -n str, --name str    name prefix for output files [sample]
  -m str, --motif str   telomeric motif sequence [TTAGGG]
  --oligoscore float    minimum alignment score fraction required
                        for capture oligo sequence match [1]
  --barscore float      minimum alignment score fraction required
                        for barcode sequence match [1]
  -t int, --threads int
                        specify number of threads [1]
  -v, --version         print version
  -q, --quiet           hide verbose
  -h, --help            show this help message and exit
```

### Generating capture oligo and barcode FASTA files
The capture oligo sequence refers to the sequence after the barcode on the telobait, while the barcode sequence refers to the start of the telobait until the end of the barcode.

See example below:
```
# Telobait 
      5'-ATCGANNNNNNNNGATGCCAGATGCACGGAGCA...-Biotin-3'
         |||||||||||||||||||||||||||||||||
3'-CCAATCTAGCTNNNNNNNNCTACGGTCTACGTGCCTCGT...-5'

where NNNNNNNN is the sample barcode

# The capture_oligo.fa should contain the following sequence
>capture_oligo
GATGCCAGATGCACGGAGCA

# The barcodes.fa should contain the following sequences
>sample1
ATCGACGGTTCAA
>sample2
ATCGAGCTGGATT
>sample3
ATCGATAACTCGG
```

### Output

| Output | Comment |
| :--- | :--- |
| sample.telomap.tsv | Main output table containing per read information |
| sample.telomap.anchors.tsv | Table containing details of subtelomeric anchor sequences for chromosome-end mapping |
| plots | Directory containg analytic plots |
| plots_QC | Directory containing QC plots |

#### Main output table information (sample.telomap.tsv)

| # | Column name | Comment |
| :--- | :--- | :--- |
| 1 | RNAME | Read ID |
| 2 | OLIGO | Name of capture oligo detected (will be NA for WGS mode) |
| 3 | BARCODE | Detected sample barcode (will be NA for WGS mode) |
| 4 | STRAND | Strand of read with reference to motif |
| 5 | MOTIF | Telomeric motif detected |
| 6 | RLEN | Read length |
| 7 | RAW_TELOMERE_LEN | Length of telomeric region in bases demarcated by TELOMERE_START and TELOMERE_END (See below) |
| 8 | TELOMERE_LEN | Length of telomeric region in bases after excluding all non-telomeric sequences within RAW_TELOMERE_LEN (See below)|
| 9 | TELOMERE_END_MOTIF | Last six 3' telomeric bases within telomeric region on read |
| 10 | TRF_COUNT | Number of TRF binding motifs detected (5'-TTAGGGTTA-3') |
| 11 | CHROM | Detected chromosomal end or unmappable sub-telomeric cluster (i.e. begins with "U") |
| 12 | TELOMERE_START | Refers to the 1-based position on a read at the 5′ end of the telomeric region as defined by containing at least two consecutive telomeric repeats (5′-TTAGGGTTAGGG-3′) |
| 13 | TELOMERE_END | Refers to the 1-based position on a read at the 3' end of the telomeric region as defined by being adjacent to the telobait sequence in telobait mode, or the 3' most complete/partial telomeric sequence in WGS mode |
| 14 | NUM_PASS | Number of full-length subreads (only for pacbio BAM data) |
| 15 | RQUAL | Read quality (only for pacbio BAM data) |

#### Telomere length calcution example:

Example read: 5'-ATAGGCATGC TTAGGGTTAGGG TTAGGG TTAGGG TTAGGG TG TTAGGG G TTAGGG TTGGG TTAGGG TTAGGGTTAG ATACAG-3'

| Term | Value | Sequence | Comment |
| :--- | :--- | :--- | :--- |
| TELOMERE_START | 11 | TTAGGGTTAGGG | Positions 11-22 |
| TELOMERE_END | 76 | TTAGGGTTAG | Positions 67-76 |
| RAW_TELOMERE_LEN | 66 | TTAGGGTTAGGG TTAGGG TTAGGG TTAGGG TG TTAGGG G TTAGGG TTGGG TTAGGG TTAGGGTTAG | TELOMERE_START(76)-TELOMERE_END(11)+1=66 |
| TELOMERE_LEN | 54 | TTAGGGTTAGGG TTAGGG TTAGGG TTAGGG TTAGGG TTAGGG TTAGGG TTAGGG | Motif count(9)*motif length(6)=54 |
| TELOMERE_END_MOTIF | GGTTAG | GGTTAG | - |

#### Subtelomeric anchor cluster table (sample.telomap.anchors.tsv)
| # | Column name | Comment |
| :--- | :--- | :--- |
| 1 | BARCODE | Detected sample barcode (will be NA for WGS mode) |
| 2 | ANCHOR_READ | Read ID of representative read for the anchor cluster |
| 3 | ANCHOR_SEQ | Consensus subtelomeric anchor sequence |
| 4 | READ_SUPPORT | Number of reads supporting subtelomeric anchor sequence |
| 5 | CHROM | Chromosomal end mapped by anchor sequence or unmappable anchor ID (i.e. begins with "U") |

### Operating system

* Linux (x86_64 architecture, tested in Ubuntu 22.04.1)

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
