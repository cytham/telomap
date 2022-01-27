#!/usr/bin/env python3

# Telomap: A bioinformatics tool for analyzing telobait-captured long-read sequencing data

__version__ = '0.0.1'
__author__ = 'CY THAM'
import os
import re
import sys
from sys import argv
import datetime
from datetime import date
import pysam
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import Counter


def main():
    start_time = datetime.datetime.now()
    print('Telovar started')
    print(start_time)
    if len(argv) != 6:
        sys.exit("Input Error - usage: python telovar.py pb-reads.bam capture_oligo.fa barcodes.fa telo_motifs.fa working_directory")

    data = argv[1]
    oligo_path = argv[2]
    barcode_path = argv[3]
    motif_path = argv[4]
    wk_dir = argv[5]
    # Fixed parameters (Some may be modified according to user preference)
    data_type = 'pacbio'  # Type of long-read data (Cannot be altered)
    file_type = 'bam'  # Type of input file (Cannot be altered)
    oligo_loc = '3prime'  # Location of capture oligo sequence on telomere (3prime/5prime)
    oligo_win = 150  # Length of window (bp) to search for capture oligo sequence
    oligo_align_thres = 1  # Minimum alignment score percentage required for capture oligo sequence match [100%]
    bar_align_thres = 1  # Minimum alignment score percentage required for barcode sequence match [100%]
    motif_win = 50  # Length of window (bp) to search for telomeric motif adjacent to capture oligo
    motif_counts = 2  # Minimum number of telomeric motifs required (adjacent to capture oligo) to be considered as telomeric read
    motif_direct = 'upstream'  # Location of telomeric motif with respect to capture oligo (upstream/downstream)
    motif_end_len = 6  # Length of sequence directly adjacent to capture oligo to consider for telomere end sequence analysis
    # max_motif_gap = 200  # Not in used
    max_gap_size = 6  # Maximum gap size to be considered for gap sequence analysis

    oligos = prep_oligo(oligo_path)
    barcodes = prep_oligo(barcode_path)
    motifs = prep_motif(motif_path)

    if data_type == 'pacbio':
        if file_type == 'bam':
            input_name = os.path.basename(data).rsplit('.bam', 1)[0]
            output_path = os.path.join(wk_dir, input_name)
            save = pysam.set_verbosity(0)  # Suppress BAM index missing warning
            sam = pysam.AlignmentFile(data, 'rb', check_sq=False)
            pysam.set_verbosity(save)  # Revert verbosity level
            total = 0
            capture = 0
            telomere = 0
            alignments = {}
            motif_out = {}
            telo_end = {}
            info_dict = {'rname': [], 'read_len': [], 'num_pass': [], 'read_qual': [], 'strand': [], 'oligo': [],
                         'barcode': [], 'oscore': [], 'bscore': [], 'junct': [], 'motifs': [], 'telo_end': [], 'telo_len': [],
                         'telo_len_str': [], 'gap_size': [], 'gap_loc': [], 'gap_dist': [], 'gap_seq': [], 'telo_len_wgap': [],
                         'chrom': [], 's_junct': [], 'trf1_motif': []}
            end_motifs = ['TTAGGG', 'GTTAGG', 'GGTTAG', 'GGGTTA', 'AGGGTT', 'TAGGGT']
            chr_motif_dict = {
                'chr2S': 'GCCTCGCTTTGGGACCCCTCGGGGCCGCATCGACGGTGAATAAAATACTTCCTGTTTGCAGCACTGAATACTCAAAGTCAGACACGAGTTAGGCTTCGGC',
                'chr3S': 'CCTGCGCCGGCGCCGCGCGCCTCTCTGCGCCTGCGCCGGCGCCGCGCGCCTCTCTGCGCCTGCGCCGGCGCCGCGCGCCTCTCTGCGCCTGCAGAGAGGG',
                'chr4SE': 'GCCTCTCCACAGACGTTGGGGGCACTGCTTCATTTTGGGACAACTGGGGGCCACATCCACGGTAAAGATCCTTCCTCTTTGCAGCCGAGAATAATGAGGG',
                'chr5S': 'TGCAGCCTTCCTAATGCACACGTAACACCCAAAATATGATATCCACATTGCTCATGTAACAAGCACCTGTATGCTAATGCACTCCCTCAATACAAAATTG',
                'chr7E': 'TGCCTCAAAACAAAATATTAATATAAGATCGGCATTCCCCACACTGCCATGCAGTGCTAAGACAGCAATGAGAACAGTCAACATAATAACCCTAACAGTG',
                'chr9S': 'GAGGGTGAGGATGGGGGTGGGGGTGGGGGTTGGGGTTGGGGTTGGGGTTAGGGTTCGGGTTCGGGTTCGGGTTCGGGTTCGGGTTCGGGTTCGGGTTGGG',
                'chr11E': 'CTTGGAGGCACGGCCTTCGTTTCGGACAATTCGGGGCCGCATCGACGGTGAATAAAATCCTTCCTCTTTGCAGCCCTGAATAATCAGGGTCAGAGATCAG',
                'chr12E': 'TGCCTCAAAACAAAATATTAATATAAGATCGGCAATCCGCACACTGCCGTGCAGTGCTAAGACAGCAATGAAAATAGTCAACATAATAACCCTAATAGTG',
                'chr13E': 'ACACCCTAAATATAACAGGCATATTACTCATGGAGGGTTAGGGTTCAGGTTCGGGTTCGGGTTCGGGTTCGGGTTCGGGGTTCGGGGTTTGGGGTTCGGG',
                'chr15E': 'GGGCTGCATTAAAGGGTCCAGTTGCAGCATCGGAACGCAAATGCAGCATTCCTAATGCACACATGATACCCAAAATATAACACCCACATTGCTCATGTGG',
                'chr18E': 'GCGCGGGCGCGCCGCGCCTCTCTGCGCCTGCGCGGGCGCGCCGCGCCTCTCTGCGCCTGCGCGGGCGCGCCGCCTTTGCGAGGATGGAGTTGCGTTCTCC',
                'chr19E': 'GCGAAGGCGGAGCACAGTTCTTCTCAGGTCAGACCCGGGTGGGCGGGCTGAGGGCACTGCGAGGGCGGAGCTGCGTTCTGTTCAGCACAGACCTGGGGGG'
            }
            # end_add = {'TTAGGG': 0, 'GTTAGG': 5, 'GGTTAG': 4, 'GGGTTA': 3, 'AGGGTT': 2, 'TAGGGT': 1}
            for seg in sam:
                total += 1
                if total % 500000 == 0:
                    print('Analyzed %i sequences' % total)
                    now_time = datetime.datetime.now()
                    print('Time elasped:', now_time - start_time)
                # if total == 100000:
                #     break
                qname = seg.query_name
                fasta = seg.query_sequence
                info_dict['rname'].append(qname)
                info_dict['read_len'].append(len(fasta))
                info_dict['num_pass'].append(seg.get_tag('np'))
                info_dict['read_qual'].append(round(seg.get_tag('rq'), 3))
                alignments[qname] = oligo_align(oligos, barcodes, fasta, oligo_window=oligo_win, direction=oligo_loc,
                                                othreshold=oligo_align_thres, bthreshold=bar_align_thres)
                if alignments[qname]:
                    if len(alignments[qname]) > 1:  # If read contains oligo
                        capture += 1
                        junct = int(alignments[qname][5])
                        # Record data
                        info_dict['oligo'].append(alignments[qname][0][1:])  # remove '>'
                        info_dict['barcode'].append(alignments[qname][1][1:])  # remove '>'
                        info_dict['oscore'].append(str(int(alignments[qname][3])))
                        info_dict['bscore'].append(str(int(alignments[qname][4])))
                        info_dict['junct'].append(str(junct))
                        if alignments[qname][2] == 'forward':
                            info_dict['strand'].append('+')
                        else:  # Reverse
                            info_dict['strand'].append('-')
                            fa = Seq(fasta)
                            fasta = str(fa.reverse_complement())
                        motif_out[qname] = motif_finder(motifs, fasta,
                                                        junction=junct, motif_window=motif_win,
                                                        min_counts=motif_counts, direction=motif_direct)
                        if motif_out[qname]:  # If read contains telomeric motifs
                            # Retrieve telomere end sequence
                            m = motif_out[qname][0]
                            motif_indexes = [(i.start(0), i.end(0)) for i in re.finditer(m, fasta)]

                            blocks, tandem = linker_tandem(motif_indexes)
                            if tandem:  # if two motif in tandem exist
                                telomere += 1
                                info_dict['motifs'].append(m)
                                telo_end[qname] = adjacent_motif(fasta, junct, direction=motif_direct, motif_len=motif_end_len)
                                trf1_motif_indexes = [(i.start(0), i.end(0)) for i in re.finditer('TTAGGGTTA', fasta)]
                                info_dict['telo_end'].append(telo_end[qname])
                                block_limit = blocks[-1][0]
                                count = 0
                                telo_motif_indexes = []
                                for k in motif_indexes:
                                    if k[0] >= block_limit:
                                        count += 1
                                        telo_motif_indexes.append(k)
                                chroms = sub_telo_label(fasta, block_limit, chr_motif_dict)
                                info_dict['telo_len_wgap'].append(int(blocks[-1][1] - blocks[-1][0]))  # + add_len)
                                info_dict['telo_len'].append(count * 6)  # + add_len)
                                info_dict['telo_len_str'].append(str(count * 6))  # + add_len))
                                gap_sizes, gap_locs, gap_dists = gaps(telo_motif_indexes, min_gap=1)
                                info_dict['gap_size'].append(gap_sizes)
                                info_dict['gap_loc'].append(gap_locs)
                                info_dict['gap_dist'].append(gap_dists)
                                temp_seq = []
                                for i in range(len(gap_sizes)):
                                    if gap_sizes[i] <= max_gap_size:
                                        temp_seq.append(fasta[gap_locs[i][0]: gap_locs[i][1]])
                                    else:
                                        temp_seq.append(None)
                                info_dict['gap_seq'].append(temp_seq)
                                info_dict['chrom'].append(chroms)
                                info_dict['s_junct'].append(str(int(block_limit)))
                                info_dict['trf1_motif'].append(trf1_motif_indexes)
                            else:
                                info_dict['motifs'].append(None)
                                info_dict['telo_end'].append(None)
                                info_dict['telo_len'].append(None)
                                info_dict['telo_len_str'].append(None)
                                info_dict['telo_len_wgap'].append(None)
                                info_dict['gap_size'].append(None)
                                info_dict['gap_loc'].append(None)
                                info_dict['gap_dist'].append(None)
                                info_dict['gap_seq'].append(None)
                                info_dict['chrom'].append(None)
                                info_dict['s_junct'].append(None)
                                info_dict['trf1_motif'].append(None)
                        else:
                            info_dict['motifs'].append(None)
                            info_dict['telo_end'].append(None)
                            info_dict['telo_len'].append(None)
                            info_dict['telo_len_str'].append(None)
                            info_dict['telo_len_wgap'].append(None)
                            info_dict['gap_size'].append(None)
                            info_dict['gap_loc'].append(None)
                            info_dict['gap_dist'].append(None)
                            info_dict['gap_seq'].append(None)
                            info_dict['chrom'].append(None)
                            info_dict['s_junct'].append(None)
                            info_dict['trf1_motif'].append(None)
                    elif alignments[qname][0] == 'multi':
                        # Record data
                        info_dict['strand'].append(None)
                        info_dict['oligo'].append('Multi')
                        info_dict['barcode'].append(None)
                        info_dict['oscore'].append(None)
                        info_dict['bscore'].append(None)
                        info_dict['junct'].append(None)
                        info_dict['motifs'].append(None)
                        info_dict['telo_end'].append(None)
                        info_dict['telo_len'].append(None)
                        info_dict['telo_len_str'].append(None)
                        info_dict['telo_len_wgap'].append(None)
                        info_dict['gap_size'].append(None)
                        info_dict['gap_loc'].append(None)
                        info_dict['gap_dist'].append(None)
                        info_dict['gap_seq'].append(None)
                        info_dict['chrom'].append(None)
                        info_dict['s_junct'].append(None)
                        info_dict['trf1_motif'].append(None)
                        motif_out[qname] = []
                    else:
                        raise Exception('Error: Unrecognised alignment object %s' % str(alignments))
                else:
                    # Record data
                    info_dict['strand'].append(None)
                    info_dict['oligo'].append(None)
                    info_dict['barcode'].append(None)
                    info_dict['oscore'].append(None)
                    info_dict['bscore'].append(None)
                    info_dict['junct'].append(None)
                    info_dict['motifs'].append(None)
                    info_dict['telo_end'].append(None)
                    info_dict['telo_len'].append(None)
                    info_dict['telo_len_str'].append(None)
                    info_dict['telo_len_wgap'].append(None)
                    info_dict['gap_size'].append(None)
                    info_dict['gap_loc'].append(None)
                    info_dict['gap_dist'].append(None)
                    info_dict['gap_seq'].append(None)
                    info_dict['chrom'].append(None)
                    info_dict['s_junct'].append(None)
                    info_dict['trf1_motif'].append(None)
                    motif_out[qname] = []
            df = pd.DataFrame.from_dict(info_dict)
            with open(output_path + ".pkl", "wb") as f:
                pickle.dump(df, f)
            plot_len(df, output_path)
            plot_pass_qual_scatter(df, output_path)
            plot_passes(df, output_path)
            plot_qual(df, output_path)
            plot_pie(total, capture, telomere, output_path)
            plot_sample_read(df, output_path, barcodes)
            plot_telo_end(df, output_path, barcodes, end_motifs)
            plot_telo_len(df, output_path)
            plot_telo_violin(df, output_path, barcodes)
            plot_telo_gap(df, output_path, barcodes)
            out_path = os.path.join(wk_dir, '%s.telovar.tsv' % input_name)
            tsv = open(out_path, 'w')
            write_metadata(tsv, data, file_type, data_type, oligo_path, barcode_path, motif_path, wk_dir, oligo_loc, oligo_win,
                           oligo_align_thres, bar_align_thres, motifs, motif_win, motif_counts, motif_direct, motif_end_len,
                           max_motif_gap, total, capture, telomere)
            tsv.close()
            df = df[['rname', 'read_len', 'num_pass', 'read_qual', 'strand', 'oligo', 'oscore', 'barcode', 'bscore', 'junct',
                     'motifs', 'telo_end', 'telo_len_str', 'chrom', 's_junct']]
            df = df.sort_values(by=['barcode'])
            df.to_csv(out_path, mode='a', header=False, index=False, sep='\t')
            end_time = datetime.datetime.now()
            print('Telovar ended')
            print(end_time)
            print('Time elasped:', end_time - start_time)
        else:
            raise Exception('Error: Input data type needs to be bam, please check the script.')
    else:
        raise Exception('Error: Input data type needs to be pacbio, please check the script.')


# Write metadata into output TSV file
def write_metadata(tsv, data_path, file_type, data_type, oligo_path, barcode_path, motif_path, wk_dir, oloc, owin, othres,
                   bthres, motifs, mwin, mcount, mdir, melen, mgap, treads, oreads, mreads):
    today = date.today()
    today_date = today.strftime("%d-%m-%Y")
    tsv.write('##date=%s\n' % today_date)
    tsv.write('##tool=Telovar-%s\n' % __version__)
    tsv.write('##source_reads=%s\n' % data_path)
    tsv.write('##file_type=%s\n' % file_type)
    tsv.write('##data_type=%s\n' % data_type)
    tsv.write('##oligo_sequence=%s\n' % oligo_path)
    tsv.write('##barcode_sequence=%s\n' % barcode_path)
    tsv.write('##motif_sequence=%s\n' % motif_path)
    tsv.write('##output_directory=%s\n' % wk_dir)
    tsv.write('##oligo_location=%s\n' % oloc)
    tsv.write('##oligo_window=%s\n' % owin)
    tsv.write('##oligo_alignment_threshold=%s\n' % othres)
    tsv.write('##barcode_alignment_threshold=%s\n' % bthres)
    tsv.write('##motifs=%s\n' % motifs)
    tsv.write('##motif_window=%s\n' % mwin)
    tsv.write('##minimum_motif_count=%s\n' % mcount)
    tsv.write('##motif_direction=%s\n' % mdir)
    tsv.write('##motif_end_length=%s\n' % melen)
    tsv.write('##motif_max_gap=%s\n' % mgap)
    tsv.write('##total_reads=%s\n' % treads)
    tsv.write('##captured_reads=%s\n' % oreads)
    tsv.write('##captured_reads_with_motif=%s\n' % mreads)
    tsv.write('#RNAME\tRLEN\tNUM_PASS\tRQUAL\tSTRAND\tOLIGO\tOLIGO_SCORE\tBARCODE\tBARCODE_SCORE\tBARCODE_JUNC\tMOTIF'
              '\tTELOMERE_END\tTELOMERE_LEN\tChrom\tS_Junct\n')


# Plot read length histogram
def plot_len(df, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.read_len.hist(bins=np.logspace(start=np.log10(10), stop=np.log10(1000000), num=200))
    ax.xaxis.grid(False)
    ax.set_xscale('log')
    plt.ylabel('Number of reads')
    plt.xlabel('Read length (bp)')
    plt.title('Distribution of Read Length')
    plt.savefig(output_path + '.read-len.png', dpi=100, edgecolor='none')


# Plot read passes histogram
def plot_passes(df, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.num_pass.hist(bins=200, range=(0, 250))
    ax.xaxis.grid(False)
    # ax.set_xscale('log')
    plt.ylabel('Number of reads')
    plt.xlabel('Number of passes')
    plt.title('Distribution of Passes')
    plt.savefig(output_path + '.num-pass.png', dpi=100, edgecolor='none')


# Plot read quality histogram
def plot_qual(df, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.read_qual.hist(bins=50)
    ax.xaxis.grid(False)
    plt.ylabel('Number of reads')
    plt.xlabel('Quality score')
    plt.title('Distribution of Read Quality')
    plt.savefig(output_path + '.read-qual.png', dpi=100, edgecolor='none')


# Plot read quality vs passes scatter plot
def plot_pass_qual_scatter(df, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.plot.scatter("num_pass", "read_qual", alpha=0.6, figsize=(8, 6))
    ax.xaxis.grid(False)
    plt.ylabel('Quality score')
    plt.xlabel('Number of passes')
    plt.title('Read Quality score vs Number of Passes')
    plt.savefig(output_path + '.pass-qual-scatter.png', dpi=100, edgecolor='none')


# Plot number of reads for each sample
def plot_sample_read(df, output_path, barcodes):
    samples = []
    sample_counts = []
    sample_telo = []
    for b in barcodes:
        samples.append(b[1:])
        sample_counts.append(df[df['barcode'] == b[1:]].shape[0])
        sample_telo.append(df[(df['barcode'] == b[1:]) & (df['motifs'].notnull())].shape[0])
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.bar(samples, sample_counts)
    ax.xaxis.grid(False)
    plt.ylabel('Number of reads with barcode')
    plt.xlabel('Samples')
    plt.title('Number of Barcoded Reads Per Sample')
    plt.savefig(output_path + '.sample-reads.png', dpi=100, edgecolor='none')
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.bar(samples, sample_telo)
    ax.xaxis.grid(False)
    plt.ylabel('Number of barcoded reads with motif')
    plt.xlabel('Samples')
    plt.title('Number of Barcoded Reads With Motif Per Sample')
    plt.savefig(output_path + '.sample-reads-motif.png', dpi=100, edgecolor='none')


# Pie chart of telomeric read capture
def plot_pie(total, capture, telomere, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    labels1 = ['With capture oligo', 'Without capture oligo']
    data1 = [capture, total - capture]
    labels2 = ['Telomeric reads', 'Non-telomeric reads']
    data2 = [telomere, capture - telomere]

    def func(pct, allvals):
        absolute = int(round(pct / 100. * np.sum(allvals)))
        return "{:.1f}%\n({:d})".format(pct, absolute)

    wedges1, texts, autotexts1 = ax1.pie(data1, autopct=lambda pct: func(pct, data1), startangle=60,
                                         colors=['#5dbd6e', '#bd615d'])
    for autotext in autotexts1:
        autotext.set_color('black')

    wedges2, texts, autotexts2 = ax2.pie(data2, autopct=lambda pct: func(pct, data2), startangle=60,
                                         colors=['#655dbd', '#e59400'])

    for autotext in autotexts2:
        autotext.set_color('black')

    plt.legend(wedges1 + wedges2, labels1 + labels2,
               loc='center',
               bbox_to_anchor=(-0.2, 0))

    plt.tight_layout()
    plt.suptitle('Read Capture Proportions')
    plt.savefig(output_path + '.capture-pie.png', dpi=100, edgecolor='none')


# Plot telomere end motif
def plot_telo_end(df, output_path, barcodes, ends):
    ends_dict = {x: [] for x in ends}
    ends_dict['Others'] = []
    for b in barcodes:
        total = df[(df['barcode'] == b[1:]) & (df['telo_end'].notnull())].shape[0]
        for e in ends:
            count = df[(df['barcode'] == b[1:]) & (df['telo_end'] == e)].shape[0]
            if count > 0:
                ends_dict[e].append(count / total * 100)
            else:
                ends_dict[e].append(0)
        count = df[(df['barcode'] == b[1:]) & (~df['telo_end'].isin(ends)) & (df['telo_end'].notnull())].shape[0]
        if count > 0:
            ends_dict['Others'].append(count / total * 100)
        else:
            ends_dict['Others'].append(0)
    df2 = pd.DataFrame(ends_dict, index=[x[1:] for x in barcodes])
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    df2.plot.bar(rot=0, width=0.9, figsize=(10, 6))
    ax.xaxis.grid(False)
    plt.ylim(0, 100)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('Percentage of telomeric reads')
    plt.xlabel('Samples')
    plt.title('Analysis of Telomere Ends Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-end-samples.png', dpi=100, edgecolor='none')


# Plot read length histogram
def plot_telo_len(df, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.telo_len.hist(bins=100)
    ax.xaxis.grid(False)
    plt.ylabel('Number of reads')
    plt.xlabel('Telomere length (bp)')
    plt.title('Distribution of Overall Telomere Length')
    plt.savefig(output_path + '.telomere-len-total.png', dpi=100, edgecolor='none')


# Plot telomere length violin plot
def plot_telo_violin(df, output_path, barcodes):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    _ = sns.violinplot(x="barcode", y="telo_len", data=df, order=[x[1:] for x in barcodes], cut=0)
    ax.xaxis.grid(False)
    plt.ylabel('Telomere length (bp)')
    plt.xlabel('Samples')
    plt.title('Distribution of Telomere Length Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-violin.png', dpi=100, edgecolor='none')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    props = {
        'boxprops': {'facecolor': 'none', 'edgecolor': 'black'},
        'medianprops': {'color': 'red'},
        'whiskerprops': {'color': 'black'},
        'capprops': {'color': 'black'}
    }
    _ = sns.boxplot(x="barcode", y="telo_len", data=df, order=[x[1:] for x in barcodes], color='w', **props)
    ax.xaxis.grid(False)
    plt.ylabel('Telomere length (bp)')
    plt.xlabel('Samples')
    plt.title('Distribution of Telomere Length Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-boxplot.png', dpi=100, edgecolor='none')
    ###
    info_tab = [[], [], [], [], []]
    for b in barcodes:
        df2 = df[(df['barcode'] == b[1:]) & (df['telo_len'].notnull())]['telo_len']
        if df2.empty:
            info_tab[0].append('-')
            info_tab[1].append('-')
            info_tab[2].append('-')
            info_tab[3].append('-')
            info_tab[4].append('-')
        else:
            info_tab[0].append(int(df2.mean()))
            info_tab[1].append(int(df2.median()))
            info_tab[2].append(str(int(df2.quantile(0.25))) + '-' + str(int(df2.quantile(0.75))))
            info_tab[3].append(int(df2.min()))
            info_tab[4].append(int(df2.max()))
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(111)
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=info_tab, rowLabels=['Mean (bp)', 'Median (bp)', 'IQR (bp)', 'Min (bp)', 'Max (bp)'],
             colLabels=[x[1:] for x in barcodes], cellLoc='center', rowLoc='center', colLoc='center', loc="center")
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-table.png', dpi=300, edgecolor='none')
    ###
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    _ = sns.violinplot(x="barcode", y="telo_len_wgap", data=df, order=[x[1:] for x in barcodes], cut=0)
    ax.xaxis.grid(False)
    plt.ylabel('Telomere length (bp)')
    plt.xlabel('Samples')
    plt.title('Distribution of Telomere Length Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-violin-wgap.png', dpi=100, edgecolor='none')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    props = {
        'boxprops': {'facecolor': 'none', 'edgecolor': 'black'},
        'medianprops': {'color': 'red'},
        'whiskerprops': {'color': 'black'},
        'capprops': {'color': 'black'}
    }
    _ = sns.boxplot(x="barcode", y="telo_len_wgap", data=df, order=[x[1:] for x in barcodes], color='w', **props)
    ax.xaxis.grid(False)
    plt.ylabel('Telomere length (bp)')
    plt.xlabel('Samples')
    plt.title('Distribution of Telomere Length Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-boxplot-wgap.png', dpi=100, edgecolor='none')
    ###
    info_tab = [[], [], [], [], []]
    for b in barcodes:
        df2 = df[(df['barcode'] == b[1:]) & (df['telo_len_wgap'].notnull())]['telo_len_wgap']
        if df2.empty:
            info_tab[0].append('-')
            info_tab[1].append('-')
            info_tab[2].append('-')
            info_tab[3].append('-')
            info_tab[4].append('-')
        else:
            info_tab[0].append(int(df2.mean()))
            info_tab[1].append(int(df2.median()))
            info_tab[2].append(str(int(df2.quantile(0.25))) + '-' + str(int(df2.quantile(0.75))))
            info_tab[3].append(int(df2.min()))
            info_tab[4].append(int(df2.max()))
    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(111)
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=info_tab, rowLabels=['Mean (bp)', 'Median (bp)', 'IQR (bp)', 'Min (bp)', 'Max (bp)'],
             colLabels=[x[1:] for x in barcodes], cellLoc='center', rowLoc='center', colLoc='center', loc="center")
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-len-table-wgap.png', dpi=300, edgecolor='none')


# Plot gap analysis
def plot_telo_gap(df, output_path, barcodes):
    # Gap size histogram
    total_dict = {'barcode': [], 'gap_size': [], 'gap_dist': [], 'gap_seq': []}
    for b in barcodes:
        listoflist = df[(df['barcode'] == b[1:]) & (df['telo_len'].notnull())]['gap_size'].tolist()
        for i in listoflist:
            for j in i:
                total_dict['barcode'].append(b[1:])
                total_dict['gap_size'].append(j)
        listoflist = df[(df['barcode'] == b[1:]) & (df['telo_len'].notnull())]['gap_dist'].tolist()
        for i in listoflist:
            for j in i:
                total_dict['gap_dist'].append(j)
        listoflist = df[(df['barcode'] == b[1:]) & (df['telo_len'].notnull())]['gap_seq'].tolist()
        for i in listoflist:
            for j in i:
                total_dict['gap_seq'].append(j)
    df2 = pd.DataFrame.from_dict(total_dict)
    all_gaps = df2['gap_size'].tolist()
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.hist(all_gaps, bins=5000, weights=np.ones(len(all_gaps)) / len(all_gaps))
    plt.xlim(0, 300)
    plt.ylim(0, 1)
    ax.xaxis.grid(False)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1))
    plt.ylabel('Percentage of gaps')
    plt.xlabel('Gap size (bp)')
    plt.title('Analysis of Telomeric Gap Sizes')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap.png', dpi=100, edgecolor='none')
    # Gap ranges
    ranges = []
    for i in all_gaps:
        if i <= 10:
            ranges.append('1-10')
        elif i <= 100:
            ranges.append('11-100')
        elif i <= 200:
            ranges.append('101-200')
        elif i <= 300:
            ranges.append('201-300')
        else:
            ranges.append('>300')
    df2['ranges'] = ranges
    x = Counter(ranges)
    temp_num = []
    temp_perc = []
    for i in ['1-10', '11-100', '101-200', '201-300', '>300']:
        temp_num.append(x[i])
        temp_perc.append(str(round(x[i] / len(all_gaps) * 100, 1)) + '%')
    tab_list = [temp_num, temp_perc]
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=tab_list, rowLabels=['Number', 'Percentage'],
             colLabels=['1-10 bp', '11-100 bp', '101-200 bp', '201-300 bp', '>300 bp'], cellLoc='center', rowLoc='center',
             colLoc='center', loc="center")
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap-table.png', dpi=300, edgecolor='none')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.bar(['1-10 bp', '11-100 bp', '101-200 bp', '201-300 bp', '>300 bp'], [float(x.strip('%')) for x in temp_perc])
    plt.ylim(0, 100)
    ax.xaxis.grid(False)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())
    plt.ylabel('Percentage of gaps')
    plt.title('Analysis of Telomeric Gap Sizes')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap-range.png', dpi=100, edgecolor='none')
    # Gap size range across samples
    range_dict = {x: [] for x in ['1-10 bp', '11-100 bp', '101-200 bp', '201-300 bp', '>300 bp']}
    for b in barcodes:
        total = df2[df2['barcode'] == b[1:]].shape[0]
        for r in ['1-10 bp', '11-100 bp', '101-200 bp', '201-300 bp', '>300 bp']:
            count = df2[(df2['barcode'] == b[1:]) & (df2['ranges'] == r.strip(' bp'))].shape[0]
            if count > 0:
                range_dict[r].append(count / total * 100)
            else:
                range_dict[r].append(0)
    df3 = pd.DataFrame(range_dict, index=[x[1:] for x in barcodes])
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    df3.plot.bar(rot=0, width=0.9, figsize=(10, 6))
    ax.xaxis.grid(False)
    plt.ylim(0, 100)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('Percentage of gaps')
    plt.xlabel('Samples')
    plt.title('Analysis of Telomeric Gap Sizes Across Samples')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap-samples.png', dpi=100, edgecolor='none')

    # Gap dist

    def combine(m, v):
        return str(m) + '-' + str(v)

    size_list = df2['gap_size'].tolist()
    dist_list = df2['gap_dist'].tolist()
    k = list(map(combine, size_list, dist_list))
    k = Counter(k)
    x = []
    y = []
    z = []
    for i in k:
        x.append(int(i.split('-')[1]))
        y.append(int(i.split('-')[0]))
        z.append(k[i])
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.scatter(x, y, s=z, alpha=0.5)
    plt.ylim(0, 600)
    plt.ylabel('Gap size (bp)')
    plt.xlabel('Gap distance from capture oligo (bp)')
    plt.title('Analysis of Telomeric Gap Size and distance')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap-dist.png', dpi=100, edgecolor='none')
    # Plot gap dist histogram
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.hist(x, bins=1000)
    plt.xlim(0, 20000)
    plt.ylabel('Number of gaps')
    plt.xlabel('Gap distance from capture oligo (bp)')
    plt.title('Analysis of Telomeric Gap distance')
    plt.tight_layout()
    plt.savefig(output_path + '.telomere-gap-dist2.png', dpi=100, edgecolor='none')
    # Gap seq
    total_seq = df2[df2['gap_seq'].notnull()]['gap_seq'].tolist()
    seq_dict = {'one': [], 'two': [], 'three': [], 'four': [], 'five': [], 'six': []}
    for i in total_seq:
        if len(i) == 1:
            seq_dict['one'].append(i)
        elif len(i) == 2:
            seq_dict['two'].append(i)
        elif len(i) == 3:
            seq_dict['three'].append(i)
        elif len(i) == 4:
            seq_dict['four'].append(i)
        elif len(i) == 5:
            seq_dict['five'].append(i)
        elif len(i) == 6:
            seq_dict['six'].append(i)
        else:
            raise Exception('Error: Number of gap bases exceeds six.')
    n = 1
    for i in seq_dict:
        labels = []
        ratios = []
        total_counts = 0
        x = Counter(seq_dict[i])
        for k in x:
            total_counts += x[k]
        for k in x:
            ratio = x[k] / total_counts
            if ratio >= 0.05:  # Threshold of 5% to list gap sequence in figure
                labels.append(k)
                ratios.append(x[k] / total_counts)
        if n == 3:
            fig = plt.figure(figsize=(8, 6))  # fig = plt.figure(figsize=(18, 6))
        else:
            fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        plt.bar(labels, ratios)
        plt.ylabel('Proportion')
        plt.xlabel('Gap sequence')
        plt.title('Analysis of Telomeric Gap sequence')
        plt.tight_layout()
        plt.savefig(output_path + '.telomere-gap-seq' + str(n) + '.png', dpi=100, edgecolor='none')
        n += 1
    # Gap ratio across sample
    g = 6  # Plot in group of 6
    barcode_list = [b[1:] for b in barcodes]
    groups = [barcode_list[i * g:(i + 1) * g] for i in range((len(barcode_list) + g - 1) // g)]
    n = 1
    size_limit = 1000
    for j in groups:
        fig, axs = plt.subplots(len(j), figsize=(10, 8))
        if len(j) > 1:
            for w, b in enumerate(j):
                data_list = df2[(df2['barcode'] == b) & (df2['gap_size'] <= size_limit)]['gap_size'].tolist()
                axs[w].hist(data_list, bins=400, density=True)
                axs[w].set_title(b)
                axs[w].set_xlim(0, size_limit)
                axs[w].set_yscale('log')
        else:
            for w, b in enumerate(j):
                data_list = df2[(df2['barcode'] == b) & (df2['gap_size'] <= size_limit)]['gap_size'].tolist()
                axs.hist(data_list, bins=400, density=True)
                axs.set_title(b)
                axs.set_xlim(0, size_limit)
                axs.set_yscale('log')
        # add big axis for labels
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Gap size (bp)')
        plt.ylabel('Gap size density', labelpad=20)
        plt.tight_layout()
        plt.savefig(output_path + '.telomere-gap-size-ratio' + str(n) + '.png', dpi=100, edgecolor='none')
        n += 1


# Prepare capture oligos for alignment
def prep_oligo(oligo_fasta):
    oligos = {}
    # words = {}
    with open(oligo_fasta) as f:
        for i in f:
            oligo_name = i.strip()
            oligo = next(f).strip()
            fa = Seq(oligo)
            rc_oligo = str(fa.reverse_complement())
            oligos[oligo_name] = [oligo, rc_oligo]
    return oligos


# Prepare motifs for identification
def prep_motif(motif_fasta):
    motifs = []
    with open(motif_fasta) as f:
        for _ in f:
            motifs.append(next(f).strip())
    return motifs


# Align oligos and capture reads
def oligo_align(oligos, barcodes, fasta, oligo_window=100, direction='5prime', othreshold=0.80, bthreshold=0.90,
                othres_multi=0.80):
    if direction == '3prime':
        forward_region = fasta[-1 * oligo_window:]
        reverse_region = fasta[:oligo_window]
    else:
        forward_region = fasta[:oligo_window]
        reverse_region = fasta[-1 * oligo_window:]
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    for oligo in oligos:
        aligner.open_gap_score = -100
        aligner.extend_gap_score = -100
        max_score = len(oligos[oligo][0])
        alignments = aligner.align(forward_region, oligos[oligo][0])
        score_o = alignments[0].score
        start = alignments[0].aligned[0][0][0]
        if score_o / max_score >= othreshold:
            if len(alignments) > 1:  # Omit read if more than one capture oligo is found in oligo window
                return ['multi']
            # Check additional oligo presence in read
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -1
            if len(fasta) >= oligo_window:
                forward_rest = fasta[:len(fasta) - oligo_window + start]
            else:
                forward_rest = fasta[:start]
            for _oligo in oligos:
                if len(forward_rest) > 0:
                    _alignments = aligner.align(forward_rest, oligos[_oligo][0])
                    if len(_alignments) > 0:
                        if _alignments[0].score / max_score >= othres_multi:
                            return ['multi']
                _alignments = aligner.align(fasta, oligos[_oligo][1])  # Try reverse complement oligo
                if _alignments[0].score / max_score >= othres_multi:
                    return ['multi']
            aligner.open_gap_score = -100
            aligner.extend_gap_score = -100
            # Screen through all barcodes and break if two hits  ## Report only first barcode that passed threshold
            hit = 0
            hit_barcode = ''
            hit_score_b = 0
            junction = 0
            for barcode in barcodes:
                max_score = len(barcodes[barcode][0])
                alignments = aligner.align(forward_region, barcodes[barcode][0])
                score_b = alignments[0].score
                if score_b / max_score >= bthreshold:
                    hit += 1
                    if hit > 1:
                        return ['multi']
                    else:
                        start = alignments[0].aligned[0][0][0]
                        if len(fasta) >= oligo_window:
                            junction = len(fasta) - oligo_window + start
                        else:
                            junction = start
                        hit_barcode = barcode
                        hit_score_b = score_b
                        if fasta[junction:junction + 13] != barcodes[barcode][0]:
                            print('forward')
                            print(fasta)
                            print(barcodes[barcode][0])
                            print(barcode)
                            raise Exception('Barcode not the same')
            if hit == 1:
                return [oligo, hit_barcode, 'forward', score_o, hit_score_b, junction]
        else:
            # Try reverse complement oligo
            alignments = aligner.align(reverse_region, oligos[oligo][1])
            score_o = alignments[0].score
            end = alignments[0].aligned[0][-1][1]
            if score_o / max_score >= othreshold:
                if len(alignments) > 1:  # Omit read if more than one capture oligo is found in oligo window
                    return ['multi']
                # Check additional oligo presence in read
                aligner.open_gap_score = -1
                aligner.extend_gap_score = -1
                reverse_rest = fasta[end:]
                for _oligo in oligos:
                    if len(reverse_rest) > 0:
                        _alignments = aligner.align(reverse_rest, oligos[_oligo][1])
                        if len(_alignments) > 0:
                            if _alignments[0].score / max_score >= othres_multi:
                                return ['multi']
                    _alignments = aligner.align(fasta, oligos[_oligo][0])  # Try un-reverse complement oligo
                    if _alignments[0].score / max_score >= othres_multi:
                        return ['multi']
                aligner.open_gap_score = -100
                aligner.extend_gap_score = -100
                # Screen through all barcodes and break if two hits  ## Report only first barcode that passed threshold
                hit = 0
                hit_barcode = ''
                hit_score_b = 0
                junction = 0
                for barcode in barcodes:
                    max_score = len(barcodes[barcode][0])
                    alignments = aligner.align(reverse_region, barcodes[barcode][1])
                    score_b = alignments[0].score
                    if score_b / max_score >= bthreshold:
                        hit += 1
                        if hit > 1:
                            return ['multi']
                        else:
                            end = alignments[0].aligned[0][-1][1]
                            junction = len(fasta) - end
                            hit_barcode = barcode
                            hit_score_b = score_b
                            if fasta[end - 13:end] != barcodes[barcode][1]:
                                print('reverse')
                                print(fasta)
                                print(barcodes[barcode][1])
                                print(barcode)
                                raise Exception('Barcode not the same')
                if hit == 1:
                    return [oligo, hit_barcode, 'reverse', score_o, hit_score_b, junction]


# Find motifs in sequence
def motif_finder(motifs, fasta, junction=0, motif_window=-1, min_counts=1, direction='downstream'):
    # Only consider motif in forward strand
    if motif_window == -1:  # If search_region is not use, search entire FASTA
        motif_window = len(fasta)
    if direction == 'upstream':
        search_region = fasta[junction - motif_window:junction]
    else:
        search_region = fasta[junction:junction + motif_window]
    output = []
    for motif in motifs:
        matches = [(i.start(0), i.end(0)) for i in re.finditer(motif, search_region)]
        if len(matches) >= min_counts:
            output = [motif, len(matches)]
            break  # Only choose 1st motif
    return output


# Find motif adjacent of reference point
def adjacent_motif(fasta, junction, direction='downstream', motif_len=6):
    if direction == 'upstream':
        return fasta[junction - motif_len:junction]
    else:
        return fasta[junction:junction + motif_len]


# Find the block of motif sequences by maximum gap size
def linker(indexes, max_gap):
    blocks = []
    start = indexes[0][0]
    for i in range(len(indexes) - 1):
        if indexes[i + 1][0] - indexes[i][1] <= max_gap:
            pass
        else:
            break_ = indexes[i][1]
            blocks.append((start, break_))
            start = indexes[i + 1][0]
    blocks.append((start, indexes[-1][1]))
    return blocks


# Find the block of motif sequences by recognizing the first two tandem telomere motif
def linker_tandem(indexes):
    tandem = False
    blocks = []
    start = indexes[0][0]
    for i in range(len(indexes) - 1):
        if indexes[i + 1][0] - indexes[i][1] == 0:  # if motifs in tandem
            blocks.append((start, indexes[-1][1]))
            tandem = True
            break
        else:  # if motifs are not in tandem, containing gaps between them
            blocks.append((start, indexes[i][1]))
            start = indexes[i + 1][0]
    return blocks, tandem


# Analyze gaps in telomeric read
def gaps(indexes, min_gap=1):
    gap_sizes = []
    gap_loc = []
    for i in range(len(indexes) - 1):
        gap_size = indexes[i + 1][0] - indexes[i][1]
        if gap_size >= min_gap:
            gap_sizes.append(gap_size)
            gap_loc.append([indexes[i][1], indexes[i + 1][0]])
    # last_telo = indexes[-1][1]
    start_telo = indexes[0][0]  # Start of tandem repeat
    # inv_gap_dist = [last_telo-x[1]-2 for x in gap_loc]
    gap_dist = [x[1] - start_telo for x in gap_loc]
    return gap_sizes, gap_loc, gap_dist


# Label chromosomes by subtelomeric region
def sub_telo_label(fasta, block_limit, chr_motif_dict, motif_len=100, min_score=97):
    sub_telo_region = fasta[:block_limit]
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    chroms = []
    if len(sub_telo_region) >= motif_len:
        for chromo in chr_motif_dict:
            alignments = aligner.align(sub_telo_region, chr_motif_dict[chromo])
            if alignments[0].score > min_score:
                chroms.append(chromo)
    return chroms


if __name__ == '__main__':
    main()
