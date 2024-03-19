#!/usr/bin/env python3

# Telomap: A tool for analyzing telomeres in telobait-captured or WGS long-read sequencing data

__author__ = 'CY Tham'

import os
import pickle
from datetime import datetime
from telomap import get_args
from telomap import TeloMap
#import telomap
from telomap.plots import *

def main():
    args = get_args()

    # Observe verbosity
    if args.quiet:
        sys.stdout = open(os.devnull, 'w')

    now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
    print(now + ' - Telomap started')

    # Run Telomap
    out = TeloMap(args.mode, args.reads, args.capoligo, args.barcodes, args.threads, args.name, args.motif, args.oligoscore, args.barscore, tsv_header=True)

    # Prepare directories
    os.makedirs(args.dir, exist_ok=True)
    os.makedirs(os.path.join(args.dir, 'plots'), exist_ok=True)
    plot_dir = os.path.join(args.dir, 'plots', out.input_name)
    os.makedirs(os.path.join(args.dir, 'plots_QC'), exist_ok=True)
    qc_dir = os.path.join(args.dir, 'plots_QC', out.input_name)
    # os.makedirs(os.path.join(args.dir, 'plots_trf'), exist_ok=True)
    # trf_dir = os.path.join(args.dir, 'plots_trf', out.input_name)
    # os.makedirs(os.path.join(args.dir, 'plots_tvs'), exist_ok=True)
    # tsv_dir = os.path.join(args.dir, 'plots_tvs', out.input_name)

    # Pickle data
    if args.pickle:
        with open(os.path.join(args.dir, out.input_name + ".telomap.pkl"), "wb") as f:
            pickle.dump(out, f)

    # Generate TSV files
    create_tsv(out, args.dir)

    # Generate plots
    # Plot overview figures
    # Plot number of reads for each sample
    plot_sample_read(out.df, plot_dir, out.barcodes)
    # Pie chart of telomeric read capture
    plot_pie(out.counts, plot_dir)
    # Plot telomere length histogram
    plot_telo_len(out.df, plot_dir)
    # Plot telomere length per sample (Violin)
    plot_telo_violin(out.df, plot_dir, out.barcodes)
    # Plot telomere end motif
    plot_telo_end(out.df, plot_dir, out.barcodes)
    # Plot barplot for number of reads mapping to each chromosome end
    plot_chrm_bar(out.df_anchors, plot_dir)
    # Plot gap analysis
    # plot_telo_gap(out.df, plot_dir, out.barcodes)  # Disabled due to memory consumption

    # Plot QC figures
    # Plot read length histogram
    plot_len(out.df, qc_dir)
    # Plot read passes histogram
    plot_passes(out.df, qc_dir)
    # Plot read quality histogram
    plot_qual(out.df, qc_dir)
    # Plot read quality vs passes scatter plot
    plot_pass_qual_scatter(out.df, qc_dir)

    # Plot TRF figures
    # Plot TRF1/TRF2 binding motif boxplot for each chromosome end
    # plot_trf_boxplot(out.barcode_reads, out.df, trf_dir)  # Disabled due to memory consumption

    # Plot TVS figures
    # Plot TVS signature figures
    # plot_tvs_sig(out.tvs_arr, out.tvs_read_counts, tsv_dir, telo_len=1000)  # Disabled due to memory consumption

    now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
    print(now + ' - Telomap ended')


# Create TSV file from dataframe
def create_tsv(out, wk_dir):
    df_path = os.path.join(wk_dir, out.input_name + ".telomap.tsv")
    dfa_path = os.path.join(wk_dir, out.input_name + ".telomap.anchors.tsv")
    df_tsv = open(df_path, 'w')
    dfa_tsv = open(dfa_path, 'w')
    df_tsv.write(''.join(out.header))
    dfa_tsv.write("BARCODE\tANCHOR_READ\tANCHOR_SEQ\tREAD_SUPPORT\tCHROM\n")
    df = out.df[['rname', 'oligo', 'barcode', 'strand', 'motifs', 'read_len', 'telo_len_wgap', 'telo_len', 'telo_end', 'trf_count', 'chrom', 
                 's_junct', 'e_junct', 'num_pass', 'read_qual']]
    dfa = out.df_anchors
    df = df.sort_values(by=['barcode'])
    df.to_csv(df_tsv, mode='a', header=False, index=False, sep='\t')
    dfa.to_csv(dfa_tsv, header=False, index=False, sep='\t')
    df_tsv.close()
    dfa_tsv.close()


if __name__ == '__main__':
    main()
