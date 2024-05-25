# Built-in functions to create plots

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from natsort import natsorted, ns
from collections import Counter


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
def plot_pie(counts, output_path):
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)
    labels1 = ['Telomeric reads', 'Non-telomeric reads', 'Multiple capture oligo', 'No capture oligo']
    telomere = counts[3]
    multi = counts[2]
    only_oligo = counts[1] - counts[3]
    no_oligo = counts[0] - counts[1] - counts[2]
    data = [telomere, only_oligo, multi, no_oligo]

    def func(pct, allvals):
        absolute = int(round(pct / 100. * np.sum(allvals)))
        return "{:.1f}%\n({:d})".format(pct, absolute)

    wedges1, texts, autotexts1 = ax1.pie(data, autopct=lambda pct: func(pct, data), startangle=60,
                                         colors=['#5dbd6e', '#655dbd', '#e38421', '#bd615d'],
                                         textprops={'fontsize': 14})
    for autotext in autotexts1:
        autotext.set_color('black')

    plt.legend(wedges1, labels1,
               loc='center',
               bbox_to_anchor=(1, 1),
               prop={"size": 14})

    plt.tight_layout()
    plt.savefig(output_path + '.capture-pie.png', dpi=100, edgecolor='none')


# Plot telomere end motif
def plot_telo_end(df, output_path, barcodes, ends=None):
    if ends is None:
        ends = ['TTAGGG', 'GTTAGG', 'GGTTAG', 'GGGTTA', 'AGGGTT', 'TAGGGT']
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
    df.telo_len_wgap.hist(bins=100)
    ax.xaxis.grid(False)
    plt.ylabel('Number of reads')
    plt.xlabel('Telomere length (bp)')
    plt.title('Distribution of Overall Telomere Length')
    plt.savefig(output_path + '.telomere-len-total.png', dpi=100, edgecolor='none')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    df.telo_len.hist(bins=100)
    ax.xaxis.grid(False)
    plt.ylabel('Number of reads')
    plt.xlabel('Telomere length (bp)')
    plt.title('Distribution of Overall Telomere Motif Length')
    plt.savefig(output_path + '.telomere-mlen-total.png', dpi=100, edgecolor='none')


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


# Plot barplot for number of reads mapping to each chromosome end
def plot_chrm_bar(dfa, output_path):
    samples = dfa['barcode'].unique().tolist()
    for s in samples:
        sort_chrm = natsorted(dfa[dfa['barcode'] == s]['chrom'].tolist(), alg=ns.LOWERCASEFIRST)
        l = [x.replace('chr', '') for x in sort_chrm]
        fig = plt.figure(figsize=(16, 6))
        ax = fig.add_subplot(111)
        bar_plot = sns.barplot(x='chrom', y='anchor_hits', data=dfa[dfa['barcode'] == s], order=sort_chrm)
        bar_plot.bar_label(bar_plot.containers[0])
        _ = bar_plot.set_xticks(bar_plot.get_xticks())
        _ = bar_plot.set_xticklabels(l)
        ax.xaxis.grid(False)
        _ = plt.xticks(rotation=45, ha="right")
        plt.margins(x=0.02, tight=True)
        plt.ylabel('Number of telomeric reads')
        plt.xlabel('Chromosome end clusters')
        plt.title(s.upper())
        plt.tight_layout()
        plt.savefig(output_path + '.chrm-end-count.' + s + '.png', dpi=100, edgecolor='none')
        plt.close()


# Plot TRF1/TRF2 binding motif boxplot for each chromosome end
def plot_trf_boxplot(barcode_reads, df, output_path):
    props = {
        'medianprops': {'color': 'red'}
    }
    kb = '1kb'
    for b in barcode_reads:
        df2 = df[(df['barcode'] == b) & df['trf_count'].notnull()][['trf_count', 'chrom']].copy()
        df2['chrom'] = df2['chrom'].replace(np.nan, 'ZNA')
        fig = plt.figure(figsize=(16, 6))
        ax = fig.add_subplot(111)
        l = natsorted(df2['chrom'].unique().tolist(), alg=ns.LOWERCASEFIRST)
        ll = [x.replace('ZNA', 'NA') for x in l]
        ll = [x.replace('chr', '') for x in ll]
        box_plot = sns.boxplot(x="chrom", y="trf_count", data=df2, order=l, width=0.6, **props)
        _ = box_plot.set_xticks(box_plot.get_xticks())
        _ = box_plot.set_xticklabels(ll)
        ax.xaxis.grid(False)
        _ = plt.xticks(rotation=45, ha="right")
        plt.margins(x=0.02, tight=True)
        _ = ax.set(xlabel=None)
        _ = plt.ylabel('Frequency of TRF1/2 binding motif\nwithin first %s of canonical telomeric sequences' % kb, fontsize=12)
        _ = plt.title(b)
        plt.tight_layout()
        plt.savefig(output_path + '.trf-boxplot.' + b + '.png', dpi=100, edgecolor='none')
        # plt.show()
        plt.close()


# Plot TVS signature
def plot_tvs_sig(tvs_arr, tvs_read_counts, output_path, telo_len=1000):
    for b in tvs_arr:
        for c in tvs_arr[b]:
            fig = plt.figure(figsize=(6, 8))
            ax = fig.add_subplot(111)
            ax.set_xlabel("Telomere position starting from sub-telomeric junction (5'-3')")
            ax.set_ylabel('TVS percentage (%)')
            ax.title.set_text(b + ' | ' + c + ' | ' + str(tvs_read_counts[b][c]) + ' reads')
            # ax.spines['top'].set_color('none')
            # ax.spines['bottom'].set_color('none')
            # ax.spines['left'].set_color('none')
            # ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
            ax.plot(tvs_arr[b][c][:telo_len])
            plt.tight_layout()
            chrom = c.replace('/', '-')
            plt.savefig(output_path + '.tvs-sig.' + b + '.' + chrom + '.png', dpi=100, edgecolor='none')
            plt.close()
