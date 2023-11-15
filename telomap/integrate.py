# Main pipeline of Telomap

import os
from Bio.Seq import Seq
from datetime import date
from .capture import TeloCapture
from .cluster import SubTeloClust
from .tvs import tvs_analyzer
from .version import __version__


class TeloMap:

    def __init__(self, read_path: str, oligo_path: str, barcode_path: str, data_type: str, cores: int, tsv_header=False):
        self.chm13_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ref/t2t-chm13-subtelo-1000-60.fa')
        self.data_type = data_type  # fasta, fastq, bam, pacbio-bam
        self.cores = cores
        self.tsv_header = tsv_header
        # Parse capture oligo sequences
        self.oligos = self.parse_fasta(oligo_path)
        # Parse barcode sequences
        self.barcodes = self.parse_fasta(barcode_path)
        self.df, self.read_fasta, self.barcode_reads, self.counts, self.input_name, self.header = \
            self.capture_telomeres(read_path)
        read_to_cluster, self.df_anchors = self.cluster_telomeres()
        self.df['chrom'] = self.df['rname'].map(read_to_cluster)
        self.tvs_arr, self.tvs_read_counts = self.telo_variant_seq_analysis()

    def capture_telomeres(self, read_path):
        cap = TeloCapture(read_path, self.oligos, self.barcodes, self.data_type)
        if self.tsv_header:
            header = self.create_tvs_header(read_path, cap)
        else:
            header = None
        return cap.df, cap.read_fasta, cap.barcode_reads, cap.counts, cap.input_name, header

    def cluster_telomeres(self):
        clust = SubTeloClust(self.read_fasta, self.barcode_reads, self.chm13_path, self.cores)
        return clust.read_to_clust, clust.dfa

    def telo_variant_seq_analysis(self):
        tvs_arr, tvs_read_counts = tvs_analyzer(self.barcode_reads, self.df)
        return tvs_arr, tvs_read_counts

    # Prepare parse fasta file
    @staticmethod
    def parse_fasta(fasta):
        seq_dict = {}
        with open(fasta) as f:
            for i in f:
                if i.strip():
                    seq_name = i.strip()
                    seq = next(f).strip()
                    fa = Seq(seq)
                    seq_rc = str(fa.reverse_complement())
                    seq_dict[seq_name] = [seq, seq_rc]
        return seq_dict

    # Create TVS header
    def create_tvs_header(self, read_path, cap):
        h = []
        today = date.today()
        today_date = today.strftime("%d-%m-%Y")
        h.append('##date=%s\n' % today_date)
        h.append('##tool=Telomap-%s\n' % __version__)
        h.append('##source_reads=%s\n' % read_path)
        h.append('##data_type=%s\n' % self.data_type)
        h.append('##oligo_sequence=%s\n' % self.oligos)
        h.append('##barcode_sequence=%s\n' % self.barcodes)
        h.append('##oligo_location=%s\n' % cap.oligo_loc)
        h.append('##oligo_window=%s\n' % cap.oligo_window)
        h.append('##oligo_alignment_threshold=%s\n' % cap.oligo_align_threshold)
        h.append('##multi_oligo_alignment_threshold=%s\n' % cap.multi_oligo_align_threshold)
        h.append('##barcode_alignment_threshold=%s\n' % cap.barcode_align_threshold)
        h.append('##motifs=%s\n' % ','.join(cap.motifs))
        h.append('##motif_window=%s\n' % cap.motif_window)
        h.append('##minimum_motif_count=%s\n' % cap.motif_counts)
        h.append('##motif_direction=%s\n' % cap.motif_direction)
        h.append('##motif_end_length=%s\n' % cap.motif_end_len)
        h.append('##trf_binding_motif=%s\n' % cap.trf_motif)
        h.append('##trf_canonical_limit=%s\n' % cap.trf_canonical_limit)
        h.append('##minimum_telomere_length=%s\n' % cap.telo_min_len)
        h.append('##total_reads=%s\n' % cap.counts[0])
        h.append('##captured_reads=%s\n' % cap.counts[1])
        h.append('##captured_reads_with_motif=%s\n' % cap.counts[2])
        h.append('#RNAME\tBARCODE\tSTRAND\tCHROM\tRLEN\tTELOMERE_LEN\tTELOMERE_M_LEN\tTRF_COUNT\tTELOMERE_END\tMOTIF'
                 '\tSTART_JUNCT\tBARCODE_JUNC\tOLIGO\tOLIGO_SCORE\tBARCODE_SCORE\tNUM_PASS\tRQUAL\n')
        return h
