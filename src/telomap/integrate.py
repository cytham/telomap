# Main pipeline of Telomap

import os
from Bio.Seq import Seq
from datetime import date
from datetime import datetime
from .capture import TeloCapture
from .cluster import SubTeloClust
from .tvs import tvs_analyzer
from .version import __version__


class TeloMap:

    def __init__(self, mode:str, read_path: str, oligo_path: str, barcode_path: str, cores: int, sample_name: str, motif: str, oligoscore: float, barscore: float, tsv_header=False):
        self.chm13_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ref/t2t-chm13-subtelo-1000-60.fa')
        # self.data_type = data_type  # fasta, fastq, bam, pacbio-bam
        self.mode = mode
        self.cores = cores
        self.tsv_header = tsv_header
        self.input_name = sample_name
        self.motif = motif
        self.oligoscore = oligoscore
        self.barscore = barscore
        # Parse capture oligo sequences
        # self.oligos = self.parse_fasta(oligo_path)
        # Parse barcode sequences
        # self.barcodes = self.parse_fasta(barcode_path)
        # Check mode and parse oligo and barcode sequences
        self.oligos, self.barcodes = self.check_mode(oligo_path, barcode_path)
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        print(now + ' - Capturing telomeric reads')
        self.df, self.read_fasta, self.barcode_reads, self.counts, self.header = \
            self.capture_telomeres(read_path)
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        print(now + ' - Clustering telomeric reads')
        read_to_cluster, self.df_anchors = self.cluster_telomeres()
        self.df['chrom'] = self.df['rname'].map(read_to_cluster)
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        print(now + ' - Analyzing TVS')
        self.tvs_arr, self.tvs_read_counts = self.telo_variant_seq_analysis()

    def capture_telomeres(self, read_path):
        cap = TeloCapture(read_path, self.oligos, self.barcodes, self.input_name, self.mode, self.motif, self.oligoscore, self.barscore)
        if self.tsv_header:
            header = self.create_tvs_header(read_path, cap)
        else:
            header = None
        return cap.df, cap.read_fasta, cap.barcode_reads, cap.counts, header

    def cluster_telomeres(self):
        clust = SubTeloClust(self.read_fasta, self.barcode_reads, self.chm13_path, self.cores)
        return clust.read_to_clust, clust.dfa

    def telo_variant_seq_analysis(self):
        tvs_arr, tvs_read_counts = tvs_analyzer(self.barcode_reads, self.df)
        return tvs_arr, tvs_read_counts

    # Check mode and barcodes
    def check_mode(self, oligo_path, barcode_path):
        if self.mode == 'wgs':
            oligos = {'>NA': []}
            barcodes = {'>NA': []}
        elif self.mode == 'telobait':
            oligos = self.parse_fasta(oligo_path)
            barcodes = self.parse_fasta(barcode_path)
            if not oligos:
                raise Exception('Error: Oligos FASTA input missing')
            if not barcodes:
                raise Exception('Error: Barcodes FASTA input missing')   
        return oligos, barcodes
    
    # Prepare parse fasta file
    @staticmethod
    def parse_fasta(fasta):
        seq_dict = {}
        if not fasta:
            return None
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
        h.append('##run_mode=%s\n' % self.mode)
        h.append('##oligo_sequence=%s\n' % self.oligos)
        h.append('##barcode_sequence=%s\n' % self.barcodes)
        h.append('##oligo_location=%s\n' % cap.oligo_loc)
        h.append('##oligo_window=%s\n' % cap.oligo_window)
        h.append('##oligo_alignment_threshold=%s\n' % cap.oligo_align_threshold)
        h.append('##multi_oligo_alignment_threshold=%s\n' % cap.multi_oligo_align_threshold)
        h.append('##barcode_alignment_threshold=%s\n' % cap.barcode_align_threshold)
        h.append('##motifs=%s\n' % ','.join(cap.motifs))
        h.append('##motif_window=%s\n' % cap.motif_window)
        h.append('##motif_window_wgs=%s\n' % cap.motif_window_wgs)
        h.append('##minimum_motif_count=%s\n' % cap.motif_counts)
        h.append('##motif_direction=%s\n' % cap.motif_direction)
        h.append('##motif_end_length=%s\n' % cap.motif_end_len)
        h.append('##trf_binding_motif=%s\n' % cap.trf_motif)
        h.append('##trf_canonical_limit=%s\n' % cap.trf_canonical_limit)
        h.append('##minimum_telomere_length=%s\n' % cap.telo_min_len)
        h.append('##total_reads=%s\n' % cap.counts[0])
        h.append('##captured_reads=%s\n' % cap.counts[1])
        h.append('##captured_reads_with_motif=%s\n' % cap.counts[3])
        h.append('#RNAME\tOLIGO\tBARCODE\tSTRAND\tMOTIF\tRLEN\tRAW_TELOMERE_LEN\tTELOMERE_LEN\tTELOMERE_END_MOTIF\tTRF_COUNT'
                 '\tCHROM\tTELOMERE_START\tTELOMERE_END\tNUM_PASS\tRQUAL\n')
        return h
