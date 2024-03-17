# Capture telomeric reads

import os
import re
import pysam
from datetime import datetime
import pandas as pd
from Bio.Seq import Seq
from Bio import Align


class TeloCapture:

    def __init__(self, read_path, oligos, barcodes, sample_name, mode, motif, oligoscore, barscore):
        self.read_path = read_path
        self.oligos = oligos
        self.barcodes = barcodes
        self.analysis_mode = mode    
        self.data_type = self.detect_data_type()
        self.input_name = sample_name
        # self.pacbio_hifi_rq = 0.99  # Q20
        self.oligo_loc = '3prime'  # Location of capture oligo sequence on telomere (3prime/5prime)
        self.oligo_window = 150  # Length of window (bp) to search for capture oligo sequence
        self.oligo_align_threshold = oligoscore  # Minimum alignment score percentage required for capture oligo sequence match [100%]
        self.multi_oligo_align_threshold = 0.8  # Minimum alignment score fraction for additional capture oligo seq match [80%]
        self.barcode_align_threshold = barscore  # Minimum alignment score percentage required for barcode sequence match [100%]
        self.motif_window = 50  # Length of window (bp) to search for telomeric motif adjacent to capture oligo, or end of read for WGS
        self.motif_window_wgs = 150  # Length of window (bp) to search for telomeric motif for WGS mode
        self.motif_counts = 2  # Minimum number of telomeric motifs required to pass (adjacent to capture oligo)
        self.motif_direction = 'upstream'  # Location of telomeric motif with respect to capture oligo (upstream/downstream)
        self.telo_start_sequence = 'TTAGGGTTAGGG'  # Defines the 5 prime start of telomere
        self.motif_end_len = 6  # Length of telomere end sequence motif
        self.max_gap_size = 6  # Maximum gap size to be considered for gap sequence analysis
        self.motifs = [motif]  # Specify telomeric motifs
        self.telo_min_len = 0
        self.trf_motif = 'TTAGGGTTA'  # Specify TRF1/2 binding motif
        self.trf_canonical_limit = 1000
        self.end_motifs = ['TTAGGG', 'GTTAGG', 'GGTTAG', 'GGGTTA', 'AGGGTT', 'TAGGGT']  # Specify telomere end motifs
        self.aligner1 = Align.PairwiseAligner()
        self.aligner1.mode = 'local'
        self.aligner1.match_score = 1
        self.aligner1.mismatch_score = 0
        self.aligner1.open_gap_score = -100
        self.aligner1.extend_gap_score = -100
        self.aligner2 = Align.PairwiseAligner()
        self.aligner2.mode = 'local'
        self.aligner2.match_score = 1
        self.aligner2.mismatch_score = 0
        self.aligner2.open_gap_score = -1
        self.aligner2.extend_gap_score = -1
        # self.input_name = self.get_input_name()
        self.data = self.read_input()
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        if self.analysis_mode == 'telobait':
            print(now, '- Running telobait capture mode')
            self.df, self.read_fasta, self.barcode_reads, self.counts = self.telo_capture()
        elif self.analysis_mode == 'wgs':
            print(now, '- Running WGS capture mode')
            self.df, self.read_fasta, self.barcode_reads, self.counts = self.telo_capture_wgs()

    def telo_capture_wgs(self):
        total_no = 0
        capture_no = 0
        multi_no = 0
        telomere_no = 0
        read_fasta = {}
        barcode_reads = {x.strip('>'): [] for x in self.barcodes}
        barcode_name = list(self.barcodes.keys())[0].strip('>')
        oligo_name = list(self.oligos.keys())[0].strip('>')
        df_dict = {'rname': [], 'read_len': [], 'num_pass': [], 'read_qual': [], 'strand': [], 'oligo': [],
                   'barcode': [], 'oscore': [], 'bscore': [], 'e_junct': [], 'motifs': [], 'telo_end': [], 'telo_len': [],
                   'gap_size': [], 'gap_loc': [], 'gap_dist': [], 'gap_seq': [], 'telo_len_wgap': [],
                   's_junct': [], 'telo_motif': [], 'trf_motif': [], 'trf_count': []}
        for seg in self.data:
            strand = oligo = barcode = oligo_score = bar_score = junct = motif = telo_end = telo_len_no_gap = telo_len_wgap = \
                gap_sizes = gap_locs = gap_dists = gap_seq = telo_start_index = telo_end_index = telo_motif_indexes = trf_motif_indexes = trf_count = \
                None
            total_no += 1
            capture_no += 1
            if total_no % 100000 == 0:
                now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
                print(now, ' - Analyzed %i sequences' % total_no)
            # if total_no > 100000:
            #     break
            qname = seg.query_name
            fasta = seg.query_sequence
            if self.data_type == 'bam':
                try:  # If pacbio-bam
                    np = seg.get_tag('np')
                    rq = round(seg.get_tag('rq'), 4)
                    # if rq <= self.pacbio_hifi_rq:
                    #     continue
                except KeyError:
                    np = None
                    rq = None
            else:
                np = None
                rq = None
            output = self.motif_finder_wgs(fasta)
            if output:  # If read contains telomeric motifs
                strand = '+'
            else:
                fa = Seq(fasta)
                fasta = str(fa.reverse_complement())
                output = self.motif_finder_wgs(fasta)
                if output:  # If read contains telomeric motifs
                    strand = '-'
            if output:
                motif = output[0]
                telo_start_index = fasta.find(self.telo_start_sequence)
                if telo_start_index > -1:  # If telomere start sequence is found
                    # Retrieve last motif sequence
                    last_motifs_fa = fasta[output[1][-1][0]:output[1][-1][0] + len(motif)*2]
                    telo_end = self.wgs_get_end_motif(last_motifs_fa, motif)
                    telo_end_index = output[1][-1][0] + len(telo_end)
                    if telo_end_index - telo_start_index >= self.telo_min_len:  # Minimum telomere length filter
                        telomere_no += 1
                        barcode_reads[barcode_name].append(qname)
                        telo_len_wgap = telo_end_index - telo_start_index
                        telo_motif_indexes, telo_len_no_gap = self.get_telo_len_nogap(motif, output[1], telo_start_index, telo_end_index)
                        # Analyze telomeric gaps
                        gap_sizes, gap_locs, gap_dists = self.identify_gaps(telo_motif_indexes)
                        gap_seq = self.get_gap_seq(fasta, gap_sizes, gap_locs)
                        # Analyze TRF1/TRF2 binding motif
                        trf_motif_indexes = [(_.start(), _.end()) for _ in re.finditer(self.trf_motif, fasta)]
                        if telo_len_no_gap <= self.trf_canonical_limit:
                            pass
                        else:
                            trf_count = self.count_trf(len(motif), telo_start_index, telo_motif_indexes, trf_motif_indexes)
                        telo_start_index += 1  # make 1-based
            # Record data
            read_fasta[qname] = fasta  # Record FASTA
            df_dict['rname'].append(qname)
            df_dict['read_len'].append(len(fasta))
            df_dict['num_pass'].append(np)
            df_dict['read_qual'].append(rq)
            df_dict['oligo'].append(oligo_name)
            df_dict['barcode'].append(barcode_name)
            df_dict['oscore'].append(oligo_score)
            df_dict['bscore'].append(bar_score)
            df_dict['strand'].append(strand)
            df_dict['motifs'].append(motif)
            df_dict['telo_end'].append(telo_end)
            df_dict['telo_len'].append(telo_len_no_gap)
            df_dict['telo_len_wgap'].append(telo_len_wgap)
            df_dict['gap_size'].append(gap_sizes)
            df_dict['gap_loc'].append(gap_locs)
            df_dict['gap_dist'].append(gap_dists)
            df_dict['gap_seq'].append(gap_seq)
            df_dict['s_junct'].append(telo_start_index)
            df_dict['e_junct'].append(telo_end_index)
            df_dict['telo_motif'].append(telo_motif_indexes)
            df_dict['trf_motif'].append(trf_motif_indexes)
            df_dict['trf_count'].append(trf_count)
        # Create dataframe
        df = pd.DataFrame.from_dict(df_dict)
        # Convert some columns from float to int64 to avoid decimal
        for col in ['s_junct', 'e_junct', 'telo_len', 'telo_len_wgap', 'trf_count']:
            # df[col] = df[col].fillna(-1).astype('int64').replace(-1, None)
            df[col] = df[col].astype(float).fillna(-1).astype('int64').replace(-1, None)
        return df, read_fasta, barcode_reads, (total_no, capture_no, multi_no, telomere_no)
    
    def telo_capture(self):
        total_no = 0
        capture_no = 0
        multi_no = 0
        telomere_no = 0
        read_fasta = {}
        df_dict = {'rname': [], 'read_len': [], 'num_pass': [], 'read_qual': [], 'strand': [], 'oligo': [],
                   'barcode': [], 'oscore': [], 'bscore': [], 'e_junct': [], 'motifs': [], 'telo_end': [], 'telo_len': [],
                   'gap_size': [], 'gap_loc': [], 'gap_dist': [], 'gap_seq': [], 'telo_len_wgap': [],
                   's_junct': [], 'telo_motif': [], 'trf_motif': [], 'trf_count': []}
        barcode_reads = {x.strip('>'): [] for x in self.barcodes}
        for seg in self.data:
            strand = oligo = barcode = oligo_score = bar_score = junct = motif = telo_end = telo_len_no_gap = telo_len_wgap = \
                gap_sizes = gap_locs = gap_dists = gap_seq = telo_start_index = telo_motif_indexes = trf_motif_indexes = trf_count = \
                None
            total_no += 1
            if total_no % 100000 == 0:
                now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
                print(now, ' - Analyzed %i sequences' % total_no)
            # if total_no > 100000:
            #     break
            qname = seg.query_name
            fasta = seg.query_sequence
            if self.data_type == 'bam':
                try:  # If pacbio-bam
                    np = seg.get_tag('np')
                    rq = round(seg.get_tag('rq'), 4)
                    # if rq <= self.pacbio_hifi_rq:
                    #     continue
                except KeyError:
                    np = None
                    rq = None
            else:
                np = None
                rq = None
            res = self.oligo_bar_align(fasta)
            if res:  # If oligo and barcode successfully aligns
                if res != 'multi':  # If sequence aligns to a single oligo and barcode
                    capture_no += 1
                    oligo = res[0].strip('>')
                    barcode = res[1].strip('>')
                    strand = res[2]
                    oligo_score = str(int(res[3]))
                    bar_score = str(int(res[4]))
                    junct = int(res[5])
                    fasta = self.get_strand(fasta, strand)
                    motif_out = self.motif_finder(fasta, junct)
                    if motif_out:  # If read contains telomeric motifs
                        # Retrieve telomere end sequence
                        motif = motif_out[0]
                        telo_start_index = fasta.find(self.telo_start_sequence)
                        if telo_start_index > -1:  # If telomere start sequence is found
                            if junct - telo_start_index >= self.telo_min_len:  # Minimum telomere length filter
                                telomere_no += 1
                                # Analyze capture sequence adjacent to capture oligo
                                telo_end = self.get_end_motif(fasta, junct)
                                barcode_reads[barcode].append(qname)
                                telo_len_wgap = junct - telo_start_index
                                telo_motif_indexes, telo_len_no_gap = self.get_motif_index(motif, fasta, telo_start_index, junct)
                                # Analyze telomeric gaps
                                gap_sizes, gap_locs, gap_dists = self.identify_gaps(telo_motif_indexes)
                                gap_seq = self.get_gap_seq(fasta, gap_sizes, gap_locs)
                                # Analyze telomeric end sequence
                                # telo_end = self.get_end_motif(fasta, junct)
                                # Analyze TRF1/TRF2 binding motif
                                trf_motif_indexes = [(_.start(), _.end()) for _ in re.finditer(self.trf_motif, fasta)]
                                if telo_len_no_gap <= self.trf_canonical_limit:
                                    pass
                                else:
                                    trf_count = self.count_trf(len(motif), telo_start_index, telo_motif_indexes, trf_motif_indexes)
                                telo_start_index += 1  # make 1-based
                elif res == 'multi':  # If sequence aligns to multiple oligos or barcodes
                    oligo = 'Multi'
                    multi_no += 1
                else:
                    raise Exception('Error: Unrecognised alignment result %s' % str(res))
            # Record data
            read_fasta[qname] = fasta  # Record FASTA
            df_dict['rname'].append(qname)
            df_dict['read_len'].append(len(fasta))
            df_dict['num_pass'].append(np)
            df_dict['read_qual'].append(rq)
            df_dict['oligo'].append(oligo)
            df_dict['barcode'].append(barcode)
            df_dict['oscore'].append(oligo_score)
            df_dict['bscore'].append(bar_score)
            df_dict['strand'].append(strand)
            df_dict['motifs'].append(motif)
            df_dict['telo_end'].append(telo_end)
            df_dict['telo_len'].append(telo_len_no_gap)
            df_dict['telo_len_wgap'].append(telo_len_wgap)
            df_dict['gap_size'].append(gap_sizes)
            df_dict['gap_loc'].append(gap_locs)
            df_dict['gap_dist'].append(gap_dists)
            df_dict['gap_seq'].append(gap_seq)
            df_dict['s_junct'].append(telo_start_index)
            df_dict['e_junct'].append(junct)
            df_dict['telo_motif'].append(telo_motif_indexes)
            df_dict['trf_motif'].append(trf_motif_indexes)
            df_dict['trf_count'].append(trf_count)
        # Create dataframe
        df = pd.DataFrame.from_dict(df_dict)
        # Convert some columns from float to int64 to avoid decimal
        for col in ['s_junct', 'e_junct', 'telo_len', 'telo_len_wgap', 'trf_count']:
            # df[col] = df[col].fillna(-1).astype('int64').replace(-1, None)
            df[col] = df[col].astype(float).fillna(-1).astype('int64').replace(-1, None)
        return df, read_fasta, barcode_reads, (total_no, capture_no, multi_no, telomere_no)
    
    # Detect data type
    def detect_data_type(self):
        if any(self.read_path.lower().endswith(s) for s in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']):
            return 'fastx'
        elif self.read_path.lower().endswith('.bam'):
            return 'bam'
        else:
            raise Exception('Error: Input data type is not accepted. Only "fasta", "fastq", or "bam" data types are accepted')

    # Read input file
    def read_input(self):
        if self.data_type == 'fastx':
            data = pysam.FastxFile(self.read_path)
            data = self.fastatx_generator(data)
        elif self.data_type == 'bam':
            save = pysam.set_verbosity(0)  # Suppress BAM index missing warning
            data = pysam.AlignmentFile(self.read_path, 'rb', check_sq=False)
            pysam.set_verbosity(save)
        return data

    # Generator for Fastatx
    def fastatx_generator(self, data):
        for d in data:
            out = self.Fastatx(d)
            yield out

    # Create class to match fasta object to bam
    class Fastatx:

        # Constructor method of inner class
        def __init__(self, d):
            self.query_name = d.name
            self.query_sequence = d.sequence

    # Align oligos and barcodes to reads
    def oligo_bar_align(self, fasta):
        prime3_region, prime5_region = self.get_oligo_region(fasta)
        hit_barcode = hit_score_b = junction = None
        for oligo in self.oligos:
            max_score = len(self.oligos[oligo][0])
            alignments = self.aligner1.align(prime3_region, self.oligos[oligo][0])
            score_o = alignments[0].score
            start = alignments[0].aligned[0][0][0]
            if score_o / max_score >= self.oligo_align_threshold:
                if len(alignments) > 1:  # Omit read if more than one capture oligo is found in oligo window
                    return 'multi'
                # Check additional oligo presence in read
                if len(fasta) >= self.oligo_window:
                    remaining_seq = fasta[:len(fasta) - self.oligo_window + start]
                else:
                    remaining_seq = fasta[:start]
                for _oligo in self.oligos:
                    if self.check_add_oligo(fasta, remaining_seq, self.oligos[_oligo][0], self.oligos[_oligo][1]):
                        return 'multi'
                # Screen through all barcodes and break if two hits
                hit = 0
                for barcode in self.barcodes:
                    max_score = len(self.barcodes[barcode][0])
                    alignments = self.aligner1.align(prime3_region, self.barcodes[barcode][0])
                    score_b = alignments[0].score
                    if score_b / max_score >= self.barcode_align_threshold:
                        hit += 1
                        if hit > 1:
                            return 'multi'  # More than one barcode aligns
                        else:
                            start = alignments[0].aligned[0][0][0]
                            if len(fasta) >= self.oligo_window:
                                junction = len(fasta) - self.oligo_window + start
                            else:
                                junction = start
                            hit_barcode = barcode
                            hit_score_b = score_b
                            # if fasta[junction:junction + max_score] != self.barcodes[barcode][0]:
                            #     # Barcode sequence check failed
                            #     raise Exception('Error: Barcode sequence %s does not match' % barcode)
                if hit == 1:
                    return [oligo, hit_barcode, '+', score_o, hit_score_b, junction]
            else:
                # Try reverse complement oligo
                alignments = self.aligner1.align(prime5_region, self.oligos[oligo][1])
                score_o = alignments[0].score
                end = alignments[0].aligned[0][-1][1]
                if score_o / max_score >= self.oligo_align_threshold:
                    if len(alignments) > 1:  # Omit read if more than one capture oligo is found in oligo window
                        return 'multi'
                    # Check additional oligo presence in read
                    remaining_seq = fasta[end:]
                    for _oligo in self.oligos:
                        if self.check_add_oligo(fasta, remaining_seq, self.oligos[_oligo][1], self.oligos[_oligo][0]):
                            return 'multi'
                    # Screen through all barcodes and break if two hits
                    hit = 0
                    for barcode in self.barcodes:
                        max_score = len(self.barcodes[barcode][0])
                        alignments = self.aligner1.align(prime5_region, self.barcodes[barcode][1])
                        score_b = alignments[0].score
                        if score_b / max_score >= self.barcode_align_threshold:
                            hit += 1
                            if hit > 1:  # More than one barcode aligns
                                return 'multi'
                            else:
                                end = alignments[0].aligned[0][-1][1]
                                junction = len(fasta) - end
                                hit_barcode = barcode
                                hit_score_b = score_b
                                # if fasta[end - max_score:end] != self.barcodes[barcode][1]:
                                #     # Barcode sequence check failed
                                #     raise Exception('Error: Barcode sequence %s does not match' % barcode)
                    if hit == 1:
                        return [oligo, hit_barcode, '-', score_o, hit_score_b, junction]

    # Get the regions of the sequence where the capture oligo would be found
    def get_oligo_region(self, fasta):
        if self.oligo_loc == '3prime':
            forward_region = fasta[-1 * self.oligo_window:]
            reverse_region = fasta[:self.oligo_window]
        else:
            forward_region = fasta[:self.oligo_window]
            reverse_region = fasta[-1 * self.oligo_window:]
        return forward_region, reverse_region

    # Check if sequence contain additional capture oligos
    def check_add_oligo(self, fasta, region, oligo1, oligo2):
        max_score = len(oligo1)
        if len(region) > 0:
            _alignments = self.aligner2.align(region, oligo1)
            if len(_alignments) > 0:
                if _alignments[0].score / max_score >= self.multi_oligo_align_threshold:
                    return True
        _alignments = self.aligner2.align(fasta, oligo2)  # Try un-reverse complement oligo
        if _alignments[0].score / max_score >= self.multi_oligo_align_threshold:
            return True

    # Get correct strand of sequence
    @staticmethod
    def get_strand(fasta, strandness):
        if strandness == '+':
            pass
        else:  # '-' Reverse complement
            fa = Seq(fasta)
            fasta = str(fa.reverse_complement())
        return fasta

    # Find telomeric motifs in sequence
    def motif_finder(self, fasta, junction):
        # Only consider motif in forward strand
        if self.motif_window == -1:  # If search_region is not use, search entire FASTA
            motif_window = len(fasta)
        else:
            motif_window = self.motif_window
        if self.motif_direction == 'upstream':  # Telomere motif upstream of barcode
            search_region = fasta[junction - motif_window:junction]
        else:  # Telomere motif downstream of barcode
            search_region = fasta[junction:junction + motif_window]
        output = []
        for motif in self.motifs:
            matches = [_.start() for _ in re.finditer(motif, search_region)]
            if len(matches) >= self.motif_counts:
                output = [motif, len(matches)]
                break  # Only choose 1st motif
        return output
    
    # Find telomeric motifs in sequence (WGS)
    def motif_finder_wgs(self, fasta, telo_loc='3prime'):
        if telo_loc == '3prime':  # Telomere at 3prime end
            search_region = fasta[len(fasta)-self.motif_window_wgs:]
        elif telo_loc == '5prime':
            search_region = fasta[:self.motif_window_wgs]
            raise Exception('Error: 5prime setting currently disabled.')
        output = []
        for motif in self.motifs:
            matches = [_.start() for _ in re.finditer(motif, search_region)]  # Quick search
            if len(matches) >= self.motif_counts:
                motif_indexes = [(_.start(), _.end()) for _ in re.finditer(motif, fasta)]  # Elaborate search
                output = [motif, motif_indexes]
                break  # Only choose 1st motif
        return output

    # Get end motif (WGS)
    @staticmethod
    def wgs_get_end_motif(fasta, motif):
        for i in reversed(range(len(motif))):
            end_motif = motif+motif[0:i]
            if fasta.startswith(end_motif):
                return end_motif[-len(motif):]

    # Get telomere length (WGS)
    @staticmethod
    def get_telo_len_nogap(motif, motif_indexes, telo_start, telo_end_index):
        telo_motif_indexes = []
        count = 0
        for k in motif_indexes:
            if k[0] >= telo_start and k[1] <= telo_end_index:
                count += 1
                telo_motif_indexes.append(k)
        telo_len_no_gap = count * len(motif)
        return telo_motif_indexes, telo_len_no_gap
    
    # Get indexes of telomeric motifs and calculate telomere length without gaps
    @staticmethod
    def get_motif_index(motif, fasta, telo_start, telo_end):
        telo_motif_indexes = []
        motif_indexes = [(_.start(), _.end()) for _ in re.finditer(motif, fasta)]
        count = 0
        for k in motif_indexes:
            if k[0] >= telo_start and k[1] <= telo_end:
                count += 1
                telo_motif_indexes.append(k)
        telo_len_no_gap = count * len(motif)
        return telo_motif_indexes, telo_len_no_gap

    # Analyze gaps/TVS in telomeric read
    @staticmethod
    def identify_gaps(indexes, min_gap=1):
        gap_sizes = []
        gap_loc = []
        for i in range(len(indexes) - 1):
            gap_size = indexes[i + 1][0] - indexes[i][1]
            if gap_size >= min_gap:
                gap_sizes.append(gap_size)
                gap_loc.append([indexes[i][1], indexes[i + 1][0]])
        start_telo = indexes[0][0]  # Start of telomere
        gap_dist = [x[1] - start_telo for x in gap_loc]
        return gap_sizes, gap_loc, gap_dist

    # Get gap/TVS sequence
    def get_gap_seq(self, fasta, gap_sizes, gap_locs):
        gap_seq = []
        for i in range(len(gap_sizes)):
            if gap_sizes[i] <= self.max_gap_size:
                gap_seq.append(fasta[gap_locs[i][0]: gap_locs[i][1]])
            else:
                gap_seq.append(None)
        return gap_seq

    # Find end motif adjacent of capture oligo
    def get_end_motif(self, fasta, junction):
        if self.motif_direction == 'upstream':  # Telomere motif upstream of barcode
            return fasta[junction - self.motif_end_len:junction]
        else:  # Telomere motif downstream of barcode
            return fasta[junction:junction + self.motif_end_len]

    # Count number of TRF1/TRF2 binding motifs for each read
    def count_trf(self, motif_len, telo_start, telo_indexes, trf1_indexes):
        telo_end = 0
        no_of_motifs = self.trf_canonical_limit // motif_len
        for j in range(len(telo_indexes)):
            if telo_indexes[j][0] == telo_start:
                telo_end = telo_indexes[j + no_of_motifs - 1][1]
                break
        if telo_end == 0:
            raise Exception('Error: Start motif not found')
        trf_count = 0
        for j in trf1_indexes:
            if j[0] >= telo_start and j[1] <= telo_end:
                trf_count += 1
        return trf_count
