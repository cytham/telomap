# Cluster telomeric reads by subtelomeric regions
from Bio import SeqIO
from Bio import Align
import random
import pandas as pd
from collections import Counter
from datetime import datetime
from multiprocessing import Lock, Process, Queue, current_process
import time
import queue


class SubTeloClust:
    def __init__(self, read_fasta: dict, barcode_reads: dict, ref_path: str, cores: int, df):
        # self.begin_time = datetime.now()
        self.read_fasta = read_fasta
        self.barcode_reads = barcode_reads
        barcodes = [b for b in self.barcode_reads]
        self.cores = int(cores)
        self.df = df
        self.ref_path = ref_path
        self.all_motif_edit = {}
        self.all_anchor_groups = {}
        self.anchor_label = {}
        self.all_cons_seq = {}
        self.all_cons_grp = {}
        # self.telo_motif = 'TTAGGGTTAGGG'
        self.telo_buf = 12
        self.min_len = 100
        self.motif_len = 100  # Anchor length
        self.subtelo_len = 150
        self.min_hit = 10
        self.t = 0.95  # Alignment threshold
        self.seed = 7
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5
        self.cluster_max_ovlp = 0.1
        self.multi_process_cluster(barcodes)
        self.generate_consensus()
        self.map_anchor()
        dfa = self.create_anchor_df()
        self.dfa = self.unique_clust_label(dfa)
        self.read_to_clust = self.read_to_cluster()
        # self.end_time = datetime.now()
        # time_dif = self.end_time - self.begin_time
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        print(now, '- Clustering finished')
        # print(self.end_time - self.begin_time)

    # Multi-process clustering step
    def multi_process_cluster(self, tasks):
        number_of_processes = self.cores
        tasks_to_accomplish = Queue()
        tasks_that_are_done = Queue()
        processes = []
        for i in tasks:
            tasks_to_accomplish.put(i)
        # Creating processes
        for w in range(number_of_processes):
            p = Process(target=self.do_job, args=(tasks_to_accomplish, tasks_that_are_done))
            processes.append(p)
            p.start()
        while True:
            running = any(p.is_alive() for p in processes)
            while not tasks_that_are_done.empty():
                t = tasks_that_are_done.get()
                self.all_motif_edit[t[1]] = t[2]
                self.all_anchor_groups[t[1]] = t[3]
                # print(t[0])
            if not running:
                break
        # Completing process
        for p in processes:
            p.join()
        return True
    
    # Carry out queued job
    def do_job(self, tasks_to_accomplish, tasks_that_are_done):
        while True:
            try:
                '''
                    try to get task from the queue. get_nowait() function will 
                    raise queue.Empty exception if the queue is empty. 
                    queue(False) function would do the same task also.
                '''
                task = tasks_to_accomplish.get_nowait()
            except queue.Empty:
                break
            else:
                '''
                    if no exception has been raised, add the task completion 
                    message to task_that_are_done queue
                '''
                # print('Doing ', task)
                motif_edit, master_hit_grp = self.subtelo_cluster(task)
                tasks_that_are_done.put((task + ' is done by ' + current_process().name, task, motif_edit, master_hit_grp))
                time.sleep(.5)
        return True
    
    # Cluster telomeres by their subtelomere sequence
    def subtelo_cluster(self, sample):
        # reads = self.barcode_reads[sample]
        # print(len(reads))
        reads = self.df[(self.df['barcode'] == sample) & (self.df['motifs'].notnull())]['rname'].tolist()
        # print(len(reads))
        # print(sample)
        # Identify subtelo-telo junction by first occurence of "telo_motif" and extract anchor and subtelomere region
        motif_dict, subtelo_dict, fail_len_counts = self.subtelo_prep(reads)
        # Begin cluster algorithm
        hit_reads, motif_edit, master_hit_grp, readnames = self.clustering(motif_dict, subtelo_dict)
        self.all_motif_edit[sample] = motif_edit
        self.all_anchor_groups[sample] = master_hit_grp
        # p = 0
        # for i in master_hit_grp:
        #     for _ in master_hit_grp[i]:
        #         p += 1
        # fail_cluster_reads[sample] = readnames
        # print('%i reads left unclustered' % len(readnames))
        # print('Ended with %i reads clustered' % p)
        return motif_edit, master_hit_grp

    # Extract anchor and subtelomere region
    def subtelo_prep(self, reads):
        subtelo_dict = {}
        motif_dict = {}
        fail_len_counts = 0
        for r in reads:
            index = self.df[self.df['rname'] == r]['s_junct'].values[0]  # Get telomere start index
            # index = self.read_fasta[r].find(self.telo_motif)  # Get index of subtelo-telo junction according to telo_motif
            subtelo_dict[r] = self.read_fasta[r][:index + self.telo_buf]
            if len(subtelo_dict[r]) < self.min_len:  # Require subtelo region to be more than or equal to "min_len"
                _ = subtelo_dict.pop(r)
                fail_len_counts += 1
            else:
                subtelo_dict[r] = subtelo_dict[r][-self.subtelo_len:]
                motif_dict[r] = subtelo_dict[r][-self.motif_len:]
        return motif_dict, subtelo_dict, fail_len_counts

    # Perform clustering to discover anchors
    def clustering(self, motif_dict, subtelo_dict):
        random.seed(self.seed)
        readnames = sorted(list(motif_dict.keys()))
        random.shuffle(readnames)
        n = 1
        hit_reads = []
        prev_readnames = []
        motif_edit = {}
        master_hit_grp = {}
        while len(prev_readnames) != len(readnames):
            # print('Iteration %i' % n)
            # print('Using %i reads' % len(readnames))
            # Perform anchor sequence enrichment
            hit_dict = self.motif_enrich(readnames, motif_dict, subtelo_dict)
            hit_telo_dict, no_hits = self.remove_low_hits(hit_dict, motif_dict)
            # print(' Found %i potential anchors' % len(hit_telo_dict))
            # print('  %i reads with hits less than threshold' % len(no_hits))
            if not hit_telo_dict:  # If empty
                # print('No hits found, left %i anchors unclustered' % len(hit_dict))
                break
            subtelo_dict2 = {}
            for i in readnames:
                if i not in no_hits:
                    subtelo_dict2[i] = subtelo_dict[i]
            # Align reads to all anchors
            _ = []
            hit_out = self.motif_align(hit_telo_dict, subtelo_dict2)
            for i in hit_out:
                _.extend(list(hit_out[i]))
            # print('%i unique reads in clusters' % len(set(_)))
            # Cluster intersection analysis to promote cluster exclusivity
            ovlp_ind, hit_ind, read_order = self.cluster_intersect(hit_telo_dict, hit_out)
            hits = [read_order[x] for x in hit_ind]
            hit_reads.extend(hits)  # Record read names of exclusive anchor sequence
            # Calibrate motif sequences according to anchor sequence
            _ = []
            for i in hits:
                motif_edit[i] = motif_dict[i]
                master_hit_grp[i] = [i]
                _.extend(list(hit_out[i]))
                for j in hit_dict[i]:
                    rid = j[0]
                    motif_edit[rid] = self.calibrate_motif_index(j)
                    master_hit_grp[i].append(rid)
            # print('%i unique reads pass clusterings' % len(set(_)))
            prev_readnames = readnames
            readnames = []
            for i in ovlp_ind:
                read = read_order[i]
                readnames.extend(hit_out[read])
            readnames.extend(no_hits)
            readnames = sorted(list(set(readnames)))
            random.shuffle(readnames)
            n += 1
            if len(prev_readnames) < len(readnames):
                raise Exception('Error3')
        return hit_reads, motif_edit, master_hit_grp, readnames
    
    # To enrich for cluster anchor sequences from subtelomeres
    def motif_enrich(self, readnames, motif_dict, subtelo_dict):
        hit_dict = {}
        skip_dict = {}
        for f in readnames:
            if f not in skip_dict:
                hit_dict[f] = []
                skip_dict[f] = 0
                for t in readnames:
                    if t not in skip_dict:
                        a = self.aligner.align(motif_dict[f], subtelo_dict[t])
                        if len(a) > 0:  # If there's an alignment
                            if a[0].score >= self.motif_len * self.t:  # Check if alignment is good
                                hit_dict[f].append((t, a))
                                skip_dict[t] = 0
        return hit_dict

    # Remove clusters with low hits
    def remove_low_hits(self, hit_dict, motif_dict):
        hit_telo_dict = {}
        no_hits = []
        for i in hit_dict:
            if len(hit_dict[i]) >= self.min_hit:
                hit_telo_dict[i] = motif_dict[i]
            else:
                no_hits.append(i)
                no_hits.extend([x[0] for x in hit_dict[i]])
        return hit_telo_dict, no_hits

    # Clustering telomeric reads by aligning to anchor sequences with/without replacement
    def motif_align(self, motif_dict, subtelo_dict, replace=True):
        hit_out = {}
        if replace:
            for i in motif_dict:
                hit_out[i] = []
                for f in subtelo_dict:
                    a = self.aligner.align(motif_dict[i], subtelo_dict[f])
                    if a[0].score >= self.motif_len * self.t:
                        hit_out[i].append(f)
                hit_out[i] = set(hit_out[i])  # For intersection later
        else:
            skip_dict = {}
            for i in motif_dict:
                hit_out[i] = []
                for f in subtelo_dict:
                    if f not in skip_dict:
                        a = self.aligner.align(motif_dict[i], subtelo_dict[f])
                        if a[0].score >= self.motif_len * self.t:
                            hit_out[i].append(f)
                            skip_dict[f] = 0
                hit_out[i] = set(hit_out[i])  # For intersection later
        return hit_out

    # Check for anchor groups overlap
    def cluster_intersect(self, hit_telo_dict, hit_out):
        data_dict = {}
        for i in hit_telo_dict:
            data_dict[i] = []
            grp1 = hit_out[i]
            for j in hit_telo_dict:
                grp2 = hit_out[j]
                data_dict[i].append(round(len(grp1.intersection(grp2)) / len(grp1), 3))
        ind_no = list(range(len(hit_telo_dict)))
        dfi = pd.DataFrame(data_dict, index=ind_no)
        read_order = dfi.columns.tolist()
        dfi.columns = ind_no
        ovlp_ind, hit_ind = self.get_intersect_indexes(dfi, ind_no)  # Get index of exclusive or overlaping clusters
        return ovlp_ind, hit_ind, read_order

    # To intersect read clusters to find out cluster exclusivity
    def get_intersect_indexes(self, df, ind_no):
        ovlp_indexes = set()
        hit_indexes = []
        for i in ind_no:
            for j in ind_no:
                if i != j:
                    if df.loc[j, i] > self.cluster_max_ovlp:
                        ovlp_indexes.add(i)
                        ovlp_indexes.add(j)
        for i in ind_no:
            if i not in ovlp_indexes:
                hit_indexes.append(i)
        return ovlp_indexes, hit_indexes

    # Calibrate motif sequences according to anchor sequence
    def calibrate_motif_index(self, j):
        t_index = j[1][0].aligned[0]
        q_index = j[1][0].aligned[1]
        start_fill = t_index[0][0]
        end_fill = self.motif_len - t_index[-1][1]
        start = q_index[0][0]
        end = q_index[-1][1]
        que = j[1][0].query
        if len(que) >= end + end_fill:
            end_seq = que[end:end + end_fill]
        else:
            end_seq = 'N' * end_fill
        if start - start_fill >= 0:
            start_seq = que[start - start_fill:start]
        else:
            start_seq = 'N' * start_fill
        mid = ''
        for i in range(len(q_index) - 1):
            mid_n = 'N' * (t_index[i + 1][0] - t_index[i][1])
            mid = mid + que[q_index[i][0]:q_index[i][1]] + mid_n
        mid_end = que[q_index[-1][0]:q_index[-1][1]]
        new = start_seq + mid + mid_end + end_seq
        return new
        
    # Consensus generation and cluster merging
    def generate_consensus(self):
        for sample in self.all_anchor_groups:
            cons_dict = {}
            for i in self.all_anchor_groups[sample]:
                seqs = [str(self.all_motif_edit[sample][i])]
                for j in self.all_anchor_groups[sample][i]:
                    seqs.append(str(self.all_motif_edit[sample][j]))
                cons_dict[i], score = self.consensus(seqs)
            skip_dict = {}
            new_cons_seq = {}
            new_cons_grp = {}
            for i in cons_dict:
                if i not in skip_dict:
                    new_cons_seq[i] = cons_dict[i]
                    new_cons_grp[i] = list(self.all_anchor_groups[sample][i])
                    skip_dict[i] = 1
                    for j in cons_dict:
                        if i != j and j not in skip_dict:
                            a = self.aligner.align(cons_dict[i], cons_dict[j])
                            if a[0].score >= self.motif_len * self.t:
                                # Merge clusters
                                new_cons_grp[i].extend(list(self.all_anchor_groups[sample][j]))
                                skip_dict[j] = 1
                    new_cons_grp[i] = set(new_cons_grp[i])
            self.all_cons_seq[sample] = new_cons_seq
            self.all_cons_grp[sample] = new_cons_grp

    # Generate consensus anchor sequence for a cluster of reads
    @staticmethod
    def consensus(seqs):
        seq_len = len(seqs[0])
        n = len(seqs)
        data_dict = {}
        cons = ''
        score = []
        for i in range(seq_len):
            data_dict[i] = []
        for s in seqs:
            for i, j in enumerate(s):
                data_dict[i].append(j.upper())
        for i in range(seq_len):
            common = Counter(data_dict[i]).most_common(1)[0]
            cons += common[0]
            score.append(common[1] / n)
        return cons, score
    
    # Map anchor sequences to T2T-CHM13
    def map_anchor(self):
        ref_dict = {}
        for seq_record in SeqIO.parse(self.ref_path, "fasta"):
            ref_dict[seq_record.id] = seq_record.seq.upper()
        for sample in self.all_cons_seq:
            self.anchor_label[sample] = {}
            for anchor in self.all_cons_seq[sample]:
                self.anchor_label[sample][anchor] = []
                for ref in ref_dict:
                    if self.quick_align(self.all_cons_seq[sample][anchor], ref_dict[ref]):
                        self.anchor_label[sample][anchor].append(ref)
            a = 1
            for anchor in self.anchor_label[sample]:
                if len(self.anchor_label[sample][anchor]) == 1:
                    self.anchor_label[sample][anchor] = self.anchor_label[sample][anchor][0]
                elif 1 < len(self.anchor_label[sample][anchor]) <= 3:
                    self.anchor_label[sample][anchor] = '/'.join(self.anchor_label[sample][anchor][0:3])
                else:
                    self.anchor_label[sample][anchor] = 'U' + str(a)
                    a += 1
    
    # Find out if two sequences align well
    def quick_align(self, seq1, seq2):
        min_len = min(len(seq1), len(seq2))
        a = self.aligner.align(seq1, seq2)
        if a[0].score >= min_len * self.t:
            return True
        else:
            return False
    
    # Create anchor dataframe
    def create_anchor_df(self):
        data_dict = {'sample': [], 'anchor_read': [], 'anchor_seq': [], 'anchor_hits': [], 'chrom': []}
        for sample in self.barcode_reads:
            if sample in self.all_cons_seq:
                for anchor in self.all_cons_seq[sample]:
                    data_dict['sample'].append(sample)
                    data_dict['anchor_read'].append(anchor)
                    data_dict['anchor_seq'].append(str(self.all_cons_seq[sample][anchor]))
                    data_dict['anchor_hits'].append(len(self.all_cons_grp[sample][anchor]))
                    data_dict['chrom'].append(self.anchor_label[sample][anchor])
        dfa = pd.DataFrame.from_dict(data_dict)
        return dfa

    # Make cluster label unique
    def unique_clust_label(self, dfa):
        for sample in self.all_cons_seq:
            chroms = dfa[dfa['sample'] == sample]['chrom'].tolist()
            dups = [item for item, count in Counter(chroms).items() if count > 1]
            for d in dups:
                c = 1
                reads = dfa[(dfa['sample'] == sample) & (dfa['chrom'] == d)]['anchor_read'].tolist()
                for a in reads:
                    dfa.loc[dfa['anchor_read'] == a, 'chrom'] = d + '-' + str(c)
                    self.anchor_label[sample][a] = d + '-' + str(c)
                    c += 1
        return dfa

    # Make read to cluster dictionary
    def read_to_cluster(self):
        read_to_clust = {}
        for s in self.anchor_label:
            for a in self.anchor_label[s]:
                for r in self.all_cons_grp[s][a]:
                    read_to_clust[r] = self.anchor_label[s][a]
        return read_to_clust

'''
class Test:
    def __init__(self, input):
        self.out = {}
        self.multi_process_cluster(input, 10)
    # Test function
    def test_func(self, task):
        test_dict = {}
        for i in range(10000):
            _ = task + '-' + str(i)
            if i % 5000 == 0:
                test_dict[i] = _
        return test_dict
    # Carry out queued job
    def do_job(self, tasks_to_accomplish, tasks_that_are_done):
        while True:
            try:
                task = tasks_to_accomplish.get_nowait()
            except queue.Empty:
                break
            else:
                print('Doing ', task)
                test_dict = self.test_func(task)
                tasks_that_are_done.put((task + ' is done by ' + current_process().name, task, test_dict))
                time.sleep(.5)
        return True
    # Multi-process clustering step
    def multi_process_cluster(self, tasks, no_proc):
        number_of_processes = no_proc
        tasks_to_accomplish = Queue()
        tasks_that_are_done = Queue()
        processes = []
        for i in tasks:
            tasks_to_accomplish.put(i)
        # creating processes
        for w in range(number_of_processes):
            p = Process(target=self.do_job, args=(tasks_to_accomplish, tasks_that_are_done))
            processes.append(p)
            p.start()
            print('after start')
        # print the output
        # while not tasks_that_are_done.empty():
        #     t = tasks_that_are_done.get()
        #     self.out[t[1]] = t[2]
        #     print(t[0])
        while True:
            running = any(p.is_alive() for p in processes)
            while not tasks_that_are_done.empty():
                t = tasks_that_are_done.get()
                self.out[t[1]] = t[2]
                print(t[0])
            if not running:
                break
        # completing process
        for p in processes:
            p.join()
        return True

t = Test(['A', 'B', 'C', 'D'])

'''
