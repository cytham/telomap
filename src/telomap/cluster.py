# Cluster telomeric reads by subtelomeric regions
from Bio import SeqIO
from Bio import Align
import time
import random
import queue
import pandas as pd
import numpy as np
from collections import Counter
from datetime import datetime
from multiprocessing import Lock, Process, Queue, current_process
from sklearn.cluster import KMeans
from kneed import KneeLocator

class SubTeloClustWGS:
    def __init__(self, read_fasta: dict, barcodes: dict, ref_path: str, cores: int, df):
        self.read_fasta = read_fasta
        self.barcodes = barcodes
        self.ref_path = ref_path
        self.cores = int(cores)
        self.df = df
        self.telo_buf = 12
        self.min_len = 100
        self.motif_len = 100  # Anchor length
        self.subtelo_len = 150
        self.min_hit = 10
        self.t = 0.95
        self.align_thres = self.t * self.motif_len  # Alignment threshold
        self.min_read_progress = 50
        self.seed = 7
        # Set up aligner
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5
        self.output = {}
        self.tmp = {}
        for b in self.barcodes:
            self.tmp = {}  # Reset tmp
            read_ids = df[(df['barcode'] == b) & (df['motifs'].notnull())]['rname'].tolist()  # Get barcoded telomeric reads
            if not read_ids:
                print("Warning: barcode {} has no reads".format(b))
                continue
            self.tmp['motifs'], self.tmp['subtelo'], self.tmp['fail_len_counts'] = self.subtelo_prep(b)  # Prepare filtered motif and subtelo region
            data_list, data_array = self.create_array(self.tmp['motifs'])  # Create array for kmeans
            cluster_grps, self.tmp['n_clust'] = self.kmeans_clustering(data_list, data_array)  # kmeans clustering
            self.tmp['clustered'] = {}
            self.tmp['unclustered'] = {}
            self.multi_process(cluster_grps, self.clustering, self.cluster_combine)  # Multiprocess alignment clustering
            uncluster_grps = np.array_split(np.array([x for x in self.tmp['unclustered']]), self.cores)  # Create array for unclustered motifs
            self.multi_process(uncluster_grps, self.reclustering, self.recluster_combine)  # Multiprocess alignment reclustering
            self.calibrate()
            self.output[b] = self.tmp
        self.generate_consensus()
        self.map_anchor()
        dfa = self.create_anchor_df()
        self.dfa = self.unique_clust_label(dfa)
        self.read_to_clust = self.read_to_cluster()
        now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
        print(now, '- Clustering finished')
    
    # Prepare subtelomeric sequences
    def subtelo_prep(self, b):
        _df = self.df[(self.df['barcode'] == b) & (self.df['motifs'].notnull())][['rname', 's_junct']]
        _df['subtelo_index'] = _df['s_junct'] - 1 + self.telo_buf
        _df['subtelo_seq'] = [self.read_fasta[x][:y] for x, y in zip(_df.rname, _df.subtelo_index)]
        _df['subtelo'] = _df['subtelo_seq'].str[-self.subtelo_len:]
        _df['motif'] = _df['subtelo'].str[-self.motif_len:]
        _df2 = _df[_df['subtelo_seq'].str.len() >= self.min_len]
        fail_len_counts = _df.shape[0] - _df2.shape[0]
        subtelo_dict = pd.Series(_df2.subtelo.values,index=_df2.rname).to_dict()
        motif_dict = pd.Series(_df2.motif.values,index=_df2.rname).to_dict()
        return motif_dict, subtelo_dict, fail_len_counts
    
    # Create motif array
    @staticmethod
    def create_array(data):
        data_list = []
        data_array = []
        for i in data:
            data_list.append(i)
            data_array.append(data[i])
        return data_list, data_array
    
    # Perform initial kmean clustering
    def kmeans_clustering(self, data_list, data_array):
        matrix = np.asarray([np.frombuffer(s.encode(), dtype=np.uint8) for s in data_array])
        n_clusters, clusterid = self.elbow(matrix)
        # print('Using {} clusters'.format(n_clusters))
        # kmeans = KMeans(init="k-means++", 
        #                 n_clusters=n_clusts,
        #                 random_state=self.seed)
        # kmeans.fit(matrix)
        # clusterid = kmeans.labels_
        cluster_grps = []
        cluster_dict = {}
        cluster_dict2 = {}
        for i in range(len(clusterid)):
            try:
                cluster_dict[clusterid[i]].append(data_list[i])
                cluster_dict2[clusterid[i]].append(data_array[i])
            except KeyError:
                cluster_dict[clusterid[i]] = [data_list[i]]
                cluster_dict2[clusterid[i]] = [data_array[i]]
        for i in cluster_dict:
            cluster_grps.append(cluster_dict[i])
        return cluster_grps, n_clusters
    
    # Elbow curve to estimate k
    def elbow(self, data):
        n_clust_list = []
        sse = []
        cluster_labels = {}
        max_clust = 100
        min_clust = 10
        step_clust = 10
        if len(data) >= 100:
            for k in range(min_clust, max_clust + step_clust, step_clust):
                kmeans = KMeans(init="k-means++", n_clusters=k, max_iter=300, random_state=self.seed)
                kmeans.fit(data)
                n_clust_list.append(k)
                sse.append(kmeans.inertia_)
                cluster_labels[k] = kmeans.labels_
            kn = KneeLocator(n_clust_list, sse, curve='convex', direction='decreasing')
            n_clusters = kn.knee
            labels = cluster_labels[kn.knee]
        else:  # For small datasets, no clustering done
            n_clusters = 1
            kmeans = KMeans(init="k-means++", n_clusters=n_clusters, max_iter=300, random_state=self.seed)
            kmeans.fit(data)
            labels = kmeans.labels_
        return n_clusters, labels
    
    # Multi-process clustering step
    def multi_process(self, tasks, job, combine_out):
        number_of_processes = self.cores
        tasks_to_accomplish = Queue()
        tasks_that_are_done = Queue()
        processes = []
        for i in tasks:
            tasks_to_accomplish.put(i)
        # Add ending task
        tasks_to_accomplish.put(None)
        # Creating processes
        for w in range(number_of_processes):
            p = Process(target=self.do_job, args=(job, tasks_to_accomplish, tasks_that_are_done))
            processes.append(p)
            p.start()
        while True:
            running = any(p.is_alive() for p in processes)
            while not tasks_that_are_done.empty():
                out = tasks_that_are_done.get()
                combine_out(out)
            if not running:
                break
        # Completing process
        for p in processes:
            p.join()
        return True
    
    # Do queued job
    @staticmethod
    def do_job(job, tasks_to_accomplish, tasks_that_are_done):
        while True:
            task = tasks_to_accomplish.get(True)  # Remove and return an item from the queue
            if task is None:  # If queue ended
                # Add back ending task
                tasks_to_accomplish.put(None)
                break  # Finish process
            else:
                out = job(task)
                # print('Task done')
                tasks_that_are_done.put(out)  # Finished task and save output
                time.sleep(.5)       
        return True
    
    # Clustering job
    def clustering(self, readnames):
        # print('{} reads clustering'.format(len(readnames)))
        # records = len(readnames)
        random.seed(self.seed)
        readnames = sorted(list(readnames))
        random.shuffle(readnames)
        n = 1
        prev_readnames = []
        clustered = {}
        unclustered = []
        while len(prev_readnames) != len(readnames):
            # print('{} Iteration {}'.format(records, n))
            # print('{} Using {} reads'.format(records, len(readnames)))
            # Perform anchor sequence enrichment
            hit_dict = self.motif_enrich(readnames)
            # print(' Found %i initial anchors' % len(hit_dict))
            # hit_dicts.append(hit_dict)
            hit_telo_dict, no_hits = self.remove_low_hits(hit_dict)
            # print(' Found %i potential anchors' % len(hit_telo_dict))
            # print('  %i reads with hits less than threshold' % len(no_hits))
            for i in hit_telo_dict:
                clustered[i] = hit_dict[i]
            if not hit_telo_dict:  # If empty
                # print('No hits found, left %i anchors unclustered, break' % len(hit_dict))
                unclustered.extend(readnames)
                break
            prev_readnames = readnames
            readnames = no_hits
            if len(readnames) == 0:
                # print('0 unclustered reads, break')
                break
            if len(prev_readnames) - len(readnames) < self.min_read_progress:
                unclustered.extend(readnames)
                break
            readnames = sorted(list(set(readnames)))
            random.shuffle(readnames)
            n += 1
            if len(prev_readnames) < len(readnames):
                raise Exception('Error3')
        # print('Finished {} reads clustering'.format(records))
        return clustered, unclustered
    
    # To enrich for cluster anchor sequences from subtelomeres
    def motif_enrich(self, readnames):
        hit_dict = {}
        skip = set()
        for f in readnames:
            if f not in skip:
                hit_dict[f] = []
                skip.add(f)
                for t in readnames:
                    if t not in skip:
                        a = self.aligner.align(self.tmp['motifs'][f], self.tmp['subtelo'][t])
                        if len(a) > 0:  # If there's an alignment
                            if a[0].score >= self.align_thres:  # Check if alignment is good
                                hit_dict[f].append((t, a[0].aligned, a[0].query, a[0].score))
                                skip.add(t)
        return hit_dict
    
    # Remove clusters with low hits
    def remove_low_hits(self, hit_dict):
        hit_telo_dict = {}
        no_hits = []
        for i in hit_dict:
            if len(hit_dict[i]) >= self.min_hit:
                hit_telo_dict[i] = self.tmp['motifs'][i]
            else:
                no_hits.append(i)
                no_hits.extend([x[0] for x in hit_dict[i]])
        return hit_telo_dict, no_hits
    
    # Combine output from each process
    def cluster_combine(self, out):
        clustered = out[0]
        unclustered = out[1]
        for i in clustered:
            self.tmp['clustered'][i] = clustered[i]
            try:
                del self.tmp['unclustered'][i]
            except KeyError:
                pass
            for j in clustered[i]:
                try:
                    del self.tmp['unclustered'][j[0]]
                except KeyError:
                    pass
        for i in unclustered:
            self.tmp['unclustered'][i] = 0
    
    # Reclustering for unclustered reads
    def reclustering(self, readnames):
        hit_dict = {}
        for i in readnames:
            temp = []
            for j in self.tmp['clustered']:
                a = self.aligner.align(self.tmp['motifs'][j], self.tmp['subtelo'][i])
                if len(a) > 0:  # If there's an alignment
                    if a[0].score >= self.align_thres:  # Check if alignment is good
                        temp.append((j, a[0].score, a[0].aligned, a[0].query, a[0].score))
            if temp:
                top_hit = sorted(temp, key=lambda x:x[1])[-1]
                try:
                    hit_dict[top_hit[0]].append((i, top_hit[2], top_hit[3], top_hit[4]))
                except KeyError:
                    hit_dict[top_hit[0]] = [(i, top_hit[2], top_hit[3], top_hit[4])]
        return hit_dict
    
    # Combine output from multiple processes after reclustering
    def recluster_combine(self, out):
        for i in out:
            self.tmp['clustered'][i].extend(out[i])
            for j in out[i]:
                del self.tmp['unclustered'][j[0]]
    
    # Calibrate motif wrt master read
    def calibrate(self):
        self.tmp['groups'] = {}
        self.tmp['motif_cal'] = {}
        for i in self.tmp['clustered']:
            self.tmp['motif_cal'][i] = self.tmp['motifs'][i]
            self.tmp['groups'][i] = [i]
            for j in self.tmp['clustered'][i]:
                rid = j[0]
                self.tmp['motif_cal'][rid] = self.calibrate_motif_index(j)  # Calibrate reads to lead
                self.tmp['groups'][i].append(rid)
    
    # Calibrate motif sequences according to anchor sequence
    def calibrate_motif_index(self, j):
        t_index = j[1][0]
        q_index = j[1][1]
        start_fill = t_index[0][0]
        end_fill = self.motif_len - t_index[-1][1]
        start = q_index[0][0]
        end = q_index[-1][1]
        que = j[2]
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
        for b in self.output:
            # self.output[b]['cons']
            cons_dict = {}
            for i in self.output[b]['groups']:
                # seqs = [self.output[b]['motif_cal'][i]]
                seqs = []
                for j in self.output[b]['groups'][i]:
                    seqs.append(self.output[b]['motif_cal'][j])
                cons_dict[i], score = self.consensus(seqs)
            skip_dict = {}
            new_cons_seq = {}
            new_cons_grp = {}
            for i in cons_dict:
                if i not in skip_dict:
                    new_cons_seq[i] = cons_dict[i]
                    new_cons_grp[i] = self.output[b]['groups'][i]
                    skip_dict[i] = 0
                    for j in cons_dict:
                        if i != j and j not in skip_dict:
                            a = self.aligner.align(cons_dict[i], cons_dict[j])
                            if a[0].score >= self.align_thres:
                                # Merge clusters
                                # print("merging clusters {} ({}) and {} ({}) from barcode {}".format(i, len(self.output[b]['groups'][i]), j, len(self.output[b]['groups'][j]), b))
                                new_cons_grp[i].extend(self.output[b]['groups'][j])
                                skip_dict[j] = 0
                    new_cons_grp[i] = set(new_cons_grp[i])
            self.output[b]['cons'] = new_cons_seq
            self.output[b]['cons_grp'] = new_cons_grp
    
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
        for b in self.output:
            self.output[b]['anchor_label'] = {}
            for anchor in self.output[b]['cons']:
                self.output[b]['anchor_label'][anchor] = []
                for ref in ref_dict:
                    if self.quick_align(self.output[b]['cons'][anchor], ref_dict[ref]):
                        self.output[b]['anchor_label'][anchor].append(ref)
            a = 1
            for anchor in self.output[b]['anchor_label']:
                if len(self.output[b]['anchor_label'][anchor]) == 1:
                    self.output[b]['anchor_label'][anchor] = self.output[b]['anchor_label'][anchor][0]
                elif 1 < len(self.output[b]['anchor_label'][anchor]) <= 3:
                    self.output[b]['anchor_label'][anchor] = '/'.join(self.output[b]['anchor_label'][anchor][0:3])
                else:
                    self.output[b]['anchor_label'][anchor] = 'U' + str(a)
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
        data_dict = {'barcode': [], 'anchor_read': [], 'anchor_seq': [], 'anchor_hits': [], 'chrom': []}
        for b in self.output:
            # if b in self.output:
            for anchor in self.output[b]['cons']:
                data_dict['barcode'].append(b)
                data_dict['anchor_read'].append(anchor)
                data_dict['anchor_seq'].append(str(self.output[b]['cons'][anchor]))
                data_dict['anchor_hits'].append(len(self.output[b]['cons_grp'][anchor]))
                data_dict['chrom'].append(self.output[b]['anchor_label'][anchor])
        dfa = pd.DataFrame.from_dict(data_dict)
        return dfa
    
    # Make cluster label unique
    def unique_clust_label(self, dfa):
        for b in self.output:
            chroms = dfa[dfa['barcode'] == b]['chrom'].tolist()
            dups = [item for item, count in Counter(chroms).items() if count > 1]
            # print([(item, count) for item, count in Counter(chroms).items() if count > 1])
            # if dups:
                # print('We have dups from barcode {}. Dups={}'.format(b, dups))
            for d in dups:
                c = 1
                reads = dfa[(dfa['barcode'] == b) & (dfa['chrom'] == d)]['anchor_read'].tolist()
                for anchor in reads:
                    dfa.loc[dfa['anchor_read'] == anchor, 'chrom'] = d + '-' + str(c)
                    self.output[b]['anchor_label'][anchor] = d + '-' + str(c)
                    c += 1
        return dfa
    
    # Make read to cluster dictionary
    def read_to_cluster(self):
        read_to_clust = {}
        for b in self.output:
            for a in self.output[b]['anchor_label']:
                for r in self.output[b]['cons_grp'][a]:
                    read_to_clust[r] = self.output[b]['anchor_label'][a]
        return read_to_clust
