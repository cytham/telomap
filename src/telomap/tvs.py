# Analyze telomere variant sequence (TVS) per chromosome end

import pandas as pd
import numpy as np
from datetime import datetime
from sklearn.cluster import KMeans
from scipy.spatial import distance

def tvs_analyzer(barcodes, df, telo_len=3000):
    now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
    print(now, ' - TVS analysis start')
    tvs_arr = {}
    tvs_read_counts = {}
    for b in barcodes:
        # print('barcode: {}'.format(b))
        tvs_arr[b] = {}
        tvs_read_counts[b] = {}
        chroms = df[(df['barcode'] == b) & (df['telo_len'].notnull())]['chrom'].unique().tolist()
        chroms = [x for x in chroms if pd.notnull(x)]  # Remove nan
        for c in chroms:
            n = 0
            # print('chrom: {}'.format(c))
            start_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['s_junct'].tolist()
            loc_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['gap_loc'].tolist()
            len_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['telo_len_wgap'].tolist()
            longest = int(max(len_list))
            arr = np.empty([len(start_list), longest], dtype=np.int8)
            arr.fill(0)
            # print('Length: {}'.format(len(start_list)))
            for i in range(len(start_list)):
                n += 1
                start = int(start_list[i]) - 1 
                loc = loc_list[i]
                # length = int(len_list[i])
                loc_range = [[x-start, y-start] for x, y in loc]
                loc_array = np.array([], dtype=int)
                for l in loc_range:
                    loc_array = np.append(loc_array, np.arange(l[0], l[1]))
                arr[i].put(loc_array, 1)
            # print('chrom done')
            kmeans = KMeans(n_clusters=2, random_state=0).fit(arr[:, :telo_len])
            labels = kmeans.labels_
            n1 = arr[labels == 0].shape[0]
            n2 = arr[labels == 1].shape[0]
            if min(n1, n2)/(n1 + n2) >= 0.25:
                n = min(n1, n2)
                arr_dict = {}
                for clust in range(kmeans.n_clusters):
                    closest_pt_dict = {}
                    cluster_arr = arr[labels == clust]
                    cluster_cen = kmeans.cluster_centers_[clust]
                    for p in range(len(cluster_arr)):
                        closest_pt_dict[p] = distance.euclidean(cluster_arr[p][:telo_len], cluster_cen)
                    ordered_idx = sorted(closest_pt_dict.items(), key=lambda item: item[1])[:n]
                    a = [i[0] for i in ordered_idx]
                    arr_dict[clust] = cluster_arr[a, :]
                arr1 = arr_dict[0]
                arr2 = arr_dict[1]
                arr3 = np.concatenate((arr1, arr2), axis=0)
                arr_sum = (np.sum(arr3, axis=0) / len(arr3)) * 100
                n = n*2
            else:
                arr_sum = (np.sum(arr, axis=0) / len(arr)) * 100
            tvs_arr[b][c] = arr_sum
            tvs_read_counts[b][c] = n
    now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
    print(now, ' - TVS analysis end')
    return tvs_arr, tvs_read_counts
