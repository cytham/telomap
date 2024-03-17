# Analyze telomere variant sequence (TVS) per chromosome end

import pandas as pd
import numpy as np
from datetime import datetime
from sklearn.cluster import KMeans
from scipy.spatial import distance

def tvs_analyzer(barcodes, df, telo_len=3000):
    tvs_arr = {}
    tvs_read_counts = {}
    for b in barcodes:
        tvs_arr[b] = {}
        tvs_read_counts[b] = {}
        chroms = df[(df['barcode'] == b) & (df['telo_len'].notnull())]['chrom'].unique().tolist()
        chroms = [x for x in chroms if pd.notnull(x)]  # Remove nan
        for c in chroms:
            start_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['s_junct'].tolist()
            loc_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['gap_loc'].tolist()
            len_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['telo_len_wgap'].tolist()
            longest = int(max(len_list))
            n = 0
            arr = np.array([])
            for i in range(len(start_list)):
                n += 1
                bin_data = []
                start = int(start_list[i])
                loc = loc_list[i]
                length = int(len_list[i])
                if loc:
                    for j in range(loc[0][0] - start):  # Assign values from start to first gap
                        bin_data.append(0)
                    for k in range(len(loc) - 1):
                        for h in range(loc[k][1] - loc[k][0]):  # Assign values for all gaps in-between
                            bin_data.append(1)
                        for g in range(loc[k + 1][0] - loc[k][1]):  # Assign values for all non-gaps in-between
                            bin_data.append(0)
                    for f in range(loc[-1][1] - loc[-1][0]):  # Assign values for last gap
                        bin_data.append(1)
                for d in range(length - len(bin_data)):  # Assign values from last gap to the read end
                    bin_data.append(0)
                for d in range(longest - len(bin_data)):  # Add filler to match longest read
                    bin_data.append(0)
                if arr.size > 0:
                    arr = np.vstack([arr, bin_data])
                else:
                    arr = np.array(bin_data)
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
    print(now, '- TVS analysis finished')
    return tvs_arr, tvs_read_counts

# def tvs_analyzer(barcode_reads, df):
#     tvs_arr = {}
#     tvs_read_counts = {}
#     for b in barcode_reads:
#         tvs_arr[b] = {}
#         tvs_read_counts[b] = {}
#         chroms = df[(df['barcode'] == b) & (df['telo_len'].notnull())]['chrom'].unique().tolist()
#         chroms = [x for x in chroms if pd.notnull(x)]  # Remove nan
#         for c in chroms:
#             start_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['s_junct'].tolist()
#             loc_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['gap_loc'].tolist()
#             len_list = df[(df['barcode'] == b) & (df['chrom'] == c)]['telo_len_wgap'].tolist()
#             longest = int(max(len_list))
#             n = 0
#             arr = np.array([])
#             for i in range(len(start_list)):
#                 n += 1
#                 bin_data = []
#                 start = int(start_list[i])
#                 loc = loc_list[i]
#                 length = int(len_list[i])
#                 if loc:
#                     for j in range(loc[0][0] - start):  # Assign values from start to first gap
#                         bin_data.append(0)
#                     for k in range(len(loc) - 1):
#                         for h in range(loc[k][1] - loc[k][0]):  # Assign values for all gaps in-between
#                             bin_data.append(1)
#                         for g in range(loc[k + 1][0] - loc[k][1]):  # Assign values for all non-gaps in-between
#                             bin_data.append(0)
#                     for f in range(loc[-1][1] - loc[-1][0]):  # Assign values for last gap
#                         bin_data.append(1)
#                 for d in range(length - len(bin_data)):  # Assign values from last gap to the read end
#                     bin_data.append(0)
#                 for d in range(longest - len(bin_data)):  # Add filler to match longest read
#                     bin_data.append(0)
#                 if arr.size > 0:
#                     arr = np.vstack([arr, bin_data])
#                 else:
#                     arr = np.array(bin_data)
#             arr_sum = (np.sum(arr, axis=0) / arr.shape[0]) * 100
#             tvs_arr[b][c] = arr_sum
#             tvs_read_counts[b][c] = n
#     now = datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")
#     print(now, ' - TVS analysis finished')
#     return tvs_arr, tvs_read_counts
