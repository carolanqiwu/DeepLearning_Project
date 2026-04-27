import os
import sys

import numpy as np
import pandas as pd
import scipy
import torch
from sklearn import metrics

import pickle
import seaborn as sns

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import warnings
import itertools
from tqdm import tqdm

OUTDIR = "/data1/lesliec/carolw/projects/chromafold/vis/BRCA_C9C8D426_8k"
os.makedirs(OUTDIR, exist_ok=True)

def get_starts(chrom_list, step = 5e4):
    startl = [] #hold all the start positions across all chromosomes
    chroml = [] #hold the chromosome label corresponding to each start
    for chrom in chrom_list:
        chrom = 'chr{}'.format(chrom)
        cur_starts = list(np.arange(start_dict[chrom],end_dict[chrom]-5000000, step).astype(int))
        startl = startl + cur_starts
        chroml = chroml + list(np.repeat(chrom, len(cur_starts)))
    return startl, chroml

def get_visualization(y_hat_list,y_true_list,avg_stripes = False):
    y_hat_list_reshaped = np.concatenate([x.reshape(1,-1) for x in y_hat_list[start_ind_origami:start_ind_origami+plot_len]], axis = 0)
    y_true_list_reshaped = np.concatenate([x.reshape(1,-1) for x in y_true_list[start_ind_origami:start_ind_origami+plot_len]], axis = 0)
    mat = []
    mat_ytrue = []
    for i in range(plot_len):
        mat.append(np.insert(np.zeros(plot_len), i, y_hat_list_reshaped[i]))
        mat_ytrue.append(np.insert(np.zeros(plot_len), i, y_true_list_reshaped[i].clip(-16,16)))
    if avg_stripes:
        summed = (pd.DataFrame(np.array(mat)).reindex(np.arange(-1*pred_len,plot_len+1,1)).fillna(0).iloc[0:plot_len+pred_len+1,0:plot_len+pred_len+1].values+pd.DataFrame(np.array(mat)).reindex(
        np.arange(-1*pred_len,plot_len+1,1)).fillna(0).T.iloc[0:plot_len+pred_len+1,0:plot_len+pred_len+1].values)/2
    else:
        summed = pd.DataFrame(np.array(mat)).reindex(np.arange(-1*pred_len,plot_len+1,1)).fillna(0).T.iloc[0:plot_len+pred_len+1,0:plot_len+pred_len+1].values
                  
    comb = np.zeros((summed.shape[0],summed.shape[0]))
    # keep upper trangular part for vis
    comb[np.triu_indices(comb.shape[0])[0],np.triu_indices(comb.shape[0])[1]] = (summed[np.triu_indices(comb.shape[0])[0],
                                                    np.triu_indices(comb.shape[0])[1]])
    # bottom triangular part is 0
    comb[np.tril_indices(comb.shape[0],k = -1)[0],np.tril_indices(comb.shape[0],k = -1)[1]] = 0

    return comb

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols
    
def get_combined_yhat(y_hat_list, start_ind, end_ind, offset = 200, avg_stripe = False): 
    pred_len = 200
    y_hat_list_reshaped = np.concatenate([x.reshape(1,-1) for x in y_hat_list[start_ind:end_ind]], axis = 0)
    chrom_length = y_hat_list_reshaped.shape[0]
    mat = []

    for i in tqdm(range(chrom_length)):
        mat.append(np.insert(np.zeros(chrom_length+offset+1), i, np.insert(y_hat_list_reshaped[i],pred_len,0)))
    summed = pd.DataFrame(
    np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
        ).fillna(0).T.iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values
    
    if avg_stripe:
        summed = (pd.DataFrame(np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
                ).fillna(0).iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values+pd.DataFrame(
            np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
                ).fillna(0).T.iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values)/2
    # if in inference use offset -2MB    
    #summed = summed[200:-200,200:-200] # remove padded region
    # if in inference use offset 0
    summed = summed[:-200,:-200] # remove padded region

    return summed

def get_metacell_profile(tile_dict, nbrs):
    metacell_tile_dict = {}
    metacell = nbrs
    for chrom in list(tile_dict.keys()):
        metacell_tile_dict[chrom] = (scipy.sparse.csr_matrix(metacell) * tile_dict[chrom])
    return metacell_tile_dict 

def cpu_jaccard_vstripe(x):
    size = x.shape[1]
    eps=1e-8
    i = 2

    x = torch.where(x>0.0, torch.tensor([1.0]), torch.tensor([0.0]))
    num = torch.mm(x, x.transpose(0,1))
    
    x = torch.where(x==0.0, torch.tensor([1.0]), torch.tensor([0.0]))
    denom = torch.mm(x, x.transpose(0,1))
    denom = size - denom

    num = torch.div(num, torch.max(denom, eps * torch.ones_like(denom)))
    
    return num

def cpu_batch_corcoeff_vstripe(x):
    c = cpu_jaccard_vstripe(x.permute(1,0))
    c[c != c] = 0
    return c

def get_preds(chroms, path):

    y_z_hat_list = []
    tmp = []
    for chrom in chroms:
        tmp.append(np.load(path + '{}.npz'.format(chrom))['arr_0'])
    tmp = np.concatenate([y for y in tmp], axis = 0)
    y_z_hat_list.append(tmp)
    y_z_hat_list = np.concatenate([np.expand_dims(y,1) for y in y_z_hat_list], axis = 1)

    return y_z_hat_list

def pcolormesh_45deg(C,vmax,vmin):
    n = C.shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    # plot
    img = plt.pcolormesh(A[:,1].reshape(n+1,n+1),A[:,0].reshape(n+1,n+1),np.flipud(C),
                  cmap = 'RdBu_r',
                         vmax = vmax, vmin = vmin
                        )
    return img

ct = 'BRCA_C9C8D426'
ct2 = 'scATAC_BRCA_C9C8D426_A3FD_4455_89A9_768BC01D66A9_X009_S02_B1_T1.fragments'
genome = 'hg38'
mod = 'chromafold_CTCFmotif'

# Load Hi-C data for evaluation
hicdc = pickle.load(open("/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA_C9C8D426_8k/hichip/zvalue/{}_hichip_zscore_dict.p".format(ct), "rb"))
pval = pickle.load(open("/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA_C9C8D426_8k/hichip/qvalue/{}_hichip_qvalue_dict.p".format(ct), "rb"))

# Load input data for visualization
ctcf_motif = pickle.load(open("/data1/lesliec/carolw/repos/ChromaFold/{}_ctcf_motif_score.p".format(genome), 'rb'))
atac = pickle.load(open("/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA_C9C8D426_8k/atac/{}_tile_pbulk_50bp_dict.p".format(ct2), 'rb'))
scatac = pickle.load(open("/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA_C9C8D426_8k/atac/{}_tile_500bp_dict.p".format(ct2), 'rb'))
metacell_path = pd.read_csv('/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA_C9C8D426_8k/atac/{}_metacell_mask.csv'.format(ct2), index_col= 0).values
scatac = get_metacell_profile(scatac, metacell_path)

chrom = list(scatac.keys())[0]
print("Tile matrix shape:", scatac[chrom].shape)

# Evaluation
chrom = '6'
hicdc_mat = hicdc['chr{}'.format(chrom)].toarray()
hicdc_pval_mat = pval['chr{}'.format(chrom)].toarray()

# Get predictions
y_hat = get_preds([chrom], '/data1/lesliec/carolw/projects/chromafold/predictions/BRCA_C9C8D426_8k/prediction_{}_chr'.format(ct2))

# Evaluation starts here
pearson_list = {}
spearman_list = {}
roc_list = {}
prc_list = {}

# Get ground-truth Hi-C data
# Removes last 501 bins from both dimensions, clip extreme zscore
true_mat = hicdc_mat[:-501,:-501].clip(-16,16)
print("true_mat:",true_mat.shape)
min_len = len(np.diag(true_mat,199)) # Compute shortest diagonal length
print("min_len:",min_len)

# Significant interaction analysis
percentile_cutoff = 90 # Zscore percentile cutoff for significant interactions
bin_true = np.concatenate([np.diag(true_mat,i)[0:min_len] for i in range(1,200)]).reshape(-1, min_len)
bin_mask = (bin_true.sum(0) >= np.percentile(bin_true.sum(0),1)) #mask for low-mappability regions
print("bin_mask length:", len(bin_mask))
all_zval = np.concatenate([np.diag(true_mat,i)[0:min_len] for i in range(1,200)])
zval_cutoff = np.percentile(all_zval, percentile_cutoff)
bin_zval = all_zval > zval_cutoff
print('Positive class proportion: {}'.format(bin_zval.sum()/bin_zval.shape[0]))
print('Zval cutoff: {}'.format(zval_cutoff))

# Get processed predictions
pred_mat = get_combined_yhat(y_hat[:,0,:], start_ind = 0, end_ind = y_hat.shape[0], avg_stripe=True)
print("pred_mat:", pred_mat.shape)
print("diag lengths:")
for i in range(1,5):
    print(i, len(np.diag(true_mat,i)), len(np.diag(pred_mat,i)))
    
pearson_list[mod] = []
spearman_list[mod] = []
roc_list[mod] = []
prc_list[mod] = []

# fixed to aovid min_len exceed pred_mat index
for i in range(1,200):
    true_diag = np.diag(true_mat, i)
    pred_diag = np.diag(pred_mat, i)
    diag_len = min(len(true_diag), len(pred_diag), len(bin_mask))
    a = true_diag[:diag_len][bin_mask[:diag_len]]
    b = pred_diag[:diag_len][bin_mask[:diag_len]]
    p = a > zval_cutoff
    
    #a = np.diag(true_mat,i)[0:min_len][bin_mask]
    #b = np.diag(pred_mat,i)[0:min_len][bin_mask]
    #p = np.diag(true_mat,i)[0:min_len][bin_mask] > zval_cutoff
    
    fpr, tpr, thresholds = metrics.roc_curve(p.astype(int), b)
    precision, recall, thresholds = metrics.precision_recall_curve(p.astype(int), b)
    pearson_list[mod].append(scipy.stats.pearsonr(a,b)[0])    
    spearman_list[mod].append(scipy.stats.spearmanr(a,b)[0])    
    roc_list[mod].append(metrics.auc(fpr, tpr))
    prc_list[mod].append(metrics.auc(recall, precision))

# eval plots
plt.rcParams['figure.figsize'] = 6,6
plt.plot(pearson_list[mod], label = mod)
# plt.legend()
plt.ylabel('Pearson Correlation', fontsize = 18)
plt.xlabel('Genomic Distance (10kb)', fontsize = 18)
plt.ylim(0,1)
plt.savefig(os.path.join(OUTDIR, "pearson_vs_distance.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.rcParams['figure.figsize'] = 6,6
plt.plot(spearman_list[mod], label = mod)
# plt.legend()
plt.ylabel('Spearman Correlation', fontsize = 18)
plt.xlabel('Genomic Distance (10kb)', fontsize = 18)
plt.ylim(0,1)
plt.savefig(os.path.join(OUTDIR, "spearman_vs_distance.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.rcParams['figure.figsize'] = 6,6
plt.plot(roc_list[mod], label = mod)
plt.legend()
plt.ylabel('AUROC', fontsize = 18)
plt.xlabel('Genomic Distance (10kb)', fontsize = 18)
plt.ylim(0,1)
plt.savefig(os.path.join(OUTDIR, "auroc_vs_distance.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.rcParams['figure.figsize'] = 6,6
plt.plot(prc_list[mod], label = mod)
plt.legend()
plt.ylabel('AUPRC', fontsize = 18)
plt.xlabel('Genomic Distance (10kb)', fontsize = 18)
plt.ylim(0,1)
plt.savefig(os.path.join(OUTDIR, "auprc_vs_distance.png"), dpi=300, bbox_inches="tight")
plt.close()

# Visualization starts here
pred_mat = get_combined_yhat(y_hat[:,0,:], start_ind = 0, end_ind = y_hat.shape[0], avg_stripe=True)

# predicted plot
# Added code
resolution = 10_000
start_bp = 43_565_123
end_bp   = 44_259_021

zoom_start_bin = start_bp // resolution
zoom_end_bin   = end_bp // resolution

#define larger plotting region
pad_bp = 1_000_000
plot_start_bp = start_bp - pad_bp
plot_end_bp   = end_bp + pad_bp

plot_start_bin = plot_start_bp // resolution
plot_end_bin   = plot_end_bp // resolution

#compute zoom coordinates relative to 
zoom_start = zoom_start_bin - plot_start_bin
zoom_end   = zoom_end_bin   - plot_start_bin

# Number of ticks you want
n_ticks = 5

# Tick positions in matrix coordinates
ticks = np.linspace(zoom_start, zoom_end, n_ticks)

# Convert bins → genomic coordinates
tick_labels = [
    f"{(plot_start_bp + int(t)*resolution)/1e6:.2f} Mb"
    for t in ticks
]

#extract larger hic region
region = pred_mat[plot_start_bin:plot_end_bin,
                  plot_start_bin:plot_end_bin]

region = np.triu(region)


plt.rcParams['figure.figsize'] = (8,3)
img = pcolormesh_45deg(region, 4, -1)
# zoom in
plt.xlim(zoom_start, zoom_end)
plt.ylim(0, zoom_end - zoom_start)
plt.xticks(ticks, tick_labels)
plt.title("Predicted VEGFA genomic region: chr6:43,565,123-44,259,021")
plt.savefig(os.path.join(OUTDIR, "predicted_hic_VEGFA.png"), dpi=300, bbox_inches="tight")
plt.close()

# Ground truth plot
#extract larger hic region
region = hicdc_mat[plot_start_bin:plot_end_bin,
                  plot_start_bin:plot_end_bin]

region = np.triu(region)


plt.rcParams['figure.figsize'] = (8,3)
img = pcolormesh_45deg(region, 4, -1)
# zoom in
plt.xlim(zoom_start, zoom_end)
plt.ylim(0, zoom_end - zoom_start)
plt.xticks(ticks, tick_labels)
plt.title("Ground Truth VEGFA genomic region: chr6:43,565,123-44,259,021")
plt.savefig(os.path.join(OUTDIR, "groundtruth_hic_VEGFA.png"), dpi=300, bbox_inches="tight")
plt.close()

# Co-accessbility plot
#VEGFA genomic region: chr6:43,565,123-44,259,021
#resolution = 10_000
#start_bp = 42_565_123 
#end_bp   = 45_259_021 
#ext_start_bp = start_bp - 1_000_000 # extend 1MB
#ext_end_bp   = start_bp + 1_000_000 
#bin in HiC resolution
#start_bin = ext_start_bp // resolution #4257
#end_bin   = ext_end_bp   = start_bp + 1_000_000  // resolution #4523

start_bin = 4256
end_bin = 4525

#convert to scATAC resolution
tmp = cpu_batch_corcoeff_vstripe(torch.tensor(scatac['chr{}'.format(chrom)][:, (start_bin)*20 : (end_bin)*20].toarray()))
#convert back to HIC resolution
tmp = tmp.reshape(tmp.shape[0]//20, 20, -1).mean(axis=1).reshape(-1, tmp.shape[1]//20, 20).mean(axis=2)

n_ticks = 5
#tick_start = 1_000_000 // resolution #100
#tick_end = (end_bin - start_bin) - 1_000_000 // resolution #166
ticks = np.linspace(100, 169, n_ticks)  # bin indices along matrix
tick_labels = [f"{(42_565_123 + int(t)*10_000)/1e6:.2f} Mb" for t in ticks]

plt.rcParams['figure.figsize'] = 8, 3
img = pcolormesh_45deg(tmp, vmax = 0.5, vmin = 0)
plt.xlim(100,169)
plt.ylim(0,zoom_end - zoom_start)
plt.xticks(ticks, tick_labels)
plt.title("scATAC vertical-stripe correlation: VEGFA locus")
plt.savefig(os.path.join(OUTDIR, "scATAC_vstripe_VEGFA.png"), dpi=300, bbox_inches="tight")
plt.close()

# CTCF plot
resolution =10_000

start_bin = 43_565_123 // resolution
end_bin   = 44_259_021 //resolution

n_ticks = 6
num_bin = (end_bin-start_bin)*200
ticks = np.linspace(0, num_bin, n_ticks)  # bin indices along matrix
tick_labels = [
    f"{(43_565_123 + (t / 200) * resolution) / 1e6:.2f} Mb"
    for t in ticks
]

plt.rcParams['figure.figsize'] = 8, 1.5
plt.plot(ctcf_motif['chr{}'.format(chrom)].toarray()[0][(start_bin)*200 : (end_bin)*200])
plt.xticks(ticks, tick_labels)
plt.ylim(2.8,3.5)
plt.ylabel('CTCF motif score')
plt.savefig(os.path.join(OUTDIR, "CTCF_VEGFA.png"), dpi=300, bbox_inches="tight")
plt.close()

#ATAC plot
plt.plot(atac['chr{}'.format(chrom)][(start_bin)*200 : (end_bin)*200])
plt.xticks(ticks, tick_labels)
#plt.xlim(0, 50000)
# plt.ylim(2.8,3.5)
plt.ylabel('ATAC-seq signal')
plt.savefig(os.path.join(OUTDIR, "ATAC_VEGFA.png"), dpi=300, bbox_inches="tight")
plt.close()
