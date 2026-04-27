import numpy as np
import pandas as pd
import scipy
import torch
from sklearn import metrics
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import itertools
from tqdm import tqdm

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
    summed = summed[200:-200,200:-200] # remove padded region
    # if in inference use offset 0
    #summed = summed[:-200,:-200] # remove padded region

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


#added
def read_and_parse_coordinates(csv_file):
    df = pd.read_csv(csv_file)

    chr_range = df["coordinate"].str.split(":", expand=True)
    df["chr"] = chr_range[0]

    start_end = chr_range[1].str.split("-", expand=True)
    df["start"] = start_end[0].astype(int)
    df["end"] = start_end[1].astype(int)

    return df

def evaluate_predictions_by_chrom(
    chrom_list,
    hicdc,
    pval,
    ct2,
    mod,
    pred_base_path,
    percentile_cutoff=90,
    max_diag=200,
    plot=True
):
    """
    Loop through unique chromosomes and evaluate predictions
    """
    yhat_all = {}
    pearson_all = {}
    spearman_all = {}
    roc_all = {}
    prc_all = {}
    valid_diags_all = {}  # ← track valid diagonals per chrom
    
    chroms = sorted([str(c).replace("chr", "") for c in chrom_list])

    #chroms = sorted(
    #    hub_coordinates_df["chr"]
    #    .str.replace("chr", "")
    #   .unique()
    #)

    for chrom in chroms:
        print(f"\n===== Chromosome {chrom} =====")

        # Load HiCDC matrices
        hicdc_mat = hicdc[f"chr{chrom}"].toarray()
        hicdc_pval_mat = pval[f"chr{chrom}"].toarray()

        # Load predictions
        try:
            y_hat = get_preds(
                [chrom],
                f"{pred_base_path}/prediction_{ct2}_chr"
            )
        except (FileNotFoundError, KeyError, Exception) as e:
            print(f"chr{chrom} not found for {ct2}, skipping... ({e})")
            continue
        yhat_all[chrom] = y_hat

        # Ground truth processing
        true_mat = hicdc_mat[:-501, :-501].clip(-16, 16)

        # Predictions
        pred_mat = get_combined_yhat(
            y_hat[:, 0, :],
            start_ind=0,
            end_ind=y_hat.shape[0],
            avg_stripe=True
        )
        
        min_len = min(
            len(np.diag(true_mat, max_diag - 1)),
            len(np.diag(pred_mat, max_diag - 1))
        )

        # Significant interaction thresholding
        bin_true = np.concatenate(
            [np.diag(true_mat, i)[:min_len] for i in range(1, max_diag)]
        ).reshape(-1, min_len)

        bin_mask = (
            bin_true.sum(0)
            >= np.percentile(bin_true.sum(0), 1)
        )

        all_zval = np.concatenate(
            [np.diag(true_mat, i)[:min_len] for i in range(1, max_diag)]
        )
        zval_cutoff = np.percentile(all_zval, percentile_cutoff)

        pearson_all[chrom] = []
        spearman_all[chrom] = []
        roc_all[chrom] = []
        prc_all[chrom] = []
        valid_diags_all[chrom] = []  # ← init per chrom

        # Per-diagonal evaluation
        for i in range(1, max_diag):
            a = np.diag(true_mat, i)[:min_len][bin_mask]
            b = np.diag(pred_mat, i)[:min_len][bin_mask]
            p = a > zval_cutoff

            if np.all(b == b[0]):  # avoid constant predictions
                continue # skipped - not recorded in valid diags

            fpr, tpr, _ = metrics.roc_curve(p.astype(int), b)
            precision, recall, _ = metrics.precision_recall_curve(p.astype(int), b)

            pearson_all[chrom].append(scipy.stats.pearsonr(a, b)[0])
            spearman_all[chrom].append(scipy.stats.spearmanr(a, b)[0])
            roc_all[chrom].append(metrics.auc(fpr, tpr))
            prc_all[chrom].append(metrics.auc(recall, precision))
            valid_diags_all[chrom].append(i)  # ← only appended when not skipped

        # Plot per chromosome
        if plot:
            plt.figure(figsize=(6, 6))
            plt.plot(pearson_all[chrom], label=f"chr{chrom}")
            plt.ylabel("Pearson Correlation", fontsize=14)
            plt.xlabel("Genomic Distance (10kb)", fontsize=14)
            plt.ylim(0, 1)
            plt.legend()
            plt.show()

            plt.figure(figsize=(6, 6))
            plt.plot(spearman_all[chrom], label=f"chr{chrom}")
            plt.ylabel("Spearman Correlation", fontsize=14)
            plt.xlabel("Genomic Distance (10kb)", fontsize=14)
            plt.ylim(0, 1)
            plt.legend()
            plt.show()

            plt.figure(figsize=(6, 6))
            plt.plot(roc_all[chrom], label=f"chr{chrom}")
            plt.ylabel("AUROC", fontsize=14)
            plt.xlabel("Genomic Distance (10kb)", fontsize=14)
            plt.ylim(0, 1)
            plt.legend()
            plt.show()

            plt.figure(figsize=(6, 6))
            plt.plot(prc_all[chrom], label=f"chr{chrom}")
            plt.ylabel("AUPRC", fontsize=14)
            plt.xlabel("Genomic Distance (10kb)", fontsize=14)
            plt.ylim(0, 1)
            plt.legend()
            plt.show()

    return {
        "pearson": pearson_all,
        "spearman": spearman_all,
        "roc": roc_all,
        "prc": prc_all,
        "yhat": yhat_all,
        "hicdc_mat":hicdc_mat,
        "valid_diags": valid_diags_all,  # ← added
    }



def generate_plot_regions(
    hub_coordinates_df,
    chrom_sizes,
    flank=1_000_000,
    end_buffer=5_000_000
):
    """
    Generate non-overlapping plot regions from hub coordinates.

    Returns a DataFrame with:
    chr, ext_start, ext_end, hub_plotted_coordinates
    """

    plot_regions = []

    # keep track of already-used regions per chromosome
    used_regions = {}

    for _, row in hub_coordinates_df.iterrows():
        chrom = row["chr"]
        start = int(row["start"])
        end = int(row["end"])
        hub = row["hub"]

        chr_len = chrom_sizes[chrom]

        # 1) check if end is within last 5Mb
        if end > chr_len - end_buffer:
            print(f"Skipping {hub} ({chrom}:{start}-{end}) — too close to chromosome end")
            continue

        # 2) midpoint and extension
        mid = (start + end) // 2
        ext_start = max(0, mid - flank)
        ext_end = min(chr_len, mid + flank)

        # initialize list for chromosome
        if chrom not in used_regions:
            used_regions[chrom] = []

        # 3) check overlap with existing regions
        already_included = False
        for s, e in used_regions[chrom]:
            if not (ext_end <= s or ext_start >= e):
                already_included = True
                print(
                    f"Region {chrom}:{ext_start}-{ext_end} "
                    f"(hub {hub}) already included in "
                    f"{chrom}:{s}-{e}"
                )
                break

        if already_included:
            continue

        # record region
        used_regions[chrom].append((ext_start, ext_end))

        plot_regions.append({
            "chr": chrom,
            "ext_start": ext_start,
            "ext_end": ext_end,
            "hub_plotted_coordinates": f"{hub}:{chrom}:{start}-{end}"
        })

    return pd.DataFrame(plot_regions)


def plot_hic_region(
    mat,
    chrom,
    ext_start,
    ext_end,
    hub_plotted_coordinates,
    resolution=10_000,
    pad_bp=1_000_000,
    title_prefix="Predicted"
):
    """
    Plot predicted Hi-C region with zoomed hub coordinates.
    """

    # hub coordinates
    start_bp = ext_start
    end_bp = ext_end

    # convert to bins
    zoom_start_bin = start_bp // resolution
    zoom_end_bin   = end_bp   // resolution

    plot_start_bp  = start_bp - pad_bp
    plot_end_bp    = end_bp   + pad_bp

    plot_start_bin = plot_start_bp // resolution
    plot_end_bin   = plot_end_bp   // resolution

    # relative zoom coords
    zoom_start = zoom_start_bin - plot_start_bin
    zoom_end   = zoom_end_bin   - plot_start_bin

    # ticks
    n_ticks = 5
    ticks = np.linspace(zoom_start, zoom_end, n_ticks)
    tick_labels = [
        f"{(plot_start_bp + int(t)*resolution)/1e6:.2f} Mb"
        for t in ticks
    ]

    # extract region
    region = mat[
        plot_start_bin:plot_end_bin,
        plot_start_bin:plot_end_bin
    ]
    region = np.triu(region)

    # plot
    plt.rcParams['figure.figsize'] = (8, 3)
    img = pcolormesh_45deg(region, 4, -1)
    plt.colorbar(img, ax=plt.gca(), fraction=0.02, pad=0.04)  # add colorbar

    plt.xlim(zoom_start, zoom_end)
    plt.ylim(0, 65)
    plt.xticks(ticks, tick_labels)

    plt.title(
        f"{title_prefix}_{hub_plotted_coordinates} : "
        f"{chrom}:{start_bp:,}-{end_bp:,}"
    )

    plt.show()


def plot_scatac_coaccessibility(
    scatac,
    chrom,
    ext_start,
    ext_end,
    hub_plotted_coordinates,
    resolution_hic=10_000,
    resolution_scatac=500,   # 500bp per scATAC bin, adjust if different
    bin_factor=20,            # how many scATAC bins per Hi-C bin
    n_ticks=5,
    vmax=0.5,
    vmin=0,
    pad_bp=1_000_000
):
    """
    Plot scATAC co-accessibility (vertical-stripe correlation) for a genomic hub.
    """
    # hub coordinates
    start_bp = ext_start
    end_bp = ext_end

    # convert to bins
    zoom_start_bin = start_bp // resolution_hic
    zoom_end_bin   = end_bp   // resolution_hic

    plot_start_bp  = start_bp - pad_bp
    plot_end_bp    = end_bp   + pad_bp

    plot_start_bin = plot_start_bp // resolution_hic
    plot_end_bin   = plot_end_bp   // resolution_hic

    # relative zoom coords
    zoom_start = zoom_start_bin - plot_start_bin
    zoom_end   = zoom_end_bin   - plot_start_bin

    # subset the matrix and compute co-accessibility
    tmp = cpu_batch_corcoeff_vstripe(
        torch.tensor(scatac['chr{}'.format(chrom)][:, plot_start_bin*bin_factor : plot_end_bin*bin_factor].toarray())
    )

    # aggregate back to Hi-C resolution
    tmp = tmp.reshape(tmp.shape[0]//bin_factor, bin_factor, -1).mean(axis=1)
    tmp = tmp.reshape(tmp.shape[0], tmp.shape[1]//bin_factor, bin_factor).mean(axis=2)

    # --- define tick labels ---
    ticks = np.linspace(zoom_start, zoom_end, n_ticks)
    tick_labels = [
        f"{(plot_start_bp + int(t)*resolution_hic)/1e6:.2f} Mb"
        for t in ticks
    ]

    # --- plot ---
    plt.rcParams['figure.figsize'] = (8, 3)
    img = pcolormesh_45deg(tmp, vmax=vmax, vmin=vmin)
    plt.colorbar(img, ax=plt.gca(), fraction=0.02, pad=0.04)  # add colorbar
    plt.xlim(zoom_start, zoom_end)
    plt.ylim(0, 65)
    plt.xticks(ticks, tick_labels)
    plt.title(f"scATAC_coaccessibility_{hub_plotted_coordinates}", fontsize=10)
    plt.xlabel(chrom)
    plt.tight_layout()
    plt.show()

def plot_atac_ctcf(
    ctcf_motif,
    atac,
    chrom,
    ext_start,
    ext_end,
    hub_plotted_coordinates,
    resolution_hic=10_000,
    bin_factor=200,  # how many scATAC bins per Hi-C bin
    n_ticks=5
):
    """
    Plot CTCF motif and ATAC signal exactly over ext_start to ext_end (no extra padding).
    """

    chrom = str(chrom).replace("chr","")  # consistent dictionary key

    # Convert genomic coordinates to Hi-C bins
    start_bin_hic = ext_start // resolution_hic
    end_bin_hic   = ext_end   // resolution_hic

    # Convert Hi-C bins to scATAC bins
    start_bin_atac = start_bin_hic * bin_factor
    end_bin_atac   = end_bin_hic   * bin_factor

    # Define ticks in scATAC bins
    num_bins = end_bin_atac - start_bin_atac
    ticks = np.linspace(0, num_bins, n_ticks)
    tick_labels = [
        f"{(ext_start + (int(t)/bin_factor)*resolution_hic)/1e6:.2f} Mb"
        for t in ticks
    ]

    # --- plot CTCF ---
    plt.figure(figsize=(8, 1.5))
    ctcf_signal = ctcf_motif[f'chr{chrom}'].toarray()[0][start_bin_atac:end_bin_atac]
    plt.plot(ctcf_signal)
    plt.xticks(ticks, tick_labels)
    plt.ylim(2.8, 3.5)
    plt.title(f"CTCF_motif_score_{hub_plotted_coordinates}", fontsize=10)
    plt.tight_layout()
    plt.show()

    # --- plot ATAC ---
    fig, ax = plt.subplots(figsize=(8,1.5))
    atac_signal = atac[f'chr{chrom}'][start_bin_atac:end_bin_atac]
    ax.plot(atac_signal)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_title(f"ATACseq_signal_{hub_plotted_coordinates}", fontsize=10) 
    plt.tight_layout()
    ax.set_ylim(0,500)
    plt.show()

def evaluate_predictions_by_region(
    eval_output,
    hicdc,
    chrom,
    start,
    end,
    bin_size=10000,
    percentile_cutoff=90,
    max_diag=200,
    plot=True
):
    chrom = str(chrom).replace("chr", "")
    chrom_key = f"chr{chrom}"

    start_bin = start // bin_size
    end_bin   = end   // bin_size
    print(f"Region: {chrom_key}:{start}-{end} → bins [{start_bin}, {end_bin}]")

    # Ground truth — same processing as original, no slicing
    hicdc_mat = hicdc[chrom_key].toarray()
    true_mat  = hicdc_mat[:-501, :-501].clip(-16, 16)

    # Predictions — same as original, no slicing
    y_hat    = eval_output["yhat"][chrom]
    pred_mat = get_combined_yhat(
        y_hat[:, 0, :],
        start_ind=0,
        end_ind=y_hat.shape[0],
        avg_stripe=True
    )

    min_len = min(
        len(np.diag(true_mat, max_diag - 1)),
        len(np.diag(pred_mat, max_diag - 1))
    )

    # For each diagonal offset i, the k-th element corresponds to bins (k, k+i).
    # We only keep entries where both bins fall within [start_bin, end_bin].
    bin_true = np.concatenate(
        [np.diag(true_mat, i)[:min_len] for i in range(1, max_diag)]
    ).reshape(-1, min_len)
    bin_mask = bin_true.sum(0) >= np.percentile(bin_true.sum(0), 1)
    all_zval = np.concatenate(
        [np.diag(true_mat, i)[:min_len] for i in range(1, max_diag)]
    )
    zval_cutoff = np.percentile(all_zval, percentile_cutoff)

    pearson_r, spearman_r, roc_r, prc_r = [], [], [], []

    for i in range(1, max_diag):
        a_full = np.diag(true_mat, i)[:min_len]
        b_full = np.diag(pred_mat, i)[:min_len]

        # Region mask: diagonal offset i → element k spans bins (k, k+i)
        # Keep k where start_bin <= k and k+i <= end_bin; Both bins in region
        #region_mask = (np.arange(min_len) >= start_bin) & (np.arange(min_len) + i <= end_bin)
        # Anchor bins in the region
        region_mask = (np.arange(min_len) >= start_bin) & (np.arange(min_len) <= end_bin)

        combined_mask = bin_mask & region_mask
        a = a_full[combined_mask]
        b = b_full[combined_mask]

        if len(a) < 2 or np.all(b == b[0]):
            continue
        p = a > zval_cutoff
        if len(np.unique(p)) < 2:
            continue

        fpr, tpr, _          = metrics.roc_curve(p.astype(int), b)
        precision, recall, _ = metrics.precision_recall_curve(p.astype(int), b)
        pearson_r.append(scipy.stats.pearsonr(a, b)[0])
        spearman_r.append(scipy.stats.spearmanr(a, b)[0])
        roc_r.append(metrics.auc(fpr, tpr))
        prc_r.append(metrics.auc(recall, precision))

    print(f"Valid diagonals: {len(pearson_r)}/{max_diag - 1}")
    print(f"Mean Pearson:  {np.nanmean(pearson_r):.4f}")
    print(f"Mean Spearman: {np.nanmean(spearman_r):.4f}")
    print(f"Mean AUROC:    {np.nanmean(roc_r):.4f}")
    print(f"Mean AUPRC:    {np.nanmean(prc_r):.4f}")

    if plot:
        fig, axes = plt.subplots(1, 4, figsize=(20, 4))
        for ax, (data, label) in zip(axes, [
            (pearson_r,  "Pearson Correlation"),
            (spearman_r, "Spearman Correlation"),
            (roc_r,      "AUROC"),
            (prc_r,      "AUPRC"),
        ]):
            ax.plot(data, label=f"{chrom_key}:{start}-{end}")
            ax.set_ylabel(label, fontsize=12)
            ax.set_xlabel("Genomic Distance (10kb)", fontsize=12)
            ax.set_ylim(0, 1)
            ax.legend(fontsize=8)
        plt.suptitle(f"Region: {chrom_key}:{start:,}-{end:,}", fontsize=13)
        plt.tight_layout()
        plt.show()

    return {
        "pearson":  pearson_r,
        "spearman": spearman_r,
        "roc":      roc_r,
        "prc":      prc_r,
    }

def aggregate_by_distance_range(eval_result, distance_ranges, chrom_list):
    """
    distance_ranges: list of (low_bin, high_bin) tuples — already in diagonal/bin units
    e.g. (0,17) = 0–170kb at 10kb resolution
    """
    metrics_keys = ["pearson", "spearman", "roc", "prc"]
    chroms = sorted([str(c).replace("chr", "") for c in chrom_list])

    results = {m: {dr: {} for dr in distance_ranges} for m in metrics_keys}

    for chrom in chroms:
        if chrom not in eval_result["valid_diags"]:
            continue

        valid_diags = eval_result["valid_diags"][chrom]

        for m in metrics_keys:
            vals = eval_result[m][chrom]

            for dr in distance_ranges:
                low_bin, high_bin = dr
                range_vals = [
                    v for diag, v in zip(valid_diags, vals)
                    if low_bin <= diag < high_bin
                ]
                results[m][dr][chrom] = np.nanmean(range_vals) if range_vals else np.nan

    return results