import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import mannwhitneyu

def get_sample_cnv_segments(ascat_file, case_id):

    df = pd.read_csv(ascat_file, sep="\t")

    sample_df = df[df["case_submitter_id"] == case_id][
        ["Chromosome", "Start", "End", "Relative_Copy_Number"]
    ].copy()

    return sample_df.reset_index(drop=True)

def annotate_bins_with_cnv(cnv_df, chrom_sizes_df, bin_size=10000):
    """
    Create genomic bins from a chrom_sizes DataFrame and annotate with CNV (Relative_Copy_Number).

    Parameters
    ----------
    cnv_df : pandas.DataFrame
        CNV calls with columns: ['Chromosome','Start','End','Relative_Copy_Number']

    chrom_sizes_df : pandas.DataFrame
        Chromosome sizes with columns: ['chr', 'size']

    bin_size : int
        Bin size (default = 10 kb)

    Returns
    -------
    pandas.DataFrame
        Columns: ['chrom','start','end','copy_number','CNA_status']
    """

    bins = []

    # create genomic bins from chrom_sizes DataFrame
    for _, row in chrom_sizes_df.iterrows():
        chrom = str(row["chr"])
        length = int(row["size"])
        # skip non-standard chromosomes if needed
        if chrom in ["ALL", "M"]:
            continue
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            bins.append((chrom, start, end))

    bins_df = pd.DataFrame(bins, columns=["chrom", "start", "end"])
    # default copy number = 1 (neutral)
    bins_df["copy_number"] = 1.0
    bins_df["chrom"] = "chr" + bins_df["chrom"].astype(str)

    # annotate bins
    for _, row in cnv_df.iterrows():
        chrom = row["Chromosome"]
        seg_start = int(row["Start"])
        seg_end = int(row["End"])
        cn = float(row["Relative_Copy_Number"])

        mask = (bins_df["chrom"] == chrom) & (bins_df["start"] < seg_end) & (bins_df["end"] > seg_start)
        overlapping_bins = bins_df[mask].copy()

        for idx, bin_row in overlapping_bins.iterrows():
            # compute fraction of bin covered by CNV
            overlap_start = max(bin_row["start"], seg_start)
            overlap_end = min(bin_row["end"], seg_end)
            frac_overlap = (overlap_end - overlap_start) / bin_size

            if frac_overlap >= 1.0:  # fully covered
                bins_df.at[idx, "copy_number"] = cn
            else:  # partially covered
                bins_df.at[idx, "copy_number"] = cn
                bins_df.at[idx, "CNA_status"] = "mixed"

    # classify CNA status
    def classify_cn(cn):
        if cn >= 2:
            return "amplification"
        elif cn < 0.75:
            return "deletion"
        else:
            return "neutral"

    bins_df["CNA_status"] = bins_df["copy_number"].apply(classify_cn)

    return bins_df

def create_bin_diag_accuracy_all(results, bins_df, metrics=('roc','pearson','prc')):
    """
    Create a DataFrame with one row per genomic bin storing diagonal accuracy lists
    for multiple metrics.

    Parameters
    ----------
    results : dict
        Output from evaluate_predictions_by_chrom
        Must contain results[metric] per chromosome

    bins_df : pd.DataFrame
        Annotated bins DataFrame with columns ['chrom', 'start', 'end']

    metrics : tuple
        Metrics to extract (default: roc, pearson, prc)

    Returns
    -------
    bin_acc_df : pd.DataFrame
        Columns:
        ['chrom','start','end','roc_list','pearson_list','prc_list']
    """

    bin_acc_records = []

    chroms = sorted(['chr' + str(c) for c in results[metrics[0]].keys()])

    for chrom in chroms:

        chrom_bins = (
            bins_df[bins_df['chrom'] == chrom]
            .sort_values('start')
            .reset_index(drop=True)
        )

        n_bins_total = chrom_bins.shape[0]
        n_bins = n_bins_total - 501
        chrom_bins_valid = chrom_bins.iloc[:n_bins].reset_index(drop=True)

        # store diagonal arrays for each metric
        chrom_diag_dict = {
            metric: results[metric][chrom.replace('chr','')]
            for metric in metrics
        }

        max_diag = len(next(iter(chrom_diag_dict.values())))

        for i in range(n_bins):

            record = {
                'chrom': chrom,
                'start': chrom_bins_valid.loc[i,'start'],
                'end': chrom_bins_valid.loc[i,'end']
            }

            for metric in metrics:

                diag_list = chrom_diag_dict[metric]

                diag_values = [
                    diag_list[k] if i + k < n_bins else np.nan
                    for k in range(max_diag)
                ]

                record[f"{metric}_list"] = diag_values

            bin_acc_records.append(record)

    bin_acc_df = pd.DataFrame(bin_acc_records)

    return bin_acc_df
    
def create_bin_diag_accuracy(results, bins_df, metric='roc'):
    """
    Create a DataFrame with one row per genomic bin, storing 
    the list of diagonal accuracies for each bin.

    Parameters
    ----------
    results : dict
        Output from evaluate_predictions_by_chrom
        Must contain results[metric] per chromosome

    bins_df : pd.DataFrame
        Annotated bins DataFrame with columns ['chrom', 'start', 'end', ...]
        Must contain bins in the same chromosomes as results

    metric : str
        Which metric to use: 'roc', 'prc', 'pearson', 'spearman'

    Returns
    -------
    bin_acc_df : pd.DataFrame
        Columns:
        ['chrom', 'start', 'end', 'diag_list'] where diag_list is a list of
        accuracy per diagonal (length = max_diag)
    """
    bin_acc_records = []
    
   # Loop over chromosomes, prepend 'chr' to match bins_df
    chroms = sorted(['chr' + str(c) for c in results[metric].keys()])
    for chrom in chroms:
        # Clip last bins to match evaluate_predictions_by_chrom
        chrom_bins = bins_df[bins_df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        n_bins_total = chrom_bins.shape[0]
        n_bins = n_bins_total - 501

        chrom_bins_valid = chrom_bins.iloc[:n_bins].reset_index(drop=True)

        # Get diagonal accuracies for this chromosome
        # remove 'chr' to match results keys
        chrom_diag_list = results[metric][chrom.replace('chr', '')]  # list of length max_diag
        max_diag = len(chrom_diag_list)

        # Assign per-bin diagonal accuracy
        for i in range(n_bins):
            diag_values = []
            for k in range(max_diag):
                if i + k < n_bins:
                    diag_values.append(chrom_diag_list[k])
                else:
                    diag_values.append(np.nan)
            bin_acc_records.append({
                'chrom': chrom,
                'start': chrom_bins_valid.loc[i, 'start'],
                'end': chrom_bins_valid.loc[i, 'end'],
                f'{metric}_list': diag_values
            })

    bin_acc_df = pd.DataFrame(bin_acc_records)
    return bin_acc_df

def filter_active_bins(bin_acc_df, hicdc, percentile=90, clip_last=501):
    """
    Filter genomic bins for "active" bins based on ground truth HiChIP matrices for all chromosomes.

    Parameters
    ----------
    bin_acc_df : pd.DataFrame
        DataFrame of genomic bins with columns ['chrom', 'start', 'end', ...].
        Must include all chromosomes present in hicdc.

    hicdc : dict
        Dictionary of sparse HiChIP matrices per chromosome, keys like 'chr1', 'chr2', etc.

    percentile : float
        Percentile cutoff to define active bins (default=90 for top 10%).

    clip_last : int
        Number of bins to remove at the end of each chromosome (default=501)

    Returns
    -------
    pd.DataFrame
        Subset of bin_acc_df containing only active bins across all chromosomes.
    """
    active_bins_list = []

    for chrom in sorted(hicdc.keys(), key=lambda x: int(x.replace('chr','')) if x[3:].isdigit() else 23):
        hicdc_mat = hicdc[chrom].toarray()

        if hicdc_mat.shape[0] <= clip_last:
            print(f"{chrom} has <= {clip_last} bins, skipping")
            continue

        # Clip last bins
        true_mat = hicdc_mat[:-clip_last, :-clip_last].clip(-16, 16)

        # Subset bin_acc_df for this chromosome
        chrom_bins = bin_acc_df[bin_acc_df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        n_bins = min(chrom_bins.shape[0], true_mat.shape[0])
        chrom_bins = chrom_bins.iloc[:n_bins]

        # Sum interactions per bin
        bin_sums = true_mat.sum(axis=1)

        # Determine cutoff for active bins
        threshold = np.percentile(bin_sums, percentile)

        # Select active bins
        active_idx = np.where(bin_sums >= threshold)[0]
        active_bins_list.append(chrom_bins.iloc[active_idx])

    # Concatenate all chromosomes
    active_bins_df = pd.concat(active_bins_list).reset_index(drop=True)

    return active_bins_df

def plot_metric_ranges(all_samples_df, metric, distance_ranges):
    """
    Plot boxplots of a metric across all distance ranges.

    Parameters
    ----------
    all_samples_df : pd.DataFrame
        Combined dataframe containing all samples.
        Must include columns:
        ['sample', 'CNA_status', f'{metric}_{range_name}']

    metric : str
        Metric prefix ('roc', 'prc', 'pearson')

    distance_ranges : dict
        Dictionary of distance ranges (keys used for column suffix).
    """

    for range_name in distance_ranges.keys():

        col = f"{metric}_{range_name}"

        if col not in all_samples_df.columns:
            print(f"{col} not found, skipping")
            continue

        plt.figure(figsize=(8,6))

        sns.boxplot(
            data=all_samples_df,
            x="sample",
            y=col,
            hue="CNA_status"
        )

        plt.title(f"{metric.upper()} for {range_name} range")
        plt.xlabel("Sample")
        plt.ylabel(metric.upper())

        plt.legend(title="Region")
        plt.tight_layout()
        plt.show()


def cliffs_delta(x, y):
    """Compute Cliff's delta effect size."""
    x = np.asarray(x)
    y = np.asarray(y)

    greater = np.sum(x[:, None] > y) #counts how many x values beat y values
    less = np.sum(x[:, None] < y)

    return (greater - less) / (len(x) * len(y)) #normalized by total comparisons


def plot_metric_ranges_sig(all_samples_df, metric, distance_ranges):

    for range_name in distance_ranges.keys():

        col = f"{metric}_{range_name}"

        plt.figure(figsize=(8,6))

        ax = sns.boxplot(
            data=all_samples_df,
            x="sample",
            y=col,
            hue="CNA_status",
            hue_order=["neutral","amplification"]
        )

        samples = all_samples_df["sample"].unique()

        pairs = [((s,"neutral"),(s,"amplification")) for s in samples]

        pvalues = []
        deltas = []

        # ---- compute p-values + effect sizes ----
        for s in samples:

            sub = all_samples_df[all_samples_df["sample"]==s]

            neutral = sub[sub["CNA_status"]=="neutral"][col].dropna()
            amp = sub[sub["CNA_status"]=="amplification"][col].dropna()

            if len(neutral)==0 or len(amp)==0:
                pvalues.append(np.nan)
                deltas.append(np.nan)
                continue

            stat,p = mannwhitneyu(
                neutral,
                amp,
                alternative="greater"   # neutral > amplification
            )

            pvalues.append(p)
            deltas.append(cliffs_delta(neutral.values, amp.values))


        annotator = Annotator(
            ax,
            pairs,
            data=all_samples_df,
            x="sample",
            y=col,
            hue="CNA_status",
            hue_order=["neutral","amplification"]
        )

        annotator.configure(
            text_format="star",
            loc="outside"
        )

        annotator.set_pvalues(pvalues)
        annotator.annotate()

        # ---- add effect size labels ----
        for i,s in enumerate(samples):

            sub = all_samples_df[all_samples_df["sample"]==s]
            ymax = sub[col].max()

            ax.text(
                i,
                ymax*1.02,
                f"δ={deltas[i]:.2f}",
                ha="center",
                fontsize=9
            )

        plt.title(f"{metric.upper()} for {range_name}")
        plt.tight_layout()
        plt.show()