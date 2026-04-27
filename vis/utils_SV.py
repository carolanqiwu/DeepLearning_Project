import pandas as pd
import numpy as np
from pyliftover import LiftOver
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import mannwhitneyu

def liftover_sv_to_hg38(sv_df, chain_file="/data1/lesliec/carolw/genome/hg19ToHg38.over.chain.gz"):
    """
    LiftOver SV breakpoints from GRCh37 (hg19) to hg38.

    Parameters
    ----------
    sv_df : DataFrame
        SV table containing:
        Site1_Chromosome, Site1_Position,
        Site2_Chromosome, Site2_Position

    chain_file : str
        Path to hg19ToHg38 chain file.

    Returns
    -------
    DataFrame with hg38 coordinates added.
    """

    lo = LiftOver(chain_file)

    def convert(chrom, pos):
        chrom = f"chr{chrom}"
        result = lo.convert_coordinate(chrom, pos)

        if len(result) == 0:
            return None, None

        new_chrom, new_pos, _, _ = result[0]
        return new_chrom, int(new_pos)

    sv_df = sv_df.copy()

    sv_df[["Site1_chr_hg38","Site1_pos_hg38"]] = sv_df.apply(
        lambda r: pd.Series(convert(r["Site1_Chromosome"], r["Site1_Position"])),
        axis=1
    )

    sv_df[["Site2_chr_hg38","Site2_pos_hg38"]] = sv_df.apply(
        lambda r: pd.Series(convert(r["Site2_Chromosome"], r["Site2_Position"])),
        axis=1
    )

    return sv_df

def create_bins_from_chrom_sizes(chrom_sizes_df, bin_size=10000):
    """
    Create genomic bins for all chromosomes.

    Parameters
    ----------
    chrom_sizes_df : pd.DataFrame
        DataFrame with columns ['chr', 'size']
    bin_size : int
        Size of each bin (default=10kb)

    Returns
    -------
    bins_df : pd.DataFrame
        Columns: ['chrom', 'start', 'end']
    """
    bins = []
    for _, row in chrom_sizes_df.iterrows():
        # skip non-standard chromosomes if needed
        if str(row["chr"]) in ["ALL", "M"]:
            continue
        chrom = 'chr' + str(row['chr']) if not str(row['chr']).startswith('chr') else str(row['chr'])
        chrom_len = int(row['size'])
        for start in range(0, chrom_len, bin_size):
            end = min(start + bin_size, chrom_len)
            bins.append((chrom, start, end))
    
    bins_df = pd.DataFrame(bins, columns=['chrom', 'start', 'end'])
    return bins_df

def annotate_sv_bins(bins_df, data_sv, flank_size=50000):
    """
    Annotate bins as SV_flank, SV_pair_region (intra only), or background.

    Parameters
    ----------
    bins_df : pd.DataFrame
        Columns: ['chrom', 'start', 'end']
    data_sv : pd.DataFrame
        Columns: ['Site1_chr_hg38','Site1_pos_hg38','Site2_chr_hg38','Site2_pos_hg38']
    flank_size : int
        Size of flanking bins to annotate around breakpoints (default 50 kb)

    Returns
    -------
    bins_df : pd.DataFrame
        Adds column 'SV_category' with values: 'SV_flank', 'SV_pair_region', 'background'
    """
    bins_df = bins_df.copy()
    bins_df['SV_category'] = 'background'

    for _, sv in data_sv.iterrows():
        chrom1, chrom2 = sv['Site1_chr_hg38'], sv['Site2_chr_hg38']
        pos1, pos2 = sv['Site1_pos_hg38'], sv['Site2_pos_hg38']

        # --- Annotate flanking bins for both intra and inter ---
        # Flank around Site1
        mask_flank1 = (
            (bins_df['chrom'] == chrom1) &
            (bins_df['end'] >= pos1 - flank_size) &
            (bins_df['start'] <= pos1 + flank_size)
        )
        bins_df.loc[mask_flank1, 'SV_category'] = 'SV_flank'

        # Flank around Site2
        mask_flank2 = (
            (bins_df['chrom'] == chrom2) &
            (bins_df['end'] >= pos2 - flank_size) &
            (bins_df['start'] <= pos2 + flank_size)
        )
        bins_df.loc[mask_flank2, 'SV_category'] = 'SV_flank'

        # --- Annotate SV_pair_region for intra-chromosomal SVs only ---
        if chrom1 == chrom2:
            sv_start, sv_end = sorted([pos1, pos2])
            mask_pair = (
                (bins_df['chrom'] == chrom1) &
                (bins_df['start'] > sv_start + flank_size) &  # exclude flanks
                (bins_df['end'] < sv_end - flank_size)
            )
            bins_df.loc[mask_pair, 'SV_category'] = 'SV_pair_region'

    return bins_df

def cliffs_delta(x, y):
    """Compute Cliff's delta effect size."""
    x = np.asarray(x)
    y = np.asarray(y)

    greater = np.sum(x[:, None] > y) #counts how many x values beat y values
    less = np.sum(x[:, None] < y)

    return (greater - less) / (len(x) * len(y)) #normalized by total comparisons


def plot_SV_metric_ranges_sig(all_samples_df, metric, distance_ranges):

    for range_name in distance_ranges.keys():

        col = f"{metric}_{range_name}"

        plt.figure(figsize=(8,6))

        ax = sns.boxplot(
            data=all_samples_df,
            x="sample",
            y=col,
            hue="SV_category",
            hue_order=["background","SV_region"]
        )

        samples = all_samples_df["sample"].unique()

        pairs = [((s,"background"),(s,"SV_region")) for s in samples]

        pvalues = []
        deltas = []

        # ---- compute p-values + effect sizes ----
        for s in samples:

            sub = all_samples_df[all_samples_df["sample"]==s]

            neutral = sub[sub["SV_category"]=="background"][col].dropna()
            amp = sub[sub["SV_category"]=="SV_region"][col].dropna()

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
            hue="SV_category",
            hue_order=["background","SV_region"]
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