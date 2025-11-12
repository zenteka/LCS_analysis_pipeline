"""
LCS main plotting functions.

"""
import os
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from statannot import add_stat_annotation
from scipy.stats import ranksums, wilcoxon, mannwhitneyu

# %% Definitions
def plot_inter_chromosomal_map_chrom_level_oe(sample_dump, sample_name, LONG_CHROMS):
    '''
        Plot the inter-chromosomal map of the Hi-C data
    '''     
    M = pd.pivot_table(sample_dump, index='chr1', columns='chr2', values='balanced_rescaled',
            observed=True).sort_index(axis=0).sort_index(axis=1)
    n_chroms = M.shape[0] - 1
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    fig.suptitle(f"{sample_name}", fontsize=16)
    sns.heatmap(M, cmap='bwr',
                square=True,
                center=1/n_chroms,
                vmin=0.5/n_chroms, vmax=1.5/n_chroms,
                linewidth=0.1,
                cbar_kws={'label': 'balanced counts', 'shrink': 0.5, 'aspect': 8},
                ax=ax)
    patch = plt.Polygon([[0, 1.1], [0, 1.05], [1, 1.05]], transform=ax.transAxes, color='grey')
    fig.add_artist(patch)
    ax.text(0.5, 1.1, "Chromosome size", transform=ax.transAxes, ha='center')
    ax.plot([0, M.shape[0]], [len(LONG_CHROMS), len(LONG_CHROMS)], color='black', linewidth=2)
    ax.plot([len(LONG_CHROMS), len(LONG_CHROMS)], [0, M.shape[0]], color='black', linewidth=2)
    ax.set_xlabel("Chromosome 2")
    ax.set_ylabel("Chromosome 1")
    return fig, ax

def plot_inter_chromosomal_map_chrom_level_foldchange(ratio, LONG_CHROMS):
    M = pd.pivot_table(ratio[ratio.chr1 != ratio.chr2],
                       index='chr1', 
                       columns='chr2', 
                       values='foldchange', 
                       observed=True)\
          .sort_index(axis=0)\
          .sort_index(axis=1)
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    sns.heatmap(M.applymap(lambda x: np.log2(x)), 
            vmin=-0.5, vmax=0.5, 
                cmap='RdBu_r',
                linewidth=0.1,
                cbar_kws={'label': '$log_2$(fold-change)\nof balanced counts', 'shrink': 0.5, 'aspect': 8}, 
                square=True,
                ax = ax)
    patch = plt.Polygon([[0, 1.1], [0, 1.05], [1, 1.05]], transform=ax.transAxes, color='grey')
    fig.add_artist(patch)
    ax.text(0.5, 1.1, "Chromosome size", transform=ax.transAxes, ha='center')
    ax.plot([0, M.shape[0]], [len(LONG_CHROMS), len(LONG_CHROMS)], color='black', linewidth=2)
    ax.plot([len(LONG_CHROMS), len(LONG_CHROMS)], [0, M.shape[0]], color='black', linewidth=2)
    ax.set_xlabel("Chromosome 2")
    ax.set_ylabel("Chromosome 1")

    ratioVals = ratio[ratio.chr1 < ratio.chr2]
    pvalString = ""
    for class1, class2 in [('LL', 'LS'), ("LL", "SS"), ('LS', "SS")]:
        pval = mannwhitneyu(ratioVals.loc[ratioVals.pair_length_type == class1, "foldchange"].values, 
                            ratioVals.loc[ratioVals.pair_length_type == class2, "foldchange"].values)\
                            .pvalue
        pvalString += "$P_{" + class1 + " - " + class2 + "}$ = " + f"{pval:.2E}\n"
    ax.text(1.05, 1.05, pvalString, transform=ax.transAxes, ha='left')
    return fig, ax


def plot_compartment_scores(ax, comp_stats, COMPARTMENTS_8_LABELS):
    # sns.pointplot(data = comp_stats, 
    #               x = "compartment", 
    #               y = "ratio", 
    #               order = COMPARTMENTS_8_LABELS,
    #               linestyle='none',
    #               color='black',
    #               scale = 1.5,
    #               markers="o",
    #               ax=ax)
    # Fix plot
    sns.pointplot(
        data=comp_stats,
        x ="compartment",
        y ="ratio",
        order=COMPARTMENTS_8_LABELS,
        errorbar=None,  # or keep ci="sd", etc.
        markers="o",
        color="black",
        linestyles='None',
        markersize=5,
        ax=ax,
    )

    ax.set_ylabel("Loss of compartment segregation\n($\Delta$C-score < 0 / $\Delta$C-score > 0)")
    ax.set_xlabel("")
    xticks = ax.get_xticks()
    #ax.set_xticks(xticks, rotation_mode='anchor', horizontalalignment='right')
    sns.despine(trim=True, ax=ax)


def plot_pairwise_compartment_scores(ax, pairwise_compartment_segregation_changes, s1, s2, COMPARTMENTS_8_LABELS):
        
    ab_aabb_values = pairwise_compartment_segregation_changes.assign(
                                                h1 = lambda x: x.compartment1.map(lambda y: y.split(".")[0]),
                                                h2 = lambda x: x.compartment2.map(lambda y: y.split(".")[0]),
                                                ll = lambda x: x.apply(lambda y: "AB" if y.h1 != y.h2 else "AA-BB", axis=1))\
                            .groupby("ll")['ratio'].agg(list)\
                            .to_dict()
    statistic, pvalue = ranksums(ab_aabb_values['AA-BB'], ab_aabb_values['AB'])

    stats = pd.Series({
            'wilcoxon_statistic': statistic,
            'wilcoxon_pvalue': pvalue
        })

    pairwise_stats = pd.pivot_table(pairwise_compartment_segregation_changes, 
                   index = 'compartment1', 
                   columns = 'compartment2', 
                   values = "negLog2Ratio")\
        .reindex(index=COMPARTMENTS_8_LABELS, columns=COMPARTMENTS_8_LABELS)\
        .T
    
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    sns.heatmap(pairwise_stats, 
            cmap='RdBu_r',vmin=-4, vmax=4, 
            center=0, linewidth=1, square=True,
                cbar_kws=dict(shrink=0.5, 
                          aspect=10, 
                          location='bottom',
                          label="Gain vs. Loss of\ncompartment segregation\nfold-change (log2)"),
            ax = ax)
    yticklabels = ax.get_yticklabels()
    xticklabels = ax.get_xticklabels()
    ax.set_yticklabels(yticklabels)
    ax.set_xticklabels(xticklabels) #rotation=0)
    ax.set_title(f"{s2} vs {s1}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.text(x = 0.45, y = .8, s = f"P-value(AB vs AA+BB):\n{stats.wilcoxon_pvalue:.2e}", transform=ax.transAxes, size=10)
    return fig, ax

def plot_shared_boundaries_insulation(ax,
                                      shared_boundaries, 
                                      topk_thresh,
                                      TOP_K,
                                      INSULATION_CHANGE_TO_COLOR,
                                      xname='sample1', 
                                      yname="sample2"):
    sns.scatterplot(data=shared_boundaries, 
                x='control_ins', 
                y='wgd_ins', 
                hue='insulation_change', 
                palette = INSULATION_CHANGE_TO_COLOR,
                alpha=0.4,
                rasterized=False,
                ax = ax)
    xmin, xmax = min(shared_boundaries.control_ins.min(), shared_boundaries.wgd_ins.min()), max(shared_boundaries.control_ins.max(), shared_boundaries.wgd_ins.max())
    ax.plot([xmin, xmax], [xmin, xmax], linestyle = '--', color = 'black', linewidth=2)

    df = shared_boundaries[(shared_boundaries.control_ins <= topk_thresh) & (shared_boundaries.wgd_ins <= topk_thresh)]    # find top-right corner
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.vlines([topk_thresh], ymin, topk_thresh, color = 'black', linestyle = '-', linewidth=1)
    ax.hlines([topk_thresh], xmin, topk_thresh, color = 'black', linestyle = '-', linewidth=1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.text(x=topk_thresh - 0.8, y = ymin + 0.4, s = f"Top-{TOP_K}\nboundaries", ha='center')
    
    pvalue = wilcoxon(shared_boundaries.control_ins.values, shared_boundaries.wgd_ins.values).pvalue
    ax.text(1, 1.05, f'Wilcoxon p = {pvalue:.2E}',
             horizontalalignment='right',
             verticalalignment='center',
             transform = ax.transAxes)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    sns.despine(ax=ax)

def plot_shared_TOPK_boundaries_insulation(ax,
                                           shared_boundaries, 
                                           topk_thresh,
                                           INSULATION_CHANGE_TO_COLOR,
                                           xname='sample1', 
                                           yname="sample2",
                                           ):
    topk_boundaries = shared_boundaries[(shared_boundaries.control_ins <= topk_thresh) & (shared_boundaries.wgd_ins <= topk_thresh)]
    sns.scatterplot(data=topk_boundaries, 
                x='control_ins',
                y='wgd_ins',
                hue='insulation_change', 
                palette = INSULATION_CHANGE_TO_COLOR,
                alpha=0.4,
                rasterized=False,
            ax = ax)
    xmin, xmax = min(topk_boundaries.control_ins.min(), topk_boundaries.wgd_ins.min()), max(topk_boundaries.control_ins.max(), topk_boundaries.wgd_ins.max())
    ax.plot([xmin, xmax], [xmin, xmax], linestyle = '--', color = 'black', linewidth=2)
    # add number of gains of the linestyle
    ax.text(x = 0.5, y = 0.1, s = f"n = {topk_boundaries[topk_boundaries.insulation_change == 'gain'].shape[0]}", transform=ax.transAxes, ha='center')
    # add number of losses below line
    ax.text(x = 0.5, y = 0.9, s = f"n = {topk_boundaries[topk_boundaries.insulation_change == 'loss'].shape[0]}", transform=ax.transAxes, ha='center')
    
    
    pvalue = wilcoxon(topk_boundaries.control_ins.values, topk_boundaries.wgd_ins.values).pvalue
    ax.text(1, 1.05, f'Wilcoxon p = {pvalue:.2E}',
             horizontalalignment='right',
             verticalalignment='center',
             transform = ax.transAxes)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    sns.despine(ax=ax)

def plot_intra_chromosomal_OE_contacts(ax, nonzero_shared_interactions_mindist, sample1="sample1", sample2="sample2"):
    sns.scatterplot(data=nonzero_shared_interactions_mindist, 
                    x='log2_control', 
                    y= 'log2_wgd', 
                    hue='change_type', 
                    palette={'gain of signal': 'red', 'loss of signal': 'blue'}, 
                    alpha=0.5,
                    s=3,
                    rasterized=True,
                    ax=ax)
    ax.set_xlim(-7, 3)
    ax.set_ylim(-7, 3)
    ax.plot([-7, 3], [-7, 3], color='black', linestyle='--', linewidth=1)
    ax.axhline(0, color = 'black', linestyle='--', linewidth=1)
    ax.axvline(0, color = 'black', linestyle='--', linewidth=1)
    ax.set_xlabel(f"{sample1} O/E (log2)")
    ax.set_ylabel(f"{sample2} O/E (log2)")
    sns.despine(ax = ax)

