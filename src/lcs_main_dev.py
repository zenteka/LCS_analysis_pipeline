'''
Loss of Chromatin segregation analysis pipeline main script.
#: -------------------------------------------------------------------------
# AUTHOR: Karol Piera

# Clean version of Luca Nanni's main LCS_main.py script.
# -------------------------------------------------------------------------

'''
import logging
import os
import argparse
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt

import pyBigWig
from pybedtools.bedtool import BedTool
from cooler import Cooler
from cooltools.lib.numutils import observed_over_expected
from itertools import chain
logging.disable(logging.INFO)

from lcs_plt import (
    plot_inter_chromosomal_map_chrom_level_oe,
    plot_inter_chromosomal_map_chrom_level_foldchange,
    plot_compartment_scores,
    plot_pairwise_compartment_scores,
    plot_shared_boundaries_insulation,
    plot_shared_TOPK_boundaries_insulation,
    plot_intra_chromosomal_OE_contacts,
)

# Reqiete parsing into a function to be able to call it from other scripts
def process_arguments():
    parser = argparse.ArgumentParser(description='Run LCS analysis')
    parser.add_argument("-c", "--comparison", type=str, help='Path to the tsv file containing the desired comparisons')
    parser.add_argument('-d','--intra_chr_res', type=int, help='Intra-chromosomal resolution', default=10000000)
    parser.add_argument('-D','--INTER_CHR_RES', type=int, help='Inter-chromosomal resolution', default=1000000000)
    parser.add_argument('-g','--genome', type=str, help='Genome assembly')
    parser.add_argument('-i', '--INPUT_DIR', type=str, help='Input data directory.')
    parser.add_argument('-o','--output', type=str, help='Path to the output directory')
    parser.add_argument('-r','--HicResolution', type=int, help='Hi-C resolution of cscores', default=10000000)
    parser.add_argument('-w','--INS_W_SIZE', type=int, help='Insulation window size from the main script', default=1000000)
    return parser.parse_args()

# %% Define functions
def load_hic_interchromosomal_dump(spath, chrom_order):
    s_ints = pd.read_csv(spath, sep="\t")
    s_ints['bin_id1'] = s_ints.chr1 + ":" + s_ints.start1.astype(str) + "-" + s_ints.end1.astype(str)
    s_ints['bin_id2'] = s_ints.chr2 + ":" + s_ints.start2.astype(str) + "-" + s_ints.end2.astype(str)
    s_ints_other = s_ints[s_ints.bin_id1 != s_ints.bin_id2].copy()
    x1 = s_ints_other[['chr1', 'start1', 'end1', 'bin_id1']].copy()
    x2 = s_ints_other[['chr2', 'start2', 'end2', 'bin_id2']].copy()
    s_ints_other.loc[:, ['chr1', 'start1', 'end1', 'bin_id1']] = x2.values
    s_ints_other.loc[:, ['chr2', 'start2', 'end2', 'bin_id2']] = x1.values
    s_ints_both = pd.concat([s_ints, s_ints_other], axis=0, ignore_index=True)
    s_ints_both['chr1'] = pd.Categorical(s_ints_both['chr1'].map(lambda x: f"chr{x}" if not x.startswith("chr") else x), categories=chrom_order, ordered=True)
    s_ints_both['chr2'] = pd.Categorical(s_ints_both['chr2'].map(lambda x: f"chr{x}" if not x.startswith("chr") else x), categories=chrom_order, ordered=True)
    return s_ints_both

def get_shared_bins(s1_compartment_scores, s2_compartment_scores):
    # Getting samples' compartment annotations and add a key for each bin
    s1_cscores = s1_compartment_scores\
        .rename(columns = {'cscore': 'cscore1'})\
        [['key', 'cscore1']]

    s2_cscores = s2_compartment_scores\
        .rename(columns = {'cscore': 'cscore2'})\
        [['key', 'cscore2']]
    # Getting only the bins having the same compartment in both samples
    shared_bin_scores = s1_cscores.merge(s2_cscores, on = 'key')    
    # Getting the bin-wise difference of C-SCOREs between sample1 and sample2
    shared_bin_scores['cscoreDiff'] = shared_bin_scores.cscore2 - shared_bin_scores.cscore1
    # If the C-SCORE goes up we have a 'gain', otherwise we have a 'loss'
    shared_bin_scores['cscoreChange'] = shared_bin_scores.cscoreDiff.map(lambda x: "gain" if x > 0 else "loss")
    # changetype = [ 'gain', 'loss']
    # shared_bin_scores['cscoreChange'] = pd.Categorical(shared_bin_scores.cscoreDiff.map(lambda x: changetype[0] if x > 0 else changetype[1]),
    # categories = changetype)
    # We recover the compartment of that bin from the key
    shared_bin_scores['compartment'] = shared_bin_scores.key.map(lambda x: x.split("_")[1])
    shared_bin_num = shared_bin_scores.shape[0]
    return shared_bin_scores, shared_bin_num
    
def get_compartment_stats(shared_bin_scores):
    """
    For each compartment level, we compute the ratio of the number of bins
    loosing segregation VS the number of bins gaining segregation
    """
    comp_stats = (
        # count caterogical values
        shared_bin_scores.groupby('compartment')['cscoreChange'].value_counts()
        .to_frame("n_bins")
        .reset_index()
        .pivot(index='compartment', columns='cscoreChange', values='n_bins')  # Fill missing values with 
        .assign(ratio=lambda x: x.loss / x.gain)
        .reset_index()
    )
    return comp_stats

def get_pairwise_compartment_segregation_changes(s1_compartment_scores, s2_compartment_scores, bins, COMPARTMENTS_8_LABELS):
    # We align the two samples so that bins are matched
    s1_compartment_scores_aligned = s1_compartment_scores.set_index("key").reindex(index = bins)
    s2_compartment_scores_aligned = s2_compartment_scores.set_index("key").reindex(index = bins)
    assert s1_compartment_scores_aligned.index.tolist() == s2_compartment_scores_aligned.index.tolist()

    # We now compute a score for each pair of compartment levels
    pairwise_compartment_segregation_changes = []
    # For each compartment level
    for comp1 in COMPARTMENTS_8_LABELS:
        # We take all the bins belonging to that compartment
        s1_comp1 = s1_compartment_scores_aligned[s1_compartment_scores_aligned.compartment == comp1]
        s2_comp1 = s2_compartment_scores_aligned[s2_compartment_scores_aligned.compartment == comp1]
        assert s1_comp1.index.tolist() == s2_comp1.index.tolist()      
        # For each "other" compartment level
        for comp2 in COMPARTMENTS_8_LABELS:
            # We compute the ratio of the number of bins 
            s12_comp12_diff = s2_comp1[comp2] - s1_comp1[comp2]
            s12_comp12_ratio = s12_comp12_diff[s12_comp12_diff < 0].shape[0] / s12_comp12_diff[s12_comp12_diff > 0].shape[0]
            pairwise_compartment_segregation_changes.append({'compartment1': comp1, 'compartment2': comp2, "ratio": s12_comp12_ratio})
    
    pairwise_compartment_segregation_changes = pd.DataFrame.from_dict(pairwise_compartment_segregation_changes)
    pairwise_compartment_segregation_changes['key'] = pairwise_compartment_segregation_changes.apply(lambda x: " - ".join(sorted([x.compartment1, x.compartment2])), axis=1)
    pairwise_compartment_segregation_changes = pairwise_compartment_segregation_changes.groupby('key')['ratio'].mean().to_frame("ratio").reset_index()
    pairwise_compartment_segregation_changes['compartment1'] = pairwise_compartment_segregation_changes.key.map(lambda x: x.split(" - ")[0])
    pairwise_compartment_segregation_changes['compartment2'] = pairwise_compartment_segregation_changes.key.map(lambda x: x.split(" - ")[1])
    pairwise_compartment_segregation_changes['negLog2Ratio'] = pairwise_compartment_segregation_changes['ratio'].map(lambda x: -1*np.log2(x))
    return pairwise_compartment_segregation_changes[['compartment1', 'compartment2', 'ratio', 'negLog2Ratio']]

def load_bigwig(path, CHROMSIZES):
    bw = pyBigWig.open(path)

    intervals = []
    for chrom in bw.chroms():
        chrom_intervals = bw.intervals(chrom)
        chrom_intervals = pd.DataFrame.from_records(chrom_intervals, columns = ["start", "end", "score"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        chrom_intervals['chr'] = chrom
        chrom_intervals = chrom_intervals[['chr', 'start', 'end', 'score']]
        intervals.append(chrom_intervals)
    intervals = pd.concat(intervals, axis=0, ignore_index=True)
    intervals = intervals[intervals.chr.isin(CHROMSIZES.chr.tolist())].reset_index(drop=True)
    return intervals

def load_boundaries(path, CHROMSIZES):
    sample_boundaries = pd.read_csv(path, sep="\t", header=None, names = ['chr', 'start', 'end', 'name', 'score', 'strand'])
    sample_boundaries['chr'] = sample_boundaries.chr.astype(str)
    if not sample_boundaries.chr.iloc[0].startswith("chr"):
        sample_boundaries['chr'] = 'chr' + sample_boundaries.chr.astype(str)
    sample_boundaries = sample_boundaries[sample_boundaries.chr.isin(CHROMSIZES.chr.tolist())]
    sample_boundaries['start'] -= 1
    sample_boundaries = sample_boundaries[['chr', 'start', 'end', 'score']]
    return sample_boundaries

def get_shared_boundaries(control_boundaries, control_insulation, wgd_boundaries, wgd_insulation):
    shared = BedTool.from_dataframe(control_boundaries[['chr', 'start', 'end', 'score']].sort_values(['chr', 'start', 'end']))\
                       .closest(BedTool.from_dataframe(wgd_boundaries[['chr', 'start', 'end', 'score']].sort_values(['chr', 'start', 'end'])), d=True)\
                       .to_dataframe(names = ['chr', 'start', 'end', 'control_score', "wgd_chr", "wgd_start", "wgd_end", "wgd_score", "distance"])
    shared = shared[shared.distance <= 1]
    shared = shared.merge(control_insulation.rename(columns = {'score': 'control_ins'}), on=['chr', 'start', 'end'])\
                     .merge(wgd_insulation.rename(columns = {'score': 'wgd_ins'}), on=['chr', 'start', 'end'])
    shared['insulation_change'] = shared.apply(lambda y: "loss" if y.wgd_ins > y.control_ins else "gain", axis=1)
    shared['control_ins_rank'] = shared.control_ins.rank()
    shared['wgd_ins_rank'] = shared.wgd_ins.rank()
    return shared

def get_intra_matrices(cool_path, resolution, chrom, balance='KR'):
    """Given a cooler file, it returns the intra-chromosomal interaction matrix of a specified chromosome
    both as counts and as o/e"""
    H = Cooler(cool_path + f"::/resolutions/{resolution}")
    if chrom in H.chromnames:
        H = H.matrix(balance=balance).fetch(chrom)
    else:
        H = H.matrix(balance=balance).fetch(chrom[3:])
    H[~np.isfinite(H)] = 0
    mask = H.sum(axis=0) > 0
    OE, _, _, _ = observed_over_expected(H, mask)
    OE[~mask, :] = 1
    OE[:, ~mask] = 1
    return H, OE, mask

def get_shared_interactions(control_OE, control_mask, wgd_OE, wgd_mask, resolution):
    """Given two intra-chromosomal Hi-C matrices and their good bins masks, it produces
    a Pandas dataframe with one row for each interaction"""

    # get the intersection of good bins
    shared_mask = control_mask & wgd_mask
    good_bins = np.where(shared_mask)[0]

    row, col = np.triu_indices_from(control_OE)
    control_values = control_OE[row, col]
    wgd_values = wgd_OE[row, col]

    shared_interactions = pd.DataFrame({
            'bin1_id': row,
            'bin2_id': col,
            'control': control_values,
            'wgd': wgd_values
        })
    shared_interactions['distance'] = (shared_interactions.bin2_id - shared_interactions.bin1_id)*resolution
    shared_interactions = shared_interactions[shared_interactions.bin1_id.isin(good_bins) & shared_interactions.bin2_id.isin(good_bins)]
    # shared_interactions = shared_interactions[(shared_interactions.control > 0) & (shared_interactions.wgd > 0)]
    return shared_interactions


def get_all_shared_intra_interactions(control_cool_path, wgd_cool_path, resolution, exclude_chroms=None):
    """ Given two Hi-C maps, it returns the aligned O/E interactions for each chromosome at the given resolution
    """
    all_shared_interactions = []

    exclude_chroms = [] if exclude_chroms is None else exclude_chroms

    for chrom in [c for c in CHROMSIZES.chr.tolist() if c not in exclude_chroms]:
        control_H, control_OE, control_mask = get_intra_matrices(control_cool_path, resolution, chrom)
        wgd_H, wgd_OE, wgd_mask = get_intra_matrices(wgd_cool_path, resolution, chrom)

        # some trivial checks
        assert control_H.shape == wgd_H.shape, "WGD and Control matrices have different shapes"
        assert control_OE.shape == wgd_OE.shape, "WGD and Control matrices have different shapes"

        shared_interactions = get_shared_interactions(control_OE, control_mask, wgd_OE, wgd_mask, resolution)
        shared_interactions['chr'] = chrom
        # print(f"{chrom} - {shared_interactions.shape[0]} interactions")
        all_shared_interactions.append(shared_interactions)
    all_shared_interactions = pd.concat(all_shared_interactions, axis=0, ignore_index=True)
    all_shared_interactions = all_shared_interactions[["chr"] + [c for c in all_shared_interactions.columns if c != 'chr']]
    return all_shared_interactions

def get_interaction_type(x):
    if x == "MATERNAL - PATERNAL":
        return "TRANS"
    elif x in ["MATERNAL - MATERNAL", "PATERNAL - PATERNAL"]:
        return "CIS"
    else:
        return "OTHER"

# %% Run main
def main(args):
    '''
    Main function to run the LCS analysis.
    '''
    # Define constants
    TOP_K=300
    MAX_DISTANCE = 20000000
    OUTPUT_PATH = args.output                                                          # Save allthe main plots
    INTER_CHROMOSOMAL_OE_OUTPUT = f"{OUTPUT_PATH}"+"/Interchromosomal_OE"                             # Save the inter-chromosomal Observed/Expected plots for each sample
    os.makedirs(INTER_CHROMOSOMAL_OE_OUTPUT, exist_ok=True)
    CHROMSIZES = pd.read_csv(args.genome, sep="\t", header=None, names = ['chr', 'length'])
    CHROMS_BY_SIZE = CHROMSIZES.sort_values('length', ascending=False).chr.tolist()
    LONG_CHROMS = CHROMS_BY_SIZE[:15]
    SHORT_CHROMS = CHROMS_BY_SIZE[15:]
    
    COMPARTMENTS_8_LABELS = ["{}.{}.{}".format(l1, l2, l3) for l1 in ['A', 'B'] for l2 in [1, 2] for l3 in [1, 2]]
    INSULATION_CHANGE_TO_COLOR = {
		'gain': 'red',
		'loss': 'blue'
	}

    # INPUT
    DATA_PATH = args.INPUT_DIR
    HIC_MAPS_PATH = f"{DATA_PATH}"+"/{hic_sample}/aligned/inter_30.mcool"
    INTER_CHROMOSOMAL_DUMPS_PATH = os.path.join(DATA_PATH, "{hic_sample}", "2_trans", f"binsize_{args.INTER_CHR_RES}.tsv")                                             # DUMP_COOL.py output - inter-chromosomal dumps
    COMPARTMENT_SEGREGATION_PATH = os.path.join(DATA_PATH, "{hic_sample}", "3_saddle", "compartment_scores.tsv")
    INSULATION_PATH = os.path.join(DATA_PATH, "{hic_sample}", "4_insulation", f"insulation_scores_r50000_KR_w{args.INS_W_SIZE}.bigWig")                                          # Insulation boundaries - output of run_module_insulation() function in ./Run_LCS_piepline.sh
    BOUNDARIES_PATH = os.path.join(DATA_PATH, "{hic_sample}", "4_insulation", f"insulation_scores_r50000_KR_w{args.INS_W_SIZE}_boundaries.bed")                                # Insulation score - output of run_module_insulation() function in ./Run_LCS_piepline.sh
    SHARED_BINS_BY_COMPARISON = {}
    
    # Get comparison file from commonly used formats
    if args.comparison.endswith('.tsv'):
        LCS_COMPARISONS = pd.read_csv(args.comparison, sep='\t',
                                      names=['sample1', 'sample2'])
        print(LCS_COMPARISONS)
    elif args.comparison.endswith('.csv'):
        LCS_COMPARISONS = pd.read_csv(args.comparison,
                                      names=['sample1', 'sample2'])
    elif args.comparison.endswith('.xlsx'):
        LCS_COMPARISONS = pd.read_excel(args.comparison,
                                        names = ['sample1', 'sample2'])
    else:
        print("Input file format not recognized. Please provide a valid tsv, csv or xlsx file")
        raise ValueError("Input file format not recognized. Please provide a valid tsv, csv or xlsx file")
    REPORT = {}
    OVERWRITE = True
    for t in chain(LCS_COMPARISONS.itertuples()):
        print("---------------------------------")
        print(f"{t.sample1} VS {t.sample2}")
        print("---------------------------------")
        REPORT[f"{t.sample1} VS {t.sample2}"] = []
        
        print("------ Inter-chromosomal contacts ------")        
        sample_dump1 = load_hic_interchromosomal_dump(os.path.join(INTER_CHROMOSOMAL_DUMPS_PATH.format(hic_sample = t.sample1)), chrom_order=CHROMS_BY_SIZE)
        name = os.path.join(INTER_CHROMOSOMAL_OE_OUTPUT, f"{t.sample1}_InterChromosomal_OE.pdf")

        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plot_inter_chromosomal_map_chrom_level_oe(sample_dump1, t.sample1, LONG_CHROMS)
            fig.savefig(name)
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
        
            sample_dump2 = load_hic_interchromosomal_dump(os.path.join(INTER_CHROMOSOMAL_DUMPS_PATH.format(hic_sample = t.sample2)), 
                                                          chrom_order=CHROMS_BY_SIZE)
            name = os.path.join(INTER_CHROMOSOMAL_OE_OUTPUT, f"{t.sample2}_InterChromosomal_OE.pdf")
            if OVERWRITE or (not os.path.exists(name)):
                fig, ax = plot_inter_chromosomal_map_chrom_level_oe(sample_dump2, t.sample2, LONG_CHROMS)
                fig.savefig(name)
                plt.close(fig)
                REPORT[f"{t.sample1} VS {t.sample2}"].append(name)

        s12_dump = sample_dump1[['chr1', 'chr2', 'balanced_rescaled']].rename(columns = {'balanced_rescaled': 'score1'})\
                                                                      .merge(sample_dump2[['chr1', 'chr2', 'balanced_rescaled']].rename(columns = {'balanced_rescaled': 'score2'}), 
                                                                             on = ['chr1', 'chr2'])\
                                                                      .assign(foldchange = lambda x: x.score2 / x.score1,
                                                                              chr1_length_type=lambda x: x.chr1.map(lambda y: 'L' if y in LONG_CHROMS else 'S'),
                                                                              chr2_length_type=lambda x: x.chr2.map(lambda y: 'L' if y in LONG_CHROMS else 'S'),
                                                                              pair_length_type=lambda x: x.apply(lambda y: "".join(sorted([y.chr1_length_type, y.chr2_length_type])), axis=1))

        comparison_path=os.path.join(OUTPUT_PATH, f"{t.sample1}_VS_{t.sample2}")
        os.makedirs(comparison_path, exist_ok=True)

        name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_InterChromosomal_Ratio.pdf")
        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plot_inter_chromosomal_map_chrom_level_foldchange(s12_dump, LONG_CHROMS)
            fig.savefig(name)
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
            
        print("------ Compartment scores ------")
        compartment_scores1 = pd.read_csv(COMPARTMENT_SEGREGATION_PATH.format(hic_sample=t.sample1), sep="\t")\
                                .assign(key = lambda x: x.chr + ":" + \
                                        x.start.astype(str) + "-" + \
                                        x.end.astype(str) + "_" + \
                                        x.compartment)
        compartment_scores2 = pd.read_csv(COMPARTMENT_SEGREGATION_PATH.format(hic_sample=t.sample2), sep="\t")\
                                .assign(key = lambda x: x.chr + ":" + \
                                        x.start.astype(str) + "-" + \
                                        x.end.astype(str) + "_" + \
                                        x.compartment)
        shared_bin_scores, bin_num = get_shared_bins(compartment_scores1, compartment_scores2)
        SHARED_BINS_BY_COMPARISON[f"{t.sample1}_VS_{t.sample2}"] = bin_num
        comp_stats = get_compartment_stats(shared_bin_scores)
        
        name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_Compartment_Scores.pdf")
        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            plot_compartment_scores(ax, comp_stats, COMPARTMENTS_8_LABELS)
            fig.savefig(name, bbox_inches='tight')
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
        
        pairwise_compartment_segregation_changes = get_pairwise_compartment_segregation_changes(
                compartment_scores1, 
                compartment_scores2, 
                shared_bin_scores.key.values,
                COMPARTMENTS_8_LABELS
            )
    
        name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_Compartment_PairedScores.pdf")
        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plot_pairwise_compartment_scores(ax, pairwise_compartment_segregation_changes, t.sample1, t.sample2, COMPARTMENTS_8_LABELS)
            fig.savefig(name, bbox_inches='tight')
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
            
        print("----- Boundary LCS -----")
        sample1_insulation = load_bigwig(INSULATION_PATH.format(hic_sample=t.sample1), CHROMSIZES)
        sample1_boundaries = load_boundaries(BOUNDARIES_PATH.format(hic_sample=t.sample1), CHROMSIZES)
        sample2_insulation = load_bigwig(INSULATION_PATH.format(hic_sample=t.sample2), CHROMSIZES)
        sample2_boundaries = load_boundaries(BOUNDARIES_PATH.format(hic_sample=t.sample2), CHROMSIZES)
    
        shared_boundaries = get_shared_boundaries(sample1_boundaries, sample1_insulation, sample2_boundaries, sample2_insulation)
    
        with pd.option_context('mode.use_inf_as_na', True):
            shared_boundaries = shared_boundaries.dropna(subset=['control_ins', 'wgd_ins'])
        
        # finding the threshold of top-k boundaries
        control_thresh = shared_boundaries[shared_boundaries.control_ins_rank <= TOP_K].control_ins.max()
        wgd_thresh = shared_boundaries[shared_boundaries.wgd_ins_rank <= TOP_K].wgd_ins.max()
        topk_thresh = max(control_thresh, wgd_thresh)
    
        name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_Insulation_AllBoundaries.pdf")
        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            plot_shared_boundaries_insulation(ax, shared_boundaries, topk_thresh, TOP_K, INSULATION_CHANGE_TO_COLOR,
                                              xname=t.sample1, yname=t.sample2)
            fig.savefig(name)
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
        
        name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_Insulation_Top{TOP_K}.pdf")
        if OVERWRITE or (not os.path.exists(name)):
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            plot_shared_TOPK_boundaries_insulation(ax, shared_boundaries, topk_thresh, INSULATION_CHANGE_TO_COLOR,
                                                   xname=t.sample1, yname=t.sample2)
            fig.savefig(name)
            plt.close(fig)
            REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
            
        print("-----Skipping Intra-chromosomal contacts ------")    
        # name = os.path.join(comparison_path, f"{t.sample1}_VS_{t.sample2}_IntraChromosomal_OE.pdf")
        # if OVERWRITE or (not os.path.exists(name)):
        #     all_shared_interactions = get_all_shared_intra_interactions(HIC_MAPS_PATH.format(hic_sample=t.sample1), 
        #                                     HIC_MAPS_PATH.format(hic_sample=t.sample2),
        #                                                                     parser.parse_args().intra_chr_res,
        #                             ['chrY', 'chrM'])
        
        #     nonzero_shared_interactions = all_shared_interactions[(all_shared_interactions.control > 0) & (all_shared_interactions.wgd > 0)].copy()
        #     nonzero_shared_interactions['fc'] = nonzero_shared_interactions.wgd / nonzero_shared_interactions.control
        #     nonzero_shared_interactions['log2_control'] = nonzero_shared_interactions['control'].map(np.log2)
        #     nonzero_shared_interactions['log2_wgd'] = nonzero_shared_interactions['wgd'].map(np.log2)
        #     nonzero_shared_interactions['log2_fc'] = nonzero_shared_interactions.fc.map(np.log2)
        #     nonzero_shared_interactions['fc_of_logs'] = nonzero_shared_interactions.log2_wgd / nonzero_shared_interactions.log2_control
        #     nonzero_shared_interactions['change_type'] = nonzero_shared_interactions.apply(lambda x: "loss of signal" if abs(x.log2_control) > abs(x.log2_wgd) else "gain of signal", axis=1)
        #     nonzero_shared_interactions_mindist = nonzero_shared_interactions[(nonzero_shared_interactions.distance <= MAX_DISTANCE)]
        
        #     fig, ax = plt.subplots(figsize=(5, 5))
        #     plot_intra_chromosomal_OE_contacts(ax, nonzero_shared_interactions_mindist, sample1 = t.sample1, sample2 = t.sample2)                
        #     fig.savefig(name)
        #     plt.close(fig)                
        # REPORT[f"{t.sample1} VS {t.sample2}"].append(name)
    # Make SHARED_BINS_DF from dict
    SHARED_BINS_DF = pd.DataFrame.from_dict(SHARED_BINS_BY_COMPARISON, orient='index').reset_index().rename(columns = {'index': 'Comparison', 0: 'Number_of_Shared_Bins'})
    SHARED_BINS_DF.to_csv(os.path.join(OUTPUT_PATH, "CompartmentLCS_Shared_Bins_By_Comparison.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    args = process_arguments()
    main(args)

        
