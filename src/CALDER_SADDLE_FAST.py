'''
Author: Karol Piera and Luca Nanni. 

Saddle analysis adapted from Luca Nanni.

Major modifications:
- Removed statannot dependecy (incompatible with seaborn), replaced with statannotations.
- Major data pandas manipulations are now done using polars
- redundant intermediate variables removed
- preprcoes functions are now in utils.py module
- Plot functions are now in saddle_plt.py module

# TODO:
- refactor quantile computation in polars (see githhub issue https://github.com/pola-rs/polars/issues/18236)
- see if fanc dependency can be removed
'''
import os
import subprocess
import logging
import argparse

import fanc
import numpy as np
import pandas as pd
import polars as pl

from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from pybedtools.bedtool import BedTool

# Get saddle modules
import saddle_plt
import saddle_utils 

import matplotlib.pyplot as plt
import seaborn as sns

CALDER_HEADER=['chr', 'start', 'end',
               'compartment_label', 'domain_rank_8_levels',
               'strand', 'tickStart', 'tickEnd', 'color']
COMPARTMENTS_8_LABELS = ["{}.{}.{}".format(l1, l2, l3) for l1 in ['A', 'B'] for l2 in [1, 2] for l3 in [1, 2]]

def set_logger():
    logger = logging.getLogger("CALDERSaddle")
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(name)-12s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger

def process_arguments():
    parser = argparse.ArgumentParser(description='Analyse compartment segregation in a cool file (Saddle plots)')
    parser.add_argument("hic_path", type=str, help='Path to the mcool file')
    parser.add_argument("compartments_path", type=str, help='Path to the CALDER bed file')
    parser.add_argument('output_path', type=str, help='Output path where to store the dumped result')
    parser.add_argument('--hic_resolution', type=int, const=10000000, nargs='?', help='Resolution of the Hi-C data to use')
    parser.add_argument('--compartments_resolution', type=int, const=50000, nargs='?', help='Resolution of the compartment calls')
    parser.add_argument('--n_tiles', type=int, const=16, nargs='?', help='Number of percentiles to use')
    parser.add_argument("--excludeChroms", type=str, default=None, help="Chromosomes to exclude")
    args = parser.parse_args()
    return args

def main():

    logger = set_logger()
    
    logger.info("Binning compartments")
    args = process_arguments()
    hic_path=args.hic_path
    hic_resolution=args.hic_resolution
    compartments_path=args.compartments_path
    compartments_resolution=args.compartments_resolution
    excluded_chroms = args.excludeChroms.split(",") if args.excludeChroms is not None else []
    n_tiles = args.n_tiles
    output_path = args.output_path
    
    logger.info("Compartment segregation analysis")
    logger.info("Parameters:")
    logger.info(f"- Hi-C: {hic_path}")
    logger.info(f"- Hi-C resolution: {hic_resolution}")
    logger.info(f"- Compartments: {compartments_path}")
    logger.info(f"- Compartments resolution: {compartments_resolution}")
    # Default list is long and pollutes the log
    #logger.info(f"- Chromosomes to exclude: " + ", ".join(excluded_chroms))  
    logger.info(f"- N. tiles: {n_tiles}")
    logger.info(f"- Output path: {output_path}")

    os.makedirs(output_path, exist_ok=True)    
    PLT_PATH = os.path.join(output_path, "plots")
    os.makedirs(PLT_PATH, exist_ok=True)

    # Init
    all_tile_level_path = os.path.join(output_path, f"tiles_{n_tiles}.tsv")
    all_comp_level_path = os.path.join(output_path, f"compartments.tsv")
    all_domain_level_path = os.path.join(output_path, f"domains.tsv")
    all_chrom_compartment_scores_path = os.path.join(output_path, f"compartment_scores.tsv")
    
    all_tile_level = []
    all_comp_level = []
    all_domain_level = []
    all_chrom_compartment_scores = []
    
    hic = fanc.load(hic_path + f"::resolutions/{hic_resolution}")

    compartments_bed = BedTool.from_dataframe(
        saddle_utils.load_compartments(compartments_path, CALDER_HEADER)
    ).sort()
    
    bins = BedTool().window_maker(b = compartments_bed, w = compartments_resolution)
    compartments_binned = bins.map(compartments_bed,
                                   c = [5, 6, 7],
                                   o = ["distinct", 'mean', 'distinct'])\
                                  .to_dataframe(names =\
                                                ['chr', 'start', 'end', "compartment", 'domain_rank', 'compartment_id'])
    
    for chrom in filter(lambda x: x not in excluded_chroms, hic.chromnames):
        print(chrom)
        logger.info(f"Analysing {chrom}")
        logger.info("- Extracting O/E values")
        c_pixels = pl.DataFrame(saddle_utils.get_oe_values(hic, chrom))
        
        # Get the binned compartments for the chromosome
        # - Extract the Hi-C bin (in the resolution of the Hi-C data) from the start of the compartment bin
        # - Assign a tile to each bin, based on its domain rank, using percentiles    
        # Comment about tiles - tiles are created based on the domain rank of the compartments. 
        logger.info("- Extracting compartments for the chromosome")
    
        chrom = "chr" + str(chrom) if not str(chrom).startswith("chr") else chrom
        # TODO: cut and quct issue. 
        # Polars cut and qcut are not consistent with pandas.
        # Find a faster work-around if possible. For now, stick with Luca's implementation.
        chrom_comps = compartments_binned[compartments_binned.chr == chrom]\
            .assign(bin = lambda x: x.start//hic_resolution,
                    tile = lambda x: pd.cut(x.domain_rank, bins = n_tiles, labels = list(range(n_tiles))))\
                    [['bin', 'tile', 'start', "compartment", 'domain_rank', 'compartment_id']]

        chrom_comps = pl.DataFrame(chrom_comps)                  # Convert to polars
        # Associate to each pixel in the Hi-C matrix the compartment information of its two bins    
        logger.info("- Associating compartments to Hi-C interactions")
    
        # Perform the joins and remove duplicates
        chrom_pixels_with_comp = c_pixels \
            .join(chrom_comps.rename({col: f"{col}1" for col in chrom_comps.columns}), on="bin1", how="inner") \
            .join(chrom_comps.rename({col: f"{col}2" for col in chrom_comps.columns}), on="bin2", how="inner") \
            .unique().filter(
                (pl.col("bin1") <= pl.col("bin2")) & (pl.col("start1") <= pl.col("start2"))
            )
    
        assert chrom_comps.shape[0]*(chrom_comps.shape[0] - 1)/2 + chrom_comps.shape[0] == chrom_pixels_with_comp.shape[0]
    
        logger.info("- Commputing segregation scores")
        #################### Segregation scores #######################
        local_scores_pl = chrom_pixels_with_comp.with_columns([
            (1 - (pl.col("domain_rank2") - pl.col("domain_rank1")).abs()).alias("compSim"),
        (pl.col("oe") * (1 - (pl.col("domain_rank2") - pl.col("domain_rank1")).abs())).alias("weight")
        ]).select([
            "start1", "start2", "oe", "domain_rank1", "domain_rank2", "compartment1", "compartment2", "compSim", "weight"
        ]).sort(
            by=["start1", "start2"]
        )
    
        lower_local_scores_pl = local_scores_pl.filter(pl.col("start1") < pl.col("start2")).sort(
            by=["start1", "start2"]
        )
        both_local_scores = pl.concat([
            local_scores_pl,
            lower_local_scores_pl.rename(
                {"start1": "start2", "start2": "start1",  "domain_rank1": "domain_rank2", "domain_rank2": "domain_rank1",
                 "compartment1": "compartment2", "compartment2": "compartment1"}
            ).select(
                ["start1", "start2", "oe", "domain_rank1", "domain_rank2", "compartment1", "compartment2", "compSim", "weight"]
            )], how="vertical" )
        
        # Computing the correlation between O/E values and compartment similarity
        scores = both_local_scores.group_by("start1").agg(
            pl.corr("compSim", "oe", method="spearman").alias("oe_sim_corr")
        ).rename({"start1": "start"})
    
        comp_scores = both_local_scores.group_by(['start1', 'compartment1', 'compartment2']).agg(
            pl.sum("oe").alias("avg_oe")).pivot(
                index = ['start1', 'compartment1'], on = 'compartment2', values = 'avg_oe'
            ).sort('start1').select(
                ["start1", "compartment1"] + COMPARTMENTS_8_LABELS
            )
        total = len(comp_scores)
        _ = chrom_comps.group_by(['compartment']).agg(
            pl.len().alias("count")
        ).with_columns(
            (pl.col("count") / total).alias("expected")
        ).drop("count")
        expected = _.drop('compartment').transpose(
            include_header=False,
            column_names=_['compartment'].to_list()
        )
    
        # Divdide by expected for a given compartment domain
        comp_scores_norm_pl = comp_scores.with_columns(
            (pl.col(c) / expected[c]).alias(c)
        for c in COMPARTMENTS_8_LABELS
        ).with_columns(
            row_sum=pl.sum_horizontal(COMPARTMENTS_8_LABELS)
        ).with_columns(
            # normalize by row_sum division
        (pl.col(c) / pl.col("row_sum")).alias(c)
        for c in COMPARTMENTS_8_LABELS
        ).drop("row_sum")
        
        # segregation scores
        # TODO: This is suboptimal. Figure out how to do it in polars.
        comp_scores = comp_scores.to_pandas().set_index(['start1', 'compartment1'])
        comp_scores_norm = comp_scores_norm_pl.to_pandas().set_index(['start1', 'compartment1'])    
        comp_scores['cscore'] = comp_scores.apply(lambda x: x[x.name[1]], axis=1)
        comp_scores_norm['cscore'] = comp_scores_norm.apply(lambda x: x[x.name[1]], axis=1)
        # Move back to polars
        comp_scores = pl.DataFrame(comp_scores.reset_index())
        comp_scores_norm_pl = pl.DataFrame(comp_scores_norm.reset_index())
    
        # Pull results together and add missing columns
        chrom_compartment_scores = chrom_comps.join(
            (both_local_scores.group_by('start1').agg(
            # Coverage
            pl.len().alias("n_interactions"),
            pl.sum("oe").alias("oe_sum"),
            pl.mean("oe").alias("oe_avg")
        ).rename({"start1": "start"})).sort("start"), on = 'start', how="inner").with_columns(
            # Add missing columns
            pl.lit(chrom).alias("chr"),
            (pl.col("start") + compartments_resolution).alias("end"),
        ).join(
            # Append the c-score segregation scores
            scores.rename({"start": "start"}), on = 'start', how="inner").select(
                # re-order to match bedformat
            ["chr", "start", "end", "compartment_id", "compartment", "tile", "domain_rank", "oe_sim_corr", "n_interactions"]
            ).join(
                # Append norm. segregation scores
            comp_scores_norm_pl.rename({"start1": "start",
                                        "compartment1": "compartment"}),
            on = ['start', 'compartment'], how="inner")
        
        chrom_compartment_scores = chrom_compartment_scores.to_pandas()
        all_chrom_compartment_scores.append(chrom_compartment_scores)
        
        # Plot
        SEG_PLT_PATH = os.path.join(PLT_PATH, "1_chrom_compartment_segregation")
        os.makedirs(SEG_PLT_PATH, exist_ok=True)
        saddle_plt.plt_comp_segreg_distribution(chrom_compartment_scores,
                                                chrom,
                                                output_path=SEG_PLT_PATH)

        #################### Tile-level analysis #######################a
        # -------------------
        # The Hi-C pixels are upper-triangular, which means that when aggregating counts at the level of 
        # tiles (tile1 - tile2 - count) we have to duplicate the interactions also for the lower triangular
        # part. This can be skipped by just ordering the tile numbers in the dataframe such that 
        # there is no difference between a pair and its reversed version.
        # Example: 
        #       tile1 = 7, tile2 = 3, count = 100
        #       tile1 = 2, tile2 = 5, count = 50
        #       tile1 = 3, tile2 = 7, count = 200
        # can be converted in the following
        #       tile1 = 3, tile2 = 7, count = 100
        #       tile1 = 2, tile2 = 5, count = 50
        #       tile1 = 3, tile2 = 7, count = 200
        # Plot distribution of segregation scores
        logger.info("- Tile-level analysis")
        
        DIST_PLT_PATH = os.path.join(PLT_PATH, "2_chrom_16-tiles_segreg_distribution")
        os.makedirs(DIST_PLT_PATH, exist_ok=True)
        saddle_plt.plt_tile_segreg_distribution(chrom_compartment_scores, chrom, output_path=DIST_PLT_PATH)
        
        tile_level = chrom_pixels_with_comp.group_by(['tile1','tile2']).agg(
            pl.sum("oe").alias("sum"),
            pl.len().alias("count")
        ).with_columns(
            tile1 = pl.when(pl.col("tile1") < pl.col("tile2")).then(pl.col("tile1")).otherwise(pl.col("tile2")),
            tile2 = pl.when(pl.col("tile1") < pl.col("tile2")).then(pl.col("tile2")).otherwise(pl.col("tile1"))
        ).group_by(['tile1','tile2']).agg(
            pl.sum("sum").alias("sum"),
            pl.sum("count").alias("count"),
            pl.lit(chrom).alias("chr")
        ).with_columns(
            (pl.col('sum') / pl.col('count')).alias("avg_oe")).select(
                ["chr", "tile1", "tile2", "avg_oe"]
            ).sort(by=["tile1", "tile2"])
        assert tile_level.shape[0] == n_tiles*(n_tiles - 1)/2 + n_tiles

        tile_level = tile_level.to_pandas()
        all_tile_level.append(tile_level)
    
        # Plot
        TILE_PLT_PATH = os.path.join(PLT_PATH, "3_tile_analysis")
        os.makedirs(TILE_PLT_PATH, exist_ok=True)
        saddle_plt.plt_tile_analysis(tile_level,
                                     chrom_comps.to_pandas(),
                                     chrom,
                                     output_path=TILE_PLT_PATH)

        #################### Subcompartment-level analysis #######################
        logger.info("- Sub-compartment-level analysis")
    
        comp_level = chrom_pixels_with_comp.group_by(['compartment1', 'compartment2']).agg(
            pl.sum("oe").alias("sum"),
            pl.len().alias("count")
        ).with_columns(
            pl.when(pl.col("compartment1") < pl.col("compartment2")).then(pl.col("compartment1")).otherwise(pl.col("compartment2")).alias("compartment1"),
            pl.when(pl.col("compartment1") < pl.col("compartment2")).then(pl.col("compartment2")).otherwise(pl.col("compartment1")).alias("compartment2")
        ).group_by(['compartment1','compartment2']).agg(
            pl.sum("sum").alias("sum"),
            pl.sum("count").alias("count"),
            pl.lit(chrom).alias("chr")
        ).with_columns(
            (pl.col('sum') / pl.col('count')).alias("avg_oe")
        ).select(
            ["chr", "compartment1", "compartment2", "avg_oe"]
        ).sort(by=["compartment1", "compartment2"]
               )
        
        n_comp_levels = len(set(comp_level['compartment1']).intersection(set(comp_level['compartment2'])))
        assert comp_level.shape[0] == n_comp_levels*(n_comp_levels - 1)/2 + n_comp_levels
        # convert to pandas
        comp_level = comp_level.to_pandas()
        all_comp_level.append(comp_level)
        # Plot
        SUBCOMP_PLT_PATH = os.path.join(PLT_PATH, "4_chrom_comp_saddle")
        os.makedirs(SUBCOMP_PLT_PATH, exist_ok=True)
        saddle_plt.plt_chrom_comp_saddle(comp_level, chrom, output_path=SUBCOMP_PLT_PATH)
        RANK_PLT_PATH = os.path.join(PLT_PATH, "5_chrom_comp_rank_distribution")
        os.makedirs(RANK_PLT_PATH, exist_ok=True)
        saddle_plt.plt_comp_rank_distribution(chrom_comps.to_pandas(), chrom, output_path=RANK_PLT_PATH)
        
        #################### Compartment domain rank level analysis #######################
        logger.info("- Domain-level analysis")
        domain_level = chrom_pixels_with_comp.group_by(['compartment_id1', 'compartment_id2',
                                                        'domain_rank1', 'domain_rank2']).agg(
                                                            pl.sum("oe").alias("sum"),
                                                            pl.len().alias("count")
                                                        ).with_columns(
                                                            pl.when(pl.col("domain_rank1") <= pl.col("domain_rank2")).then(pl.col("compartment_id1")).otherwise(pl.col("compartment_id2")).alias("compartment_id1"),
                                                            pl.when(pl.col("domain_rank1") <= pl.col("domain_rank2")).then(pl.col("compartment_id2")).otherwise(pl.col("compartment_id1")).alias("compartment_id2"),
                                                            pl.when(pl.col("domain_rank1") <= pl.col("domain_rank2")).then(pl.col("domain_rank1")).otherwise(pl.col("domain_rank2")).alias("domain_rank1"),
                                                            pl.when(pl.col("domain_rank1") <= pl.col("domain_rank2")).then(pl.col("domain_rank2")).otherwise(pl.col("domain_rank1")).alias("domain_rank2")
                                                        ).group_by(['compartment_id1', 'compartment_id2', 'domain_rank1','domain_rank2']).agg(
                                                            pl.sum("sum").alias("sum"),
                                                            pl.sum("count").alias("count"),
                                                            pl.lit(chrom).alias("chr")
                                                        ).with_columns(
                                                            (pl.col('sum') / pl.col('count')).alias("avg_oe")
                                                        ).select(
                                                            ["chr", "compartment_id1", "domain_rank1", "compartment_id2", "domain_rank2", "avg_oe"]
                                                        ).sort(by=["compartment_id1", "compartment_id2"])
        
        n_domains = len(chrom_comps['domain_rank'].unique())
        assert domain_level.shape[0] == n_domains*(n_domains - 1)/2 + n_domains
        
        domain_level = domain_level.to_pandas()
        all_domain_level.append(domain_level)
    
        # Plot
        DOM_PATH = os.path.join(PLT_PATH, "6_domain_analysis")
        os.makedirs(DOM_PATH, exist_ok=True)
        saddle_plt.plt_chrom_domainRank_saddle(domain_level, chrom, output_path=DOM_PATH)
    
    # Save
    logger.info("Saving data")
    all_tile_level = pd.concat(all_tile_level, axis=0, ignore_index=True)
    all_tile_level.to_csv(all_tile_level_path, sep="\t", index=False, header=True)
    all_comp_level = pd.concat(all_comp_level, axis=0, ignore_index=True)
    all_comp_level.to_csv(all_comp_level_path, sep="\t", index=False, header=True)
    all_domain_level = pd.concat(all_domain_level, axis=0, ignore_index=True)
    all_domain_level.to_csv(all_domain_level_path, sep="\t", index=False, header=True)
    all_chrom_compartment_scores = pd.concat(all_chrom_compartment_scores, axis=0, ignore_index=True)
    all_chrom_compartment_scores.to_csv(all_chrom_compartment_scores_path, sep="\t", index=False, header=True)

    chrom_sizes = hic.chromsizes.rename(index = lambda x: "chr" + str(x) if not str(x).startswith("chr") else str(x))

    # Converting compartment segregation scores into bigWig and bedGraph
    if not os.path.isfile(os.path.join(output_path, "compartment_segregation.bigWig")):
        logger.info("Converting compartment segregation to bigWig and bedGraph")
        chrom_sizes.drop(excluded_chroms, axis=0, errors = 'ignore')\
                   .to_frame('size')\
                   .reset_index()\
                   .to_csv(os.path.join(output_path, "chrom.sizes"), sep="\t", index=False, header=False)

        all_chrom_compartment_scores[['chr', 'start', 'end', 'oe_sim_corr']]\
            .dropna()\
            .assign(end = lambda x: x.apply(lambda y: min(y.end, chrom_sizes[y.chr]), axis=1))\
            .sort_values(['chr', 'start', 'end'])\
            .to_csv(os.path.join(output_path, "compartment_segregation.bedGraph"), sep="\t", index=False, header=False)
        subprocess.run("bedGraphToBigWig " + \
                       os.path.join(output_path, "compartment_segregation.bedGraph") + \
                       " " + \
                       os.path.join(output_path, "chrom.sizes") + \
                       " " + \
                       os.path.join(output_path, "compartment_segregation.bigWig"), shell=True)

    # Converting compartment ranks into bigWig and bedGraph
    if not os.path.isfile(os.path.join(output_path, "compartment_rank.bigWig")):
        logger.info("Converting compartment ranks to bigWig and bedGraph")
    
        all_chrom_compartment_scores[['chr', 'start', 'end', 'domain_rank']]\
            .dropna()\
            .assign(end = lambda x: x.apply(lambda y: min(y.end, chrom_sizes[y.chr]), axis=1))\
            .sort_values(['chr', 'start', 'end'])\
            .to_csv(os.path.join(output_path, "compartment_rank.bedGraph"), sep="\t", index=False, header=False)
        
        subprocess.run("bedGraphToBigWig " + \
                       os.path.join(output_path, "compartment_rank.bedGraph") + \
                       " " + \
                       os.path.join(output_path, "chrom.sizes") + \
                       " " + \
                       os.path.join(output_path, "compartment_rank.bigWig"), shell=True)
        
    # Converting compartment cscores into bigWig and bedGraph
    if not os.path.isfile(os.path.join(output_path, "compartment_cscore.bigWig")):
        logger.info("Converting compartment C-SCORE to bigWig and bedGraph")
        all_chrom_compartment_scores[['chr', 'start', 'end', 'cscore']]\
            .dropna()\
            .assign(end = lambda x: x.apply(lambda y: min(y.end, chrom_sizes[y.chr]), axis=1))\
            .sort_values(['chr', 'start', 'end'])\
            .to_csv(os.path.join(output_path, "compartment_cscore.bedGraph"), sep="\t", index=False, header=False)
    
        subprocess.run("bedGraphToBigWig " + \
                       os.path.join(output_path, "compartment_cscore.bedGraph") + \
                       " " + \
                       os.path.join(output_path, "chrom.sizes") + \
                       " " + \
                       os.path.join(output_path, "compartment_cscore.bigWig"), shell=True)

    # Deriving genome-wide statistcs
    # ------------------------------
    # - Distribution of compartment segregation both genomewide and chromosome level
    # - Saddle plot at compartment level genomewide
    # - Saddle plot at tile level genomewide
    # def __plot_mean(data, column=None, ypos=0.035, **kwrgs):
    #       avg = data[column].mean()
    #       plt.axvline(avg, color = 'black', linestyle='--', linewidth=1)
    #       plt.text(x = avg-0.01, y = ypos, s = f"Avg. {avg:.3f}", horizontalalignment='right')

    logger.info("Computing compartment segregation distributions")
    # Genome-wide plotsu

    logger.info("Computing genome-wide saddle plots")
    saddle_plt.plt_genomewide_saddle(all_comp_level, output_path=output_path)

    logger.info("Computing genome-wide tile saddle plots")
    saddle_plt.plt_genomewide_tile_saddle(all_tile_level, output_path=output_path)

    logger.info("Statistical validation of segregation")        
    logger.info("Computing segregation of A-A and B-B paris vs A-B pairs ")

    all_comp_2level = saddle_plt.plt_AB_segregation(all_comp_level, output_path=output_path)

    logger.info('Saving genomewide segregation scores')
    seg_scores = pl.concat([
        pl.DataFrame(all_comp_2level).group_by("compartment_pair_diff_same").agg(
            pl.col('avg_oe').median().alias("median_oe"),
            pl.lit("genomewide").alias("level")
        ).select(['level', 'compartment_pair_diff_same', 'median_oe']),
        pl.DataFrame(all_comp_2level).group_by(['chr', "compartment_pair_diff_same"]).agg(
            pl.col('avg_oe').median().alias("median_oe")
        ).rename({"chr": "level"}),
    ])
    seg_scores.filter(pl.col('compartment_pair_diff_same') == "AB").rename(
        {"median_oe": "AB"}
    ).drop("compartment_pair_diff_same").join(
        seg_scores.filter(pl.col('compartment_pair_diff_same') == "AA-BB").rename(
            {"median_oe": "AA-BB"}
        ).drop("compartment_pair_diff_same"), on="level", how="inner").with_columns(
            segregation_score = pl.col('AA-BB') / pl.col('AB')
        ).write_csv(
            os.path.join(output_path, "segregation_stats.tsv"),
            separator = "\t"
    )
    
if __name__ == '__main__':
    main()

