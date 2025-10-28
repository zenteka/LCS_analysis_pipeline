import os
import subprocess
import logging
import argparse
import fanc
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
import numpy as np
from pybedtools.bedtool import BedTool

import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation


CALDER_HEADER=['chr', 'start', 'end', 'compartment_label', 'domain_rank_8_levels', 
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

def load_compartments(calder_compartments_path):
	compartments = pd.read_csv(calder_compartments_path, sep="\t", header=None, names = CALDER_HEADER)
	compartments = compartments.groupby(['chr', 'compartment_label'])\
							   .agg({'start': 'min', 'end': 'max'})\
							   .reset_index()[['chr', 'start', 'end', 'compartment_label']]\
							   .sort_values(['chr', 'start', 'end'])\
							   .reset_index(drop=True)

	assert BedTool.from_dataframe(compartments).merge(d = -1).to_dataframe().shape[0] == compartments.shape[0]

	compartments['id'] = compartments.chr + ":" + compartments.start.astype(str) + "-" + compartments.end.astype(str) + "_" + compartments.compartment_label
	compartments['compartment_label_8'] = compartments.compartment_label.map(lambda x: x[:5])

	domain_rank = compartments.groupby("chr").compartment_label.rank(ascending=False, pct=True)
	compartments['domain_rank'] = domain_rank

	min_max_domain_rank = compartments.groupby('chr')['domain_rank']\
								 .min()\
								 .to_frame("min_domain_rank")\
								 .reset_index()\
								 .merge(compartments.groupby('chr')['domain_rank'].max().to_frame("max_domain_rank").reset_index())

	compartments = compartments.merge(min_max_domain_rank)
	compartments['domain_rank'] = (compartments.domain_rank - compartments.min_domain_rank)/(compartments.max_domain_rank - compartments.min_domain_rank)
	return compartments[["chr", 'start', 'end', 'compartment_label', 'compartment_label_8', 'domain_rank', 'id']]


def get_oe_values(hic, chrom):
	M = hic.matrix((chrom, chrom), oe=True)
	row, col = np.triu_indices_from(M)
	vals = M[row, col]

	pixels = pd.DataFrame({
			"bin1": row,
			'bin2': col,
			"oe": vals.data,
			"mask": vals.mask
		})
	pixels = pixels[~pixels["mask"]].drop("mask", axis=1)
	return pixels


def main():

	logger = set_logger()

	args = process_arguments()
	hic_path=args.hic_path
	hic_resolution=args.hic_resolution
	compartments_path=args.compartments_path
	compartments_resolution=args.compartments_resolution
	excluded_chroms = args.excludeChroms.split(",") if args.excludeChroms is not None else []
	n_tiles = args.n_tiles
	output_path = args.output_path

	# hic_path="data/database/hic/CP-A/WGD/CPA_WGD_run2_mega/aligned/inter_30.mcool"
	# hic_resolution = 1000000
	# compartments_path="data/processed/hic/compartments/compartments_CALDER/IGV/CPA_WGD_run2_mega.bed"
	# compartments_resolution=50000
	# excluded_chroms = ['chrY', 'chrM', 'chrMT', 'Y', 'M', "MT"]
	# n_tiles = 16
	# output_path = "data/processed/hic/compartments/calder_saddle/CPA_WGD_run2_mega"

	logger.info("Compartment segregation analysis")
	logger.info("Parameters:")
	logger.info(f"- Hi-C: {hic_path}")
	logger.info(f"- Hi-C resolution: {hic_resolution}")
	logger.info(f"- Compartments: {compartments_path}")
	logger.info(f"- Compartments resolution: {compartments_resolution}")
	logger.info(f"- Chromosomes to exclude: " + ", ".join(excluded_chroms))
	logger.info(f"- N. tiles: {n_tiles}")
	logger.info(f"- Output path: {output_path}")

	# Creating the output path
	os.makedirs(output_path, exist_ok=True)
	os.makedirs(os.path.join(output_path, "plots"), exist_ok=True)

	logger.info("Loading compartments")
	compartments = load_compartments(compartments_path)

	# Add 'chr' to chromosome names if necessary
	compartments['chr'] = compartments.chr.map(lambda x: "chr" + str(x) if not str(x).startswith("chr") else str(x))


	# Binning compartments
	logger.info("Binning compartments")
	compartments_bed = BedTool.from_dataframe(compartments).sort()
	bins = BedTool().window_maker(b = compartments_bed, w = compartments_resolution)
	compartments_binned = bins.map(compartments_bed, c = [5, 6, 7], o = ["distinct", 'mean', 'distinct'])\
							  .to_dataframe(names = ['chr', 'start', 'end', "compartment", 'domain_rank', 'compartment_id'])

	# loading Hi-C
	logger.info("Loading Hi-C")
	hic = fanc.load(hic_path + f"@{hic_resolution}")

	all_tile_level_path = os.path.join(output_path, f"tiles_{n_tiles}.tsv")
	all_comp_level_path = os.path.join(output_path, f"compartments.tsv")
	all_domain_level_path = os.path.join(output_path, f"domains.tsv")
	all_chrom_compartment_scores_path = os.path.join(output_path, f"compartment_scores.tsv")
	if not all([os.path.isfile(all_tile_level_path), 
				os.path.isfile(all_comp_level_path), 
				os.path.isfile(all_domain_level_path),
				os.path.isfile(all_chrom_compartment_scores_path)]):
		all_tile_level = []
		all_comp_level = []
		all_domain_level = []
		all_chrom_compartment_scores = []
		# We iterate for each chromosome (which is not excluded)
		for chrom in filter(lambda x: x not in excluded_chroms, hic.chromnames):
			logger.info(f"Analysing {chrom}")

			# Get the OE values for the chromosome
			logger.info("- Extracting O/E values")
			chrom_pixels = get_oe_values(hic, chrom)
			# Get the binned compartments for the chromosome
			# - Extract the Hi-C bin (in the resolution of the Hi-C data) from the start of the compartment bin
			# - Assign a tile to each bin, based on its domain rank, using percentiles
			logger.info("- Extracting compartments for the chromosome")
			chrom = "chr" + str(chrom) if not str(chrom).startswith("chr") else chrom
			chrom_comps = compartments_binned[compartments_binned.chr == chrom]\
								.assign(bin = lambda x: x.start//hic_resolution,
										tile = lambda x: pd.cut(x.domain_rank, bins = n_tiles, labels = list(range(n_tiles))))\
								[['bin', 'tile', 'start', "compartment", 'domain_rank', 'compartment_id']]
			
			# Associate to each pixel in the Hi-C matrix the compartment information
			# of its two bins
			logger.info("- Associating compartments to Hi-C interactions")
			chrom_pixels_with_comp = chrom_pixels\
						.merge(chrom_comps.rename(columns = lambda x: x + "1"), on = 'bin1')\
						.merge(chrom_comps.rename(columns = lambda x: x + "2"), on = 'bin2')\
						.drop_duplicates()
			chrom_pixels_with_comp = chrom_pixels_with_comp[(chrom_pixels_with_comp.bin1<=chrom_pixels_with_comp.bin2) & \
															(chrom_pixels_with_comp.start1 <= chrom_pixels_with_comp.start2)].copy()
			assert chrom_comps.shape[0]*(chrom_comps.shape[0] - 1)/2 + chrom_comps.shape[0] == chrom_pixels_with_comp.shape[0]

			# Computing local compartmentalization scores
			# -------------------------------------------
			# Adding the lower triangular matrix to the pixels
			logger.info("- Commputing segregation scores")
			local_scores = chrom_pixels_with_comp\
									.assign(compSim = lambda x: 1 - (x.domain_rank2 - x.domain_rank1).abs(),
										    weight = lambda x: x.compSim * x.oe)\
									[['start1', 'start2', 
									  'oe', 
									  'domain_rank1', 'domain_rank2', 
									  'compartment1', 'compartment2',
									  'compSim', 'weight']]
			lower_local_scores = local_scores[local_scores.start1 < local_scores.start2].copy()
			start1 = lower_local_scores.start1.copy().values
			start2 = lower_local_scores.start2.copy().values
			dom1 = lower_local_scores.domain_rank1.copy().values
			dom2 = lower_local_scores.domain_rank2.copy().values
			comp1 = lower_local_scores.compartment1.copy().values
			comp2 = lower_local_scores.compartment2.copy().values
			lower_local_scores['start1'] = start2
			lower_local_scores['start2'] = start1
			lower_local_scores['domain_rank1'] = dom2
			lower_local_scores['domain_rank2'] = dom1
			lower_local_scores['compartment1'] = comp2
			lower_local_scores['compartment2'] = comp1

			both_local_scores = pd.concat([local_scores, lower_local_scores], axis=0, ignore_index=True)

			# Computing the correlation between O/E values and compartment similarity
			scores = both_local_scores.groupby('start1')\
									  .apply(lambda x: spearmanr(x.compSim.values, x.oe.values)[0])\
									  .to_frame("oe_sim_corr")\
									  .reset_index()\
									  .rename(columns = {'start1': 'start'})

			# Computing compartment scores at the 8-levels
			# For each bin, 
			# - We sum all the O/E interactions between the bin and all the 8 compartments
			# - We divide those values by the total number of interactions, thus obtaining a percentage
			# - We divide those values by the expected percentage of interactions for each compartment level,
			#	which is obtained as the percentage of bins covered by each level in that chromosome
			# - We divide those values by the total number of interactions, thus obtaining a percentage
			comp_scores = both_local_scores.groupby(['start1', 'compartment1', 'compartment2'])\
											['oe']\
											.sum()\
											.to_frame("avg_oe")\
											.reset_index()\
											.pivot(index = ['start1', 'compartment1'], columns = 'compartment2', values = 'avg_oe')\
											.reindex(columns = COMPARTMENTS_8_LABELS)

			comp_scores = pd.DataFrame(comp_scores.values / comp_scores.sum(1).values[:, np.newaxis],
									   index = comp_scores.index, columns = comp_scores.columns)

			comp_expected_scores = chrom_comps.compartment.value_counts(normalize=True).reindex(index = COMPARTMENTS_8_LABELS)

			comp_scores_norm = pd.DataFrame(comp_scores.values / comp_expected_scores.values, 
											index = comp_scores.index, columns = comp_scores.columns)
			comp_scores_norm = pd.DataFrame(comp_scores_norm.values / comp_scores_norm.sum(1).values[:, np.newaxis],
											index = comp_scores_norm.index, columns = comp_scores_norm.columns)

			comp_scores['cscore'] = comp_scores.apply(lambda x: x[x.name[1]], axis=1)
			comp_scores_norm['cscore'] = comp_scores_norm.apply(lambda x: x[x.name[1]], axis=1)


			# Adding coverage info
			coverage = both_local_scores.groupby('start1')\
											['oe']\
											.agg(['count', 'sum', 'mean'])\
											.reset_index()\
											.rename(columns = {'start1': 'start',
															   'count': 'n_interactions',
															   'sum': 'oe_sum',
															   'mean': 'oe_avg'})
			# Compiling all compartment segration data
			chrom_compartment_scores = chrom_comps.merge(coverage, on = 'start').merge(scores, on = 'start')
			chrom_compartment_scores['chr'] = chrom
			chrom_compartment_scores['end'] = chrom_compartment_scores.start + compartments_resolution
			chrom_compartment_scores = chrom_compartment_scores[['chr', 'start', 'end', 
																 'compartment_id', 'compartment', 
																 'tile', 'domain_rank',
																 'oe_sim_corr', 'n_interactions']]

			chrom_compartment_scores = chrom_compartment_scores.merge(comp_scores_norm.reset_index()\
																.rename(columns = {'start1': 'start', 'compartment1': 'compartment'}), 
																		on = ['start', 'compartment'])
			all_chrom_compartment_scores.append(chrom_compartment_scores)


			# fig = plt.figure()
			# sns.boxplot(data=chrom_compartment_scores, x='compartment', order = reversed(COMPARTMENTS_8_LABELS), 
			# 			y='oe_sim_corr')
			# plt.xlabel("Compartment")
			# plt.ylabel("Segregation score")
			# plt.title(chrom)
			# sns.despine()
			# fig.savefig(os.path.join(output_path, "plots", f"{chrom}_compartments_segregation_distribution.pdf"))
			# plt.close(fig)

			# fig = plt.figure()
			# sns.boxplot(data=chrom_compartment_scores, x = 'tile', 
			# 			y='oe_sim_corr', palette = 'bwr')
			# plt.xlabel("Compartment percentile")
			# plt.ylabel("Segregation score")
			# plt.title(chrom)
			# sns.despine()
			# fig.savefig(os.path.join(output_path, "plots", f"{chrom}_{n_tiles}_tiles_segregation_distribution.pdf"))
			# plt.close(fig)

			# fig = plt.figure()
			# sns.scatterplot(data = chrom_compartment_scores, 
			# 				x="domain_rank", y = "oe_sim_corr", hue = "compartment",
			# 				hue_order=reversed(COMPARTMENTS_8_LABELS), palette = "bwr",
			# 				rasterized=False)
			# plt.xlabel("Compartment rank")
			# plt.ylabel("Segregation score")
			# plt.title(chrom)
			# sns.despine()
			# plt.legend(bbox_to_anchor=(1, 1))
			# fig.savefig(os.path.join(output_path, "plots", f"{chrom}_domains_segregation_distribution.pdf"))
			# plt.close(fig)

			# Tile-level analysis
			# -------------------
			# The Hi-C pixels are upper-triangular, which means that when aggregating counts at the level of 
			# tiles (tile1 - tile2 - count) we have to duplicate the interactions also for the lower triangular
			# part. This can be skipped by just ordering the tile numbers in the dataframe such that 
			# there is no difference between a pair and its reversed version.
			# Example: 
			#	tile1 = 7, tile2 = 3, count = 100
			#	tile1 = 2, tile2 = 5, count = 50
			#	tile1 = 3, tile2 = 7, count = 200
			# can be converted in the following
			#	tile1 = 3, tile2 = 7, count = 100
			#	tile1 = 2, tile2 = 5, count = 50
			#	tile1 = 3, tile2 = 7, count = 200
			logger.info("- Tile-level analysis")
			tile_level = chrom_pixels_with_comp.groupby(['tile1','tile2'])["oe"].agg(['sum', 'count']).reset_index()
			tile1 = tile_level.apply(lambda x: int(min([x.tile1, x.tile2])), axis=1)
			tile2 = tile_level.apply(lambda x: int(max([x.tile1, x.tile2])), axis=1)
			tile_level['tile1'] = tile1.values
			tile_level['tile2'] = tile2.values
			tile_level = tile_level.groupby(['tile1','tile2'])[['sum', 'count']].sum().reset_index()
			tile_level['avg_oe'] = tile_level['sum'] / tile_level["count"]
			tile_level['chr'] = chrom
			tile_level = tile_level[['chr', 'tile1', 'tile2', 'avg_oe']]
			assert tile_level.shape[0] == n_tiles*(n_tiles - 1)/2 + n_tiles

			all_tile_level.append(tile_level)

			fig = plt.figure()
			H = tile_level.pivot(index = 'tile1', columns = 'tile2', values = "avg_oe")
			H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
			H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)
			sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1, 
						cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
			plt.xlabel("Compartment percentile")
			plt.ylabel("Compartment percentile")
			plt.xticks(rotation=0)
			plt.yticks(rotation=0)
			plt.title(f"{chrom}")
			fig.savefig(os.path.join(output_path, "plots", f"{chrom}_{n_tiles}_tiles_saddle.pdf"))
			plt.close(fig)


			fig = plt.figure()
			sns.boxplot(data = chrom_comps, x='tile', y = 'domain_rank', palette = 'bwr')
			plt.xlabel("Compartment percentile")
			plt.ylabel("Compartment rank")
			sns.despine()
			fig.savefig(os.path.join(output_path, "plots", f"{chrom}_{n_tiles}_tiles_rank_distribution.pdf"))
			plt.close(fig)

			# Compartment-level analysis
			# --------------------------
			# Same logic as before, but applied to 8-level compartments 
			logger.info("- Compartment-level analysis")
			comp_level = chrom_pixels_with_comp.groupby(['compartment1', 'compartment2'])['oe'].agg(['sum', 'count']).reset_index()
			comp1 = comp_level.apply(lambda x: sorted([x.compartment1, x.compartment2])[0], axis=1)
			comp2 = comp_level.apply(lambda x: sorted([x.compartment1, x.compartment2])[1], axis=1)
			comp_level['compartment1'] = comp1.values
			comp_level['compartment2'] = comp2.values
			comp_level = comp_level.groupby(['compartment1','compartment2'])[['sum', 'count']].sum().reset_index()
			comp_level['avg_oe'] = comp_level['sum'] / comp_level["count"]
			comp_level['chr'] = chrom
			comp_level = comp_level[['chr', 'compartment1', 'compartment2', 'avg_oe']]
			n_comp_levels = len(set(comp_level.compartment1).intersection(set(comp_level.compartment2)))
			assert comp_level.shape[0] == n_comp_levels*(n_comp_levels - 1)/2 + n_comp_levels

			all_comp_level.append(comp_level)

			fig = plt.figure()
			H = comp_level.pivot(index = 'compartment1', columns = 'compartment2', values = "avg_oe")
			H = H.reindex(index=COMPARTMENTS_8_LABELS, columns=COMPARTMENTS_8_LABELS)
			H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
			H = H.sort_index(ascending=True, axis=0).sort_index(ascending=False, axis=1)
			sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1, 
						cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
			plt.xlabel("Compartment")
			plt.ylabel("Compartment")
			plt.xticks(rotation=0)
			plt.yticks(rotation=0)
			plt.title(f"{chrom}")
			fig.savefig(os.path.join(output_path, "plots", f"{chrom}_compartments_saddle.pdf"))
			plt.close(fig)


			# fig = plt.figure()
			# sns.boxplot(data = chrom_comps, x='compartment', y = 'domain_rank', order = reversed(COMPARTMENTS_8_LABELS))
			# plt.xlabel("Compartment")
			# plt.ylabel("Compartment rank")
			# sns.despine()
			# fig.savefig(os.path.join(output_path, "plots", f"{chrom}_compartments_rank_distribution.pdf"))
			# plt.close(fig)

			# Compartment domain analysis
			# ---------------------------
			# Same logic as before, but applied to each compartment domain separately
			logger.info("- Domain-level analysis")
			domain_level = chrom_pixels_with_comp.groupby(['compartment_id1', 'compartment_id2', 
														   'domain_rank1', 'domain_rank2'])['oe'].agg(['sum', 'count']).reset_index()
			def __align(x):
				if x.domain_rank1 <= x.domain_rank2:
					return pd.Series([x.compartment_id1, x.compartment_id2, x.domain_rank1, x.domain_rank2])
				else:
					return pd.Series([x.compartment_id2, x.compartment_id1, x.domain_rank2, x.domain_rank1])
			aligned_domains = domain_level.apply(__align, axis=1)

			domain_level.loc[:, ['compartment_id1', 'compartment_id2', 'domain_rank1', 'domain_rank2']] = aligned_domains.values
			domain_level = domain_level.groupby(['compartment_id1', 'compartment_id2', 'domain_rank1','domain_rank2'])\
												[['sum', 'count']].sum().reset_index()
			domain_level['avg_oe'] = domain_level['sum'] / domain_level["count"]
			domain_level['chr'] = chrom
			domain_level = domain_level[['chr', 'compartment_id1', 'domain_rank1', 'compartment_id2', 'domain_rank2', 'avg_oe']]
			n_domains = chrom_comps.domain_rank.nunique()
			assert domain_level.shape[0] == n_domains*(n_domains - 1)/2 + n_domains

			all_domain_level.append(domain_level)

			fig = plt.figure()
			H = domain_level.pivot(index = 'domain_rank1', columns = 'domain_rank2', values = "avg_oe")
			H = H.sort_index(axis=0).sort_index(axis=1)
			H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
			H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)
			sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, xticklabels=False, yticklabels=False,
						cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5), rasterized=True)
			plt.xlabel("Compartment domain\n"r"(Inactive $\rightarrow$ Active)")
			plt.ylabel("Compartment domain\n"r"(Inactive $\rightarrow$ Active)")
			plt.title(f"{chrom}")
			fig.savefig(os.path.join(output_path, "plots", f"{chrom}_domains_saddle.pdf"))
			plt.close(fig)

		logger.info("Saving data")
		all_tile_level = pd.concat(all_tile_level, axis=0, ignore_index=True)
		all_tile_level.to_csv(all_tile_level_path, sep="\t", index=False, header=True)
		all_comp_level = pd.concat(all_comp_level, axis=0, ignore_index=True)
		all_comp_level.to_csv(all_comp_level_path, sep="\t", index=False, header=True)
		all_domain_level = pd.concat(all_domain_level, axis=0, ignore_index=True)
		all_domain_level.to_csv(all_domain_level_path, sep="\t", index=False, header=True)
		all_chrom_compartment_scores = pd.concat(all_chrom_compartment_scores, axis=0, ignore_index=True)
		all_chrom_compartment_scores.to_csv(all_chrom_compartment_scores_path, sep="\t", index=False, header=True)

	else:
		logger.info("Data already computed. Loading...")
		all_tile_level = pd.read_csv(all_tile_level_path, sep="\t")
		all_comp_level = pd.read_csv(all_comp_level_path, sep="\t")
		all_domain_level = pd.read_csv(all_domain_level_path, sep="\t")
		all_chrom_compartment_scores = pd.read_csv(all_chrom_compartment_scores_path, sep="\t")

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
	def __plot_mean(data, column=None, ypos=0.035, **kwrgs):
		avg = data[column].mean()
		plt.axvline(avg, color = 'black', linestyle='--', linewidth=1)
		plt.text(x = avg-0.01, y = ypos, s = f"Avg. {avg:.3f}", horizontalalignment='right')


	logger.info("Computing compartment segregation distributions")
	if not os.path.isfile(os.path.join(output_path, "plots", "genomewide_compartment_segregation_distribution.pdf")):
		g = sns.displot(data = all_chrom_compartment_scores, x = "oe_sim_corr", kind='hist', stat='probability')\
				.set(xlim=(-1, 1))\
				.set_xlabels("Compartment segregation")\
				.set_ylabels("Frequency")\
				.map_dataframe(__plot_mean, column = "oe_sim_corr")
		g.savefig(os.path.join(output_path, "plots", "genomewide_compartment_segregation_distribution.pdf"))
		plt.close(g.fig)

	if not os.path.isfile(os.path.join(output_path, "plots", "genomewide_compartment_segregation_distribution_by_chrom.pdf")):
		g = sns.displot(data = all_chrom_compartment_scores, 
						col = "chr", 
						col_order = list(filter(lambda x: x not in excluded_chroms, chrom_sizes.index.tolist())), 
						col_wrap=8,
						x = "oe_sim_corr", 
						kind='hist', 
						stat='probability',
						common_norm=False)\
				.set(xlim=(-1, 1))\
				.set_xlabels("Compartment segregation")\
				.set_ylabels("Frequency")\
				.map_dataframe(__plot_mean, column = "oe_sim_corr", ypos=0.11)\
				.set_titles(row_template="{row_name}", col_template="{col_name}",size = 20 )
		g.savefig(os.path.join(output_path, "plots", "genomewide_compartment_segregation_distribution_by_chrom.pdf"))
		plt.close(g.fig)
	logger.info("Computing genome-wide saddle plots")
	genomewide_comp_level = all_comp_level.groupby(['compartment1', 'compartment2'])['avg_oe']\
											.mean().to_frame("avg_oe").reset_index()
	
	fig = plt.figure()
	H = genomewide_comp_level.pivot(index = 'compartment1', columns = 'compartment2', values = "avg_oe")
	H = H.reindex(index=COMPARTMENTS_8_LABELS, columns=COMPARTMENTS_8_LABELS)
	H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
	H = H.sort_index(ascending=True, axis=0).sort_index(ascending=False, axis=1)
	sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1, 
				cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
	plt.xlabel("Compartment")
	plt.ylabel("Compartment")
	plt.xticks(rotation=0)
	plt.yticks(rotation=0)
	fig.savefig(os.path.join(output_path, "plots", f"genomewide_compartments_saddle.pdf"))
	plt.close(fig)


	genomewide_tile_level = all_tile_level.groupby(['tile1', 'tile2'])['avg_oe']\
											.mean().to_frame("avg_oe").reset_index()

	# fig = plt.figure()
	# H = genomewide_tile_level.pivot(index = 'tile1', columns = 'tile2', values = "avg_oe")
	# H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
	# H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)
	# sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1, 
	# 			cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
	# plt.xlabel("Compartment percentile")
	# plt.ylabel("Compartment percentile")
	# plt.xticks(rotation=0)
	# plt.yticks(rotation=0)
	# fig.savefig(os.path.join(output_path, "plots", f"genomewide_{n_tiles}_tiles_saddle.pdf"))
	# plt.close(fig)


	# Aggregating compartments at 2-level hierarchy (A/B)
	logger.info("Statistical validation of segregation")
	all_comp_2level = all_comp_level.assign(compartment1 = lambda x: x.compartment1.map(lambda y: y.split(".")[0]),
											compartment2 = lambda x: x.compartment2.map(lambda y: y.split(".")[0]),
											compartment_pair = lambda x: x.compartment1 + x.compartment2,
											compartment_pair_diff_same = lambda x: x.compartment_pair.map(lambda y: y if y == "AB" else "AA-BB"))

	comp_2leve_palette = { 
		"AA-BB": "lightblue",
		"AB": "salmon"
	}


	fig = plt.figure(figsize=(3, 5))
	sns.boxplot(data = all_comp_2level, 
				x = 'compartment_pair_diff_same', y = "avg_oe",
				order = list(comp_2leve_palette.keys()),
				palette = comp_2leve_palette,
				showfliers=False,
				boxprops={'alpha':0.3})
	ax = sns.stripplot(data = all_comp_2level, 
					x = 'compartment_pair_diff_same', y = "avg_oe",
					order = list(comp_2leve_palette.keys()),
					palette = comp_2leve_palette)
	add_stat_annotation(ax, 
						data = all_comp_2level,
			    x = 'compartment_pair_diff_same', y = "avg_oe",
						order = list(comp_2leve_palette.keys()),
						box_pairs=[
							("AA-BB", "AB")
						],
						test='Mann-Whitney', text_format='full',
						comparisons_correction=None)
	plt.ylim(0, 5)
	plt.xlabel("Compartment interaction")
	plt.ylabel("Average O/E count")
	sns.despine(trim=True)
	fig.savefig(os.path.join(output_path, "plots", f"genomewide_AB_segregation_boxplot.pdf"))
	plt.close(fig)

	logger.info('Saving genomewide segregation scores')
	pd.concat([
		all_comp_2level.groupby("compartment_pair_diff_same")\
					['avg_oe']\
					.median()\
					.to_frame("median_oe")\
					.T\
					.reset_index(drop=True)\
					.assign(segregation_score = lambda x: x['AA-BB'] / x['AB'],
							level = "genomewide")\
					[['level', 'AA-BB', 'AB', 'segregation_score']],
		all_comp_2level.groupby(['chr', "compartment_pair_diff_same"])\
					['avg_oe']\
					.median()\
					.to_frame("median_oe")\
					.reset_index()\
					.pivot(index = 'chr', columns = 'compartment_pair_diff_same', values = "median_oe")\
					.reset_index()\
					.assign(segregation_score = lambda x: x['AA-BB'] / x['AB'])\
					.rename(columns = {'chr': 'level'})
	], axis=0, ignore_index=True)\
					.to_csv(os.path.join(output_path, "segregation_stats.tsv"), sep="\t", index=False, header=True)


if __name__ == '__main__':
	main()


