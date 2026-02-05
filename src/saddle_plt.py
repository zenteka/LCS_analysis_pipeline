import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statannotations.Annotator import Annotator

import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend that doesn't rely on tkinter
import matplotlib.pyplot as plt


COMPARTMENTS_8_LABELS = ["{}.{}.{}".format(l1, l2, l3) for l1 in ['A', 'B'] for l2 in [1, 2] for l3 in [1, 2]]

def plt_comp_segreg_distribution(indf, chrom, output_path=None):
    """
    Plot the distribution of compartment segregation scores for a given chromosome.
    Plot:
    -Boxplot of segregation scores for each compartment.
    -Scatter plot of segregation scores against compartment rank.
    
    Every data point is a correlation between the observed over
    expected (O/E) values and the compartment similiarity scores
    (between two bins).
    
    Parameters:
    - indf: DataFrame containing compartment segregation scores.
    - chrom: Name of the chromosome.
    - compartments_resolution: Resolution of the compartments.
    - output_path: Path to save the plots.
    """
    # Create a boxplot for compartment segregation scores
    fig = plt.figure()
    sns.boxplot(data=indf,
                x='compartment',
                y='oe_sim_corr',
                order = list(reversed(COMPARTMENTS_8_LABELS)),
                hue = 'compartment',
                hue_order = list(reversed(COMPARTMENTS_8_LABELS)),
                palette = "bwr",
                )
    plt.xlabel("Compartment")
    plt.ylabel("Segregation score")
    plt.title(chrom)
    sns.despine()
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_compartments_segregation_distribution.pdf"))
        plt.close(fig)
    else :
        plt.show()
    
    fig = plt.figure()
    sns.scatterplot(data = indf, 
                    x="domain_rank", y="oe_sim_corr",
                    hue="compartment", hue_order=list(reversed(COMPARTMENTS_8_LABELS)),
                    palette="bwr",
                    legend=False,
                    rasterized=False)
    plt.xlabel("Compartment rank")
    plt.ylabel("Segregation score")
    plt.title(chrom)
    sns.despine() 
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_domains_segregation_distribution.pdf"))
        plt.close(fig)
    else :
        plt.show()                    
    
def plt_tile_segreg_distribution(indf, chrom, output_path=None):
    """
    Plot the distribution of tile segregation scores for a given chromosome.
    Parameters:
    - indf: DataFrame containing compartment segregation scores.
    - chrom: Name of the chromosome.
    - output_path: Path to save the plots.
    """
    fig = plt.figure()
    sns.boxplot(data=indf,
                x = 'tile', 
                y='oe_sim_corr',
                hue = "tile",
                palette = 'bwr')
    plt.xlabel("Compartment percentile")
    plt.ylabel("Segregation score")
    n_tiles =  len(indf['tile'].unique())
    plt.title(chrom)
    sns.despine()
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_{n_tiles}_tiles_segregation_distribution.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_tile_analysis(tiledf, chrom_comps, chrom, output_path=None):
    """
    Plot the output of the tile analysis for a given chromosome.
    Parameters:
    - indf: DataFrame containing compartment segregation scores.
    - chrom: Name of the chromosome.
    - output_path: Path to save the plots.
    Plots:
        - Tile saddle plot of average O/E values.
        - Tile rank distribution plot.
    """
    n_tiles = len(chrom_comps['tile'].unique())
    fig = plt.figure()                                                                    
    H = tiledf.pivot(index = 'tile1', columns = 'tile2', values = "avg_oe")           
    H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T                            
    H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)          
    sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1,                  
                          cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))  
    plt.xlabel("Compartment percentile")                                                  
    plt.ylabel("Compartment percentile")                                                  
    plt.xticks(rotation=0)                                                                
    plt.yticks(rotation=0)                                                                
    plt.title(f"{chrom}")
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_{n_tiles}_tiles_saddle_plot.pdf"))
        plt.close(fig)
    else:
        plt.show()
    
    fig = plt.figure()
    sns.boxplot(data = chrom_comps,
                x='tile',
                y='domain_rank',
                hue='tile', 
                palette='bwr')
    plt.xlabel("Compartment percentile")
    plt.ylabel("Compartment rank")
    sns.despine()
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_{n_tiles}_tiles_rank_distribution.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_chrom_comp_saddle(comp_level, chrom, output_path=None):
    '''
    Plot the saddle plots for a given chromosome for 8 subcompartments.
    '''
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
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_compartments_saddle_plot.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_comp_rank_distribution(chrom_comps, chrom, output_path=None):
    '''
        Plot the compartment domain rank distribution for a given chromosome.    
    '''
    fig = plt.figure()
    sns.boxplot(data=chrom_comps,
                x='compartment', y='domain_rank',
                hue='compartment',
                order=list(reversed(COMPARTMENTS_8_LABELS)),
                hue_order=list(reversed(COMPARTMENTS_8_LABELS)),
                palette='bwr')
    plt.xlabel("Compartment")
    plt.ylabel("Compartment rank")
    sns.despine()
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_compartments_rank_distribution.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_chrom_domainRank_saddle(domain_level, chrom, output_path=None):
    '''
    Plot the saddle plots for a given chromosome for continuous domain ranks.
    '''
    fig = plt.figure()
    H = domain_level.pivot(index = 'domain_rank1', columns = 'domain_rank2', values = "avg_oe")
    H = H.sort_index(axis=0).sort_index(axis=1)
    H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
    H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)    
    sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, xticklabels=False, yticklabels=False,
                cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5), rasterized=True)
    plt.xlabel("Compartment domain\n"r"(Inactive $\rightarrow$ Active)")
    plt.ylabel("Compartment domain\n"r"(Inactive $\rightarrow$ Active)")
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.title(f"{chrom}")
    
    if output_path is not None:
        fig.savefig(os.path.join(output_path, f"{chrom}_domains_saddle_plot.pdf"))
        plt.close(fig)
    else:
        plt.show()
        
#################### Genome-wide plots ####################
def plt_genomewide_saddle(all_comp_level, output_path=None):
    '''
    Plot the genomewide saddle plots for a given chromosome for 8 subcompartments.
    '''
    GW = all_comp_level.groupby(['compartment1', 'compartment2'])['avg_oe']\
                       .mean().to_frame("avg_oe").reset_index()
    fig = plt.figure()
    H = GW.pivot(index = 'compartment1', columns = 'compartment2', values = "avg_oe")
    H = H.reindex(index=COMPARTMENTS_8_LABELS, columns=COMPARTMENTS_8_LABELS)
    H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
    H = H.sort_index(ascending=True, axis=0).sort_index(ascending=False, axis=1)
    sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1,
                cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
    plt.xlabel("Compartment")
    plt.ylabel("Compartment")
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.title("Genome-wide")
    
    if output_path is not None:
        fig.savefig(os.path.join(output_path, "plots", f"genomewide_compartments_saddle_plot.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_genomewide_tile_saddle(all_tile_level, output_path=None):
    '''
    Plot the genome-wide saddle plots for all 16 tiles.
    TODO: Convert computation to poalrs.
    '''
    n_tiles = len(all_tile_level['tile1'].unique())
    tile_GW = all_tile_level.groupby(['tile1', 'tile2'])['avg_oe']\
                                              .mean().to_frame("avg_oe").reset_index()
    
    fig = plt.figure()
    H = tile_GW.pivot(index = 'tile1', columns = 'tile2', values = "avg_oe")
    H.loc[:, :] = np.triu(H.values) + np.triu(H.values, k=1).T
    H = H.sort_index(ascending=False, axis=0).sort_index(ascending=True, axis=1)
    sns.heatmap(H, cmap='bwr', vmin=0, vmax=2, square=True, linewidth=1, 
                          cbar_kws=dict(label = "Avg. O/E counts", aspect=5, shrink=.5))
    plt.xlabel("Compartment percentile")
    plt.ylabel("Compartment percentile")
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    if output_path is not None:
        fig.savefig(os.path.join(output_path, "plots", f"genomewide_{n_tiles}_tiles_saddle_plot.pdf"))
        plt.close(fig)
    else:
        plt.show()

def plt_AB_segregation(all_comp_level,output_path=None):
    '''
    This function validates compartment segregation at level of A/B compartments.
    '''
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
                    hue = 'compartment_pair_diff_same',
                    palette = comp_2leve_palette,
                    showfliers=False,
                    boxprops={'alpha':0.3})
    ax = sns.stripplot(data = all_comp_2level, 
                       x = 'compartment_pair_diff_same', y = "avg_oe",
                       hue = 'compartment_pair_diff_same',
                       order = list(comp_2leve_palette.keys()),
                       palette = comp_2leve_palette)
    annotator = Annotator(ax, 
    pairs=[("AA-BB", "AB")],
    data=all_comp_2level,
    x='compartment_pair_diff_same', y='avg_oe',
    order=list(comp_2leve_palette.keys()))

    annotator.configure(test='Mann-Whitney', text_format='full', loc='inside')
    annotator.apply_and_annotate()
    # Isuue _BoxPlotter calss removed, stat_annot is not compatible.
    # add_stat_annotation(ax, 
    #                     data = all_comp_2level,
    #                     x = 'compartment_pair_diff_same', y = "avg_oe",
    #                     order = list(comp_2leve_palette.keys()),
    #                     box_pairs=[
    #                         ("AA-BB", "AB")
    #                     ],
    #                     test='Mann-Whitney', text_format='full',
    #                     comparisons_correction=None)
    plt.ylim(0, 5)
    plt.xlabel("Compartment interaction")
    plt.ylabel("Average O/E count")
    sns.despine(trim=True)
    if output_path is not None:
        fig.savefig(os.path.join(output_path, "plots", f"genomewide_AB_segregation_boxplot.pdf"))
        plt.close(fig)
    else:
        plt.show()
        
    return all_comp_2level
