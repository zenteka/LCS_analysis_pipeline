'''
Helper functions for the calder saddle tool.
'''

import pandas as pd
import numpy as np
import fanc
from pybedtools import BedTool


def load_compartments(calder_compartments_path, CALDER_HEADER):
    '''
        Loads the compartments from a CALDER bed file and processes them.
        Args:
                calder_compartments_path (str): Path to the CALDER bed file.
                CALDER_HEADER (list): List of column names for the compartments data.
    '''
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
    '''
    Extracts the observed/expected (O/E) values for a given chromosome from a Hi-C matrix.
    '''
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
