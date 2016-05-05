import pandas as pd
import os
from collections import Counter
import utils

# logging
import logging
logger = logging.getLogger(__name__)

def process_cgc(path):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    s = df['Mutation Types']
    is_small = s.str.contains('Mis|F|N|S').fillna(False)
    is_somatic = ~df['Tumour Types(Somatic)'].isnull()
    df = df[is_small & is_somatic].copy()

    # get gene names
    cgc_genes = df['Gene Symbol'].tolist()

    return cgc_genes


def num_cgc_overlap(signif_genes, cgc_list):
    cgc_result = {}
    for method in signif_genes:
        intersect = len(set(signif_genes[method]) & set(cgc_list))
        cgc_result[method] = intersect
    return cgc_result


def main(opts):
    logger.info('Running cgc_overlap sub-command . . .')
    config = utils.load_config(opts['config'])

    # get the significant genes for each method
    signif_dict = utils.fetch_significant_genes(opts['input_dir'],
                                                opts['qvalue'],
                                                config)
    num_methods = len(signif_dict)

    # get cgc genes
    cgc_genes = process_cgc(opts['cgc'])

    # count the overlap
    num_cgc_dict = num_cgc_overlap(signif_dict, cgc_genes)
    num_signif_dict = {k: len(signif_dict[k]) for k in signif_dict}

    # format result
    overlap_df = pd.DataFrame({'# CGC': pd.Series(num_cgc_dict),
                               '# significant': pd.Series(num_signif_dict)})
    overlap_df['Fraction overlap (CGC)'] = overlap_df['# CGC'].astype(float) / overlap_df['# significant']

    # save result
    cgc_path = os.path.join(opts['output'], 'cgc_overlap.txt')
    overlap_df.to_csv(cgc_path, sep='\t')

    if opts['plot']:
        import plot_data
        cgc_path = os.path.join(opts['output'], 'cgc_overlap.pdf')
        plot_data.cgc_overlap(overlap_df, cgc_path)

    logger.info('Finished cgc_overlap sub-command.')
