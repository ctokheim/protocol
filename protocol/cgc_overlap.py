import pandas as pd
import os
from collections import Counter
import utils

# logging
import logging
logger = logging.getLogger(__name__)

def num_cgc_overlap(signif_genes, cgc_list):
    """Intersect significant genes with CGC or other driver gene list."""
    cgc_result = {}
    for method in signif_genes:
        intersect = len(set(signif_genes[method]) & set(cgc_list))
        cgc_result[method] = intersect
    return cgc_result


def read_custom_list(path):
    """Read in a custom driver gene list."""
    with open(path) as handle:
        gene_list = [l.strip() for l in handle]
    return gene_list


def main(opts):
    logger.info('Running list_overlap sub-command . . .')
    config = utils.load_config(opts['config'])

    # get the significant genes for each method
    signif_dict = utils.fetch_significant_genes(opts['input_dir'],
                                                opts['qvalue'],
                                                config)
    num_methods = len(signif_dict)

    # read in driver genes
    if opts['cgc']:
        # get cgc genes
        driver_genes = utils.process_cgc(opts['cgc'])
        list_name = 'CGC'
    else:
        # read in custom driver gene list
        driver_genes = read_custom_list(opts['gene_list'])
        list_name = 'Custom List'

    # count the overlap
    num_cgc_dict = num_cgc_overlap(signif_dict, driver_genes)
    num_signif_dict = {k: len(signif_dict[k]) for k in signif_dict}

    # format result
    overlap_df = pd.DataFrame({'# '+list_name: pd.Series(num_cgc_dict),
                               '# significant': pd.Series(num_signif_dict)})
    overlap_df['Fraction overlap ({0})'.format(list_name)] = overlap_df['# '+list_name].astype(float) / overlap_df['# significant']

    # save result
    cgc_path = os.path.join(opts['output'], 'gene_list_overlap.txt')
    overlap_df.to_csv(cgc_path, sep='\t')

    if opts['plot']:
        import plot_data
        cgc_path = os.path.join(opts['output'], 'gene_list_overlap.pdf')
        plot_data.cgc_overlap(overlap_df, cgc_path, list_name)

    logger.info('Finished list_overlap sub-command.')
