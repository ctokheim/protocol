import pandas as pd
import os
from collections import Counter
import utils

# logging
import logging
logger = logging.getLogger(__name__)

def num_cgc_overlap(signif_genes, cancer_type, cgc_list):
    """Intersect significant genes with CGC or other driver gene list."""
    cgc_result = {}
    for method in signif_genes:
        mygenes = signif_genes[method].get(cancer_type, [])
        intersect = len(set(mygenes) & set(cgc_list))
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
                                                config)
    num_methods = len(signif_dict)

    # read in driver genes
    driver_list_names = ['cgc', 'landscapes', 'kandoth',
                         'high_confidence_list', 'gene_list']
    driver_genes = set()
    for l in driver_list_names:
        # skip if list not provided
        if opts[l] is None: continue

        # add genes to list
        if l == 'cgc':
            # get cgc genes
            driver_genes |= set(utils.process_cgc(opts['cgc']))
            list_name = 'CGC'
        else:
            # read in custom driver gene list
            driver_genes |= set(read_custom_list(opts[l]))
            list_name = 'Custom List'
    driver_genes = list(driver_genes)

    all_overlap_list = []
    cancer_types = signif_dict[signif_dict.keys()[0]].keys()
    for cancer_type in cancer_types:
        # count the overlap
        num_cgc_dict = num_cgc_overlap(signif_dict, cancer_type, driver_genes)
        num_signif_dict = {k: len(signif_dict[k].get(cancer_type, []))
                           for k in signif_dict}

        # format result
        overlap_df = pd.DataFrame({'# '+list_name: pd.Series(num_cgc_dict),
                                '# significant': pd.Series(num_signif_dict)})
        frac_ovlp_col = 'Fraction overlap ({0})'.format(list_name)
        overlap_df[frac_ovlp_col] = overlap_df['# '+list_name].astype(float) / overlap_df['# significant']
        overlap_df[frac_ovlp_col] = overlap_df[frac_ovlp_col].fillna(0)

        # append result
        overlap_df['CODE'] = cancer_type
        all_overlap_list.append(overlap_df)

        if opts['plot']:
            import plot_data
            cgc_path = os.path.join(opts['output'], cancer_type+'.driver_list_overlap.pdf')
            plot_data.cgc_overlap(overlap_df, cgc_path, list_name)

    # save result
    cgc_path = os.path.join(opts['output'], 'gene_list_overlap.txt')
    pd.concat(all_overlap_list).to_csv(cgc_path, sep='\t')

    logger.info('Finished list_overlap sub-command.')
