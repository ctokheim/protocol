import numpy as np
import pandas as pd
import eval_stats
import utils
import config as cfg
import os
import cgc_overlap

# logging
import logging
logger = logging.getLogger(__name__)

def calculate_mlfc(df, meth_name, config):
    pval_cols = ['pvalue']

    # calculate the mlfc scores for each pvalue column
    mlfc_score_list = []
    for pval_col in pval_cols:
        #pvals = df[pval_col].dropna().copy()
        pvals = df[pval_col].fillna(1.0).copy()
        tmp_score = eval_stats.mean_log_fold_change(pvals)
        mlfc_score_list.append(tmp_score)

    # calculate the mean score, if multiple values
    mean_mlfc = np.mean(mlfc_score_list)

    return mean_mlfc


def main(opts):
    logger.info('Running p-value sub-command . . .')
    # load config
    config = utils.load_config(opts['config'])

    # Read in driver lists
    cgc = utils.process_cgc(opts['cgc'])
    landscapes = cgc_overlap.read_custom_list(opts['landscapes'])
    kandoth = cgc_overlap.read_custom_list(opts['kandoth'])
    tamborero = cgc_overlap.read_custom_list(opts['high_confidence_list'])
    blacklist = set(cgc) | set(landscapes) | set(kandoth) | set(tamborero)

    # iterate over methods
    gene_methods = cfg.fetch_level_names(config, level='gene')
    mlfc_list = []
    for method_name in gene_methods:
        # get full data for each method
        full_data = utils.fetch_filtered_dataframes(opts['input_dir'],
                                                    blacklist,
                                                    method_name)

        # compute MLFC scores
        mlfc_result = {m: calculate_mlfc(full_data[m], m, config)
                       for m in full_data}
        mlfc_series = pd.Series(mlfc_result)
        mlfc_series.name = method_name
        mlfc_list.append(mlfc_series)

    # creat a dataframe
    mlfc_df = pd.concat(mlfc_list, axis=1)

    # save output
    mlfc_path = os.path.join(opts['output'], 'mlfc_scores.txt')
    mlfc_df.to_csv(mlfc_path, sep='\t')

    #if opts['plot']:
        #import plot_data
        #mlfc_path = os.path.join(opts['output'], 'mlfc_scores.pdf')
        #plot_data.mlfc_score(mlfc_series, mlfc_path)

    logger.info('Finished p-value sub-command.')
    return mlfc_df
