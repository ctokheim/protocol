import numpy as np
import pandas as pd
import eval_stats
import utils
import os

# logging
import logging
logger = logging.getLogger(__name__)

def calculate_mlfc(df, meth_name, config):
    # check whether there is custom p-value columns
    if utils.is_valid_config(config, meth_name, 'pvalue'):
        pval_cols = config[meth_name]['pvalue']
    else:
        pval_cols = ['pvalue']

    # calculate the mlfc scores for each pvalue column
    mlfc_score_list = []
    for pval_col in pval_cols:
        pvals = df[pval_col].dropna().copy()
        tmp_score = eval_stats.mean_log_fold_change(pvals)
        mlfc_score_list.append(tmp_score)

    # calculate the mean score, if multiple values
    mean_mlfc = np.mean(mlfc_score_list)

    return mean_mlfc


def main(opts):
    logger.info('Running p-value sub-command . . .')
    # load config
    config = utils.load_config(opts['config'])

    # get full data for each method
    full_data = utils.fetch_filtered_dataframes(opts['input_dir'],
                                                opts['output'],
                                                opts['min'])

    # compute MLFC scores
    mlfc_result = {m: calculate_mlfc(full_data[m], m, config)
                   for m in full_data}
    mlfc_series = pd.Series(mlfc_result)

    # save output
    mlfc_path = os.path.join(opts['output'], 'mlfc_scores.txt')
    mlfc_series.to_csv(mlfc_path, sep='\t')

    if opts['plot']:
        import plot_data
        mlfc_path = os.path.join(opts['output'], 'mlfc_scores.pdf')
        plot_data.mlfc_score(mlfc_series, mlfc_path)

    logger.info('Finished p-value sub-command.')
    return mlfc_result
