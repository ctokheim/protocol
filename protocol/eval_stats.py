import numpy as np
import scipy.stats as stats
import pandas as pd
import utils
import sys


def top_drop_overlap(s1, s2, depth):
    # get gene of interest by specified depth
    s1, s2 = s1[:depth], s2[:2*depth]

    # genes are supposed to be the index of the series
    s1_genes = set(s1.index)
    s2_genes = set(s2.index)

    # calculate jaccard index
    num_intersect = len(s1_genes & s2_genes)
    num_total = len(s1_genes)
    if num_total:
        # provided series are not empty
        ov_sim = num_intersect / float(num_total)
    else:
        # empty series case
        ov_sim = 0
    return ov_sim


def calculate_sem(wp):
    """Calculates the standard error of the mean for a pd.Panel object.

    **Note:** The pd.Panel.apply method seems to have a bug preventing
    me from using it. So instead I am using the numpy apply function
    for calculating sem.

    Parameters
    ----------
    wp : pd.Panel
        panel that stratifies samples

    Returns
    -------
    tmp_sem : pd.DataFrame
        standard error of the mean calculated along the sample axis
    """
    #tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)  # hack because pandas apply method has a bug
    #tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           #columns=wp.minor_axis,
                           #index=wp.major_axis)
    return tmp_sem


def calculate_stats(df,
                    metrics=['precision', 'recall', 'ROC AUC', 'PR AUC', 'count']):
    """Computes mean and sem of classification performance metrics.

    Parameters
    ----------
    result_dict : dict
        dictionary with the i'th sample as the key and data frames
        with "oncogene"/"tsg" (row) classification performance metrics
        (columns) as values

    Returns
    -------
    result_df : pd.DataFrame
        Data frame with mean and sem of classification performance
        metrics. (rows: "oncogene"/"tsg", columns: summarized metrics)
    """
    tmp_means = df.mean()
    tmp_sem = df.sem()
    name = df.columns[0]
    result_df = pd.DataFrame({name+' mean': tmp_means,
                              name+' sem': tmp_sem})
    result_df.rename(index=lambda x: df.index[0], inplace=True)
    return result_df


def consistency_comparison(df1, df2, mydepth,
                           method, config):
    """Function called by multiprocessing to run predictions.

    """
    # figure out which column to use
    if utils.is_valid_config(config, method, 'consistency'):
        pval_cols = config[method]['consistency']
    else:
        pval_cols = ['pvalue']

    # check to see if the columns exist
    all_exist = all([p in df1.columns for p in pval_cols])
    if not all_exist:
        print('Not all p-value columns were named correctly. Please specify '
              'correct name in config file.')
        sys.exit(1)

    results = pd.DataFrame(index=[method])

    pval1 = df1[pval_cols].min(axis=1).copy()
    pval2 = df2[pval_cols].min(axis=1).copy()
    pval1.sort_values(ascending=True, inplace=True)
    pval2.sort_values(ascending=True, inplace=True)

    # add top-drop scores
    top_drop = top_drop_overlap(pval1, pval2, depth=mydepth)
    td_scores = [top_drop]

    tmp_results = pd.DataFrame({
        'TopDrop {0} overlap'.format(mydepth): td_scores
    },
    index=[method])
    results = pd.concat([results, tmp_results], axis=1)

    return results


def mean_log_fold_change(data):
    tmp = data.copy()
    tmp.sort_values(ascending=True, inplace=True)
    tmp[tmp==0] = tmp[tmp>0].min()  # avoid infinity in log by avoiding zero pvals
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp))
    mlfc = np.mean(np.abs(np.log2(tmp/dist_quant)))
    return mlfc
