import numpy as np
import pandas as pd


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
    tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)  # hack because pandas apply method has a bug
    tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           columns=wp.minor_axis,
                           index=wp.major_axis)
    return tmp_sem


def calculate_stats(result_dict,
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
    wp = pd.Panel(result_dict)
    tmp_means = wp.mean(axis=0)
    tmp_sem = calculate_sem(wp)
    result_df = pd.merge(tmp_means, tmp_sem,
                         left_index=True, right_index=True,
                         suffixes=(' mean', ' sem'))
    return result_df


def consistency_comparison(df1, df2, mydepth):
    """Function called by multiprocessing to run predictions.

    """

    results = pd.DataFrame(index=['driver'])

    pval1 = df1['pvalue'].copy()
    pval2 = df2['pvalue'].copy()
    pval1.sort_values(ascending=True, inplace=True)
    pval2.sort_values(ascending=True, inplace=True)

    # add top-drop scores
    top_drop = sim.top_drop_overlap(pval1, pval2, depth=mydepth)
    td_scores = [top_drop]

    tmp_results = pd.DataFrame({
        'TopDrop {0} overlap'.format(mydepth): td_scores
    },
    index=['driver'])
    results = pd.concat([results, tmp_results], axis=1)

    return results


def mean_log_fold_change(data):
    tmp = data.copy()
    tmp.sort_values(ascending=True, inplace=True)
    tmp[tmp==0] = tmp[tmp>0].min()  # avoid infinity in log by avoiding zero pvals
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp))
    mlfc = np.mean(np.abs(np.log2(tmp/dist_quant)))
    return mlfc
