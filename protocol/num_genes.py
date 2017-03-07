import pandas as pd
import utils
import os
import csv

# logging
import logging
logger = logging.getLogger(__name__)


def iqr_outlier_cutoff(myseries, multiplier=1.5):
    """Returns the cutoff for an outlier based on the interquartile range."""
    # calculate IQR
    q1 = myseries.quantile(.25)
    q3 = myseries.quantile(.75)
    iqr = q3 - q1

    # get outlier cutoff
    cutoff = q3 + iqr*multiplier

    return cutoff



def main(opts):
    logger.info('Running num_signif sub-command . . .')
    config = utils.load_config(opts['config'])

    # get the significant genes for each method
    signif_dict = utils.fetch_significant(opts['input_dir'],
                                          config,
                                          level='gene')
    num_methods = len(signif_dict)

    # count the number of significant for each method
    num_signif_list = []
    outlier_dict = {}
    for m in signif_dict:
        # calculate the total number of significant
        num_signif = {c: len(signif_dict[m][c])
                      for c in signif_dict[m]
                      if c != 'PANCAN'}
        num_series = pd.Series(num_signif)
        num_series.name = m

        # figure out outliers
        cutoff = iqr_outlier_cutoff(num_series)
        outliers = num_series[num_series>cutoff]

        # append outliers
        outlier_dict[m] = outliers.index.tolist()

        # append number of significant for method to list
        num_signif_list.append(num_series)

    # save result
    if not os.path.exists(opts['output']): os.makedirs(opts['output'])
    # save number of significant
    out_path = os.path.join(opts['output'], 'num_significant.txt')
    pd.concat(num_signif_list, axis=1).to_csv(out_path, sep='\t')
    # save outliers
    out_path = os.path.join(opts['output'], 'outliers.txt')
    with open(out_path, 'w') as handle:
        mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for meth in outlier_dict:
            tmp_line = [meth, ','.join(outlier_dict[meth])]
            mywriter.writerow(tmp_line)

    logger.info('Finished num_signif sub-command.')
    return num_series
