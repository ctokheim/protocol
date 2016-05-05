#!/usr/bin/env python
import utils
import pandas as pd
import numpy as np
from scipy import stats
import datetime
import argparse
import logging
import os
import sys
import glob
import re
import eval_stats
import utils

logger = logging.getLogger(__name__)  # module logger

def main(opts):
    logger.info('Running consistency sub-command . . .')
    # read in config file
    config = utils.load_config(opts['config'])

    # read in results
    sim_results = {}
    # iterate over random splits
    for i, mydir in enumerate(os.listdir(opts['consistency_dir'])):
        # figure out the methods
        iter_dir = os.path.join(opts['consistency_dir'], mydir)
        methods = set(map(lambda x: re.split( '_[12]', os.path.basename(x))[0],
                          os.listdir(iter_dir)))

        # iterate over methods
        for method in methods:
            # prepare path
            first_outfix = '{prefix}_1.txt'.format(prefix=method)
            second_outfix = '{prefix}_2.txt'.format(prefix=method)
            first_path = os.path.join(iter_dir, first_outfix)
            second_path = os.path.join(iter_dir, second_outfix)

            # read data
            first_df = pd.read_csv(first_path, sep='\t', index_col=0)
            second_df = pd.read_csv(second_path, sep='\t', index_col=0)

            # calculate consistency
            consistency_df = eval_stats.consistency_comparison(first_df, second_df,
                                                               opts['depth'],
                                                               method,
                                                               config)
            sim_results.setdefault(method, [])
            #sim_results[method][i] = consistency_df
            sim_results[method].append(consistency_df)

    # concatenate each iteration into dataframe
    for method in sim_results:
        sim_results[method] = pd.concat(sim_results[method])

    # record result
    final_results = []
    for method in sim_results:
        tmp_results = eval_stats.calculate_stats(sim_results[method])
        final_results.append(tmp_results)
    final_df = pd.concat(final_results)
    output_path = os.path.join(opts['output'], 'consistency.txt')
    final_df.to_csv(output_path, sep='\t')

    # save plot if specified
    if opts['plot']:
        import plot_data
        output_path = os.path.join(opts['output'], 'consistency.pdf')
        plot_data.consistency(final_df, opts['depth'], output_path)

    logger.info('Finished consistency sub-command.')
    return final_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
