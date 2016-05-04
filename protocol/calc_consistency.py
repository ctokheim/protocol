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
import simulation as sim

logger = logging.getLogger(__name__)  # module logger

def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def parse_arguments():
    # make a parser
    info = 'Calculates similarity metrics between results of sampled data sets'
    parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')

    # program arguments
    help_str = 'Directory of output from run_consistency_simulation.sh'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Maximum depth of genes from top of list to consider for consistency. '
                'Only used when random half split flag is specified.')
    parser.add_argument('-depth', '--depth',
                        type=int,
                        default=100,
                        help=help_str)
    help_str = ('Prefix of filename containing p-values. The prefix should be '
                'followed by either "1.txt" or "2.txt" representing the first '
                'and second split, respectively (Default: "result").')
    parser.add_argument('-p', '--prefix',
                        type=str,
                        default='result',
                        help=help_str)
    help_str = 'Output file containing results'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging

    opts = vars(args)

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))

    return opts


def main(opts):
    # read in results
    sim_results = {}
    # iterate over methods
    for method_dir in os.listdir(opts['input']):
        full_method_dir = os.path.join(opts['input'], method_dir)
        # iterate over random splits
        for i, mydir in enumerate(full_method_dir):
            # prepare path
            first_outfix = '{prefix}1.txt'.format(prefix=opts['prefix'])
            second_outfix = '{prefix}2.txt'.format(prefix=opts['prefix'])
            first_path = os.path.join(full_method_dir, mydir, first_outfix)
            second_path = os.path.join(full_method_dir, mydir, second_outfix)

            # read data
            first_df = pd.read_csv(first_path, sep='\t', index_col=0)
            second_df = pd.read_csv(second_path, sep='\t', index_col=0)

            # calculate consistency
            consistency_df = consistency_comparison(first_df, second_df)
            sim_results.setdefault(method_dir, {})
            sim_results[method_dir][i] = consistency_df

    # record result for a specific sample rate
    final_results = {}
    for method in sim_results:
        tmp_results = calculate_stats(sim_results)
        final_results[method] = tmp_results
        #tmp_results.to_csv(opts['output'], sep='\t')
    final_df = pd.DataFrame(final_results)
    final_df.to_csv(opts['output'], sep='\t')

    return final_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
