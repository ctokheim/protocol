import pandas as pd
from scipy import stats
import numpy as np

# logging import
import warnings
import logging
import datetime
import os
import sys

logger = logging.getLogger(__name__)  # module logger

def start_logging(log_file='', log_level='INFO', verbose=False):
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

    # ignore warnings if not in debug
    if log_level.upper() != 'DEBUG':
        warnings.filterwarnings('ignore')

    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

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


def fetch_significant_genes(input_dir, qval):
    signif_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        # get the significant genes
        signif_df = df[df['qvalue']<=qval]
        signif_dict[method_name] = signif_df['gene'].tolist()
    return signif_dict


def fetch_raw_dataframes(input_dir):
    data_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        data_dict[method_name] = df

    return data_dict
