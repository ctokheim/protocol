import pandas as pd
from scipy import stats
import numpy as np
import yaml

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


def fetch_significant_genes(input_dir, qval, config):
    signif_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        # figure out which columns are being used
        if is_valid_config(config, method_name, 'qvalue'):
            qval_cols = config[method_name]['qvalue']
        else:
            qval_cols = ['qvalue']

        tmp_genes = set()
        for qval_col in qval_cols:
            # get the significant genes
            signif_df = df[df[qval_col]<=qval]
            tmp_genes |= set(signif_df['gene'].tolist())
        signif_dict[method_name] = list(tmp_genes)
    return signif_dict


def read_method_overlap_genes(path, min_methods):
    # read data
    df = pd.read_table(path)

    # add up all genes with min overlap
    gene_list = []
    for i in range(min_methods, len(df.columns)+1):
        tmp = df[str(i)].dropna().tolist()
        gene_list += tmp

    return gene_list


def fetch_raw_dataframes(input_dir):
    data_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        data_dict[method_name] = df

    return data_dict


def fetch_filtered_dataframes(input_dir, output_dir, min_methods):
    # get the raw results
    df_dict = fetch_raw_dataframes(input_dir)

    # get the overlapping genes
    overlap_path = os.path.join(output_dir, 'gene_overlap_counts.txt')
    overlap_genes = read_method_overlap_genes(overlap_path, min_methods)

    # iterate over each method
    for method in df_dict:
        # filter out agreed upon genes
        df_copy = df_dict[method].copy()
        df_copy = df_copy[~df_copy['gene'].isin(overlap_genes)]

        # update dictionary
        df_dict[method] = df_copy

    return df_dict


def load_config(path):
    """Load YAML configuration file."""
    if path is not None:
        with open(path) as handle:
            config = yaml.load(handle)
        return config
    else:
        return None


def is_valid_config(myconfig, method_name, attribute):
    # check if the above attributes are valid
    has_config = myconfig is not None
    if not has_config:
        return False
    has_method = method_name in myconfig
    if not has_method:
        return False
    has_attribute = attribute in myconfig[method_name]
    if not has_attribute:
        return False

    # no problem, return True
    return True
