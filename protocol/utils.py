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
import config as cfg

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


def fetch_significant_genes(input_dir, method_name, config):
    """Read the significant driver genes for each method."""
    signif_dict = {}
    gene = 'gene'
    meth_input_dir = os.path.join(input_dir, method_name)
    for method_file in os.listdir(meth_input_dir):
        if not method_file.endswith('.txt'): continue
        if method_file.upper().startswith('README'): continue
        cancer_type_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(meth_input_dir, method_file)
        df = pd.read_table(full_path)

        # get the treshold for significance
        thresh_col, score_val, top_direction = cfg.fetch_threshold(config, method_name)

        # figure out if a custom score or q-value is used
        if cfg.is_valid_config(config, method_name, 'threshold'):
            #score_val = float(config[method_name]['threshold']['score'])
            tmp_genes = set()
            # get the top scoring genes
            if top_direction == 'low':
                signif_df = df[df[thresh_col]<=score_val]
            else:
                signif_df = df[df[thresh_col]>=score_val]
            tmp_genes |= set(signif_df[gene].tolist())
            signif_dict[cancer_type_name] = list(tmp_genes)
        else:
            # use q-value for threshold
            tmp_genes = set()
            # get the significant genes
            signif_df = df[df[thresh_col]<=score_val]
            tmp_genes |= set(signif_df[gene].tolist())
            signif_dict[cancer_type_name] = list(tmp_genes)
    return signif_dict


def fetch_single_method_significant_genes(input_dir, method_name, config):
    """Read the significant driver genes for a single method."""
    signif_dict = {}
    gene = 'gene'
    meth_input_dir = os.path.join(input_dir, method_name)
    for method_file in os.listdir(meth_input_dir):
        if not method_file.endswith('.txt'): continue
        if method_file.upper().startswith('README'): continue
        cancer_type_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(meth_input_dir, method_file)
        df = pd.read_table(full_path)

        # get the treshold for significance
        thresh_col, score_val, top_direction = cfg.fetch_threshold(config, method_name)

        # figure out if a custom score or q-value is used
        if cfg.is_valid_config(config, method_name, 'threshold'):
            tmp_genes = set()
            # get the top scoring genes
            top_direction = config[method_name]['threshold']['top']
            if top_direction == 'low':
                signif_df = df[df[thresh_col]<=score_val]
            else:
                signif_df = df[df[thresh_col]>=score_val]
            tmp_genes |= set(signif_df[gene].tolist())
            signif_dict[cancer_type_name] = list(tmp_genes)
        else:
            # use q-value for threshold
            tmp_genes = set()
            # get the significant genes
            signif_df = df[df[thresh_col]<=score_val]
            tmp_genes |= set(signif_df[gene].tolist())
            signif_dict[cancer_type_name] = list(tmp_genes)
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


def process_cgc(path):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    s = df['Mutation Types']
    is_small = s.str.contains('Mis|F|N|S').fillna(False)
    is_somatic = ~df['Tumour Types(Somatic)'].isnull()
    df = df[is_small & is_somatic].copy()

    # get gene names
    cgc_genes = df['Gene Symbol'].tolist()

    return cgc_genes


def fetch_raw_dataframes(input_dir, method_name):
    data_dict = {}
    meth_input_dir = os.path.join(input_dir, method_name)
    for method_file in os.listdir(meth_input_dir):
        #method_name = os.path.splitext(method_file)[0]
        cancer_type = os.path.splitext(method_file)[0]

        # skip READMEs
        if method_file.upper().startswith('README'): continue

        # read in data
        full_path = os.path.join(meth_input_dir, method_file)
        df = pd.read_table(full_path)

        #data_dict[method_name] = df
        data_dict[cancer_type] = df

    return data_dict


def fetch_filtered_dataframes(input_dir, output_dir, min_methods, config, cgc_path=None):
    """Return the result files for each method as a dataframe, but with certain
    genes filtered out. This includes those agreed upon by some minimum
    number of methods or from the CGC.
    """
    # get the raw results
    df_dict = fetch_raw_dataframes(input_dir)

    # get the overlapping genes
    overlap_path = os.path.join(output_dir, 'gene_overlap_counts.txt')
    overlap_genes = read_method_overlap_genes(overlap_path, min_methods)

    # gene column name
    if cfg.is_valid_config(config, method_name, 'gene_col'):
        gene = config[method_name]['gene_col']
    else:
        gene = 'gene'

    # add in cgc genes, if available
    if cgc_path is not None:
        cgc_genes = process_cgc(cgc_path)
        overlap_genes = list(set(overlap_genes + cgc_genes))

    # iterate over each method
    for method in df_dict:
        # filter out agreed upon genes
        df_copy = df_dict[method].copy()
        df_copy = df_copy[~df_copy[gene].isin(overlap_genes)]

        # update dictionary
        df_dict[method] = df_copy

    return df_dict


def read_filtered_pvalues(input_dir, blacklist, config, method_name):
    """Get p-values without 'likely' driver genes"""
    df_dict = fetch_raw_dataframes(input_dir, method_name)

    # gene column name
    if cfg.is_valid_config(config, method_name, 'gene_col'):
        gene = config[method_name]['gene_col']
    else:
        gene = 'gene'

    # pvalue column
    if cfg.is_valid_config(config, method_name, 'pvalue'):
        pval_col = config[method_name]['pvalue'][0]
    else:
        pval_col = 'pvalue'

    # get pvalues
    pval_dict = {}
    for cancer_type in df_dict:
        tmp = df_dict[cancer_type]
        tmp = tmp[~tmp[gene].isin(blacklist)]
        pval_dict[cancer_type] = tmp[[pval_col]]

    return pval_dict


def load_config(path):
    """Load YAML configuration file."""
    if path is not None:
        with open(path) as handle:
            config = yaml.load(handle)
        return config
    else:
        return None


