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

agreed_predictions = set([
    'PLCG1', 'CRLF2', 'SMARCD1', 'SH2B3', 'STK11', 'MEN1', 'IKBKB', 'AKT1', 'B2M', 'MLH1',
    'USP28', 'TSHR', 'FGFR4', 'GPS2', 'CDC73', 'PIK3CA', 'MAP3K1', 'CACNA1D', 'FGFR3', 'TSC1',
    'ZC3H13', 'CBFB', 'ARHGAP26', 'CDH1', 'JAK2', 'ABCD2', 'NFKBIE', 'CUL1', 'GNAS', 'BRCA1',
    'ERBB3', 'BRCA2', 'DDX3X', 'PIK3R1', 'IDH2', 'IDH1', 'JAK1', 'KLF4', 'STAT5B', 'IWS1',
    'RPL10', 'SMAD2', 'CSF3R', 'KRAS', 'SETD2', 'FGFR2', 'GATA3', 'UBR5', 'ALK', 'KAT8',
    'SOCS1', 'MAX', 'PTPN11', 'ASXL1', 'POLE', 'TCF7L2', 'KIT', 'FOXA1', 'DDX5', 'FAT1',
    'RUNX1', 'CYLD', 'CD79B', 'ARID1B', 'APC', 'DAXX', 'ARID1A', 'CTCF', 'KDM5C', 'IL7R',
    'JAK3', 'CD1D', 'ACVR1B', 'CDKN1B', 'HLA-A', 'RPL5', 'COPS4', 'GGCT', 'TNFRSF14', 'FBXO11',
    'SETBP1', 'TRAF7', 'DNM2', 'CCAR1', 'ARID4B', 'XPO1', 'SPEN', 'KANSL1', 'CHD4', 'WDR47',
    'U2AF1', 'TTR', 'GK2', 'SPOP', 'PPP6C', 'NF2', 'NF1', 'PCBP1', 'AMPH', 'MYD88',
    'GRM3', 'MYH2', 'ATP1A1', 'STAG2', 'ARID2', 'RNF43', 'TUBA3C', 'TRRAP', 'PAX5', 'CUX1',
    'HRAS', 'RAD21', 'HMCN1', 'ING1', 'DNER', 'MGA', 'TP53', 'GNAQ', 'ESR1', 'FAM47B',
    'MPL', 'CBL', 'STK31', 'KRT15', 'PRDM1', 'NFE2L2', 'NSD1', 'MYCN', 'AGTR1', 'CREBBP',
    'ZRSR2', 'PDGFRA', 'GIGYF2', 'SMAD4', 'GATA1', 'ATM', 'EPHA2', 'POT1', 'SMAD3', 'SMO',
    'TBX3', 'CBLB', 'ATR', 'BIRC3', 'ABL1', 'AMER1', 'FAM46C', 'RAF1', 'FLT3', 'NCOR1',
    'CD79A', 'NCOA2', 'TP63', 'STAT3', 'H3F3B', 'MYH9', 'ECT2L', 'CDKN2A(p14)', 'RET', 'MAPK1',
    'KDM6A', 'BRWD3', 'PENK', 'CDKN1A', 'KDR', 'FBXW7', 'TGFBR2', 'FUBP1', 'PDYN', 'RQCD1',
    'GATA2', 'EGFR', 'TET2', 'ZFP36L1', 'ZFP36L2', 'MTOR', 'PRKAR1A', 'BCOR', 'ATRX', 'EP300',
    'TNFAIP3', 'DICER1', 'TBL1XR1', 'KCNJ5', 'MAP2K2', 'COL2A1', 'ALB', 'MAP2K1', 'KEAP1', 'EZH2',
    'CDK4', 'RAC1', 'PBRM1', 'CMTR2', 'BRE', 'RHEB', 'LPAR4', 'AMOT', 'CIC', 'PPP2R1A',
    'ACVR1', 'WT1', 'ZNF318', 'MSH2', 'SF3B1', 'MSH6', 'CTNNB1', 'VHL', 'USP9X', 'SOX9',
    'NOTCH2', 'MAP2K4', 'ELF3', 'SMARCA4', 'H3F3A', 'CEBPA', 'AXIN2', 'AXIN1', 'TWIST1', 'FAS',
    'NRAS', 'RB1', 'CDKN2A', 'KLF6', 'MED12', 'HNF1A', 'ETNK1', 'ATG5', 'CNOT3', 'NRG3',
    'TERT', 'AJUBA', 'NT5C2', 'BRAF', 'KMT2C', 'KMT2B', 'KMT2A', 'DNMT3A', 'SMARCB1', 'KMT2D',
    'PTEN', 'RBM10', 'CARD11', 'GNA11', 'RHOA', 'PTPRB', 'MAP3K13', 'HIST1H3B', 'PHF6', 'TP53BP1',
    'TSC2', 'SUFU', 'DACH1', 'TRIP12', 'PHOX2B', 'NPM1', 'RASA1', 'MYOD1', 'PTCH1', 'ERBB2',
    'CALR', 'SRSF2', 'DKK2', 'NOTCH1', 'CASP8', 'CDK12', 'GRIN2A', 'PTPRC', 'ARHGAP35', 'CHD8',
    'FOXL2', 'BAP1', 'BCL6', 'MET'
])

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


