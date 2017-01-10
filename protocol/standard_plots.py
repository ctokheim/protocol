import sys
import argparse
import utils
import os
import logging
import pandas as pd
import cgc_overlap
import method_overlap
import p_value
import num_genes
import consistency
import plot_data
import utils
import numpy as np

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    # make a parser
    info = 'Evaluates cancer driver gene methods'
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
                        default='stdout',
                        help='Path to log file. (accepts "stdout")')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

    # program arguments
    help_str = 'directory containing results from methods on full data'
    parser.add_argument('-i', '--input-dir',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Mutation file (MAF file)'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Method name (should be same as listed inside config.yaml))'
    parser.add_argument('-n', '--method-name',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Name of column with cancer type information in MAF file'
    parser.add_argument('-t', '--tumor-type-col',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Path to Cancer Gene Census file'
    parser.add_argument('-c', '--cgc',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Landscapes driver gene list'
    parser.add_argument('-landscapes', '--landscapes',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Configuration file (YAML format)'
    parser.add_argument('-config', '--config',
                        type=str, default=None,
                        help=help_str)
    help_str = ('Q-value threshold for significance (Default: 0.1). '
                'This option is only specified if you do not include a threshold in '
                'the config file. The value in the configuration file will take '
                'precedence.')
    parser.add_argument('-q', '--qvalue',
                        type=float, default=.1,
                        help=help_str)
    help_str = 'output directory'
    parser.add_argument('-o', '--output',
                        type=str, default=None,
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
    utils.start_logging(log_file=log_file,
                        log_level=log_level,
                        verbose=args.verbose)  # start logging

    opts = vars(args)

    # log user entered command
    #logger.info('Version: {0}'.format(prob2020.__version__))
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


# count number of samples for each tumor type
def get_num_samples(df_dict):
    num_sample_dict = dict()
    for key, df in df_list:
        num_sample_dict[key] = df['Tumor_Sample_Barcode'].nunique()
    return num_sample_dict


# get the average number of drivers per sample assuming all non-silent
# mutations "count" for significant genes
def get_avg_drivers(df_dict, sig_genes):
    num_avg_drivers_dict = dict()
    drivers_per_samp_dict = dict()
    for ttype, df in df_dict.iteritems():
        # handle samples
        all_samples = df['Tumor_Sample_Barcode'].unique()
        df = df[df['is_nonsilent']==1]  # keep nonsilent mutations
        df = df.drop_duplicates(subset=['Tumor_Sample_Barcode', 'Hugo_Symbol'])

        # count number of driver genes with non-silent mutations
        if type(sig_genes) is dict:
            df_tmp = df[df['Hugo_Symbol'].isin(sig_genes[ttype])]
        else:
            df_tmp = df[df['Hugo_Symbol'].isin(sig_genes)]
        num_drivers = df_tmp.groupby('Tumor_Sample_Barcode')['is_nonsilent'].sum()

        # handle cases with samples with no driver genes mutated
        no_drivers = list(set(all_samples) - set(num_drivers.index))
        no_drivers_series = pd.Series(np.zeros(len(no_drivers)),
                                      index=no_drivers)
        num_drivers_per_samp = pd.concat([num_drivers, no_drivers_series])

        num_avg_drivers_dict[ttype] = np.mean(num_drivers_per_samp)
        drivers_per_samp_dict[ttype] = num_drivers_per_samp
    return num_avg_drivers_dict, drivers_per_samp_dict


def main(opts):
    # make output directory if it doesn't exist
    if not os.path.exists(opts['output']):
        os.makedirs(opts['output'])

    # load config file
    config = utils.load_config(opts['config'])

    # run commands
    cgc = utils.process_cgc(opts['cgc'])
    landscapes = cgc_overlap.read_custom_list(opts['landscapes'])

    # get significant genes
    signif_dict = utils.fetch_single_method_significant_genes(opts['input_dir'],
                                                              opts['method_name'],
                                                              opts['qvalue'],
                                                              config)
    pancan_df = pd.read_table(os.path.join(opts['input_dir'], 'PANCAN.txt'))


    ###########################
    # Pan-cancer plots
    ###########################
    # overlap with gene lists
    logger.info('Overlapping genes with CGC and Cancer Genome Landscapes . . .')
    pancan_genes = set(signif_dict['PANCAN'])
    intersect_cgc = len(pancan_genes & set(cgc))
    intersect_landscapes = len(pancan_genes & set(landscapes))
    intersect_all = len(pancan_genes & (set(cgc) | set(landscapes)))
    s = pd.Series([intersect_cgc, intersect_landscapes, intersect_all, len(pancan_genes)],
                   index=['CGC', 'Landscapes', 'Either', 'Method Total'])
    out_path = os.path.join(opts['output'], 'gene_list_overlap.pdf')
    plot_data.single_method_overlap(s, out_path)
    logger.info('Finished')

    # qq-plot
    # figure out which columns are being used
    method_name = opts['method_name']
    if utils.is_valid_config(config, method_name, 'pvalue'):
        pval_col = config[method_name]['pvalue'][0]
    else:
        pval_col = 'pvalue'
    logger.info('Creating p-value QQ plot . . .')
    out_path = os.path.join(opts['output'], 'qq_plot.pdf')
    plot_data.single_method_qqplot(pancan_df[pval_col], out_path)
    logger.info('Finished.')

    ###########################
    # Cancer type specific plots
    ###########################
    del signif_dict['PANCAN']
    if not signif_dict:
        return
    # process maf file
    mut_df = pd.read_table(opts['mutations'])
    drop_variants = ["3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "Intron",]
    mut_df = mut_df[~mut_df['Variant_Classification'].isin(drop_variants)]
    mut_df['is_nonsilent'] = 0
    mut_df.loc[mut_df['Variant_Classification']!='Silent', 'is_nonsilent'] = 1

    # plot number of drivers per sample
    logger.info('Analyzing number of drivers per sample . . .')
    ttype_data = dict()
    for ttype in mut_df['Tumor_Type'].dropna().unique():
        ttype_data[ttype] = mut_df[mut_df['Tumor_Type']==ttype].copy()
    driver_mean, driver_per_sample = get_avg_drivers(ttype_data, signif_dict)
    out_path = os.path.join(opts['output'], 'cancer_type_per_sample.pdf')
    order = plot_data.single_method_driver_per_sample(driver_per_sample,
                                                      out_path)
    logger.info('Finished.')

    # plot number of driver genes
    logger.info('Analyzing number of driver in each cancer type . . .')
    signif_ct = pd.Series({k: len(signif_dict[k]) for k in signif_dict})
    out_path = os.path.join(opts['output'], 'cancer_type_num_drivers.pdf')
    plot_data.single_method_num_drivers_per_type(signif_ct, order, out_path)
    logger.info('Finished.')


def cli_main():
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()