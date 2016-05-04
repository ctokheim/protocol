import sys
import argparse
import utils
import logging

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    # make a parser
    info = 'Evaluates cancer driver gene methods'
    parent_parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parent_parser.add_argument('-ll', '--log-level',
                               type=str,
                               action='store',
                               default='',
                               help='Write a log file (--log-level=DEBUG for debug mode, '
                               '--log-level=INFO for info mode)')
    parent_parser.add_argument('-l', '--log',
                               type=str,
                               action='store',
                               default='stdout',
                               help='Path to log file. (accepts "stdout")')
    parent_parser.add_argument('-v', '--verbose',
                               action='store_true',
                               default=False,
                               help='Flag for more verbose log output')

    # add subparsers
    subparsers = parent_parser.add_subparsers(title='Sub-commands', dest='kind')
    parser_split = subparsers.add_parser('split_mutations',
                                         help='Splits mutations in a MAF-like format into two random halves',
                                         description='Splits mutations in a MAF-like format into two random halves. '
                                         'Each split maintains the proportion of samples in each cancer type.')
    parser_cgc = subparsers.add_parser('cgc_overlap',
                                       help='Evaluate the overlap of significant genes with the Cancer Gene Census (CGC)',
                                       description='Evaluate the overlap of significant genes with the Cancer Gene Census (CGC)')
    parser_ovlp = subparsers.add_parser('method_overlap',
                                        help='Counts the number of methods that find specific genes significant',
                                        description='Counts the number of methods that find specific genes significant')
    parser_pval = subparsers.add_parser('pvalue',
                                        help='Examine the p-value distribution',
                                        description='Examine the p-value distribution')
    parser_signif = subparsers.add_parser('num_signif',
                                          help='Examine the number of significant genes',
                                          description='Examine the number of significant genes')
    help_info = 'Find statistically significant Tumor Suppressor-like genes.'
    parser_tsg = subparsers.add_parser('tsg',
                                       help=help_info,
                                       description=help_info + ' Evaluates for a higher proportion '
                                       'of inactivating mutations than expected.')

    # program arguments
    for i, parser in enumerate([parser_split, parser_cgc, parser_ovlp, parser_pval,
                                parser_signif, parser_tsg]):
        # group of parameters
        major_parser = parser.add_argument_group(title='Major options')
        advance_parser = parser.add_argument_group(title='Advanced options')

        if i > 0:
            help_str = 'directory containing input files'
            major_parser.add_argument('-i', '--input-dir',
                                    type=str, default=None,
                                    help=help_str)
        else:
            help_str = 'Mutation file to split'
            major_parser.add_argument('-m', '--mutations',
                                      type=str, required=True,
                                      help=help_str)
        help_str = 'Configuration file (YAML format)'
        major_parser.add_argument('-config', '--config',
                                  type=str, default=None,
                                  help=help_str)
        help_str = 'output directory'
        major_parser.add_argument('-o', '--output',
                                  type=str, default=None,
                                  help=help_str)

        if i == 0:
            help_str = ('Column name containing sample IDs (Default: checks '
                        '"Tumor_Sample_Barcode" or "Tumor_Sample")')
            advance_parser.add_argument('-s', '--sample-col',
                                        type=str,
                                        help=help_str)
            help_str = 'Number of iterations to randomly split data (Default: 10)'
            advance_parser.add_argument('-n', '--number',
                                        type=int, default=10,
                                        help=help_str)
        elif i == 1:
            help_str = 'Path to Cancer Gene Census file'
            major_parser.add_argument('-c', '--cgc',
                                      type=str, required=True,
                                      help=help_str)
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
        elif i == 2:
            help_str = ('Compare overlap with other methods already evaluated '
                        'in Tokheim et al. 2016 (choose "pancan", "breast", '
                        '"pancreatic", "head", "lung"). Default, doesn\'t include '
                        'previously computed genes.')
            advance_parser.add_argument('-p', '--pre-computed',
                                        type=str, default=None,
                                        help=help_str)
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
        elif i == 3:
            help_str = ('Minimum number of methods finding a gene significant to '
                        'not include that gene\' p-value (Default: 3)')
            major_parser.add_argument('-m', '--min',
                                      type=int, default=3,
                                      help=help_str)
        elif i == 4:
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
    args = parent_parser.parse_args()

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


def main(opts):
    if opts['kind'] == 'split_mutations':
        import split_mutations
        split_mutations.main(opts)
    elif opts['kind'] == 'cgc_overlap':
        import cgc_overlap
        cgc_overlap.main(opts)
    elif opts['kind'] == 'method_overlap':
        import method_overlap
        method_overlap.main(opts)
    elif opts['kind'] == 'pvalue':
        import p_value
        p_value.main(opts)
    elif opts['kind'] == 'num_signif':
        import num_genes
        num_genes.main(opts)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


