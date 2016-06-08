import sys
import argparse
import utils
import os
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
    parser_pipeline = subparsers.add_parser('pipeline',
                                            help='Run all sub-commands evaluating methods',
                                            description='Run all sub-commands evaluating methods')
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
    help_info = 'Evaluate method consistency'
    parser_consis = subparsers.add_parser('consistency',
                                          help=help_info,
                                          description='Evaluate method consistency')

    # program arguments
    for i, parser in enumerate([parser_pipeline, parser_split, parser_cgc, parser_ovlp,
                                parser_pval, parser_signif, parser_consis]):
        # group of parameters
        major_parser = parser.add_argument_group(title='Major options')
        advance_parser = parser.add_argument_group(title='Advanced options')

        if i != 1 and i != 6:
            help_str = 'directory containing results from methods on full data'
            major_parser.add_argument('-i', '--input-dir',
                                    type=str, default=None,
                                    help=help_str)
        elif i == 1:
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
            help_str = 'Path to Cancer Gene Census file'
            major_parser.add_argument('-c', '--cgc',
                                      type=str, required=True,
                                      help=help_str)
            help_str = ('Minimum number of methods finding a gene significant to '
                        'not include that gene\' p-value (Default: 3)')
            major_parser.add_argument('-m', '--min',
                                      type=int, default=3,
                                      help=help_str)
            help_str = 'Directory containing the consistency results'
            major_parser.add_argument('-consis-dir', '--consistency-dir',
                                      type=str, required=True,
                                      help=help_str)
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
            help_str = 'Ranking depth to consider for consistency (Default: 100)'
            advance_parser.add_argument('-d', '--depth',
                                        type=int, default=100,
                                        help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        if i == 1:
            help_str = ('Column name containing sample IDs (Default: checks '
                        '"Tumor_Sample_Barcode" or "Tumor_Sample")')
            advance_parser.add_argument('-s', '--sample-col',
                                        type=str,
                                        help=help_str)
            help_str = 'Number of iterations to randomly split data (Default: 10)'
            advance_parser.add_argument('-n', '--number',
                                        type=int, default=10,
                                        help=help_str)
        elif i == 2:
            help_str = 'Path to Cancer Gene Census file'
            major_parser.add_argument('-c', '--cgc',
                                      type=str, required=True,
                                      help=help_str)
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 3:
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 4:
            help_str = ('Minimum number of methods finding a gene significant to '
                        'not include that gene\' p-value (Default: 3)')
            major_parser.add_argument('-m', '--min',
                                      type=int, default=3,
                                      help=help_str)
            help_str = ('Path to Cancer Gene Census file (optional). Additionaly '
                        'removes genes from the CGC.')
            major_parser.add_argument('-c', '--cgc',
                                      type=str, default=None,
                                      help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 5:
            help_str = 'Q-value threshold for significance (Default: 0.1)'
            advance_parser.add_argument('-q', '--qvalue',
                                        type=float, default=.1,
                                        help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 6:
            help_str = 'Directory containing the consistency results'
            major_parser.add_argument('-consis-dir', '--consistency-dir',
                                      type=str, required=True,
                                      help=help_str)
            help_str = 'Ranking depth to consider for consistency (Default: 100)'
            advance_parser.add_argument('-d', '--depth',
                                        type=int, default=100,
                                        help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
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
    # make output directory if it doesn't exist
    if not os.path.exists(opts['output']):
        os.makedirs(opts['output'])

    # run commands
    if opts['kind'] == 'pipeline':
        import cgc_overlap
        import method_overlap
        import p_value
        import num_genes
        import consistency
        cgc_overlap.main(opts)
        method_overlap.main(opts)
        p_value.main(opts)
        num_genes.main(opts)
        consistency.main(opts)
    elif opts['kind'] == 'split_mutations':
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
    elif opts['kind'] == 'consistency':
        import consistency
        consistency.main(opts)


def cli_main():
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()
