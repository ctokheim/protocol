import sys
import argparse
import utils
import os
import logging

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    # make a parser
    info = 'Consensus driver gene analysis'
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
    parser_cgc = subparsers.add_parser('list_overlap',
                                       help='Evaluate the overlap of significant genes with a previous driver gene list',
                                       description='Evaluate the overlap of significant genes with a previous driver gene list')
    parser_ovlp = subparsers.add_parser('method_overlap',
                                        help='Counts the number of methods that find specific genes significant',
                                        description='Counts the number of methods that find specific genes significant')
    parser_pval = subparsers.add_parser('pvalue',
                                        help='Examine the p-value distribution',
                                        description='Examine the p-value distribution')
    parser_signif = subparsers.add_parser('num_signif',
                                          help='Examine the number of significant genes',
                                          description='Examine the number of significant genes')
    help_info = 'Perform driver gene consensus analysis'
    parser_consensus = subparsers.add_parser('consensus',
                                             help=help_info,
                                             description='Perform driver gene consensus analysis')

    # program arguments
    for i, parser in enumerate([parser_cgc, parser_ovlp,
                                parser_pval, parser_signif,
                                parser_consensus]):

        # group of parameters
        major_parser = parser.add_argument_group(title='Major options')
        advance_parser = parser.add_argument_group(title='Advanced options')

        if i != 4:
            help_str = 'directory containing results from methods on full data'
            major_parser.add_argument('-i', '--input-dir',
                                    type=str, default=None,
                                    help=help_str)

        help_str = 'Configuration file (YAML format)'
        major_parser.add_argument('-config', '--config',
                                  type=str, default=None,
                                  help=help_str)
        help_str = 'output directory'
        major_parser.add_argument('-o', '--output',
                                  type=str, default=None,
                                  help=help_str)

        if i == 0 or i == 2:
            #list_parser = major_parser.add_mutually_exclusive_group(required=True)
            help_str = 'Path to Cancer Gene Census file'
            major_parser.add_argument('-c', '--cgc',
                                      type=str, default=None,
                                      help=help_str)
            help_str = 'Cancer genome landscapes driver gene list'
            major_parser.add_argument('-landscapes', '--landscapes',
                                type=str, default=None,
                                help=help_str)
            help_str = 'List of Pancan12 drivers from Kandoth et al.'
            major_parser.add_argument('-k', '--kandoth',
                                type=str, default=None,
                                help=help_str)
            help_str = 'List of High confidence Pancan12 drivers from Tamborero et al.'
            major_parser.add_argument('-hcd', '--high-confidence-list',
                                type=str, default=None,
                                help=help_str)
            help_str = 'Custom driver gene list'
            advance_parser.add_argument('-g', '--gene-list',
                                type=str, default=None,
                                help=help_str)
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 1:
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 3:
            help_str = 'Generate plots examining evaluation (Default: False)'
            advance_parser.add_argument('-p', '--plot',
                                        action='store_true', default=False,
                                        help=help_str)
        elif i == 4:
            help_str = 'Flag indicating not to weight methods'
            major_parser.add_argument('--not-weighted',
                                      action='store_true', default=False,
                                      help=help_str)
            help_str = 'File containing MLFC scores'
            major_parser.add_argument('-m', '--mlfc',
                                      type=str,
                                      help=help_str)
            help_str = 'File containing outlier cancer types'
            major_parser.add_argument('-outlier', '--outlier',
                                      type=str,
                                      help=help_str)
            help_str = 'File containing driver list overlap'
            major_parser.add_argument('-do', '--driver-overlap',
                                      type=str,
                                      help=help_str)
            help_str = 'Directory containing methods calls for each method'
            major_parser.add_argument('-k', '--cocall',
                                      type=str, required=True,
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
    elif opts['kind'] == 'list_overlap':
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
    elif opts['kind'] == 'consensus':
        import consensus
        consensus.main(opts)


def cli_main():
    opts = parse_arguments()
    main(opts)


if __name__ == '__main__':
    cli_main()
