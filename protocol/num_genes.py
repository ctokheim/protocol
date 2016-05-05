import pandas as pd
import utils
import os

# logging
import logging
logger = logging.getLogger(__name__)

def main(opts):
    logger.info('Running num_signif sub-command . . .')
    config = utils.load_config(opts['config'])

    # get the significant genes for each method
    signif_dict = utils.fetch_significant_genes(opts['input_dir'],
                                                opts['qvalue'],
                                                config)
    num_methods = len(signif_dict)

    # count the number of significant for each method
    num_signif = {m: len(signif_dict[m]) for m in signif_dict}
    num_series = pd.Series(num_signif)

    # save result
    out_path = os.path.join(opts['output'], 'num_significant.txt')
    num_series.to_csv(out_path, sep='\t')

    logger.info('Finished num_signif sub-command.')
    return num_series
