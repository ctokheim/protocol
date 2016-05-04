from random_split import RandomSplit
import pandas as pd
import sys
import argparse
import os

def mymkdir(dir_path):
    try:
        os.mkdir(dir_path)
    except OSError:
        pass


def parse_arguments():
    d = ('Generate splits of the data while respecting the proportion '
         'of each tumor type in the split.')
    parser = argparse.ArgumentParser(description=d)
    help_str = 'Mutations in correct format'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Number of times to perform random split of data'
    parser.add_argument('-n', '--number',
                        type=int, default=1,
                        help=help_str)
    help_str = ('Column name containing the sample ID. Default: automatically'
                ' checks whether "Tumor_Sample_Barcode" or "Tumor_Sample" is a column.')
    parser.add_argument('-s', '--sample-col',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Output directory for mutations split into two files'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # try to make directory
    mymkdir(opts['output'])

    # read in data
    df = pd.read_csv(opts['input'], sep='\t')

    # figure out the sample column
    if opts['sample_col'] is not None:
        samp_col = opts['sample_col']
    if 'Tumor_Sample_Barcode' in df.columns:
        samp_col = 'Tumor_Sample_Barcode'
    elif 'Tumor_Sample' in df.columns:
        samp_col = 'Tumor_Sample'
    else:
        logger.warning('Please specify the column name for the sample ID (--sample-col)')
        sys.exit(1)

    # setup random splitting object
    SAMPLE_RATE = .5  # half splits
    dfg = RandomSplit(df.copy(),
                      col_name=sample_col,
                      sub_sample=SAMPLE_RATE,
                      num_iter=opts['number'])

    # make random splits
    for i, (left_df, right_df) in enumerate(dfg.dataframe_generator()):
        output_dir = os.path.join(opts['output'], 'Iteration_{0}'.format(i))
        mymkdir(output_dir)

        # make sure data is sorted by genes, so no problem with entropy script
        # assuming sorted order
        left_df.sort(columns=['Gene'], inplace=True)
        right_df.sort(columns=['Gene'], inplace=True)

        lout_path = os.path.join(output_dir, 'first.txt')
        left_df.to_csv(lout_path,
                       sep='\t', index=False)
        rout_path = os.path.join(output_dir, 'second.txt')
        right_df.to_csv(rout_path,
                        sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
