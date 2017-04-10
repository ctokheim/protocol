"""
File: consensus.py
Author: Collin Tokheim
Email: ctokheim@jhu.edu
Github: ctokheim
Description: Creates a driver gene consensus amongst methods
"""
import pandas as pd
from collections import Counter
import csv
import os
import glob

def read_outliers(mypath):
    """Read which cancer types are outliers for each method."""
    with open(mypath) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        outliers = {line[0]: line[1].split(',') for line in myreader}
    return outliers


def main(opts):
    # read in mlfc/outliers
    if opts['mlfc']:
        mlfc_df = pd.read_table(opts['mlfc'], index_col=0)
        mlfc_med = mlfc_df.median(axis=1)
        mlfc_q1 = mlfc_df.quantile(.25, axis=1)
        #import IPython ; IPython.embed()
        #mlfc_bad = mlfc_df.apply(lambda x, m: x>=m, args=(1,), axis=1)
    #rename_dict = {'Unnamed: 0': 'method'}
    if opts['driver_overlap']:
        driver_ovlp_df = pd.read_table(opts['driver_overlap'], index_col=0) # .rename(columns=rename_dict)
        driver_meds = driver_ovlp_df.groupby('CODE')['Fraction overlap (Custom List)'].median()
        driver_q3 = driver_ovlp_df.groupby('CODE')['Fraction overlap (Custom List)'].quantile(.75)
    if opts['outlier']:
        outliers = read_outliers(opts['outlier'])
        counts = Counter([c for k in outliers for c in outliers[k]])
        # outlier if majority of methods agree
        outlier_cancer_types = [c
                                for c in counts
                                if counts[c]>=(len(outliers)/2)]

    # iterate over all cancer types/PANCAN
    score_list = []
    pattern = os.path.join(opts['cocall'], '*.overlap.txt')
    # temporary hack
    limited_methods = [#'NetBox-snv',
                       #'OncoIMPACT-snv',
                       #'DriverNet-snv',
                       'HotSpot3D',
                       ]
    bad_method_list = []
    for path in glob.glob(pattern):
        # file names are CANCER_TYPE.gene_overlap.txt
        cancer_type = os.path.basename(path).split('.')[0]

        # read in method overlap dataframe
        df = pd.read_table(path, index_col=0)

        # add weighting unless said no by user
        if not opts['not_weighted']:
            # figure out bad mlfc
            #bad_mlfc = mlfc_df.ix[cancer_type] > 0
            bad_mlfc = mlfc_df.ix[cancer_type] >= mlfc_med[cancer_type]

            # figure out if an outlier
            mydict = {}
            for m in outliers:
                mydict[m] = cancer_type in outliers[m]
            outlier_series = pd.Series(mydict)

            # figure if overlap problem
            ovlp_series = driver_ovlp_df[driver_ovlp_df['CODE']==cancer_type]['Fraction overlap (Custom List)'] <= driver_meds.ix[cancer_type]

            # combine the two criteria
            if cancer_type in outlier_cancer_types:
                bad_method = pd.Series({i: True for i in bad_mlfc.index})
            elif cancer_type == 'PANCAN':
                tmp = (bad_mlfc.astype(int) + outlier_series.astype(int) + ovlp_series.astype(int))
                bad_method = tmp > 0
            else:
                tmp = (bad_mlfc.astype(int) + outlier_series.astype(int) + ovlp_series.astype(int))
                bad_method = tmp > 1
                bad_method[limited_methods] = tmp[limited_methods] > 0
            good_methods = bad_method[~bad_method].index.tolist()
            #good_methods = [g for g in good_methods if g not in hack_remove]

            # make multiplier adjustment
            df.loc[:, good_methods] = 2 * df[good_methods]

            bad_method.name = cancer_type
            bad_method_list.append(bad_method)

        # calculate score
        score = df.sum(axis=1)
        score.name = cancer_type
        score_list.append(score)

    # save results
    full_score_df = pd.concat(score_list, axis=1)
    out_path = os.path.join(opts['output'], 'consensus.txt')
    full_score_df.to_csv(out_path, sep='\t')
    bad_meth_df = pd.concat(bad_method_list, axis=1)
    out_path = os.path.join(opts['output'], 'less_reliable.txt')
    bad_meth_df.to_csv(out_path, sep='\t')
    import IPython ; IPython.embed()
