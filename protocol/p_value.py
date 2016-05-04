import pandas as pd
import eval_stats
import utils
import os

def calculate_mlfc(df):
    pvals = df['pvalue']
    mlfc_score = eval_stats.mean_log_fold_change(pvals)
    return mlfc_score


def main(opts):
    # get full data for each method
    full_data = utils.fetch_raw_dataframes(opts['input_dir'])

    # compute MLFC scores
    mlfc_result = {m: calculate_mlfc(full_data[m])
                   for m in full_data}
    mlfc_series = pd.Series(mlfc_result)

    # save output
    mlfc_path = os.path.join(opts['output'], 'mlfc_scores.txt')
    mlfc_series.to_csv(mlfc_path, sep='\t')

    return mlfc_result
