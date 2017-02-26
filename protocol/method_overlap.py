import pandas as pd
import os
from collections import Counter
import utils
import config as cfg

# logging
import logging
logger = logging.getLogger(__name__)

def gene_overlap_count(cts, num_methods):
    # now count the number of overlaps
    ovlp_genes = {}
    for nmeth in range(1, num_methods+1):
        mygenes = [x for x in cts if cts[x]==nmeth]
        ovlp_genes[nmeth] = mygenes

    # format the results to a dataframe
    max_genes = max(map(len, ovlp_genes.values()))
    for i in range(1, num_methods+1):
        list_len = len(ovlp_genes[i])
        ovlp_genes[i] = ovlp_genes[i] + [None]*(max_genes-list_len)
    overlap_df = pd.DataFrame(ovlp_genes)

    return overlap_df


def method_overlap_count(signif_genes, config, level='gene'):
    meth_names = cfg.fetch_level_names(config, level=level, exclude=True)
    output_list = []
    for method in meth_names:
        # merge method counts for each method individually
        method_cts = [signif_genes[c][signif_genes[c][method]>0].sum(axis=1)
                      for c in signif_genes]
        myresult = pd.concat(method_cts, axis=1)
        max_method_ct = myresult.max(axis=1)

        # name series according to method
        max_method_ct.name = method

        # append counts to list
        output_list.append(max_method_ct)

    # format into dataframe
    output_df = pd.concat(output_list, axis=1)
    #output_df = pd.DataFrame(output_list, columns=header_names)

    return output_df


def main(opts):
    logger.info('Running method_overlap sub-command . . .')
    config = utils.load_config(opts['config'])

    all_signif_dict = utils.fetch_significant(opts['input_dir'], config, level='gene')
    call_dict = {}
    gene_methods = cfg.fetch_level_names(config, level='gene')
    cancer_types = all_signif_dict[all_signif_dict.keys()[0]].keys()
    for cancer_type in cancer_types:
        # setup a dataframe with the union of significant genes
        union_signif = reduce(lambda x, y: x | y,
                              (set(all_signif_dict[m].get(cancer_type, []))
                               for m in all_signif_dict))
        call_df = pd.DataFrame(index=union_signif)

        # add results from each method
        for meth in gene_methods:
            call_df[meth] = 0
            call_df.loc[all_signif_dict[meth].get(cancer_type, []), meth] = 1

        # add results to dictionary
        call_dict[cancer_type] = call_df

        # save dataframe
        tmp_path = os.path.join(opts['output'], cancer_type+'.gene_overlap.txt')
        call_df.to_csv(tmp_path, sep='\t')

    # eliminate any excluded methods
    if 'exclude' in config:
        exclude_methods = config['exclude']
        cancer_type_list = list(call_dict.keys())
        for c in cancer_type_list:
            for meth in exclude_methods:
                if meth in call_dict[c]:
                    del call_dict[c][meth]

    # Analyze method counts
    method_overlap_df = method_overlap_count(call_dict, config, level='gene')
    tmp_path = os.path.join(opts['output'], 'gene_overlap.txt')
    method_overlap_df.to_csv(tmp_path, sep='\t')
    all_cols = method_overlap_df.columns.tolist()
    meth_melt_df = pd.melt(method_overlap_df, value_vars=all_cols).dropna()
    import IPython ; IPython.embed()
    raise

    # save result
    method_path = os.path.join(opts['output'], 'method_overlap.txt')
    method_ovlp_df.to_csv(method_path, sep='\t', index=False)

    # save plot
    if opts['plot']:
        import plot_data
        plot_path = os.path.join(opts['output'], 'method_overlap.pdf')
        plot_data.method_overlap(method_ovlp_df, plot_path)

    logger.info('Finished method_overlap sub-command.')
