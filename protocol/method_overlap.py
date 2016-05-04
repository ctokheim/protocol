import pandas as pd
import os
from collections import Counter
import utils


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


def method_overlap_count(signif_genes, gene_counts):
    output_list = []
    for method in signif_genes:
        # list of num ovlps
        method_cts = [gene_counts[g]-1 for g in signif_genes[method]]

        # get overlap counts
        num_uniq = len([x for x in method_cts if x==0])
        num_one = len([x for x in method_cts if x==1])
        num_two = len([x for x in method_cts if x==2])
        num_three = len([x for x in method_cts if x==3])
        num_total = len(method_cts)

        # append result
        tmp_list = [method, num_total, num_uniq, num_one, num_two, num_three]
        output_list.append(tmp_list)

    # format into dataframe
    header_names = ['Method', 'Total', 'unique', '1', '2', '>=3']
    output_df = pd.DataFrame(output_list, columns=header_names)

    return output_df


def main(opts):
    config = utils.load_config(opts['config'])

    # get the significant genes for each method
    signif_dict = utils.fetch_significant_genes(opts['input_dir'],
                                                opts['qvalue'],
                                                config)
    num_methods = len(signif_dict)

    # count how many time each gene is significant
    gene_cts = Counter([g for method in signif_dict for g in signif_dict[method]])
    gene_overlap_df = gene_overlap_count(gene_cts, num_methods)

    # calculate the number of overlaps for each method
    method_ovlp_df = method_overlap_count(signif_dict, gene_cts)

    # save result
    gene_path = os.path.join(opts['output'], 'gene_overlap_counts.txt')
    gene_overlap_df.to_csv(gene_path, sep='\t', index=False)
    method_path = os.path.join(opts['output'], 'method_overlap_counts.txt')
    method_ovlp_df.to_csv(method_path, sep='\t', index=False)
