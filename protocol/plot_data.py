import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns


def multi_qqplot(data, max_pval=1.0):
    with sns.axes_style('ticks'), sns.plotting_context('paper', font_scale=2.5):
        # change dpi
        import matplotlib as mpl
        mpl.rc('savefig', dpi=300)

        # make qq plot for each method
        g = sns.FacetGrid(data, col="method", col_wrap=3,
                          sharey=False, aspect=1.5)
        g.map(qqplot, "p-value")
        plt.tight_layout()

        #set_axes_label(g.fig, 'Theoretical ($log_{10}(p)$)', 'Ovserved ($log_{10}(p)$)', ylab_xoffset=-0.02, xlab_yoffset=-0.02)
        set_axes_label(g.fig, 'Expected p-value', 'Ovserved p-value',
                       ylab_xoffset=-0.03, ylab_yoffset=.62, xlab_yoffset=-0.02)
        g.set_titles('{col_name}')

        # set ylabel
        #g.axes[0].set_ylabel('p-value')
        #g.axes[3].set_ylabel('p-value')
        #g.axes[6].set_ylabel('p-value')

        # set xlim
        for myax in g.axes:
            myax.set_xlim((0, max_pval))
            myax.set_ylim((0, max_pval))


def fix_formatting(ax, title, remove_xlab=True, remove_ylab=False):
    """simple function to set title and remove xlabel"""
    ax.set_title(title)
    if remove_xlab:
        ax.set_xlabel('')
    if remove_ylab:
        ax.set_ylabel('')


def set_axes_label(fig, xlab, ylab,
                   ylab_yoffset=.55, ylab_xoffset=0.04,
                   xlab_yoffset=.04, xlab_xoffset=0.5):
    txt1 = fig.text(xlab_xoffset, xlab_yoffset, xlab, ha='center', size=22)
    txt2 = fig.text(ylab_xoffset, ylab_yoffset, ylab, ha='center', size=22, rotation=90)
    return txt1, txt2


def qqplot(data,
           ax=None,
           log=False, title=None,
           use_xlabel=True, use_ylabel=True,
           **kwargs):
    """qq-plot with uniform distribution"""
    tmp = data.copy()
    tmp.sort()
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp)+1)
    if log:
        log_quant = -np.log10(dist_quant)
        if ax is None:
            plt.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            plt.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        else:
            ax.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            ax.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        # set axis labels
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical (-log10(p))')
            else: ax.set_xlabel('Theoretical (-log10(p))')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed (-log10(p))')
            else: ax.set_ylabel('Observed (-log10(p))')
    else:
        if ax is None:
            plt.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            plt.plot([0, 1], [0, 1], ls="-", color='red')
        else:
            ax.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            ax.plot([0, 1], [0, 1], ls="-", color='red')
            ax.set_ylabel('p-value')
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical p-value')
            else: ax.set_xlabel('Theoretical p-value')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed p-value')
            else: ax.set_ylabel('Observed p-value')
    if title:
        ax.set_title(title)
    sns.despine()
    return ax


def single_method_qqplot(pval_series, filepath):
    """Perform a qq plot for a single method."""
    with sns.axes_style('ticks'), sns.plotting_context('paper', font_scale=1.5):
        qqplot(pval_series, log=True)
        ax = plt.gca()
        #ax.set_xlim((0, 5)); ax.set_ylim((0, 6))
        plt.gcf().set_size_inches((3.5,2.5))
        plt.title('p-value QQ plot')
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()


def single_method_overlap(s, filepath):
    """Plot overlap with cancer gene census and landscapes list.

    Parameters
    ----------
    s : pd.Series
        contains overlap counts
    """
    with sns.axes_style('ticks'), sns.plotting_context('paper', font_scale=1.5):
        ax = sns.barplot(s.index, s, color='black')

        sns.despine()
        plt.xticks(rotation=45, ha='right', va='top')

        ax.set_ylabel('# of overlapping genes')
        plt.gcf().set_size_inches(2, 3)
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()


def single_method_driver_per_sample(driver_per_sample, filepath):
    with sns.axes_style('ticks'), sns.plotting_context('paper', font_scale=1.0):
        driver_list = [pd.DataFrame({'tumor type': ttype, 'type': 'All Driver', '# Per Sample': count})
                    for ttype, count in driver_per_sample.iteritems()]
        ttype_df = pd.concat(driver_list)
        ttype_df['# Per Sample'] = ttype_df['# Per Sample'].astype(int)
        mean_per_sample = ttype_df.groupby('tumor type')['# Per Sample'].mean()
        custom_order = mean_per_sample.sort_values().index
        sns.barplot('tumor type', '# Per Sample', data=ttype_df,
                    order=custom_order, color=sns.color_palette()[0], ci=95)
        plt.ylabel('Avg. # mutated drivers per sample')
        plt.xlabel('')
        sns.despine()
        plt.xticks(rotation=90)
        plt.gcf().set_size_inches((3.5,3.5))
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()
    return custom_order


def single_method_num_drivers_per_type(gene_ttype_ct, custom_order, filepath):
    with sns.axes_style('ticks'), sns.plotting_context('paper', font_scale=1.0):
        sns.barplot(custom_order, gene_ttype_ct[custom_order], color='black')
        plt.xticks(rotation=90)
        sns.despine()
        plt.xlabel('cancer type')
        plt.ylabel('No. of significant genes')
        plt.gcf().set_size_inches((3.5,3.5))
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()


def method_overlap(overlap_df, path):
    # calculate the fractions
    overlap_df['unique'] = 1.0
    mycols = ['1', '2', '>=3']
    overlap_df.loc[:, mycols] = overlap_df.loc[:, mycols].astype(float).div(overlap_df['Total'], axis=0)
    overlap_df['2'] = overlap_df['2'] + overlap_df['>=3']
    overlap_df['1'] = overlap_df['1'] + overlap_df['2']

    # order methods in increasing order
    mymethods = overlap_df.sort_values('1')['Method'].tolist()

    colors = ['white'] + sns.cubehelix_palette(3)[:3]
    for i, col in enumerate(['unique', '1', '2', '>=3']):
        with sns.axes_style('ticks', rc={'xtick.major.pad':-1.0}), sns.plotting_context('talk', font_scale=1.5):
            sns.barplot('Method', col, data=overlap_df,
                        color=colors[i], label=col, order=mymethods,)

            # Finishing touches
            lgd = plt.legend(bbox_to_anchor=(1, .75), loc='upper left',
                             ncol=1, title='Additional\nmethods')
            plt.ylim((0, 1))
            plt.ylabel('Fraction of predicted drivers')
            plt.gca().set_xticklabels(mymethods, rotation=45, ha='right')
            fig = plt.gcf()
            fig.set_size_inches(7, 7)

            # set bar width to 1
            for container in plt.gca().containers:
                plt.setp(container, width=1)
            # remove extra ticks
            plt.gca().get_xaxis().tick_bottom()
            plt.gca().get_yaxis().tick_left()

            # change tick padding
            plt.gca().tick_params(axis='x', which='major', pad=0)

    plt.tight_layout()
    plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def cgc_overlap(cgc_overlap_df, path, list_name='CGC'):
    # Function to label bars
    def autolabel(rects):
        # attach some text labels
        for ii, rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., height+.005, '%s' % (name[ii]),
                     ha='center', va='bottom', size=16)

    # plot barplot
    mymethods = cgc_overlap_df.sort('Fraction overlap ({0})'.format(list_name)).index.tolist()
    name = cgc_overlap_df.ix[mymethods]['# '+list_name].tolist()
    with sns.axes_style('ticks'), sns.plotting_context('talk', font_scale=1.5):
        ax = sns.barplot(cgc_overlap_df.index,
                         cgc_overlap_df['Fraction overlap ({0})'.format(list_name)],
                         order=mymethods, color='black')

        # label each bar
        autolabel(ax.patches)

        # fiddle with formatting
        ax.set_xlabel('Methods')
        ax.set_ylabel('Fraction of predicted drivers\nfound in '+list_name)
        sns.despine()
        plt.xticks(rotation=45, ha='right', va='top')
        plt.gcf().set_size_inches(7, 7)
        # change tick padding
        plt.gca().tick_params(axis='x', which='major', pad=0)

    # save plot
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def mlfc_score(mlfc, path):
    mlfc.sort_values(inplace=True, ascending=False)
    with sns.axes_style('ticks'), sns.plotting_context('talk', font_scale=1.5):
        sns.barplot(mlfc.index, mlfc, color='black')
        sns.despine()
        plt.gca().set_xticklabels(mlfc.index, rotation=45, ha='right')
        plt.ylabel('Mean Log2 Fold Change')
        plt.gcf().set_size_inches(7, 7)
        plt.gca().tick_params(axis='x', which='major', pad=0)
        plt.title('P-value divergence from uniform')

    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def num_signif(num_genes, path):
    num_genes.sort_values(inplace=True, ascending=True)
    with sns.axes_style('ticks'), sns.plotting_context('talk', font_scale=1.5):
        sns.barplot(num_genes.index, num_genes, color='black')
        sns.despine()
        plt.gca().set_xticklabels(num_genes.index, rotation=45, ha='right')
        plt.ylabel('# Significant Genes')
        plt.gcf().set_size_inches(7, 7)
        plt.gca().tick_params(axis='x', which='major', pad=0)

    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def consistency(consis_df, depth, path):
    # get col names and format order
    mean_col = 'TopDrop {0} overlap mean'.format(depth)
    sem_col = 'TopDrop {0} overlap sem'.format(depth)
    method_order = consis_df.sort(mean_col).index.tolist()

    # make plot
    with sns.axes_style('ticks'), sns.plotting_context('talk', font_scale=1.5):
            myax = sns.barplot(method_order, consis_df[mean_col].ix[method_order],
                               order=method_order, color=sns.xkcd_rgb["grey"])
            myax.errorbar(np.arange(len(method_order)), consis_df[mean_col].ix[method_order],
                          yerr=consis_df[sem_col].ix[method_order],
                          fmt=None, ecolor=sns.xkcd_rgb["black"],
                          elinewidth=2, capsize=9, capthick=2)
            myax.set_xticklabels(method_order, rotation=45, ha='right')
            myax.yaxis.set_ticks_position('left')
            myax.xaxis.set_ticks_position('bottom')
            myax.set_ylabel('TopDrop Consistency ({0})'.format(depth))
            myax.set_ylim((0, 1.0))

            sns.despine()
            ax = plt.gca()
            plt.gca().tick_params(axis='x', which='major', pad=0)
            plt.title('Consistency')

    # save output
    plt.gcf().set_size_inches(7, 7)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
