import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import IPython

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
            if ax is None: plt.xlabel('Theoretical ($-log_{10}(p)$)')
            else: ax.set_xlabel('Theoretical ($-log_{10}(p)$)')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed ($-log_{10}(p)$)')
            else: ax.set_ylabel('Observed ($-log_{10}(p)$)')
    else:
        if ax is None:
            plt.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            plt.plot([0, 1], [0, 1], ls="-", color='red')
        else:
            ax.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            ax.plot([0, 1], [0, 1], ls="-", color='red')
            ax.set_ylabel('p-value')
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical Quantile')
            else: ax.set_xlabel('Theoretical Quantile')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed Quantile')
            else: ax.set_ylabel('Observed Quantile')
    if title:
        ax.set_title(title)
    sns.despine()


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
                             ncol=1, title='Overlapping\nmethods')
            plt.ylim((0, 1))
            plt.ylabel('Fraction Overlap (methods)')
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


def cgc_overlap():
    overlaps = [ad_list, oncofm_list, oncofml_list, oncoclust_list, music_list, tuson_list, mutsig_list, twentyPlus_list]
    #name = ['{0}/{1}'.format(*l) for l in overlaps]
    name = ['{0}'.format(*l) for l in overlaps]

    # Function to label bars
    def autolabel(rects):
        # attach some text labels
        for ii, rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x()+rect.get_width()/2., height+.005, '%s' % (name[ii]),
                    ha='center', va='bottom', size=16)

    # get data ready to plot
    mymethods = ['ActiveDriver', 'OncodriveFM', 'OncodriveFML', 'OncodriveClust', 'MuSiC', 'TUSON', 'MutsigCV', '20/20+']
    plot_df = pd.DataFrame({'# CGC genes': [ad_list[0], oncofm_list[0], oncofml_list[0], oncoclust_list[0], music_list[0],
                                            tuson_list[0], mutsig_list[0], twentyPlus_list[0]],
                            'total': [ad_list[1], oncofm_list[1], oncofml_list[1], oncoclust_list[1], music_list[1],
                                    tuson_list[1], mutsig_list[1], twentyPlus_list[1]]},
                        index=mymethods)
    plot_df['Fraction Overlap'] = plot_df['# CGC genes'].astype(float) / plot_df['total']
    plot_df.to_clipboard()

    # plot barplot
    with sns.axes_style('ticks'), sns.plotting_context('talk', font_scale=1.5):
        ax = sns.barplot(plot_df.index, plot_df['Fraction Overlap'], order=mymethods, color='black')

        # label each bar
        autolabel(ax.patches)

        # fiddle with formatting
        ax.set_xlabel('Methods')
        ax.set_ylabel('Fraction Overlap (CGC)')
        #ax.set_title('Overlap = CGC Genes')
        sns.despine()
        plt.xticks(rotation=45, ha='right', va='top')
        #ax.set_aspect(11)
        plt.gcf().set_size_inches(7, 7)
        plt.ylim((0, .5))
        # change tick padding
        plt.gca().tick_params(axis='x', which='major', pad=0)
