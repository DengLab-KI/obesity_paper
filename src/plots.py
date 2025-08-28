import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def function_sets(cell_type,function_set, function_name, mode):
    ## read in the data
    # ad_zscore = sc.pp.scale(ad_clean[ad_clean.obs.final_celltypes==cell_type], max_value=10, zero_center=True, layer = 'log_norm', copy = True)
    ad = ad_clean[ad_clean.obs.final_celltypes==cell_type]
    df = pd.read_csv(f'output/DEGs/final_negbinom_all/score_tsv/{cell_type}.tsv', sep='\t')
    ## order df by score
    df.sort_values('score' , ascending=False, inplace=True)
    common_set = df.query("convergence=='high'&p<0.05")['Unnamed: 0']
    divergent_set = df.query("convergence=='low'&p<0.05")['Unnamed: 0']
    common_function = [element for element in common_set if element in function_set]
    divergent_function = [element for element in divergent_set if element in function_set]
    print(len(common_function))
    print(len(divergent_function))
    if mode == 'common':
        function2plot = common_function
    if mode == 'divergent':
        function2plot = divergent_function
    elif mode=='both':
        function2plot = common_function+divergent_function
    if len(function2plot)>0:
        print(function2plot)
        function_name = function_name.replace('/', '_')
        sc.pl.dotplot(ad, function2plot, groupby='group', dendrogram=False, layer = 'log_norm', show=False, swap_axes=True, cmap='YlOrRd', standard_scale='var', save=f'fig2_{cell_type}_{function_name}_standard.pdf', figsize=(2.5, len(function2plot)/2), use_raw=False, categories_order=['Normal_AGA', 'Obese_AGA', 'Obese_LGA'], dot_max
=1, dot_min=0.1)

## stackbar

def plot_stacked_bar_mixture_group(
    df: pd.DataFrame,
    group_order=None,
    ax=None,
    colors=None,
    savepath=None,
    figsize=(6, 4),
):
    """
    Create a stacked bar chart for STB and CTB by Mixture, ordered by Group.

    Expects columns: 'Mixture', 'Group', 'STB', 'CTB'.
    Default group order: ['Medium control', 'Adipose spheroid'].
    """
    if group_order is None:
        group_order = ['Medium control', 'Adipose spheroid']

    data = df.copy()
    data['Group'] = pd.Categorical(data['Group'], categories=group_order, ordered=True)
    data = data.sort_values(['Group', 'Mixture'])

    x_labels = data['Mixture'].tolist()
    stb_values = data['STB'].to_numpy()
    ctb_values = data['CTB'].to_numpy()
    x_positions = np.arange(len(x_labels))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if colors is None:
        colors = {'STB': '#7BCCC4', 'CTB':'#1474B2'}

    ax.bar(x_positions, stb_values, label='STB', color=colors['STB'])
    ax.bar(x_positions, ctb_values, bottom=stb_values, label='CTB', color=colors['CTB'])

    ax.set_xlabel('Mixture')
    ax.set_ylabel('Proportion')
    ax.set_ylim(0, 1)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.legend(title='Component')
    ax.margins(x=0.02)
    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    return ax
