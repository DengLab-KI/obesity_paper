# %%
import os
import pandas as pd
import scanpy as sc
from anndata import AnnData
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from matplotlib.pyplot import rc_context
import pygame
from matplotlib import cm
from slingshot import Slingshot


# %%
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_theme(style="ticks", rc=custom_params)

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
sc.logging.print_header()
sc.set_figure_params(dpi=900, color_map='viridis_r')
%matplotlib inline

os.chdir('/mnt/data/hong/2022/human_placenta/')


# %%
matplotlib.rcParams['figure.figsize'] = [5, 5]

# %%
sct= sc.read_h5ad('output/10x_h5/h5ad/stb_v8.h5ad')


# %%
sct.X = sct.layers['raw']
sc.experimental.pp.recipe_pearson_residuals(sct, n_top_genes=2000, batch_key='group')

# %%
sc.pl.pca_variance_ratio(sct)

# %%
sc.pp.neighbors(sct, n_neighbors=15, n_pcs=10, random_state=404)
sc.tl.umap(sct, min_dist=0.2, spread=0.5)
sc.tl.leiden(sct, resolution=0.3, random_state=404)

# %%
sc.pl.umap(sct,  color=['leiden', 'C_scANVI', 'group'], title='STB clusters test')

# %%
top_genes = ['ITGB1', 'COL4A2', 'COL4A1', 'FN1', 'TIMP3', 'FBN1', 'ADGRL3', 'GRB14','FUT8', 'MGAT5', 'STT3B', 'MMP11', 'THSD4']

# %%
top_genes = ['MEG3', 'H19', 'ARL15', 'FN1', 'SEMA6D', 'ADGRL3', 'GRB14','PKIB', 'RHOBTB1', 'HSD17B1', 'KLRD1', 'ADGRG6']

# %%
sc.pl.heatmap(sct, top_genes, groupby='maturation_inhouse_three', layer='log_norm', use_raw=False, swap_axes=True, standard_scale='var', cmap='Reds', save='STB_top_test.pdf')

# %%
from collections import defaultdict
maturation = {'immature': ['nascent', 'premature1-a', 'mature2-b'],
            'mature': ['premature1-b', 'mature1-a', 'mature1-b', 'mature1-c', 'mature2-a']}
annotate_inv = defaultdict()
for k, v in maturation.items():
    for vi in v:
        annotate_inv[vi] = k
sct.obs['maturation_my'] = sct.obs['C_scANVI'].map(annotate_inv)

# %%
from collections import defaultdict
maturation = {'STB-a': ['4'],
    'STB-b': ['0', '7'],
            'STB-c': ['1', '2', '3', '5', '6']}
annotate_inv = defaultdict()
for k, v in maturation.items():
    for vi in v:
        annotate_inv[vi] = k
sct.obs['maturation_inhouse_three'] = sct.obs['leiden'].map(annotate_inv)

# %%
sct_time = pd.read_csv('data/trajectory/sct_time_df.tsv', sep='\t')
sct.obs = pd.merge(sct.obs, sct_time, left_ibndex=True, right_index=True)

# %%
sct_resample = sc.pp.subsample(sct, fraction=1, n_obs=None, random_state=404, copy=True)
sc.pl.umap(sct_resample, color='maturation_inhouse_three', palette=['#ffbb78', '#98df8a', '#ff9896'],save='figure1k')

# %%
sct_resample = sc.pp.subsample(sct, fraction=1, n_obs=None, random_state=404, copy=True)
sc.pl.umap(sct_resample, color=['maturation_inhouse_three','Pseudotime_sct_x', 'group'], cmap='YlOrRd', save='STB_maturation_my_time_fig2.pdf')

# %%
sct_resample.obs['n_genes'] = sct_resample.obs['n_genes'].astype(int)
sct_resample.obs['log_n_genes'] = np.log(sct_resample.obs['n_genes'])

# %%
# import matplotlib.colors
# cvals  = [-2., 0,  2]
# colors = ["dimgray","salmon","red"]

# norm=plt.Normalize(min(cvals),max(cvals))

sc.pl.umap(sct_resample, color=['maturation_inhouse_three', 'TENM3', 'PDE4D', 'ADAMTS6'], cmap = 'YlOrRd', layer='log_norm', save='STB_maturation_my_genes_leiden_three_ngenes_fig2.pdf')

# %%
# fig1n_genes = pd.read_csv('data/Figure1n_SourceData.txt', sep='\t')['STB-a_markers_EMT']
sc.pl.dotplot(sct, fig1n_genes, groupby='maturation_inhouse_three', layer='log_norm', cmap='YlOrRd', standard_scale='var', use_raw=False, categories_order=['STB-a', 'STB-b', 'STB-c'], save='figure1n')

# %%
## violin of ADAAMTS6 and GRB14
sc.pl.violin(sct, keys=['ADAMTS6', 'GRB14'], groupby='maturation_inhouse_v7', layer='log_norm', gitter=False, stripplot=False, use_raw=False, save='STB_maturation_my_genes_violin.pdf')

# %%
sc.pl.umap(sct_resample, color=['PSG8', 'SH3TC2', 'PAPPA', 'FLT1'], cmap = 'YlOrRd', layer='log_norm', save='STB_maturation_ng_genes.pdf')

# %%
sc.pl.dotplot(sct, ['SLC2A1', 'SLC2A3', 'SLC2A8', 'SLC2A12'], groupby='group')

# %%
# from collections import Counter
go = pd.read_csv('data/n_reg_transport.tsv', sep='\t')


# %%
sct.obs['n_genes_by_counts'] = sct.obs['n_genes_by_counts'].astype('int64')
sct.obs['log1p_ngenes'] = np.log1p(sct.obs['n_genes_by_counts'])


# %%
sc.tl.score_genes(sct, go.genesymbol.tolist(), score_name='n_transport_score', use_raw=False, random_state=404)


# %%
sc.pl.violin(sct, ['n_transport_score'], stripplot=False, groupby='group', rotation=45)


# %%
sc.pp.calculate_qc_metrics(sct, layer='raw', inplace=True, percent_top=[10])


# %%
sct_full = sct[sct.obs['n_genes_by_counts']>100, :].copy()
sct_full.X = sct_full.layers['raw'].copy()

sc.pp.scale(sct_full, max_value=10)


# %%
sc.pl.umap(sct, color='SLC30A8', layer='log_norm')

# %%
sc.pl.dotplot(sct, ['CGA', 'PSG4', 'PSG11', 'CYP19A1', 'STS', 'KISS1'], use_raw=False, groupby='group', mean_only_expressed=False, standard_scale='var')

# %%
sc.tl.rank_genes_groups(sct, groupby='group', method='wilcoxon',
                        layer='log_norm', use_raw=False, reference='Normal_AGA')
with rc_context({'figure.figsize': (5, 1.5)}):
    sc.pl.rank_genes_groups_dotplot(
        sct, n_genes=5, groups=['Obese_AGA', 'Obese_LGA'], save='n_transport_score_dot_sct.pdf')


# %%
def rgba_hex_colors_to_rgb(adata):
    for key in adata.uns.keys():
        if key.endswith('colors'):
            adata.uns[key] = np.array(
                [(c if len(c) <= 7 else c[:-2]) for c in adata.uns[key]])
    return adata


def write_adata(adata, file_name):
    # adata.obs as str
    adata.obs = adata.obs.astype(str)

    adata = rgba_hex_colors_to_rgb(adata)
    adata.write(f'{file_name}.h5ad', compression='gzip')


write_adata(sct, 'output/10x_h5/h5ad/SCT_v0.5')


# %%
sct.obs["n_genes_by_counts"] = pd.to_numeric(sct.obs["n_genes_by_counts"])


# %%
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(sct, color='n_genes_by_counts', add_outline=True, size=20,
               legend_fontsize=12, legend_fontoutline=2, frameon=False,
               title='recluster SCT', save='SCT_recluster_n_genes_by_counts.pdf')



