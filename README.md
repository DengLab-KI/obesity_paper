# Notebooks replicating the figures in Obesity manuscript

## Fig.1

1. pre-process and cluster: ipynb/process_10x-cellbender.ipynb
    This is wrapped into a snakemake pipeline.

2. fig 1i, cell type markers dot plot: ipynb/marker_dotplots.ipynb

## Fig.2

1. Syncytialtrophoblast states
    downstream analysis and visualization: ipynb/stb_recluster.py

2. Merge cell type and states annotation: ipynb/merge_annotation.ipynb

3. Cell type propotional in Fig.2
    - test: ipynb/scCODA.ipynb
    - visualization: ipynb/celltype_general_stats.ipynb

3. Venn
    draw venn diagram: ipynb/venn.ipynb

## Fig.2-4

1. C-score line plot
    ipynb/score_lineplot.ipynb

2. GSEA analysis
    - GSEA: ipynb/decouple_bp_hallmark.ipynb
    - adjust p-values by BH method: Rmd/adjust_gsea.Rmd

3. Dotplots of the genes involved in given pathways
    ipynb/dotplots pathway genes.ipynb

## Pie plot in Fig.2,Fig.3, and Extended Data Fig.2

ipynb/pie.ipynb

## Fig 5, Extended Data Fig 3a

1. generate the adjency matrix for Extended Data Fig 3a: ipynb/de_liana.ipynb
2. network analysis for Fig.5: ipynb/network_rl.ipynb

## Fig.6

1. intersect of all sets: Rmd/intersect_TO_tissue.Rmd
2. Venn: draw venn diagram: ipynb/venn.ipynb
