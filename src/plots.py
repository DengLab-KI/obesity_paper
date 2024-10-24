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