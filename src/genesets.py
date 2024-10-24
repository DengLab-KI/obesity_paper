def get_set_genes(set_df, set_name):
    gene_str= set_df.loc[set_df.set_name==set_name, 'genes'].to_list()
    gene_list = gene_str[0].split(',')
    return gene_list