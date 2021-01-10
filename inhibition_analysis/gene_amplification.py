import pandas
cn = pandas.read_csv('C:/Users/james/Desktop/Sumoylation_Analysis/data/CCLE_gene_cn.csv')
info = pandas.read_csv('C:/Users/james/Desktop/Sumoylation_Analysis/data/sample_info.csv')
cn.head()
info.head()

cn = cn.rename(columns = {'Unnamed: 0':'ID'})

myc_df = cn[['ID', 'MYC (4609)']]
myc_df.head()

info = info.rename(columns = {'DepMap_ID':'ID'})
info = info[['ID', 'lineage', 'stripped_cell_line_name']]


