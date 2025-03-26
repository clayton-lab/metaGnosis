import sys
import pandas as pd
import numpy as np
import itertools as it
import collections
import pathlib

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

# Disables the stupid pandas SettingWithCopyWarning
pd.options.mode.chained_assignment = None

### Delete these after testing is done
#read_sample = ['ABX-CJ-MANG-14', 'ABX-CJ-MANG-29', 'ABX-CJ-MANG-43', 'ABX-CJ-MANG-56', 'ABX-CJ-IRIS-14', 
#               'ABX-CJ-IRIS-29', 'ABX-CJ-IRIS-43', 'ABX-CJ-IRIS-56']
#
#coverages = ['output/mapping/minimap2/coverage_tables/genes/ABX-CJ-IRIS-29_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-IRIS-56_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-IRIS-14_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-IRIS-43_gene_coverage.txt',
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-MANG-29_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-MANG-56_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-MANG-14_gene_coverage.txt', 
#             'output/mapping/minimap2/coverage_tables/genes/ABX-CJ-MANG-43_gene_coverage.txt']
#
#lengths = ['output/mapping/minimap2/lengths/genes/ABX-CJ-IRIS-29_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-IRIS-56_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-IRIS-14_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-IRIS-43_gene_lengths.txt',
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-MANG-29_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-MANG-56_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-MANG-14_gene_lengths.txt', 
#           'output/mapping/minimap2/lengths/genes/ABX-CJ-MANG-43_gene_lengths.txt']
#
#derep_clust = "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_cluster.tsv"
#annot = 'output/annotate_bins/annotate_bin_pathways/merged_annotations/annotations.tsv' 
#pathways = 'output/annotate_bins/annotate_bin_pathways/annotated_pathways/metabolism_summary.xlsx'
#gene_info_outfile = 'output/quant_bins/quantified_gene_info.txt'
#gene_quant_outfile = 'output/quant_bins/quantified_genes.txt'

# Wildcards, params, input and output files are read in
read_sample = snakemake.params.get('read_sample')
coverages = snakemake.input.get('coverages')
lengths = snakemake.input.get('lengths')
annot = snakemake.input.get('function', '')
pathways = snakemake.input.get('pathways', '')
derep_clust = snakemake.input.get('derep_clust', '')
gene_quant_outfile = snakemake.output.get('quant_genes', '')
gene_info_outfile = snakemake.output.get('gene_info', '')

mod_annot = pathlib.Path(annot).parent.joinpath('updated_annotations.tsv')
pd_args = {'sep': '\t', 'engine': 'c', 'header': None}

# Gene lengths, Bin_IDs, sample-wise coverage, ANI-cluster mapping, and gene info file(s), are aggregated.
gene_length = pd.concat([pd.read_csv(df, **pd_args, names=['ORF_ID', 'Length']) for df in lengths], axis=0, copy=False)\
            .drop_duplicates(keep='first').set_index('ORF_ID')
gene_coverage = pd.concat([pd.read_csv(df, **pd_args, index_col=0, names=['ORF_ID', read_samp]) for df, read_samp in zip(coverages, read_sample)], axis=1, copy=False)
annotations = pd.read_csv(annot, sep='\t', engine='c', index_col=0, dtype='object').rename_axis('ORF_ID').rename(columns={'fasta': 'Bin_ID', 'scaffold': 'Contig_ID'})

# Genes were clustered with Average Nucleotide Identity (ANI), so cluster representatives (Clust_ID) and 
# their constituent members (ORF_ID) are later merged with ORF-wise gene annotations
gene_clust = pd.read_csv(derep_clust, **pd_args, names=['Clust_ID', 'ORF_ID']).drop_duplicates().set_index('ORF_ID')

# A modified annotation file is created with the cluster id in case future analyses require the full annotation set
annotations.join(gene_clust).to_csv(mod_annot, index=True, sep='\t')

# The output of DRAM distill is parsed, providing pathway info about each gene
dist_dfs = [df.assign(sheet=([sname] * df.shape[0])).drop_duplicates() for sname, df in pd.read_excel(pathways, sheet_name=None).items()]
distill = pd.concat(dist_dfs, axis=0, copy=False).reset_index(drop=True)
gene_info = ['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader', 'specific_reaction', 'Notes', 'EC', 'oxygen', 'alt_ids']
gene_samps = distill.columns[~distill.columns.isin(gene_info)]

# Rows that didn't match any samples are removed
distill = distill.loc[distill[gene_samps].dropna(how='all').index, :].reset_index(drop=True)


# Commas are removed from some Gene_IDs
distill_mult_geneids = distill[distill['gene_id'].str.contains(',')]
distill_mult_geneids.loc[:, 'gene_id'] = distill_mult_geneids['gene_id'].str.replace(r',*\s*EC:(?:\d+\.*){4}\s*', '', regex=True)
distill_mult_geneids.loc[:, 'gene_id'] = distill_mult_geneids['gene_id'].str.strip(' ').str.split(',')
distill_mult_geneids['alt_ids'] = distill_mult_geneids['gene_id'].str[1:]
distill_mult_geneids.loc[:, 'gene_id'] = distill_mult_geneids['gene_id'].str[0]
distill = pd.concat([distill[~distill['gene_id'].str.contains(',')], distill_mult_geneids], axis=0, copy=False)

# The EC column is filled in for genes that have the EC in the description
#TODO: Optimize this section better to avoid processing the ECs differently for found and non-found ECs
no_ec = distill[distill['EC'].isna()]
found_ecs = no_ec['gene_description'].str.extractall(r'(?P<ECs>(?:\(|\[)?(?:EC:?)?\s?(?:(?:\d|-)+\.){3}(?:\d|-)+(?:\)|]|))')

found_ecs = found_ecs['ECs'].str.strip('[( )]').str.replace(r'^(?:EC)?\W?', 'EC:', regex=True)\
        .unstack(sort=False, fill_value='').agg(lambda row: np.unique(row[row.astype(bool)]), axis=1)
distill.loc[~distill['EC'].isna(), 'EC'] = distill.loc[~distill['EC'].isna(), 'EC'].str.replace(r',?(?:EC)\W?', ',EC:', regex=True)\
        .apply(lambda row: np.sort(row.split(',')[1:]))

distill.loc[found_ecs.index, 'EC'] = found_ecs.values
#distill.loc[:, 'EC'] = distill['EC'].str.join(',')

# The distill output is de-replicated by grouping to the level of gene_id. Descriptions/Bin_IDs/etc that differ within a gene_id are joined in a list
agg_dict = {col: tuple if col in gene_info else lambda c: np.unique(c)[0] for col in distill.columns[distill.columns != 'gene_id']}
distill.loc[:, gene_samps] = distill.loc[:, gene_samps].apply(lambda col: col.str.split(','), result_type='expand', axis=0)
distill_derep = distill.groupby('gene_id', sort=False, as_index=True).agg(agg_dict, axis=0)

## Individual ORF_IDs that mapped to each gene are compiled in a dict, allowing total copy number of each gene to be counted later
orf2gene_dict = collections.defaultdict(list)
for tup in distill_derep[gene_samps].fillna('NA').itertuples(name=None):
    for val in it.chain.from_iterable(it.filterfalse(lambda item: item == 'NA', tup[1:])):
        orf2gene_dict[val].append(tup[0])

orf2gene = pd.DataFrame.from_dict(orf2gene_dict, orient='index').rename_axis('ORF_ID').agg(lambda row: list(row.dropna()), axis=1).rename('gene_id')

# The distill pathway information is reformatted after the ORF_IDs were extracted
distill_info = distill_derep[gene_info[1:]].reset_index(names='gene_id')
distill_info = distill_info.explode(gene_info[1:]).explode('EC').explode('alt_ids').groupby('gene_id', sort=False).agg(lambda c: tuple(c.dropna().unique()))

# All the files from earlier are finally merged into a single dataframe. Similar to before, the dataframe is aggregated to
# the gene_id level, and since some genes can have multiple annotations/ORFs/Bins (and even Clust_IDs), those are merged into lists
gene_stats = pd.concat([gene_clust, annotations[['Bin_ID', 'ko_id']], orf2gene, gene_coverage, gene_length], axis=1, copy=False).reset_index()
gene_agg_dict = {col: tuple if col in ['ORF_ID', 'Bin_ID', 'ko_id', 'gene_id'] else lambda c: list(c)[0] for col in gene_stats.columns[gene_stats.columns != 'Clust_ID']}
gene_stats_agg = gene_stats.reset_index().explode('gene_id').groupby('Clust_ID', sort=False, as_index=False).agg(gene_agg_dict)

# Genes that had coverage but weren't assigned an identifier are removed
gene_filt = gene_stats_agg.explode(['gene_id']).loc[~gene_stats_agg['gene_id'].explode().isna()].index.unique()
gene_stats_filt = gene_stats_agg.loc[gene_filt, :].explode(['ORF_ID', 'Bin_ID', 'ko_id', 'gene_id'])

# Gene copy number is calculated by counting the total number of ORFs that map to each gene across samples
cpy_nums = gene_stats_filt['gene_id'].value_counts().rename('Copy_Num')

sum_cov = gene_stats_filt.groupby('gene_id', sort=False)[gene_coverage.columns.to_list() + gene_length.columns.to_list()].agg('sum').join(cpy_nums)

# If desired, the complete ORF_ID/Clust_ID list of each gene can be written to the output by including their columns here. Those were removed
# from the final output because they make the file 10x larger
final_df_list = [gene_stats_filt.loc[:, ['gene_id', 'Bin_ID', 'ko_id']].groupby('gene_id', sort=False).agg(lambda c: tuple(c.dropna().unique())),
                 distill_info, sum_cov]

# Some KOs that were incomplete are filled in from the alt_ids column, which is then removed
# Only the extra IDs from the CAMPER dataset have this problem, so overwriting valid KOs shouldn't be possible
gene_stats_final = pd.concat(final_df_list, copy=False, axis=1)
missing_ko = gene_stats_final.loc[~gene_stats_final['alt_ids'].isin(gene_stats_final['ko_id'])]
gene_stats_final.loc[missing_ko.index, 'ko_id'] = gene_stats_final.loc[missing_ko.index, 'alt_ids']
gene_stats_final = gene_stats_final.drop(columns='alt_ids').T
gene_stats_final.to_csv(gene_info_outfile, index=True, sep='\t')

# Total coverage and gene length are summed across samples, divided by copy number, and TPM is calculated
norm_cov = sum_cov.T.div(sum_cov.T.loc['Copy_Num'], axis=1)
rpk = norm_cov.div((norm_cov.loc['Length'] / 1000), axis=1).loc[read_sample]
scale_factor = rpk.sum(axis=1)/1000000
tpm = rpk.div(scale_factor, axis=0)
tpm.to_csv(gene_quant_outfile, index=True, sep='\t')

# Alternative way of doing gene normalization for a transposed dataframe
#gene_sampsT = gene_samps.T
#norm_covT = gene_sampsT.div(gene_sampsT['Copy_Num'], axis=0)
#rpkT = norm_covT.div((norm_covT['Length'] / 1000), axis=0)[read_sample]
#scale_factorT = rpkT.sum(axis=0)/1000000
#tpmT = (rpkT / scale_factorT)
