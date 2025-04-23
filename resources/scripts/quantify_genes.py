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

# Wildcards, params, input and output files are read in
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
#bin_annot = 'output/annotate_bins/annotate_bins/bin_annotations.tsv'
#contig_annot = 'output/annotate_bins/annotate_contigs/annotations.tsv'
#pathways = 'output/annotate_bins/annotate_contig_pathways/metabolism_summary.xlsx'
#gene_info_outfile = 'output/quant_bins/quantified_gene_info.txt'
#gene_quant_outfile = 'output/quant_bins/quantified_genes.txt'

# Wildcards, params, input and output files are read in
read_sample = snakemake.params.get('read_sample')
coverages = snakemake.input.get('coverages')
lengths = snakemake.input.get('lengths')
bin_annot = snakemake.input.get('bin_annotation', '')
contig_annot = snakemake.input.get('contig_annotation', '')
pathways = snakemake.input.get('pathways', '')
gene_bincount_outfile = snakemake.output.get('gene_bincounts', '')
gene_quant_outfile = snakemake.output.get('quant_genes', '')
gene_info_outfile = snakemake.output.get('gene_info', '')

# Contig-level and bin-level gene annotations files are read
bin_annotations = pd.read_csv(bin_annot, engine='c', sep='\t', dtype='object', index_col=0).rename_axis('ORF_ID').reset_index()\
        .set_index('Assembly_ORFID').rename(columns={'fasta': 'Bin_ID', 'bin_taxonomy': 'Taxonomy', 'taxid': 'Tax_ID'})
cont_annotations = pd.read_csv(contig_annot, engine='c', sep='\t', dtype='object', index_col=0)[['fasta', 'ORF_ClustID', 'ko_id']].rename_axis('ORF_ID')

# "Unclassified" is added to bacterial taxa that are unclassified at the species level
# Unfortunately this can't be done sooner in the pipeline because I think DRAM distill needs the empty species blank to work properly
bin_annotations['Taxonomy'] = bin_annotations['Taxonomy'].replace(regex={r'([a-z])__;': r'\g<1>__Unclassified;', r's__$': 's__Unclassified'})

pd_args = {'sep': '\t', 'engine': 'c', 'header': None}

# Gene lengths, Bin_IDs, sample-wise coverage, ANI-cluster mapping, and gene info file(s), are aggregated.
gene_length = pd.concat([pd.read_csv(df, **pd_args, names=['ORF_ID', 'Length']) for df in lengths], axis=0, copy=False)\
            .drop_duplicates(keep='first').set_index('ORF_ID')
gene_coverage = pd.concat([pd.read_csv(df, **pd_args, index_col=0, names=['ORF_ID', read_samp]) for df, read_samp in zip(coverages, read_sample)], axis=1, copy=False)


# The output of DRAM distill is parsed, providing pathway info about each gene
dist_dfs = [df.assign(sheet=([sname] * df.shape[0])).drop_duplicates() for sname, df in pd.read_excel(pathways, sheet_name=None).items()]
distill = pd.concat(dist_dfs, axis=0, copy=False).reset_index(drop=True)
gene_info = ['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader', 'specific_reaction', 'Notes', 'EC', 'oxygen', 'alt_ids']
gene_samps = cont_annotations['fasta'].unique().tolist()

# Rows that didn't match any samples and columns that were misnamed are removed
distill = distill.loc[distill[gene_samps].dropna(how='all').index, :].reset_index(drop=True)
distill = distill[distill.columns[distill.columns.isin(gene_info + gene_samps)]]

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
        .unstack(fill_value='').agg(lambda row: np.unique(row[row.astype(bool)]), axis=1)
distill.loc[~distill['EC'].isna(), 'EC'] = distill.loc[~distill['EC'].isna(), 'EC'].str.replace(r',?(?:EC)\W?', ',EC:', regex=True)\
        .apply(lambda row: np.sort(row.split(',')[1:]))

distill.loc[found_ecs.index, 'EC'] = found_ecs.values

# The distill output is de-replicated by grouping to the level of gene_id. Descriptions/Bin_IDs/etc that differ within a gene_id are joined in a list
agg_dict = {col: tuple if col in gene_info else lambda c: np.unique(c)[0] for col in distill.columns[distill.columns != 'gene_id']}
distill.loc[:, gene_samps] = distill.loc[:, gene_samps].apply(lambda col: col.str.split(','), result_type='expand', axis=0)
distill_derep = distill.groupby('gene_id', sort=False, as_index=True).agg(agg_dict, axis=0)

# Individual ORF_IDs that mapped to each gene are compiled in a dict, allowing total copy number of each gene to be counted later
# Even though ORFs were clustered before DRAM annotation, clustered ORFs were given the same annotations as their representative ORF (ORF_ClustID) 
# during the contig annotation script. So DRAM distill counted those to make the metabolism summary used here.
orf2gene_dict = collections.defaultdict(list)
for tup in distill_derep[gene_samps].fillna('NA').itertuples(name=None):
    for val in it.chain.from_iterable(it.filterfalse(lambda item: item == 'NA', tup[1:])):
        orf2gene_dict[val].append(tup[0])

# A few extra genes can be found that have a KO id, but not gene ID from DRAM distill. These are added to the orf2gene dataframe,
# since a KO is technically an ID
orf2gene_df = pd.DataFrame.from_dict(orf2gene_dict, orient='index').rename_axis('ORF_ID').agg(lambda row: list(row.dropna()), axis=1).rename('gene_id')
orf2gene = pd.concat([orf2gene_df, cont_annotations['ko_id'].dropna().str.split(',').str[0].groupby('ORF_ID', sort=False).agg(list)], axis=1, sort=False)
orf2gene.loc[orf2gene['gene_id'].isna(), 'gene_id'] = orf2gene.loc[orf2gene['gene_id'].isna(), 'ko_id']
orf2gene.drop(columns='ko_id', inplace=True)

# The distill pathway information is reformatted after the ORF_IDs were extracted
distill_info = distill_derep[gene_info[1:]].reset_index(names='gene_id')
distill_info = distill_info.explode(gene_info[1:]).explode('EC').explode('alt_ids').groupby('gene_id', sort=False).agg(lambda c: tuple(c.dropna().unique()))

annot2gene = bin_annotations.join(orf2gene).dropna(subset=['gene_id']).explode('gene_id')[['Bin_ID', 'Taxonomy', 'Tax_ID', 'gene_id']]\
        .reset_index()
bin2gene = annot2gene[['Bin_ID', 'gene_id']].value_counts(sort=False).rename('Count').reset_index()\
        .pivot_table(index='Bin_ID', columns='gene_id', values='Count', fill_value=0, aggfunc='sum')
bin2gene_final = annot2gene.set_index('Bin_ID')[['Taxonomy', 'Tax_ID']].drop_duplicates().join(bin2gene)

# Write bin2gene to file here
bin2gene_final.to_csv(gene_bincount_outfile, sep='\t')

## All the files from earlier are finally merged into a single dataframe. Similar to before, the dataframe is aggregated to
## the gene_id level, and since some genes can have multiple annotations/ORFs/Bins (and even Clust_IDs), those are merged into lists
gene_stats = pd.concat([cont_annotations[['ko_id', 'ORF_ClustID']], orf2gene, gene_coverage, gene_length], axis=1, copy=False).reset_index()
gene_agg_dict = {col: tuple if col in ['ORF_ID', 'ko_id', 'gene_id'] else lambda c: list(c)[0] for col in gene_stats.columns[gene_stats.columns != 'ORF_ClustID']}
gene_stats_agg = gene_stats.explode('gene_id').groupby('ORF_ClustID', sort=False, as_index=False).agg(gene_agg_dict)

## Genes that had coverage but weren't assigned an identifier, or vice versa, are removed.
# ORFs aren't exploded, so duplicate gene_ids only occur if 2 independent ORF cluster reps were given the
# same gene ID by DRAM (in which case they would also have separate coverage information).
# Cases where one ORF was given 2 different gene_ids are okay however, even though each gene would have the
# same coverage (since they came from the same ORF)
gene_filt = gene_stats_agg.explode(['gene_id']).loc[~gene_stats_agg['gene_id'].explode().isna()].index.unique()
gene_stats_filt = gene_stats_agg.loc[gene_filt, :].explode(['ko_id', 'gene_id']).dropna(subset=read_sample, how='all').drop_duplicates()

# Gene copy number is calculated by counting the total number of ORFs that map to each gene across samples
cpy_nums = gene_stats_filt.explode(['ORF_ID'])['gene_id'].value_counts().rename('Copy_Num')
sum_cov = gene_stats_filt.groupby('gene_id', sort=False)[gene_coverage.columns.to_list() + gene_length.columns.to_list()].agg('sum').join(cpy_nums)

# If desired, the complete ORF_ID/Clust_ID list of each gene can be written to the output by including their columns here. Those were removed
# from the final output because they make the file 10x larger
final_df_list = [gene_stats_filt.loc[:, ['gene_id', 'ko_id']].groupby('gene_id', sort=False).agg(lambda c: tuple(c.dropna().unique())),
                 distill_info, sum_cov]

# Some KOs that were incomplete are filled in from the alt_ids column, which is then removed
# Only the extra IDs from the CAMPER dataset have this problem, so overwriting valid KOs shouldn't be possible
all_gene_stats = pd.concat(final_df_list, copy=False, axis=1)
missing_ko = all_gene_stats.loc[~all_gene_stats['alt_ids'].isin(all_gene_stats['ko_id'])]
all_gene_stats.loc[missing_ko.index, 'ko_id'] = all_gene_stats.loc[missing_ko.index, 'alt_ids']

gene_stats_final = all_gene_stats.dropna(subset=read_sample, how='all').drop(columns='alt_ids').fillna('()').T
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
