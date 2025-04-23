import sys
import pandas as pd
import numpy as np
import pathlib
import shutil

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

#read_sample = ['ABX-CJ-IRIS-14', 'ABX-CJ-IRIS-29', 'ABX-CJ-IRIS-43', 'ABX-CJ-IRIS-56', 'ABX-CJ-MANG-14',
#               'ABX-CJ-MANG-29', 'ABX-CJ-MANG-43', 'ABX-CJ-MANG-56']
#contig_sample = ['ABX-CJ-IRIS', 'ABX-CJ-MANG', 'MANG_Assembly_Day_43']
#coverages = ['output/mapping/minimap2/coverage_tables/bins/' + f'{readSamp}_bin_coverage.txt' for readSamp in read_sample]
#countfiles = ['output/quant_bins/sample_read_counts/' + f'{readSamp}_read_counts.csv' for readSamp in read_sample]
#genome_stats="output/annotate_bins/annotate_bin_pathways/genome_stats.tsv"  
#annot =  'output/annotate_bins/annotate_bins/bin_annotations.tsv'
#bin_info_outfile = 'output/quant_bins/quantified_bin_info.tsv'
#bin_quant_outfile = "output/quant_bins/quantified_bins.txt"

# Wildcards, params, input and output files are read in
read_sample = snakemake.params.get('read_sample', '')
coverages = snakemake.input.get('coverages', '')
countfiles = snakemake.input.get('read_counts', '') if snakemake.input.get('read_counts', '') else []
genome_stats=snakemake.input.get('genome_stats', '')
annot = snakemake.input.get('bin_annotation', '')
bin_quant_outfile = snakemake.output.get('quant_bins', '')
bin_info_outfile = snakemake.output.get('bin_info', '')


tax_folder = pathlib.Path(bin_quant_outfile).parent.joinpath('bin_abundance_by_taxa_level')

pd_args = {'sep': '\t', 'engine': 'c', 'header': None}


annotations = pd.read_csv(annot, sep='\t', engine='c', index_col=0, dtype='object').rename_axis('ORF_ID').reset_index()\
        .rename(columns={'scaffold': 'Contig_ID', 'fasta': 'Bin_ID', 'bin_taxonomy': 'Taxonomy', 'taxid': 'Tax_ID', 'bin_completeness': 'Completeness', 'bin_contamination': 'Contamination'})
annotations['Seq_ID'] = annotations['Bin_ID'] + '_' + annotations['Contig_ID']
annot_subset = annotations[['Seq_ID', 'Taxonomy', 'Bin_ID', 'Tax_ID', 'Completeness', 'Contamination', 'Bin_ClustID']].drop_duplicates().set_index('Seq_ID')

## Sample-wise contig coverage, contig2bin file(s), and nonhost/host read lengths are aggregated.
coverage = pd.concat([pd.read_csv(df, **pd_args, index_col=0, names=['Seq_ID', read_samp]) for df, read_samp in zip(coverages, read_sample)], axis=1, copy=False)
bin_coverage = annot_subset.join(coverage)
genomes = pd.read_csv(genome_stats, sep='\t', engine='c', header=0).rename(columns={'genome': 'Bin_ID'}).set_index('Bin_ID')

# Bin dataframes are merged, coverage is summed across contigs into bins, and blank spots in the taxonomy are replaced with "Unclassified"
bin_stats = bin_coverage.reset_index().set_index('Bin_ID').join(genomes[['5S rRNA', '16S rRNA', '23S rRNA', 'tRNA count']]).reset_index()
cont_agg_dict = {col: (lambda c: tuple(c.unique())) if col not in read_sample else lambda c: c.dropna().sum() for col in bin_stats.columns[bin_stats.columns != 'Seq_ID']}

# After summing coverage across contigs, an aggregated dataframe at the bin level is created. Since bins were de-replicated before mapping,
# coverage is only calculated for cluster representative bins (i.e., cluster representatives). This ensures that species with multiple high quality bins
# aren't overcounted during relative abundance calculation.
bin_stats_agg = bin_stats.groupby(['Bin_ClustID', 'Bin_ID'], as_index=False, sort=False).agg(cont_agg_dict).explode('Taxonomy')
bin_stats_agg['Taxonomy'] = bin_stats_agg['Taxonomy'].replace(regex={r'([a-z])__;': r'\g<1>__Unclassified;', r's__$': 's__Unclassified'})

# Because genomes were clustered at 95% ANI similarity before taxonomy prediction, it is likely that each unclassified species is unique.
# So in cases where multiple unclassified species from the same genus were identified, they are numbered Unclassified_1, Unclassified_2, etc.
rep_taxa = bin_stats_agg[['Bin_ClustID', 'Taxonomy']].drop_duplicates().explode('Bin_ClustID')
dup_taxa = rep_taxa['Taxonomy'].duplicated(keep=False)
rep_taxa.loc[dup_taxa, 'Taxonomy'] += '_' + rep_taxa.groupby('Taxonomy').cumcount().add(1).astype(str)
rep_taxa.set_index('Bin_ClustID', inplace=True)
bin_stats_agg = bin_stats_agg.explode(bin_stats_agg.columns[~bin_stats_agg.columns.isin(['Taxonomy'] + read_sample)].to_list())
bin_stats_agg.loc[:, 'Taxonomy'] = rep_taxa.loc[bin_stats_agg['Bin_ClustID'], 'Taxonomy'].array

# Bins are again aggregated at the level of their cluster representatives
bin_agg_dict = {col: (lambda c: tuple(c.dropna())) if col not in read_sample else lambda c: c.sum() for col in bin_stats_agg.columns[~bin_stats_agg.columns.isin(['Bin_ClustID', 'Taxonomy'])]}
bin_agg_dict.update({'Taxonomy': 'unique'})
derep_bins = bin_stats_agg.groupby('Bin_ClustID', sort=False).agg(bin_agg_dict).explode('Taxonomy').set_index('Taxonomy')
derep_bins = derep_bins[derep_bins.columns[~derep_bins.columns.isin(read_sample)].to_list() + read_sample]

#print(bin_stats_agg[['Bin_ClustID', 'Taxonomy']].drop_duplicates().set_index('Bin_ClustID').value_counts())
#test_stats = bin_stats_agg.sort_values(['Taxonomy', 'Bin_ClustID'])[['Bin_ClustID', 'Bin_ID', 'Taxonomy', 'ABX-CJ-IRIS-14']]
#test_stats = bin_stats_agg[['Bin_ClustID', 'Bin_ID', 'Taxonomy', 'Tax_ID']].groupby('Bin_ClustID', sort=False).agg(lambda x: tuple(x.unique()))
#print(test_stats.shape, test_stats.explode('Taxonomy')[['Bin_ID', 'Taxonomy']])
#[print(x) for x in test_stats.explode('Taxonomy').sort_values('Taxonomy')['Taxonomy']]
#print(bin_stats_agg.explode(bin_stats_agg.columns[~bin_stats_agg.columns.isin(['5S rRNA', '16S rRNA', '23S rRNA', 'tRNA count'])].to_list()))
#print(bin_stats_agg.explode(['5S rRNA', '16S rRNA', '23S rRNA', 'tRNA count']))
#print(bin_stats_agg.groupby('Taxonomy', sort=False).agg(tuple)['Bin_ClustID'])
#print(bin_stats_agg, bin_stats_agg[['Bin_ID', 'Tax_ID', 'Completeness', 'Contamination', '5S rRNA', '16S rRNA', '23S rRNA', 'tRNA count']], bin_stats_agg.sort_values('Taxonomy')[['Taxonomy'] + read_sample], bin_stats_agg.columns)
#bin_stats = pd.concat([coverage, annot_subset, genomes[['5S rRNA', '16S rRNA', '23S rRNA', 'tRNA count']]], axis=1, copy=False)
#print(bin_coverage.groupby('Bin_ClustID').agg(tuple))

#print(bin_coverage.groupby('Bin_ClustID').agg(lambda col: tuple(col.dropna()))['ABX-CJ-MANG-56'].shape)



#bin_stats = bin_coverage.groupby('Bin_ID', sort=False).sum(numeric_only=True).join(genomes)[genomes.columns.to_list() + coverage.columns.to_list()]
#bin_stats.reset_index(inplace=True)
#taxids['Taxonomy'] = taxids['Taxonomy'].replace(regex={r'([a-z])__;': r'\g<1>__Unclassified;', r's__$': 's__Unclassified'})
#taxids = taxids.drop_duplicates().set_index('Taxonomy')

# If host filtering was applied, the microbial load is approximated as the ratio of nonhost to host-mapped reads and added to the bin_stats dataframe
if countfiles:
    read_counts = pd.concat([pd.read_csv(countdf, sep=',', index_col=0, header=0) for countdf in countfiles], axis=1)
    read_counts.loc['Microbial load'] = read_counts.loc['Nonhost readcount', :] / read_counts.loc['Host readcount', :]
    bin_stats_final = pd.concat([derep_bins, read_counts], axis=0, copy=False).T
    bin_samps = bin_stats_final.loc[read_sample]
    bin_tax = bin_stats_final.loc[bin_samps.index, ~bin_stats_final.columns.isin(['Nonhost readcount', 'Host readcount', 'Microbial load'])].astype('float')

else:
    bin_stats_final = derep_bins.T
    bin_samps = bin_stats_final.loc[read_sample]
    bin_tax = bin_stats_final.loc[bin_samps.index, :].astype('float')


# Extra rows with the filtered bin statistics are kept to allow later manual inspection.
bin_stats_final.to_csv(bin_info_outfile, index=True, sep='\t')

# Relative abundance is calculated for the bins at every taxonomic level (species, genus, etc.)
tax_split_dict = {0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6:'species'}
tax_split = bin_tax.columns.to_series().str.split(';', expand=True).rename(tax_split_dict).join(bin_tax.T)
for tax_level in range(0, 7):
    basepath = tax_folder.joinpath(f'{tax_split_dict[tax_level]}_level_bin_abundance')
    basepath.mkdir(parents=True, exist_ok=True) 

    if tax_level < 6:
        tax_agg = tax_split.groupby(tax_split.columns[tax_level], sort=False, as_index=False).agg({col: 'sum' if col in bin_samps.index else lambda c: c.iloc[0] for col in tax_split.columns})
        tax_agg = tax_agg.set_index(tax_agg[tax_split.columns[:tax_level+1]].apply(lambda c: list(c), axis=1).str[:tax_level+1].str.join(';')).rename_axis('Taxonomy')[bin_tax.index].T
    else:
        tax_agg = tax_split[bin_tax.index].T

    bin_rel_abund = tax_agg.div(tax_agg.sum(axis=1), axis=0)
    bin_rel_abund.to_csv(basepath.joinpath(f'{tax_split_dict[tax_level]}_relative_abundance.tsv'), index=True, sep='\t')

    # Microbial load-normalized abundance is also calculated if host filtering was applied 
    if countfiles:
        bin_rel_norm = bin_rel_abund.mul(bin_stats_final.loc[bin_tax.index, 'Microbial load'], axis=0)
        bin_rel_norm.to_csv(basepath.joinpath(f'{tax_split_dict[tax_level]}_relative_abundance_microbial_loadnormed.tsv'), index=True, sep='\t')

shutil.copy(tax_folder.joinpath('species_level_bin_abundance/species_relative_abundance.tsv'), bin_quant_outfile)

# This section is a working version of the code that can calculate several abundance measures including TPM, metaWRAP abundance, relative and 
# and microbial load-normalized abundance. It requires the length of each contig, similar to the gene_lengths output. It doesn't
# integrate taxonomy information, so tweaking is required to get it to work the way the above code does. It also requires the length of
# each contig, which can be generated similar to how it is for the gene quantification
#######################################################################_#############################################################
#import sys
#import pathlib
#import pandas as pd
#import numpy as np
#import itertools as it
#import functools

# Wildcards, params, input and output files are read in
#read_sample = snakemake.params.get('read_sample')
#contig_sample = snakemake.wildcards.get('contig_sample')
#coverages = snakemake.input.get('coverages')
#lengths = snakemake.input.get('lengths')
#contig2bins = snakemake.input.get('contigs2bins', '')[0]
#countfiles = snakemake.input.get('read_counts')
#outfile = snakemake.output[0]
#pd_args = {'sep': '\t', 'engine': 'c', 'header': None}

## The microbial load is approximated as the ratio of nonhost to host-mapped reads.
#read_counts = pd.concat([pd.read_csv(countdf, sep=',', index_col=0, header=0) for countdf in countfiles], axis=1)
#read_counts.loc['Microbial load'] = read_counts.loc['Nonhost readcount', :] / read_counts.loc['Host readcount', :]
#read_counts_dict = read_counts.to_dict()
#
## Contig lengths, sample-wise coverage, contig2bin file(s), and nonhost/host read lengths are aggregated.
#length = pd.concat([pd.read_csv(df, **pd_args, names=['Contig_ID', 'Length']) for df in lengths], axis=0, copy=False)\
#            .drop_duplicates(keep='first').set_index('Contig_ID')
#
#cov_dfs = [pd.read_csv(df, **pd_args, index_col=0, names=['Contig_ID', read_samp]) for df, read_samp in zip(coverages, read_sample)]
#coverage = pd.concat(cov_dfs, axis=1, copy=False)
#contig2bin = pd.read_csv(contig2bins, **pd_args, names=['Contig_ID', 'Bin_ID'], index_col=0)
#samples = {}
#
## MetaWRAP abundance is calculated. Sample-wise coverage for each contig in an assembly file is converted to kilobases (called Weight).
#bin_stats = functools.reduce(lambda left, right: pd.merge(left, right, on='Contig_ID'), [contig2bin, length, coverage])
#bin_stats.insert(2, 'Weight', bin_stats['Length'] // 1000)
#
#bin_dict = bin_stats.to_dict()
#
#bin_stats_agg = bin_stats.reset_index().groupby('Bin_ID').sum(numeric_only=True).drop(columns='Weight')
#bin_stats_agg.insert(1, 'Total_Coverage', bin_stats_agg.iloc[:, 1:].sum(axis=1))
#rpk = bin_stats_agg[read_sample].div((bin_stats_agg['Length'] / 1000), axis=0)
#scale_factor = rpk.sum(axis=0)/1000000
#tpm = (rpk / scale_factor)
#rel_abund = (bin_stats_agg[read_sample] / bin_stats_agg[read_sample].sum(axis=0))
#bin_stats_final = bin_stats_agg[['Length', 'Total_Coverage']].rename(columns={'Length': 'Total_Length'})
#
## For each contig, the coverage is duplicated n times (n=Weight), and the median of that is taken to approximate coverage at the
## bin level. This calculation is performed for each sample.
#for read_samp in read_sample:
#    samples.update({read_samp: {bin_id: [] for bin_id in contig2bin['Bin_ID'].unique()}})
#    for contig in bin_stats.index.to_list():
#        samples[read_samp][bin_dict['Bin_ID'][contig]].extend(it.repeat(bin_dict[read_samp][contig], bin_dict['Weight'][contig]))
#
#    for bin_id in contig2bin['Bin_ID'].unique():
#        samples[read_samp].update({bin_id: np.median(samples[read_samp][bin_id])})
#
#    # Relative, TPM-normalized, and metaWRAP (i.e., average) abundance are included in the final output for each sample
#    bin_stats_final[read_samp + '_Rel_Abund'] = rel_abund[read_samp]
#    bin_stats_final[read_samp + '_TPM_Abund'] = tpm[read_samp]
#    bin_stats_final[read_samp + '_Avg_Abund'] = samples[read_samp]
#
#    # Finally, the sample-wise abundance calculated from earlier steps is normalized by the microbial load for each sample,
#    # minimizing the effect of host contamination and (hopefully) bypassing the compositional nature of the sequence data
#    bin_stats_final[read_samp + '_Rel_Abund_MicrobLoadNormed'] = rel_abund[read_samp] / read_counts_dict[read_samp]['Microbial load']
#    bin_stats_final[read_samp + '_TPM_Abund_MicrobLoadNormed'] = tpm[read_samp] / read_counts_dict[read_samp]['Microbial load']
#    bin_stats[read_samp + '_Avg_Abund_MicrobLoadNormed'] = pd.Series(samples[read_samp]) / read_counts_dict[read_samp]['Microbial load']
#    
#
##bin_stats_final.columns = ['Bin_ID', 'Total_Length', 'Total_Coverage'] + [f'{sample}_Abundance' for sample in read_sample]
##bin_stats_final.columns = ['Bin_ID', 'Total_Length', 'Total_Coverage'] + [f'{sample}_Average_Abundance' for sample in read_sample]
#
##bin_stats_final.iloc[:, 3:] = bin_stats_final.iloc[:, 3:] / read_counts.loc['Microbial load'].transpose().to_frame().T.values
#bin_stats_final.to_csv(outfile, index=True, sep='\t')
