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
#contig2bins = ['output/refine_bins/DAS_Tool/minimap2/run_DAS_Tool/' + f'{contigSamp}/{contigSamp}_DASTool_scaffolds2bin.txt' for contigSamp in contig_sample]
#countfiles = ['output/quant_bins/sample_read_counts/' + f'{readSamp}_read_counts.csv' for readSamp in read_sample]
#genome_stats="output/annotate_bins/annotate_bin_pathways/annotated_pathways/genome_stats.tsv"  
#bin_info_outfile = 'output/quant_bins/quantified_bin_info.tsv'
#bin_quant_outfile = "output/quant_bins/quantified_bins.txt"

# Wildcards, params, input and output files are read in
read_sample = snakemake.params.get('read_sample', '')
coverages = snakemake.input.get('coverages', '')
contig2bins = snakemake.input.get('contigs2bins', '')
countfiles = snakemake.input.get('read_counts', '')
genome_stats=snakemake.input.get('genome_stats', '')
bin_quant_outfile = snakemake.output.get('quant_bins', '')
bin_info_outfile = snakemake.output.get('bin_info', '')

tax_folder = pathlib.Path(bin_quant_outfile).parent.joinpath('bin_abundance_by_taxa_level')

pd_args = {'sep': '\t', 'engine': 'c', 'header': None}

## The microbial load is approximated as the ratio of nonhost to host-mapped reads.
read_counts = pd.concat([pd.read_csv(countdf, sep=',', index_col=0, header=0) for countdf in countfiles], axis=1)
read_counts.loc['Microbial load'] = read_counts.loc['Nonhost readcount', :] / read_counts.loc['Host readcount', :]
read_counts_dict = read_counts.to_dict()

## Sample-wise contig coverage, contig2bin file(s), and nonhost/host read lengths are aggregated.
cov_dfs = [pd.read_csv(df, **pd_args, index_col=0, names=['Seq_ID', read_samp]) for df, read_samp in zip(coverages, read_sample)]
coverage = pd.concat(cov_dfs, axis=1, copy=False)
contig2bin = pd.concat([pd.read_csv(df, **pd_args, names=['Contig_ID', 'Bin_ID']) for df in contig2bins], axis=0, copy=False)
contig2bin['Seq_ID'] = contig2bin['Bin_ID'] + '_' + contig2bin['Contig_ID']
contig2bin.set_index('Seq_ID', inplace=True)
bin_coverage = coverage.join(contig2bin)
genomes = pd.read_csv(genome_stats, sep='\t', engine='c', header=0).rename(columns={'genome': 'Bin_ID'}).set_index('Bin_ID')

# Bin dataframes are merged, coverage is summed across contigs into bins, and blank spots in the taxonomy are replaced with "Unclassified"
bin_stats = bin_coverage.groupby('Bin_ID', sort=False).sum(numeric_only=True).join(genomes)[genomes.columns.to_list() + coverage.columns.to_list()]
bin_stats.reset_index(inplace=True)
bin_stats['taxonomy'] = bin_stats['taxonomy'].replace(regex={r'([a-z])__;': r'\g<1>__Unclassified;', r's__$': 's__Unclassified'})

# Bins with the same taxonomic classification are summed, and microbial load is included in the final output. 
# Extra rows with the filtered bin statistics are kept to allow later manual inspection.
agg_dict = {col: 'sum' if col in coverage.columns else lambda row: tuple(row.dropna()) for col in bin_stats.columns}
bin_stats_agg = bin_stats.reset_index().groupby('taxonomy', sort=False).agg(agg_dict).drop(columns='taxonomy')

bin_stats_final = pd.concat([bin_stats_agg, read_counts], axis=0, copy=False).T
bin_stats_final.to_csv(bin_info_outfile, index=True, sep='\t')

bin_samps = bin_stats_final.loc[read_sample]
bin_tax = bin_stats_final.loc[bin_samps.index, ~bin_stats_final.columns.isin(['Nonhost readcount', 'Host readcount', 'Microbial load'])].astype('float')

# Relative and microbial load normalized relative abundance are calculated for the bins at every taxonomic level (species, genus, etc.)
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
    bin_rel_norm = bin_rel_abund.mul(bin_stats_final.loc[bin_tax.index, 'Microbial load'], axis=0)
    bin_rel_abund.to_csv(basepath.joinpath(f'{tax_split_dict[tax_level]}_relative_abundance.tsv'), index=True, sep='\t')
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
