import numpy as np
import pandas as pd
import pathlib
import subprocess
import sys

#TODO: Add a thing here that chooses the best mapper (would output a file named best_mapper.txt that could be used by later rules)
if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

summary_outfile = pathlib.Path(snakemake.output.get('summary', ''))
tsv_outfile = pathlib.Path(snakemake.output.get('stats_tsv', ''))
csv_outfile = pathlib.Path(snakemake.output.get('stats_csv', ''))
bin_dfs = pd.concat(pd.read_csv(bin_df, sep='\t') for bin_df in snakemake.input.get('sum_files', '')) 
#bin_stats = bin_dfs.loc[:, ~bin_dfs.columns.isin(['Mapper', 'Contig_Sample'])]
bin_stats = bin_dfs.loc[:, ~bin_dfs.columns.isin(['Mapper'])]
bin_stats.to_csv(tsv_outfile, index=False, sep='\t')


# A .csv version of the bin stats is written for downstream rules to use
bin_stats.loc[:, 'Name'] = bin_dfs['Name'] + '.fa'
bin_stats.rename(columns={'Name': 'genome', 'Completeness': 'completeness', 'Contamination': 'contamination'}).to_csv(csv_outfile, index=False, sep=',')

with open(summary_outfile, 'w') as f:
    f.write('Mapper\tContig_Sample\tTotal_Bins\tTotal_Pass\tPercent_Pass\n')
    for mapper in bin_dfs['Mapper'].unique():
        for contig_sample in bin_dfs['Contig_Sample'].unique():
            mask = ((bin_dfs['Mapper'] == mapper) & (bin_dfs['Contig_Sample'] == contig_sample))
            bin_df = bin_dfs[mask]
            filtered_bins = bin_df['Passed_Filter'].value_counts().to_dict()
            total_bins = filtered_bins[True] + filtered_bins[False]
            total_pass = filtered_bins[True]
            percent_pass = round((total_pass / total_bins), 4) * 100
            f.write(f'{mapper}\t{contig_sample}\t{total_bins}\t{total_pass}\t{percent_pass}\n')

with open(snakemake.output.get('filt_paths', ''), 'w') as outfile:
    for pathfile in snakemake.input.get('path_files'):
        with open(pathfile, 'r') as f:
            for line in f:
                outfile.write(line)
