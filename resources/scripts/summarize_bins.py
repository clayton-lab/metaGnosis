import numpy as np
import pandas as pd
import pathlib
import subprocess
import sys

#TODO: Add a thing here that chooses the best mapper (would output a file named best_mapper.txt that could be used by later rules)
if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

outfile = pathlib.Path(snakemake.output[0])
bin_dfs = pd.concat(pd.read_csv(bin_df, sep='\t') for bin_df in snakemake.input) 
bin_dfs.loc[:, ~bin_dfs.columns.isin(['Mapper', 'Contig_Sample'])]\
    .to_csv(outfile.parent.joinpath('compiled_bin_statistics.tsv'), index=False, sep='\t')

with open(outfile, 'w') as f:
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
