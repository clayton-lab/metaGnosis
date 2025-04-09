import pandas as pd
import pathlib
import subprocess
import sys
import shutil
import concurrent

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

bin_df = pd.read_csv(snakemake.input[0], sep='\t')
percent_compl = snakemake.params.get('completion_cutoff')
percent_contam = snakemake.params.get('contam_cutoff')
mapper = snakemake.wildcards.get('mapper')
contig_sample = snakemake.wildcards.get('contig_sample')
bin_dir = pathlib.Path(snakemake.params.get('bin_dir'))
out_binpath = pathlib.Path(snakemake.output.get('filtered_dir', ''))
out_summary = pathlib.Path(snakemake.output.get('temp_summary', ''))
out_manifest = snakemake.output.get('bin_paths', '')
out_summary.parent.mkdir(parents=True, exist_ok=True)
bin_df['Mapper'] = mapper
bin_df['Contig_Sample'] = contig_sample
bin_df['Passed_Filter'] = ((bin_df['Completeness'] >= percent_compl) & (bin_df['Contamination'] <= percent_contam))
filtered_bins = bin_df.loc[bin_df['Passed_Filter'] == True]

bin_df.to_csv(out_summary, index=False, sep='\t', mode='a')
out_binpath.mkdir(parents=True, exist_ok=True)

with open(out_manifest, 'w') as manifest:
    with concurrent.futures.ThreadPoolExecutor(100) as exe:
        for bin_id in filtered_bins['Name'].to_list():
            from_file = bin_dir.joinpath(f'{bin_id}.fa')
            to_file = out_binpath.joinpath(f'{bin_id}.fa')
            _ = exe.submit(shutil.copy(from_file.as_posix(), to_file.as_posix()))
            print(f'Copied "{from_file}" to "{to_file}"', file=sys.stderr)
            manifest.write(f'{to_file}\n')
