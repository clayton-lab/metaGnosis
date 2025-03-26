import ncbi_genome_download as ngd
import pandas as pd
import pathlib
import subprocess
import sys

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

# This would be snakemake.wildcards("accn", "") and asm_name (or whatever) respectively
accn = snakemake.wildcards.get("accn", "")
file_formats = snakemake.params.get("file_formats", "")
genome_db = snakemake.params.get("genome_db", "")
metafile = snakemake.input[0]
outfile = snakemake.output[0]
genome_meta = pd.read_csv(metafile, sep='\t')
mask = genome_meta['assembly_accession'].values == accn
record = genome_meta[mask].to_dict(orient='records')[0]

conf = ngd.NgdConfig.from_kwargs(groups=record['group'], section='genbank', file_formats=file_formats, dry_run=False, flat_output=True, assembly_accessions=accn,
                                         refseq_categories = 'reference,representative', output=genome_db, use_cache=False)

job = ngd.core.create_downloadjob(record, record['group'], conf)
if job:
    download = ngd.core.worker(job[0])
