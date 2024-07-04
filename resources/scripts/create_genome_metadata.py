import ncbi_genome_download as ngd
import pandas as pd
import pathlib
import subprocess
import sys

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

taxids = snakemake.params.get("taxids", "")
accns = snakemake.params.get("accns", "")
file_formats = snakemake.params.get("file_formats", "")
outfile = snakemake.output

cwd = pathlib.Path.cwd()
setattr(ngd.core, 'CACHE_DIR', cwd.joinpath('ngd_cachefiles'))
if taxids:
    taxfile = 'taxid_list.txt'
    subprocess.run(f'gimme_taxa.py -d temp_ete_db -j -o {taxfile} {taxids}', shell=True)
    tax_conf = ngd.NgdConfig.from_kwargs(groups='all', section='genbank', file_formats=file_formats, dry_run=True, flat_output=True, species_taxids=taxfile,
                                         refseq_categories = 'reference,representative', output=cwd.as_posix(), use_cache=True)
    tax_entries = ngd.core.select_candidates(tax_conf)
else:
    tax_entries = []

if accns:
    accn_conf = ngd.NgdConfig.from_kwargs(groups='all', section='genbank', file_formats=file_formats, dry_run=True, flat_output=True, assembly_accessions=accns,
                                         refseq_categories = 'reference,representative', output=cwd.as_posix(), use_cache=True)
    accn_entries = ngd.core.select_candidates(accn_conf)
else:
    accn_entries = []

entries = tax_entries + accn_entries
pd.DataFrame([entry[0] for entry in entries]).replace({'asm_name': {'\s+': '_'}}, regex=True).to_csv(outfile[0], sep='\t', index=False)
