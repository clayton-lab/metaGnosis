import sys
import pandas as pd
import pathlib
from Bio import SeqIO

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

#derep_genes = 'miniout/dereplicated_genes/dereplicated_genes_rep_seq.fasta'
#contig_samples = ['ABX-CJ-IRIS', 'ABX-CJ-MANG']
#outdir = pathlib.Path('miniout/dereplicated_gene_files')

# Wildcards, params, input and output files are read in
contig_sample = snakemake.params.get('contig_sample', '')
derep_genes = snakemake.input[0]
outdir = pathlib.Path(snakemake.output[0])

#if outdir.exists():
#    for i in outdir.iterdir():
#        i.unlink()
#    outdir.rmdir()
outdir.mkdir(parents=True, exist_ok=True)

outfiles = {contig_samp: outdir.joinpath(f'{contig_samp}.faa') for contig_samp in contig_sample}

uniq_headers = set()
dup_list = []
records = SeqIO.index(derep_genes, 'fasta')
for contig_samp in contig_sample:
    mod_records = []
    for key in records.keys():
        if key.startswith(contig_samp):
            record = records[key]
            record.id = record.id.lstrip(contig_samp + "_")
            record.description = record.description[(len(contig_samp)+1):]
            if record.description in uniq_headers:
                dup_num = dup_list.count(record.description) + 2
                dup_list.append(record.description)
                print(record.description)
                record.description = f'{record.description};{dup_num}'
                print(record.description)
            uniq_headers.add(record.description)
            record.seq = record.seq.translate(table=11)
            mod_records.append(record)

    with open(outfiles[contig_samp], 'w') as outfile:
        SeqIO.write(mod_records, outfile, "fasta")

records.close()

