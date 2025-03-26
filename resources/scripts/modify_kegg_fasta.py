#!/usr/bin/env python

from Bio.SeqIO.FastaIO import SimpleFastaParser
import collections

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

accn = snakemake.wildcards.get("accn", "")
kegg_db = snakemake.params.get('kegg_db', '')
kegg_link = snakemake.params.get('kegg_link', '')
mod_fasta = snakemake.output[0]
gene_to_ko = collections.defaultdict(list)

with open(kegg_link, 'r') as linkfile:
    for line in linkfile:
        pair = line.strip().split('\t')
        gene_to_ko[pair[0]].append(pair[1])

with open(kegg_db, 'r') as infile:
    with open(mod_fasta, 'w') as outfile:
        for header, seq in SimpleFastaParser(infile):
            desc = header.strip().split()
            match = gene_to_ko.get(desc[0])
            if not match:
                outfile.write(f'>{desc[0]} ()  {(" ").join(desc[1:])}\n{seq}\n')
            else:
                for ko in match:
                    outfile.write(f'>{desc[0]} ({ko})  {(" ").join(desc[1:])}\n{seq}\n')
