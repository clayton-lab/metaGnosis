#!/usr/bin/env python

import argparse
import pathlib
import re
import pandas as pd
import gzip
from Bio import SeqIO

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

accn = snakemake.wildcards.get("accn", "")
# The reason this needs a [0] but the other input file doesn't is probably because the syntax
# is different for the 2 input files (i.e., with snakemake.rule.output vs an explicit filestring)
input_seqfile = snakemake.input.get("seqfile", "")[0]
metafile = snakemake.input.get("genome_meta", "")
reformatted_seqfile = snakemake.output[0]
pat = re.compile(r'^(\S+)(\s.+)$')
genome_meta = pd.read_csv(metafile, sep='\t')
mask = genome_meta['assembly_accession'].values == accn
taxid = genome_meta[mask]['taxid'].values[0]

with gzip.open(input_seqfile, 'rt') as infile:
    with open(reformatted_seqfile, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            record.id = pat.sub(r'\1|kraken:taxid|{}'.format(taxid), record.description)
            record.description = pat.sub(r'\1|kraken:taxid|{}\2'.format(taxid), record.description)
            SeqIO.write(record, outfile, 'fasta')
