import sys
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

#read_sample = ['ABX-CJ-IRIS-14', 'ABX-CJ-IRIS-29', 'ABX-CJ-IRIS-43', 'ABX-CJ-IRIS-56']
#assembly = '/mnt/nrdstor/claytonlab/jhernandez/projects/metaGnosis/output/assemble/megahit/ABX-CJ-IRIS.contigs.fasta'
#coverages = ['/mnt/nrdstor/claytonlab/jhernandez/projects/metaGnosis/output/mapping/minimap2/coverage_tables/contigs/' + f'{readSamp}_Mapped_To_ABX-CJ-IRIS_coverage.txt' for readSamp in read_sample]
#outfile = 'ABX-CJ-IRIS_Coverage_Table.tsv'

# Wildcards, params, input and output files are read in
read_sample = snakemake.params.get('read_sample', '')
assembly = snakemake.input.get('contigs', '')[0]
coverages = snakemake.input.get('coverages', '')
outfile = snakemake.output[0]

pd_args = {'sep': '\t', 'engine': 'c', 'header': None}

## Sample-wise contig coverage is aggregated.
cov_dfs = [pd.read_csv(df, **pd_args, index_col=0, names=['contigname', read_samp]) for df, read_samp in zip(coverages, read_sample)]
coverage = pd.concat(cov_dfs, axis=1, copy=False)

contig_headers = []

# The coverage file is reorganized to be in the same order as the assembly headers
with open(assembly, 'r') as infile:
    for header, seq in SimpleFastaParser(infile):
        contig_headers.append(header.split(' ')[0])

coverage.loc[contig_headers].to_csv(outfile, index=True, sep='\t')
