import sys
import pandas as pd
import pathlib
import shutil
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser
from mag_annotator import annotate_bins  # DRAM's underlying python tools are imported
from mag_annotator.utils import setup_logger

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

# Wildcards, params, input and output files are read in
#cluster_reps = 'output/refine_bins/dereplicated_bins/data_tables/Wdb.csv'
#clusters = 'output/refine_bins/dereplicated_bins/data_tables/Cdb.csv'
#annot = pathlib.Path('annotate_test/annotations.tsv')
#bin_statistics = 'output/refine_bins/summarize_bins/compiled_bin_statistics.tsv'
#tax = 'output/annotate_bins/annotate_bin_taxonomy/gtdbtk.bac120.summary.tsv'
#outdir = pathlib.Path('annot_bins_test')
#tempdir = pathlib.Path('annot_bins_temp')
#threads = 8
#contigfile = 'annotate_test/scaffolds.fna'
#contig_faa = 'annotate_test/genes.faa'
#contig_fna = 'annotate_test/genes.fna'
#contig_gff = 'annotate_test/genes.gff'
#cont2bin = ['output/refine_bins/DAS_Tool/minimap2/run_DAS_Tool/ABX-CJ-IRIS/ABX-CJ-IRIS_DASTool_scaffolds2bin.txt',
#            'output/refine_bins/DAS_Tool/minimap2/run_DAS_Tool/ABX-CJ-MANG/ABX-CJ-MANG_DASTool_scaffolds2bin.txt',
#            'output/refine_bins/DAS_Tool/minimap2/run_DAS_Tool/MANG_Assembly_Day_43/MANG_Assembly_Day_43_DASTool_scaffolds2bin.txt']

# Shouldn't need contig_samples in actual snakemake script
cont2bin = snakemake.input.get('contigs2bins', '')
cluster_reps = snakemake.input.get('clust_reps', '')
clusters = snakemake.input.get('clusters', '')
annot = pathlib.Path(snakemake.input.get('annot', ''))
bin_statistics = snakemake.input.get('bin_stats', '')
tax = snakemake.input.get('taxonomy', '')[0]
contigfile = snakemake.input.get('annot_contigs', '')
contig_faa = snakemake.input.get('contig_faa', '')
contig_fna = snakemake.input.get('contig_fna', '')
contig_gff = snakemake.input.get('contig_gff', '')
binned_contigs = snakemake.output.get('binned_contigs', '')
bin_faa = snakemake.output.get('bin_faa', '')
bin_fna = snakemake.output.get('bin_fna', '')
bin_gff = snakemake.output.get('bin_gff', '')
threads = snakemake.threads
outdir = pathlib.Path(snakemake.params.get('outdir'))

outdir.mkdir(parents=True, exist_ok=True)
log_file_path = outdir.joinpath("annotate.log")
logger = annotate_bins.logging.getLogger("annotation_log")
annotate_bins.setup_logger(logger, log_file_path)

trna_loc=snakemake.output.get('trnas', '')
rrna_loc=snakemake.output.get('rrnas', '')

# Taxonomy, bin cluster, and contig2bin files are read and concatenated together
taxonomy = pd.read_csv(tax, sep='\t', usecols=['user_genome', 'classification', 'closest_genome_reference'])\
        .rename(columns={'user_genome': 'Bin_ClustID', 'classification': 'bin_taxonomy', 'closest_genome_reference': 'taxid'}).set_index('Bin_ClustID')
contigs2bins = pd.concat([pd.read_csv(contig2bin, engine='c', header=None, sep='\t', names=['Contig_ID', 'Bin_ID']) for contig2bin in cont2bin], axis=0, copy=False)
clust_reps = pd.read_csv(cluster_reps, sep=',', usecols=['genome', 'cluster']).rename(columns={'genome': 'Bin_ClustID'}).set_index('cluster')
clusts = pd.read_csv(clusters, sep=',', usecols=['genome', 'secondary_cluster']).rename(columns={'genome': 'Bin_ID', 'secondary_cluster': 'cluster'}).groupby('cluster', sort=False).agg(list)
clust_agg = pd.concat([clust_reps, clusts], axis=1, copy=False).explode('Bin_ID').reset_index().drop(columns=['cluster'])
clust_agg.loc[:, ['Bin_ClustID', 'Bin_ID']] = clust_agg[['Bin_ClustID', 'Bin_ID']].apply(lambda col: col.str.rstrip('.fa'))
tax_agg = clust_agg.set_index('Bin_ClustID').join(taxonomy).reset_index()

# Bin completion and contamination (bin_stats) is read and only filtered bins are kept
bin_stats = pd.read_csv(bin_statistics, sep='\t', usecols=['Name', 'Completeness', 'Contamination', 'Contig_Sample', 'Passed_Filter'])\
        .rename(columns={'Name': 'Bin_ID', 'Completeness': 'bin_completeness', 'Contamination': 'bin_contamination'})
filt_bins = bin_stats[bin_stats['Passed_Filter'] == True].drop(columns=['Passed_Filter'])

# Bin_stats is merged together with the previous files and with the contig-level DRAM annotation
filt_bin_stats = pd.concat([tax_agg.set_index('Bin_ID'), filt_bins.set_index('Bin_ID')], axis=1, copy=False)
bin_agg = filt_bin_stats.join(contigs2bins.groupby('Bin_ID', sort=False).agg(tuple)).reset_index().explode('Contig_ID')
annotations = pd.read_csv(annot, engine='c', sep='\t', header=0, dtype='object', index_col = 0).reset_index(names='ORF_ID')
annotations['Assembled_Contig'] = annotations['fasta'] + '_' + annotations['scaffold']
bin_agg['Assembled_Contig'] = bin_agg['Contig_Sample'] + '_' + bin_agg['Contig_ID']
annot_contigs = annotations.groupby('Assembled_Contig', sort=False).agg(tuple)

# Contigs that weren't part of the filtered bins are removed from the annotations
annot_bins = bin_agg.set_index('Assembled_Contig').join(annot_contigs)
annot_bins = annot_bins.explode(annot_contigs.columns[annot_contigs.columns.isin(annot_bins.columns)].to_list())

# Because some cluster representative ORFs weren't binned, the assembly-level ORF_ID and ORF_ClustID are kept.
# The Assembly_ORFID refers to the ORF_ID used in the assembly. A new bin-level ORF_ID is created as well,
# which has the Bin_ID (e.g., samp1.concoct.bin.1) instead of the contig sample (e.g., samp1) for a prefix.
annot_bins.rename(columns={'ORF_ID': 'Assembly_ORFID'}, inplace=True)
annot_bins.loc[:, 'fasta'] = annot_bins['Bin_ID']
annot_bins['Binned_Contig'] = annot_bins['Bin_ID'] + '_' + annot_bins['scaffold']
annot_bins['Bin_ORFID'] = annot_bins['Binned_Contig'] + '_' + annot_bins['gene_position'].astype('object')

# Annotations are sorted
sort_cols = ['fasta', 'scaffold', 'gene_position']
annot_bins = annot_bins.set_index('Bin_ORFID').sort_values(sort_cols)

# Assembly-level annotated sequences are extracted if they were binned, and headers are renamed to match which bin they were in
orf2bin = annot_bins.reset_index(names='ORF_ID').set_index('Assembly_ORFID')['ORF_ID'].dropna().to_dict()
assem2bin = annot_bins[['fasta', 'scaffold', 'Contig_Sample']].reset_index().dropna()
assem2bin['Key'] = assem2bin['Contig_Sample'] + '_' + assem2bin['scaffold']
assem2bin['Value'] = assem2bin['fasta'] + '_' + assem2bin['scaffold']
asm2bin = assem2bin[['Key', 'Value']].drop_duplicates().set_index('Key').to_dict()['Value']

# .faa and .fna file sequences are extracted
for inputfile, outputfile in zip((contig_fna, contig_faa), (bin_fna, bin_faa)):
    with open(inputfile, 'r') as infile:
        with open(outputfile, 'w') as outfile:
            for header, seq in SimpleFastaParser(infile):
                line = header.split()
                binned = orf2bin.get(line[0])
                if binned:
                    outfile.write(f'>{binned}{header[len(line[0]):]}\n{seq}\n')
# Binned contigs are extracted
with open(contigfile, 'r') as infile:
    with open(binned_contigs, 'w') as outfile:
        for header, seq in SimpleFastaParser(infile):
            line = header.split()
            binned = asm2bin.get(line[0])
            if binned:
                outfile.write(f'>{binned}{header[len(line[0]):]}\n{seq}\n')

# GFF headers are extracted
with open(contig_gff, 'r') as infile:
    with open(bin_gff, 'w') as outfile:
        headerlist = []
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)
            elif line.startswith('#'):
                headerlist.append(line)
            else:
                old_name = line.split()[0]
                new_name = asm2bin.get(old_name)
                if new_name:
                    if headerlist:
                        for header in headerlist:
                            outfile.write(header)
                        headerlist.clear()
                    outfile.write(line.replace(old_name, new_name))
                else:
                    headerlist.clear()

# tRNAs and rRNAs are predicted for the bin-level annotations
len_dict = {
    i.metadata["id"]: len(i)
    for i in annotate_bins.read_sequence(binned_contigs, format="fasta")
}
trna_table = annotate_bins.run_trna_scan(
        binned_contigs,
        outdir,
        None,
        logger,
        threads=threads,
        verbose=True,
)

rrna_table = annotate_bins.run_barrnap(
    binned_contigs, None, logger, threads=threads, verbose=True
)

# The correct bin information (i.e. which RNA came from which bin) is imputted for the called tRNAs and rRNAs
rna_info = annot_bins.reset_index()[['Binned_Contig', 'fasta']].dropna().drop_duplicates().set_index('Binned_Contig')
trna_table.loc[:, ['Name', 'fasta']] = rna_info.loc[trna_table['Name'].str.rstrip()].rename_axis('Name').reset_index().values
trna_table.to_csv(trna_loc, sep='\t', index=False)
rrna_table.loc[:, ['scaffold', 'fasta']] = rna_info.loc[rrna_table['scaffold'].str.rstrip()].rename_axis('scaffold').reset_index().values
rrna_table.to_csv(rrna_loc, sep='\t', index=False)

# The RNA locations are also added to the bin GFF files
annotate_bins.add_intervals_to_gff(
    trna_loc, bin_gff, len_dict, annotate_bins.make_trnas_interval, "Name", logger
)

annotate_bins.add_intervals_to_gff(
    rrna_loc, bin_gff, len_dict, annotate_bins.make_rrnas_interval, "scaffold", logger
)

# The final bin annotations are written and tempfiles are removed
annot_bins = annot_bins[sort_cols + annot_bins.columns[~annot_bins.columns.isin(sort_cols)].to_list()].rename_axis(None).drop(columns=['Bin_ID', 'Binned_Contig', 'Contig_ID'])
annot_bins.to_csv(outdir.joinpath('bin_annotations.tsv'), sep='\t', index=True)

outdir.joinpath('raw_trnas.txt').unlink()
outdir.joinpath('annotate.log').unlink()
