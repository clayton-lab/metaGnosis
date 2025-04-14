import sys
import pandas as pd
import pathlib
from Bio import SeqIO
from mag_annotator import annotate_bins  # DRAM's underlying python tools are imported
from mag_annotator.utils import setup_logger
import shutil

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")


# Wildcards, params, input and output files are read in
#contig_sample = snakemake.params.get('contig_sample', '')
#derep_genes = snakemake.input[0]
#outdir = pathlib.Path(snakemake.output[0])

#annot_file = "output/annotate_bins/annotate_genes/annotations.tsv"
#concat_fna = "output/refine_bins/dereplicated_genes/concatenated_genes.fna"
#contig_samples = ['ABX-CJ-IRIS', 'ABX-CJ-MANG', 'MANG_Assembly_Day_43']
#prodigal_faas = ['output/refine_bins/predict_genes/ABX-CJ-IRIS_predicted_genes.faa', 'output/refine_bins/predict_genes/ABX-CJ-MANG_predicted_genes.faa',
#                 'output/refine_bins/predict_genes/MANG_Assembly_Day_43_predicted_genes.faa']
#prodigal_fnas = ['output/refine_bins/predict_genes/ABX-CJ-IRIS_predicted_genes.fna', 'output/refine_bins/predict_genes/ABX-CJ-MANG_predicted_genes.fna',
#                 'output/refine_bins/predict_genes/MANG_Assembly_Day_43_predicted_genes.fna']
#prodigal_gffs = ['output/refine_bins/predict_genes/ABX-CJ-IRIS_predicted_genes.gff', 'output/refine_bins/predict_genes/ABX-CJ-MANG_predicted_genes.gff',
#                 'output/refine_bins/predict_genes/MANG_Assembly_Day_43_predicted_genes.gff']
#
#contigs = ['output/assemble/megahit/ABX-CJ-IRIS.contigs.fasta', 'output/assemble/megahit/ABX-CJ-MANG.contigs.fasta', 
#           'output/assemble/megahit/MANG_Assembly_Day_43.contigs.fasta']
#outdir = pathlib.Path("annotate_test")
#temp_outdir = pathlib.Path("annotate_temp")
#threads = 1

annot_file = snakemake.input.get('annotations', '')
concat_fna = snakemake.input.get('concat_fasta', '')
contigs = snakemake.input.get('contigs', '')
prodigal_faas = snakemake.input.get('faas', '')
prodigal_fnas = snakemake.input.get('fnas', '')
prodigal_gffs = snakemake.input.get('gffs', '')
contig_samples = snakemake.params.get('contig_samples', '')
outdir = pathlib.Path(snakemake.params.get('outdir', ''))
temp_outdir = pathlib.Path(snakemake.params.get('tempdir', ''))
threads = snakemake.threads


#if outdir.exists():
#    for i in outdir.iterdir():
#        i.unlink()
#    outdir.rmdir()
#outdir.mkdir(parents=True, exist_ok=True)

#if temp_outdir.exists():
#    for i in temp_outdir.iterdir():
#        if i.is_dir():
#            for j in i.iterdir():
#                j.unlink()
#            i.rmdir()
#        else:
#            i.unlink()
#    temp_outdir.rmdir()
#temp_outdir.mkdir(parents=True, exist_ok=True)

log_file_path = temp_outdir.joinpath("annotate.log")
logger = annotate_bins.logging.getLogger("annotation_log")
annotate_bins.setup_logger(logger, log_file_path)

# Cluster representatives are in the left column, members in the right
orf_clustfile = "output/refine_bins/dereplicated_genes/dereplicated_genes_cluster.tsv"

# Maybe create a dict of the faa, fna, etc. for each contig sample?

samp_ins = {}
for contig_samp in contig_samples:
    samplist = [f for flist in [contigs, prodigal_faas, prodigal_fnas, prodigal_gffs] for f in flist if contig_samp in f]
    samp_ins.update({contig_samp: {'contigs': samplist[0], 'faa': samplist[1], 'fna': samplist[2], 'gff': samplist[3]}})

print('Reading annotations', file=sys.stderr)
# The annotations outputted by DRAM are parsed
annotations = annotations = pd.read_csv(annot_file, sep='\t', engine='c', index_col=0, dtype='object')
# The actual contig sample ID is found for each ORF. Splitting by _ isn't used in case the user wants to use underscores in their naming scheme
orf_clust = pd.read_csv(orf_clustfile, header=None, engine='c', sep='\t', names=['ORF_Clust_ID', 'ORF_ID'], dtype='object').set_index('ORF_ID')

# Smaller practice versions of the files
#annotations = annotations = pd.read_csv(annot_file, sep='\t', engine='c', index_col=0, dtype='object', nrows=1000)
#orf_clust = pd.read_csv(orf_clustfile, header=None, engine='c', sep='\t', names=['ORF_Clust_ID', 'ORF_ID'], dtype='object', nrows=1000).set_index('ORF_ID')

orf_source = {}
for orf in orf_clust.index:
    for contig_samp in contig_samples:
        if contig_samp in orf:
            orf_source.update({orf: contig_samp})
            break
orf_source_df = pd.DataFrame.from_dict(orf_source, orient='index', columns=['Contig_Sample']).rename_axis('ORF_ID')
orf_clust = pd.concat([orf_clust, orf_source_df], axis=1, sort=False).reset_index()

# The clustered ORFs (which were dereplicated by ANI similarity) are given the same annotations as their representative ORF
clust_agg = orf_clust.groupby('ORF_Clust_ID', as_index=True, sort=False).agg(tuple)
annotated_clust = clust_agg.join(annotations).explode(['ORF_ID', 'Contig_Sample']).reset_index().set_index('ORF_ID')\
        .drop(columns='fasta').rename(columns={'Contig_Sample': 'fasta'})

# Columns are reordered so the non-DRAM columns are moved to the end
annotated_clust = annotated_clust[annotated_clust.columns[annotated_clust.columns != 'ORF_Clust_ID'].to_list() + ['ORF_Clust_ID']]


print('Creating annotations for all genes', file=sys.stderr)
# Contig information including the start and end positions of each ORF are extracted from the full gene fna file
# (i.e., the file used to dereplicate genes)
full_annotations = pd.concat([annotate_bins.get_gene_data(concat_fna), annotated_clust], axis=1, sort=False)

#TODO: Use the full version of this above when the final annotations are ready
full_annotations['scaffold'].replace({f'{prefix}_': '' for prefix in full_annotations['fasta'].unique()}, regex=True, inplace=True)
full_annotations['file_header'] = full_annotations['scaffold'] + '_' + full_annotations['gene_position'].astype(str)
annotations_list = []
for contig_samp in contig_samples:
    samp_outdir = temp_outdir.joinpath(contig_samp)
    samp_outdir.mkdir(parents=True, exist_ok=True)
    annot = full_annotations[full_annotations['fasta'] == contig_samp].reset_index(names='Full_ID').set_index('file_header')
    renamed_scaffolds = samp_outdir.joinpath(f'contigs.renamed.faa').as_posix()
    annotated_faa = samp_outdir.joinpath(f'genes.annotated.faa').as_posix()
    annotated_fna = samp_outdir.joinpath(f'genes.annotated.fna').as_posix()
    renamed_gffs = samp_outdir.joinpath(f'genes.annotated.gff').as_posix()
    trna_loc = samp_outdir.joinpath(f'trnas.tsv')
    rrna_loc = samp_outdir.joinpath(f'rrnas.tsv')
    current_gbk = outdir.joinpath(f'{contig_samp}.gbk') # This never actually gets made because it takes waaaay too long. Shhh....
    
    print('Annotating input files', file=sys.stderr)
    annotate_bins.rename_fasta(samp_ins[contig_samp]['contigs'], renamed_scaffolds, prefix=contig_samp)
    annotate_bins.create_annotated_fasta(samp_ins[contig_samp]['faa'], annot, annotated_faa, name=contig_samp)
    annotate_bins.create_annotated_fasta(samp_ins[contig_samp]['fna'], annot, annotated_fna, name=contig_samp)
    annotate_bins.annotate_gff(samp_ins[contig_samp]['gff'], renamed_gffs, annot, prefix=contig_samp)

    print('Annotating tRNAs', file=sys.stderr)
    len_dict = {
        i.metadata["id"]: len(i)
        for i in annotate_bins.read_sequence(renamed_scaffolds, format="fasta")
    }
    trna_table = annotate_bins.run_trna_scan(
        renamed_scaffolds,
        samp_outdir,
        contig_samp,
        logger,
        threads=threads,
        verbose=True,
    )
    trna_table.to_csv(trna_loc, sep="\t", index=False)
    annotate_bins.add_intervals_to_gff(
        trna_loc, renamed_gffs, len_dict, annotate_bins.make_trnas_interval, "Name", logger
    )

    print('Annotating rRNAs', file=sys.stderr)
    rrna_table = annotate_bins.run_barrnap(
        renamed_scaffolds, contig_samp, logger, threads=threads, verbose=True
    )
    rrna_table.to_csv(rrna_loc, sep="\t", index=False)
    annotate_bins.add_intervals_to_gff(
        rrna_loc, renamed_gffs, len_dict, annotate_bins.make_rrnas_interval, "scaffold", logger
    )

    annotations_list.append(
        annotate_bins.Annotation(
            name=contig_samp,
            scaffolds=renamed_scaffolds,
            genes_faa=annotated_faa,
            genes_fna=annotated_fna,
            gff=renamed_gffs,
            gbk=current_gbk,
            annotations=annot,
            trnas=trna_loc,
            rrnas=rrna_loc,
        )
    )

print('Annotation complete. Merging annotations', file=sys.stderr)
all_annotations = full_annotations.sort_values(
        ["fasta", "scaffold", "gene_position"]
    )

# merge gene files
annotate_bins.merge_files(
    [i.genes_fna_loc for i in annotations_list if i is not None],
    outdir.joinpath("genes.fna").as_posix(),
)
annotate_bins.merge_files(
    [i.genes_faa_loc for i in annotations_list if i is not None],
    outdir.joinpath("genes.faa").as_posix(),
)
annotate_bins.merge_files(
    [i.scaffolds_loc for i in annotations_list if i is not None],
    outdir.joinpath("scaffolds.fna").as_posix(),
)
annotate_bins.merge_files(
    [i.gff_loc for i in annotations_list if i is not None],
    outdir.joinpath("genes.gff").as_posix(),
    True,
)
trnas_locs = [
    i.trnas_loc
    for i in annotations_list
    if i is not None
    if i.trnas_loc is not None
]
if len(trnas_locs) > 0:
    annotate_bins.merge_files(trnas_locs, outdir.joinpath("trnas.tsv").as_posix(), True)
rrnas_locs = [
    i.rrnas_loc
    for i in annotations_list
    if i is not None
    if i.rrnas_loc is not None
]
if len(rrnas_locs) > 0:
    annotate_bins.merge_files(rrnas_locs, outdir.joinpath("rrnas.tsv").as_posix(), True)

all_annotations.to_csv(outdir.joinpath("annotations.tsv"), sep="\t")
shutil.rmtree(temp_outdir)
print("Yay! It's done!", file=sys.stderr)

#full_annotations.set_index('file_header', inplace=True)
#full_annotations['rank'].fillna('F', inplace=True)
#print(full_annotations.loc['k141_66447_1'])
#
#
###########################3
#full_annotations = full_annotations[full_annotations['fasta'] == 'ABX-CJ-IRIS'].reset_index().set_index('file_header')
###########################
#print(full_annotations[full_annotations['scaffold'] == 'k141_66447'])
#print(full_annotations.loc['k141_66447_1'])

###############################33
#annotate_bins.create_annotated_fasta(prodigal_faas[0], full_annotations, annotated_faa, name='ABX-CJ-IRIS')
###################################

#print(full_annotations[full_annotations['scaffold'] == 'k141_66447_1'])
#no = []
#for orfrep in annotations.index:
#    if clust_dict.get(orfrep):
#        count += 1
#
#    #for orf in clust_dict[orfrep]:
#    #    print(
#
#    if count > 10:
#        break
##        no.append(i)
#
