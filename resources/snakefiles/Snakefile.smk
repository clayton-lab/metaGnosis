import numpy as np
import pandas as pd
from os.path import join

configfile: "config.yaml"

samples_fp = config['samples']

units_fp = config['units']

reads = config['reads']

sample_table = pd.read_csv(samples_fp,
						   sep='\t',
						   header=0,
						   index_col=0,
						   na_filter=False)

units_table = pd.read_csv(units_fp, sep='\t', header=0)
units_table.set_index(['Sample_ID', 'Unit_ID'], inplace=True)
if pd.unique(units_table.index).shape[0] != units_table.shape[0]:
	raise ValueError('Every Sample_ID+Unit_ID combination must be unique for proper sample indexing!')

if len(config['assemblers']) > 1:
    log_msg = '''WARNING: More than one assembler was specified in config.yaml. This is encouraged for comparing performance between assemblers, but is not currently supported beyond assembly qc. Subsequent steps (mapping, binning, etc.) will only be performed using contigs from the first assembler.
              '''
    print(log_msg, file=sys.stderr)
selected_assembler = config['assemblers'][0]
if len(config['mappers']) > 1:
    log_msg = '''WARNING: More than one mapper was specified in config.yaml. This is encouraged for comparing performance between mappers, but is not currently supported beyond bin refinement. Subsequent steps (bin annotation, quantification, etc.) will only be performed using mapping from the first mapper.
              '''
    print(log_msg, file=sys.stderr)

selected_mapper=config['mappers'][0]

samples = sample_table.index
units = units_table.index

# Make it so that any empty ID columns are filled in with Sample_ID (future option)
#sample_table = sample_table.apply(lambda row: row.fillna(value=row.name), axis=1)

def parse_groups(group_series):
	groups = {}
	for sample, grps in group_series.items():
		if not grps or grps == 'None':
			continue
		grp_list = grps.split(',')
		for grp in grp_list:
			if grp not in groups:
				groups.update({grp: [sample]})
			else:
				groups[grp].append(sample)
	return(groups)

def make_pairings(sample_df):
	map_groups = {}
	for row in sample_df.itertuples():
		if row.Contig_ID and row.Contig_ID != 'None':
			if not map_groups.get(row.Mapping_Group):
				map_groups.update({row.Mapping_Group: {'contig_id': [row.Contig_ID], 'read_id': [row.Index]}})
			else:
				map_groups[row.Mapping_Group]['contig_id'].append(row.Contig_ID)

		if row.Index not in map_groups[row.Mapping_Group]['read_id']:
			map_groups[row.Mapping_Group]['read_id'].append(row.Index)

	contig_pairings = {}
	for grp in map_groups.keys():
		for ctg in set(map_groups[grp]['contig_id']):
			if not contig_pairings.get(ctg):
				contig_pairings.update({ctg: list(set(map_groups[grp]['read_id']))})
			else:
				for read in set(map_groups[grp]['read_id']):
					if read not in contig_pairings[ctg]:
						contig_pairings[ctg].append(read)

	return (map_groups, contig_pairings)

read_groups = [sample for sample in sample_table.index]
contig_groups = parse_groups(sample_table['Contig_ID'])
map_groups, contig_pairings = make_pairings(sample_table)

# contig_groups is useful for tracing assembled contigs back to samples (i.e., for co-assembly),
# and contig_pairings is useful for mapping sample reads to contigs (multiple reads to multiple contigs)
print('Read samples: %s\n' % read_groups)
print('Contig samples: %s\n' % contig_groups)
print('Mapping groups: %s\n' % map_groups)
print('Contig Pairings: %s\n' % contig_pairings)

# For now, samples can only be part of a single mapping group (i.e., no duplicate sample rows). Future versions could get around
# this by re-writing the pipeline to behave like units.tsv (which can be duplicate samples), but no time for that yet.
def get_read(sample, unit, read):
	return(units_table.loc[(sample, unit), read])

def get_contig_id(sample, sample_table):
	contigs = sample_table.replace(to_replace={'^None$': np.NaN, '^$': np.NaN}, regex=True).dropna(axis='index')
	return(contigs.loc[sample, 'Contig_ID'])

#def check_config(*config_items):
#    if all([key.get() for key in config_items]):
#        return config_items
#    else:
#        return []

# To make host filter optional, need to modify: (maybe build_db, qc), assemble, prototype, profile, mapping, and quant_bins
include: "resources/snakefiles/build_db.smk"
include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/prototype_selection.smk"
include: "resources/snakefiles/profile.smk"
include: "resources/snakefiles/mapping.smk"
include: "resources/snakefiles/binning.smk"
include: "resources/snakefiles/refine_bins.smk"
include: "resources/snakefiles/annotate_bins.smk"
include: "resources/snakefiles/quant_bins.smk"

# Can trigger metaquast and metaquast assemble by specifying "output/assemble/multiqc_metaquast/multiqc.html" below.
# Should do this to make sure all of assemble.smk works before adding the bin stuff (which can also be specified here to trigger it).

# Each of the output files below are essentially a switch for their respective snakefiles to enable/disable functionality:
# multiqc                   - does everything up to filtering out host reads
# multiqc_assemble          - does assembly with either/both selected assemblers
# selected_prototypes       - looks at similarity between reads to find potential prototypes (i.e., a representative sample for contig binning)
# merged_abundance_table    - read profiling with metaphlan

# Goal: have a list in config.yaml with general modules that the user can comment/uncomment to their liking.
# Then rule all will be an expand with {output_prefix}/{module}, that automatically defines which end-files to produce
rule all:
    input:
        #expand("{output_prefix}/qc/multiqc/multiqc.html",
        #        output_prefix=config['user_paths']['output_prefix'],
        #        sample=samples)
#        expand("{output_prefix}/qc/multiqc/multiqc.html",
#               "{output_prefix}/assemble/multiqc_assemble/multiqc.html",
#                output_prefix=config['user_paths']['output_prefix'] if config['user_paths']['output_prefix'] else "output")
        "output/qc/multiqc/multiqc.html",
        "output/assemble/multiqc_assemble/multiqc.html",

        #"output/prototype_selection/sourmash_plot",
        #"output/prototype_selection/prototype_selection/selected_prototypes.yaml",
        #"output/profile/metaphlan/merged_abundance_table.txt",
        #"output/profile/kraken2/merged_kreport2mpa_table.txt",
        "output/refine_bins/summarize_bins/bin_filter_summary.txt",
        #"output/annotate_bins/annotate_genes/annotations.tsv",
        #"output/refine_bins/dereplicated_bins/dereplicated_bins.fa",
        "output/mapping_qc/multiqc/multiqc_mapping.html",
        "output/quant_bins/quantified_bin_abundance.tsv",
        "output/quant_bins/quantified_gene_abundance.tsv"

        #lambda wildcards: expand("output/refine_bins/{mapper}/run_CheckM/run_checkm/{contig_sample}",
        #                     mapper=config['mappers'],
        #                     contig_sample=contig_pairings.keys())

#output/binning/maxbin2/{mapper}/bin_fastas/{contig_sample}/
        #lambda wildcards: expand("output/binning/metabat2/{mapper}/bin_fastas/{contig_sample}/",
        #                          mapper=config['mappers'],
        #                          contig_sample=contig_pairings.keys(),
        #                          read_sample=contig_pairings.values()),
