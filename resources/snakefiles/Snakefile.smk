import numpy as np
import pandas as pd
from os.path import join


# binning_fp = config['binning']
# 
# binning_df = pd.read_csv(binning_fp,
#						   header=0,
#						   index_col=0,
#						   sep='\t',
#						   na_filter=False)
# 
# def get_read(sample, unit, read):
#	  return(units_table.loc[(sample, unit), read])
# 
# def parse_groups(group_series):
#	  groups = {}
#	  for sample, grps in group_series.items():
#		  if not grps:
#			  continue
#		  grp_list = grps.split(',')
#		  for grp in grp_list:
#			  if grp not in groups:
#				  groups[grp] = [sample]
#			  else:
#				  groups[grp].append(sample)
#	  return(groups)
# 
# def make_pairings(read_grp, ctg_grp):
#	  if read_grp.keys() != ctg_grp.keys():
#		  raise ValueError('Not all keys in both from and to groups!')
# 
#	  pairings = []
#	  contig_pairings = {}
#	  for grp in read_grp.keys():
#		  r = read_grp[grp]
#		  c = ctg_grp[grp]
# 
#		  for i in r:
#			  for  j in c:
#				  pairings.append((i, j))
#				  if j not in contig_pairings:
#					  contig_pairings[j] = [i]
#				  else:
#					  contig_pairings[j].append(i)
# 
#	  return(pairings, contig_pairings)
# 
# contig_groups = parse_groups(binning_df['Contig_Groups'])
# read_groups = parse_groups(binning_df['Read_Groups'])
# pairings, contig_pairings = make_pairings(read_groups, contig_groups)
# 
# print('Contig samples: %s' % contig_groups)
# print('Read samples: %s' % read_groups)
# print('Pairings: %s' % pairings)
# print('Contig Pairings: %s' % contig_pairings)
# 
# def get_contigs(sample, binning_df):
#	  return(binning_df.loc[sample, 'Contigs'])
# 
# include: "resources/snakefiles/qc.smk"
# include: "resources/snakefiles/assemble.smk"
# include: "resources/snakefiles/mapping.smk"
# include: "resources/snakefiles/binning.smk"
# include: "resources/snakefiles/selected_bins.smk"
# 
# 
# # This looks like the final target rule that triggers every other rule. The created .done file is empty and
# # has no point though, so it is probably just used to ensure the pipeline finishes without errors. Furthermore,
# # it's probably a sign that the pipeline isn't comple. selected_bins.smk has 2 final rules that were under
# # development when this repo was cloned, so figuring out what the goal of those rules was and finishing them
# # will be a next step to making this pipeline stable.
# rule select_bins:
#	  input:
#		  lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
#								   mapper=config['mappers'],
#								   contig_sample=contig_pairings.keys())
# 
# #
# # From what I can tell, these 2 rules never actually get ran during normal pipeline operations. They're probably
# # for debugging the mapping and binning steps, since you can call specific rules with snakemake. They can be removed
# # when we eventually merge this snakefile with the main snakefile.
# rule bin_all:
#	  input:
#		  expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
#				 mapper=config['mappers'],
#				 contig_sample=contig_pairings.keys()),
#		  expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
#				 mapper=config['mappers'],
#				 contig_sample=contig_pairings.keys()),
#		  expand("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
#				 mapper=config['mappers'],
#				 contig_sample=contig_pairings.keys())
# 
# rule map_all:
#	  input:
#		  expand("output/mapping/{mapper}/sorted_bams/{pairing[0]}_Mapped_To_{pairing[1]}.bam",
#				 mapper=config['mappers'],
#				 pairing=pairings)
# 
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

read_groups = {sample: [sample] for sample in sample_table.index}
contig_groups = parse_groups(sample_table['Contig_ID'])
map_groups, contig_pairings = make_pairings(sample_table)

# contig_groups is useful for tracing assemby contigs back to samples (i.e., for co-assembly),
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

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/prototype_selection.smk"
include: "resources/snakefiles/profile.smk"

print(get_contig_id('ABX-CJ-IRIS-14', sample_table))
# Can trigger metaquast and metaquast assemble by specifying "output/assemble/multiqc_metaquast/multiqc.html" below.
# Should do this to make sure all of assemble.smk works before adding the bin stuff (which can also be specified here to trigger it).
rule all:
	input:
		"output/qc/multiqc/multiqc.html",
		"output/assemble/multiqc_assemble/multiqc.html",
		"output/prototype_selection/sourmash_plot",
		"output/prototype_selection/prototype_selection/selected_prototypes.yaml",
		"output/profile/metaphlan/merged_abundance_table.txt",
