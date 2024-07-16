from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob

# This rule isn't run, but it can replace the x_Fasta_to_Scaffolds2Bin rules below. The binning.smk
# file has been modified to have semi-generic output directories for bins, allowing a generic
# Fasta_to_Scaffolds2Bin rule. 

#Create a rule to gather all binning fastas to get rid of wildcard in rule Fasta_to_contigs2Bin. This is because having a wildcard in the first rule of a snakefile generates a workflow error, which causes the rule to be skipped.

#Determine the bin-fastas for all files
#BINNER, MAPPER, CONTIG_SAMPLE, = glob_wildcards("output/binning/{binner}/{mapper}/bin_fastas/{contig_sample}/")
#BINNER = config['binners']
#MAPPER = config['mappers']
#ASSEMBLE = config['assemblers'] 

rule Fasta_to_Contig2Bin:
    """
    Uses Fasta_to_Contig2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins="output/binning/{binner}/{mapper}/bin_fastas/{contig_sample}"
    output:
        scaffolds="output/selected_bins/{binner}/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    params:
        extension=lambda wildcards: "fasta" if wildcards.binner == "maxbin2" else "fa"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}.log"
    shell:
       """
           Fasta_to_Scaffolds2Bin.sh \
           -i {input.bins} \
           -e {params.extension} > {output.scaffolds} \
           2> {log}
       """

# This rule still needs to be harmonized with the generic Fasta_to_Scaffolds2Bin rule. The
# input binners can be reduced to a single generic 'bins' parameter, similar to the 'contigs'
# input paramater. Then it will just be a matter of fixing the output and possibly shell command
# to handle the generic 'bins' input parameter
rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                                            assembler = config['assemblers'],
                                            contig_sample = wildcards.contig_sample),
        bins = lambda wildcards: expand("output/selected_bins/{binner}/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                                         mapper = config['mappers'],
                                         binner = config['binners'],
                                         contig_sample = wildcards.contig_sample)
                #contig_sample = contig_pairings.keys())
    output:
        out="output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.txt"
    params:
        basename = "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine'],
        binner_paths=lambda wildcards, input: ",".join(input.bins),
        binners=",".join(config['binners']),
    conda:
        "../env/selected_bins.yaml"
    threads:
        config['threads']['run_DAS_Tool']
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}.log"
    shell:
        """ 
           DAS_Tool \
           --bins {params.binner_paths} \
           --contigs {input.contigs} \
           --outputbasename {params.basename} \
           --labels {params.binners} \
           --write_bins 1 \
           --write_bin_evals 1 \
           --threads {threads} \
           --search_engine {params.search_engine} \
           2> {log} 1>&2
        """

rule consolidate_summary:
	input:
	  source = expand("output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.txt",
			mapper = config['mappers'],
			contig_sample =contig_pairings.keys())		   	
	output:
	  summary="output/selected_bins/DAS_Tool_summary.txt"
	log:
	   "output/logs/selected_bins/DAS_Tool_summary.log"
	shell:	
	   "cat {input.source} > {output}"

#rule aggregate:
#	input:
#	    lambda wildcards: expand("output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DasTool_summary.txt",
#				mapper=config['mappers'],
#				contig_sample=contig_pairings.keys())

