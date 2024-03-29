from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob

# This rule isn't run, but it can replace the x_Fasta_to_Scaffolds2Bin rules below. The binning.smk
# file has been modified to have semi-generic output directories for bins, allowing a generic
# Fasta_to_Scaffolds2Bin rule. 
rule Fasta_to_Scaffolds2Bin:
	"""
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
	input:
		bins = lambda wildcards: expand("output/binning/{binner}/{mapper}/bin_fastas/{contig_sample}/",
				mapper = config['mappers'],
				binner = config['binners'],
				contig_sample = wildcards.contig_sample)
	output:
		scaffolds2bin="output/selected_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}_scaffolds2bin.tsv"
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
            -e fa > {output.scaffolds2bin}
        """

rule metabat2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/metabat2/{mapper}/bin_fastas/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/{mapper}/scaffolds2bin/metabat2/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Scaffolds2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """


rule maxbin2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/maxbin2/{mapper}/bin_fastas/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/{mapper}/scaffolds2bin/maxbin2/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Scaffolds2Bin.sh \
            -i {input.bins} \
            -e fasta > {output.scaffolds2bin}
        """

rule concoct_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/bin_fastas/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/{mapper}/scaffolds2bin/concoct/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Scaffolds2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """

# rule concoct_Fasta_to_Scaffolds2Bin:
#     """
#     Uses perl to create a scaffolds2bin.tsv file from a clustering_merged.csv file.
#     """
#     input:
#         merged = lambda wildcards: expand("output/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_clustering_merged.csv",
#                 mapper = config['mappers'],
#                 contig_sample = wildcards.contig_sample)
#     output:
#         scaffolds2bin="output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
#     conda:
#         "../env/selected_bins.yaml"
#     benchmark:
#         "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
#     shell:
#         """
#             perl -pe "s/,/\tconcoct_bins./g;" {input.merged} > {output.scaffolds2bin}
#         """


# This rule still needs to be harmonized with the generic Fasta_to_Scaffolds2Bin rule. The
# input binners can be reduced to a single generic 'bins' parameter, similar to the 'contigs'
# input paramater. Then it will just be a matter of fixing the output and possibly shell command
# to handle the generic 'bins' input parameter
rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        metabat2 = lambda wildcards: expand("output/selected_bins/{mapper}/scaffolds2bin/metabat2/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        maxbin2 = lambda wildcards: expand("output/selected_bins/{mapper}/scaffolds2bin/maxbin2/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        concoct = lambda wildcards: expand("output/selected_bins/{mapper}/scaffolds2bin/concoct/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                    selected_assembler = selected_assembler,
                    contig_sample = wildcards.contig_sample)
    output:
        out="output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.txt"
    params:
        basename = "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine']
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
            --bins {input.metabat2},{input.maxbin2},{input.concoct} \
            --contigs {input.contigs} \
            --outputbasename {params.basename} \
            --labels metabat2,maxbin2,concoct \
            --write_bins 1 \
            --write_bin_evals 1 \
            --threads {threads} \
            --search_engine {params.search_engine}
        """


# This rule initially had output.out as the output file in the fasta_files for loop,
# but that finished with a error. The changes I made make the pipeline finish without
# errors, but don't completely finish the purpose of this rule (described in triple
# quotes below). So the last 2 rules still need to be tweaked and probably combined.
# Renaming shouldn't be a problem now, so all that's left is to clean up the output
# directory.

# One idea (that could also be implemented with binning.smk) is to have DAS_Tool_Fastas,
# and another folder for the remaining DAS_Tool output, in the root selected_bins directory
# along with scaffolds2bin (already there). Then DAS_Tool_Fastas and the other output folder
# could have sample-specific subdirectories, with the bins/output in their respective sample
# directories.

# So something like /selected_bins/DAS_Tool_Fastas/{mapper}/{Sample_1,Sample_2,etc.}/x.bin.fa
# I tried a similar scheme in binning.smk with some success, though it's not perfect.

# One other issue (a matter of personal preference, really), is if the directory order should be
# {Das_Tool_Output}/{mapper} or {mapper}/{Das_Tool_Output}. Binning.smk has {binner}/{mapper},
# so whatever the choice, it should be consistent for both.
rule consolidate_DAS_Tool_bins:
    """
    Consolidates and renames bin fastas generated by DAS_Tool into a single folder
    """
    input:
        "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.txt"
    output:
        done = touch("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done")
    log:
        "output/logs/selected_bins/{mapper}/consolidate_DAS_Tool_bins/{contig_sample}.log"
    run:
        sample = wildcards.contig_sample 
        fasta_dir = join(dirname(input[0]),
                         sample + '_DASTool_bins')
        output_dir = dirname(output.done)

        fasta_files = glob(join(fasta_dir, '*.fa'))

        for file in fasta_files:
            copyfile(file,
                     join(output_dir,
                          sample + '_' + basename(file)))

# This doesn't seem to be ran in the pipeline, and it seems to build from the input of the previous rule
rule consolidate_DAS_Tool_bins_all:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())


