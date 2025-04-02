from os.path import splitext
# The rest of the workflow breaks when multiple assemblers are specified. If multiple assemblers were specified,
# only the first assembler's contigs are used from this point onward. Future workarounds include concatenating output from
# multiple assemblers into a single .fasta file before indexing them, rewriting the rest of the workflow to
# accomodate multiple assemblers (similar to how it handles multiple mappers), or creating a rule that selects
# the best assembler automatically. A high-quality assembly would generally have higher N50, larger contigs, 
# fewer contigs, higher assembly length approaching the reference genome length, fewer misassemblies, higher 
# genome fraction coverage, lower duplication ratio, consistent GC content, and fewer ambiguous bases.
# The only challenge would be automatically selecting an assembler from multiqc_assemble and/or multiqc_metaquast
# simultaneously. Definitely doable, but would require some time to figure out.

#TODO: Reorganize directory structure so there is a contig and a gene mapping output directory
rule index_contigs_bt2:
    input:
        lambda wildcards: f"output/assemble/{selected_assembler}/{wildcards.contig_sample}.contigs.fasta"
        #contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
        #                                    assembler=config['assemblers'],
        #                                    contig_sample=wildcards.contig_sample)
        #contigs = lambda wildcards: get_contigs(wildcards.contig_sample,
        #                                        binning_df)
    output:
        temp(multiext("output/mapping/bowtie2/indexed_contigs/{contig_sample}",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/mapping/bowtie2/indexed_contigs/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mapping/bowtie2/indexed_contigs/{contig_sample}_benchmark.txt"
    conda:
        "../env/mapping.yaml"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/mapping/bowtie2/indexed_contigs/{contig_sample}"
    threads:
        config['threads']['bowtie2_build']
    shell:
        """
        {params.bt2b_command} --threads {threads} \
        {input} {params.indexbase} 2> {log} 1>&2
        """

rule map_reads_to_contigs_bt2:
    """
    Maps reads to contig files using bowtie2.
    """
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_contigs_bt2.output
    output:
        aln=temp("output/mapping/bowtie2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        ref="output/mapping/bowtie2/indexed_contigs/{contig_sample}",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['map_reads']
    benchmark:
        "output/benchmarks/mapping/bowtie2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.benchmark.txt"
    log:
        "output/logs/mapping/bowtie2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        # Map reads against reference genome
        {params.bt2_command} {params.extra} -p {threads} -x {params.ref} \
          -1 {input.reads[0]} -2 {input.reads[1]} \
          2> {log} | samtools view -bS - > {output.aln}

        """
rule index_contigs_minimap2:
    input:
        lambda wildcards: f"output/assemble/{selected_assembler}/{wildcards.contig_sample}.contigs.fasta"
    output:
        index = temp("output/mapping/minimap2/indexed_contigs/{contig_sample}.mmi")
    log:
        "output/logs/mapping/minimap2/indexed_contigs/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mapping/minimap2/indexed_contigs/{contig_sample}_benchmark.txt"
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['minimap2_index']
    shell:
        """
        minimap2 -d {output.index} {input} -t {threads} 2> {log} 1>&2

        """

rule map_reads_to_contigs_minimap2:
    """
    Maps reads to contig files using minimap2.
    """
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_contigs_minimap2.output.index
    output:
        aln=temp("output/mapping/minimap2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        x=config['params']['minimap2']['x'],
        k=config['params']['minimap2']['k']
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['minimap2_map_reads']
    benchmark:
        "output/benchmarks/mapping/minimap2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}_benchmark.txt"
    log:
        "output/logs/mapping/minimap2/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        # Map reads against contigs
        minimap2 -a {input.db} {input.reads} -x {params.x} -K {params.k} -t {threads} \
        2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_contig_index_bam:
    """
    Sorts and indexes bam file.
    """
    input:
        aln="output/mapping/{mapper}/mapped_reads/mapped_to_contigs/{read_sample}_Mapped_To_{contig_sample}.bam"
    output:
        bam=temp("output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam"),
        index=temp("output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam.bai")
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/mapping/{mapper}/sort_index_bam/contigs/{read_sample}_Mapped_To_{contig_sample}.txt"
    log:
        "output/logs/mapping/{mapper}/sort_index_bam/contigs/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}
        samtools index -b -@ {threads} {output.bam} 2>> {log}
        """

rule calculate_contig_coverage:
    """
       Commands to generate a coverage table using `samtools coverage` for input into some subsequent steps
    """
    input:
        bams="output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam",
        index="output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam.bai"
    output:
        temp("output/mapping/{mapper}/coverage_tables/contigs/{read_sample}_Mapped_To_{contig_sample}_coverage.txt"),
    conda:
        "../env/mapping.yaml"
    benchmark:
        "output/benchmarks/mapping/{mapper}/calculate_coverage/contigs/{read_sample}_Mapped_To_{contig_sample}.txt"
    log:
        "output/logs/mapping/{mapper}/calculate_coverage/contigs/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
       """
          samtools coverage {input.bams} | \
          tail -n +2 | \
          sort -k1 | \
          cut -f 1,7 > {output} 2> {log}
       """

use rule index_contigs_bt2 as index_genes_bt2 with:
    input:
        "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_rep_seq.fasta"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/mapping/bowtie2/indexed_genes/indexed_genes"
    output:
        temp(multiext("output/mapping/bowtie2/indexed_genes/indexed_genes",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/mapping/bowtie2/indexed_genes/index_genes.log"
    benchmark:
        "output/benchmarks/mapping/bowtie2/indexed_genes/index_genes_benchmark.txt"

use rule map_reads_to_contigs_bt2 as map_reads_to_genes_bt2 with:
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_genes_bt2.output
    output:
        aln=temp("output/mapping/bowtie2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.bam")
    params:
        ref="output/mapping/bowtie2/indexed_genes/indexed_genes",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    benchmark:
        "output/benchmarks/mapping/bowtie2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.benchmark.txt"
    log:
        "output/logs/mapping/bowtie2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.log"

use rule index_contigs_minimap2 as index_genes_minimap2 with:
    input:
        "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_rep_seq.fasta"
    output:
        index = temp("output/mapping/minimap2/indexed_genes/indexed_genes.mmi")
    log:
        "output/logs/mapping/minimap2/indexed_genes/index_genes.log"
    benchmark:
        "output/benchmarks/mapping/minimap2/indexed_genes/index_genes_benchmark.txt"


use rule map_reads_to_contigs_minimap2 as map_reads_to_genes_minimap2 with:
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_genes_minimap2.output.index
    output:
        aln=temp("output/mapping/minimap2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.bam")
    benchmark:
        "output/benchmarks/mapping/minimap2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.benchmark.txt"
    log:
        "output/logs/mapping/minimap2/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.log"

# TODO: Mark output as temp when finished.
use rule sort_contig_index_bam as sort_gene_index_bam with:
    input:
        aln="output/mapping/{mapper}/mapped_reads/mapped_to_genes/{read_sample}_Mapped_To_Genes.bam"
    output:
        bam=temp("output/mapping/{mapper}/sorted_bams/genes/{read_sample}_Mapped_To_Genes.bam"),
        index=temp("output/mapping/{mapper}/sorted_bams/genes/{read_sample}_Mapped_To_Genes.bam.bai")
    benchmark:
        "output/benchmarks/mapping/{mapper}/sort_bams/genes/{read_sample}_Mapped_To_Genes.txt"
    log:
        "output/logs/mapping/{mapper}/sort_bams/genes/{read_sample}_Mapped_To_Genes.log"

#TODO: Mark the output files for this as temp when you figure out what to do.
rule calculate_gene_coverage:
    input:
        bams="output/mapping/{mapper}/sorted_bams/genes/{read_sample}_Mapped_To_Genes.bam",
        index="output/mapping/{mapper}/sorted_bams/genes/{read_sample}_Mapped_To_Genes.bam.bai"
    output:
        coverage_table=temp("output/mapping/{mapper}/coverage_tables/genes/{read_sample}_gene_coverage.txt"),
        lengths=temp("output/mapping/{mapper}/lengths/genes/{read_sample}_gene_lengths.txt")
    conda:
        "../env/mapping.yaml"
    log:
        "output/logs/mapping/{mapper}/calculate_coverage/genes/{read_sample}_Mapped_To_Genes.log"
    shell:
       """
          samtools coverage {input.bams} | \
          tail -n +2 | \
          sort -k1 | \
          cut -f 1,4 > {output.coverage_table} 2> {log}

          samtools idxstats {input.bams} | \
          cut -f 1,2 > {output.lengths} 2> {log}
       """
use rule index_contigs_bt2 as index_bins_bt2 with:
    input:
        lambda wildcards: "output/annotate_bins/annotate_bin_pathways/merged_annotations/scaffolds.fna"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/mapping/bowtie2/indexed_bins/indexed_bins"

    output:
        temp(multiext("output/mapping/bowtie2/indexed_bins/indexed_bins",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/mapping/bowtie2/indexed_bins/index_bins.log"
    benchmark:
        "output/benchmarks/mapping/bowtie2/indexed_bins/index_bins_benchmark.txt"

use rule map_reads_to_contigs_bt2 as map_reads_to_bins_bt2 with:
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_bins_bt2.output
    output:
        aln=temp("output/mapping/bowtie2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.bam")
    params:
        ref="output/mapping/bowtie2/indexed_bins/indexed_bins",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    benchmark:
        "output/benchmarks/mapping/bowtie2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.benchmark.txt"
    log:
        "output/logs/mapping/bowtie2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.log"

use rule index_contigs_minimap2 as index_bins_minimap2 with:
    input:
        "output/annotate_bins/annotate_bin_pathways/merged_annotations/scaffolds.fna"
    output:
        index = temp("output/mapping/minimap2/indexed_bins/indexed_bins.mmi")
    log:
        "output/logs/mapping/minimap2/indexed_bins/index_bins.log"
    benchmark:
        "output/benchmarks/mapping/minimap2/indexed_bins/index_bins_benchmark.txt"


use rule map_reads_to_contigs_minimap2 as map_reads_to_bins_minimap2 with:
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_bins_minimap2.output.index
    output:
        aln=temp("output/mapping/minimap2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.bam")
    benchmark:
        "output/benchmarks/mapping/minimap2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.benchmark.txt"
    log:
        "output/logs/mapping/minimap2/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.log"

# TODO: Mark output as temp when finished.
use rule sort_contig_index_bam as sort_bin_index_bam with:
    input:
        aln="output/mapping/{mapper}/mapped_reads/mapped_to_bins/{read_sample}_Mapped_To_Bins.bam"
    output:
        bam=temp("output/mapping/{mapper}/sorted_bams/bins/{read_sample}_Mapped_To_Bins.bam"),
        index=temp("output/mapping/{mapper}/sorted_bams/bins/{read_sample}_Mapped_To_Bins.bam.bai")
    benchmark:
        "output/benchmarks/mapping/{mapper}/sort_bams/bins/{read_sample}_Mapped_To_Bins.txt"
    log:
        "output/logs/mapping/{mapper}/sort_bams/bins/{read_sample}_Mapped_To_Bins.log"

#TODO: Mark the output files for this as temp when you figure out what to do.
use rule calculate_contig_coverage as calculate_bin_coverage with:
    input:
        bams="output/mapping/{mapper}/sorted_bams/bins/{read_sample}_Mapped_To_Bins.bam",
        index="output/mapping/{mapper}/sorted_bams/bins/{read_sample}_Mapped_To_Bins.bam.bai"
    output:
        temp("output/mapping/{mapper}/coverage_tables/bins/{read_sample}_bin_coverage.txt"),
    log:
        "output/logs/mapping/{mapper}/calculate_coverage/bins/{read_sample}_Mapped_To_Bins.log"
    benchmark:
        "output/benchmarks/mapping/{mapper}/calculate_coverage/contigs/{read_sample}_Mapped_To_Bins.txt"

