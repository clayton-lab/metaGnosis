rule metaspades:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

    """
    input:
        fastq1_list=lambda wildcards: expand("output/qc/host_filter/nonhost/{contig_seqs}.R1.fastq.gz",
                                  contig_seqs=contig_groups[wildcards.contig_sample]),
        fastq2_list=lambda wildcards: expand("output/qc/host_filter/nonhost/{contig_seqs}.R2.fastq.gz",
                                  contig_seqs=contig_groups[wildcards.contig_sample])
    output:
        contigs="output/assemble/metaspades/{contig_sample}.contigs.fasta",
    params:
        temp_dir=directory("output/{contig_sample}_temp"),
        fastq1=lambda wildcards, input: expand("--pe-1 1 {contig_seq}",
                                                      contig_seq=input.fastq1_list),
        fastq2=lambda wildcards, input: expand("--pe-2 1 {contig_seq}",
                                                      contig_seq=input.fastq2_list)
    conda:
        "../env/assemble.yaml"
    threads:
        config['threads']['spades']
    benchmark:
        "output/benchmarks/assemble/metaspades/{contig_sample}_benchmark.txt"
    log:
        "output/logs/assemble/metaspades/{contig_sample}.log"
    resources:
        mem_mb=config['mem_mb']['spades']
    shell:
        """
        # Make temporary output directory
        mkdir -p {params.temp_dir}

        # run the metaspades assembly
        metaspades.py --threads {threads} \
          -o {params.temp_dir}/ \
          --memory $(({resources.mem_mb}/1024)) \
          {params.fastq1} \
          {params.fastq2} \
          2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {params.temp_dir}/contigs.fasta {output.contigs}
        rm -rf {params.temp_dir}
        """

rule megahit:
    """

    Performs a metagenomic assembly on a sample using MEGAHIT.

    """
    input:
        fastq1_list=lambda wildcards: expand("output/qc/host_filter/nonhost/{contig_seqs}.R1.fastq.gz",
                                  contig_seqs=contig_groups[wildcards.contig_sample]),
        fastq2_list=lambda wildcards: expand("output/qc/host_filter/nonhost/{contig_seqs}.R2.fastq.gz",
                                  contig_seqs=contig_groups[wildcards.contig_sample])

#
#        lambda wildcards: get_contig_id(wildcards.contig_sample,
#                                                sample_table),
        #fastq1=rules.host_filter.output.nonhost_R1,
        #fastq2=rules.host_filter.output.nonhost_R2,
    output:
        contigs="output/assemble/megahit/{contig_sample}.contigs.fasta"
    params:
        presets=config['params']['megahit']['presets'],
        temp_dir=directory("output/{contig_sample}_temp"),

        # If there are multiple .fastq files for a single contigs, they're joined into a comma-separated list, otherwise a single .fastq is passed
        fastq1=lambda wildcards, input: ",".join(input.fastq1_list) if len(input.fastq1_list) > 1 else input.fastq1_list,
        fastq2=lambda wildcards, input: ",".join(input.fastq2_list) if len(input.fastq2_list) > 1 else input.fastq2_list
    conda:
        "../env/assemble.yaml"
    threads:
        config['threads']['megahit']
    benchmark:
        "output/benchmarks/assemble/megahit/{contig_sample}_benchmark.txt"
    log:
        "output/logs/assemble/megahit/{contig_sample}.log"
    resources:
        mem_mb=config['mem_mb']['megahit']
    shell:
        """
        megahit -t {threads} \
                -o {params.temp_dir}/ \
                --memory $(({resources.mem_mb}*1024*1024)) \
                --presets {params.presets} \
                -1 {params.fastq1} \
                -2 {params.fastq2} \
                2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {params.temp_dir}/final.contigs.fa {output.contigs}
        rm -rf {params.temp_dir}
        """

rule quast:
    """
    Does an evaluation of assembly quality with Quast
    """
    input:
        lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                                 assembler=config['assemblers'],
                                 contig_sample=wildcards.contig_sample)
    output:
        report="output/assemble/{assembler}/quast/{contig_sample}/report.txt",
    params:
        outdir=directory("output/assemble/{assembler}/quast/{contig_sample}/")
    threads:
        1
    log:
        "output/logs/assemble/{assembler}/quast/{contig_sample}.log"
    conda:
        "../env/assemble.yaml"
    benchmark:
        "output/benchmarks/assemble/{assembler}/quast/{contig_sample}_benchmark.txt"
    shell:
        """
        quast.py \
          -o {params.outdir} \
          -t {threads} \
          {input}
          touch {output.report}
        """

rule multiqc_assemble:
    input:
        lambda wildcards: expand("output/assemble/{assembler}/quast/{contig_sample}/report.txt",
                                 assembler=config['assemblers'],
                                 contig_sample=contig_groups.keys())
    output:
        "output/assemble/multiqc_assemble/multiqc.html"
    params:
        "--dirs " + config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/assemble/multiqc_assemble/multiqc_assemble.log"
    benchmark:
        "output/benchmarks/assemble/multiqc_assemble/multiqc_assemble_benchmark.txt"
    wrapper:
        "v2.7.0/bio/multiqc"

rule metaquast:
    """
    Does an evaluation of assembly quality with Quast
    """
    input:
        lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                                 assembler=wildcards.assembler,
                                 contig_sample=wildcards.contig_sample)
    output:
        report="output/assemble/{assembler}/metaquast/{contig_sample}/report.html"
    threads:
        config['threads']['metaquast']
    log:
        "output/logs/assemble/{assembler}/metaquast/{contig_sample}.log"
    params:
        outdir=directory("output/assemble/{assembler}/metaquast/{contig_sample}"),
        refs=config['params']['metaquast']['reference_dir'],
        extra=config['params']['metaquast']
    conda:
        "../env/assemble.yaml"
    benchmark:
        "output/benchmarks/assemble/{assembler}/metaquast/{contig_sample}_benchmark.txt"
    shell:
        """
        metaquast.py \
          -r {params.refs} \
          -o {params.outdir} \
          -t {threads} \
          {params.extra} \
          {input}
        """

rule multiqc_metaquast:
    input:
        expand(rules.metaquast.output.report,
               assembler=config['assemblers'],
               contig_sample=contig_groups.keys())
    output:
        "output/assemble/multiqc_metaquast/multiqc.html"
    params:
        "--dirs " + config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/assemble/multiqc_metaquast/multiqc_metaquast.log"
    benchmark:
        "output/benchmarks/assemble/multiqc_metaquast/multiqc_metaquast_benchmark.txt"
    wrapper:
        "0.72.0/bio/multiqc"
