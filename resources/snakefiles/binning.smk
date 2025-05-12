# These work (I think) because sample is the name of the firt variable passed to the
# function. So wildcards.contig_sample (which is contig_pairings.keys()) is passed as "sample"
def get_bam_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/contigs/{contig_pairings}_Mapped_To_{sample}.bam",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

def get_index_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/contigs/{contig_pairings}_Mapped_To_{sample}.bam.bai",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

rule make_metabat2_coverage_table:
    """
    Uses jgi_summarize_bam_contig_depths to generate a depth.txt file.
    """
    input:
        bams = lambda wildcards: get_bam_list(wildcards.contig_sample, wildcards.mapper, contig_pairings)
    output:
        coverage_table="output/binning/metabat2/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}.log"
    shell:
        """
            jgi_summarize_bam_contig_depths --outputDepth {output.coverage_table} {input.bams} 2> {log}
        """

rule run_metabat2:
    """
    Runs Metabat2:
    MetaBAT2 clusters metagenomic contigs into different "bins", each of which should correspond to a putative genome.

    MetaBAT2 uses nucleotide composition information and source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs) to perform binning.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),
        coverage_table = lambda wildcards: expand("output/binning/metabat2/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt",
                mapper=config['mappers'],
                contig_sample=wildcards.contig_sample)
    output:
        bins = directory("output/binning/metabat2/{mapper}/bin_fastas/{contig_sample}/")
    params:
        basename = "output/binning/metabat2/{mapper}/bin_fastas/{contig_sample}/{contig_sample}.metabat2.bin",
        min_contig_length = config['params']['metabat2']['min_contig_length'],
        extra = config['params']['metabat2']['extra']  # optional parameters
    threads:
        config['threads']['run_metabat2']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/metabat2/{mapper}/run_metabat2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/metabat2/{mapper}/run_metabat2/{contig_sample}.log"
    shell:
        """
            metabat2 {params.extra} --numThreads {threads} \
            --inFile {input.contigs} \
            --outFile {params.basename} \
            --abdFile {input.coverage_table} \
            --minContig {params.min_contig_length} \
            2> {log} 1>&2
        """

rule make_maxbin2_abund_list:
    """
       Combines the file paths from 'make_maxbin2_coverage_table' for MaxBin2
    """
    input:
        lambda wildcards: expand("output/mapping/{mapper}/coverage_tables/contigs/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample])
    output:
        abund_list = "output/binning/maxbin2/{mapper}/abundance_lists/{contig_sample}_abund_list.txt"
    benchmark:
        "output/benchmarks/binning/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list_benchmark.txt"
    log:
        "output/logs/binning/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list.log"
    run:
        with open(output.abund_list, 'w') as f:
            for fp in input:
                f.write('%s\n' % fp)


rule run_maxbin2:
    """
    Runs MaxBin2:
    MaxBin2 clusters metagenomic contigs (assembled contiguous genome fragments) into different "bins", each of which corresponds to a putative population genome. It uses nucleotide composition information, source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs), and phylogenetic marker genes to perform binning through an Expectation-Maximization (EM) algorithm.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),
        abund_list = lambda wildcards: expand("output/binning/maxbin2/{mapper}/abundance_lists/{contig_sample}_abund_list.txt",
                mapper=config['mappers'],
                contig_sample=wildcards.contig_sample),

        coverage_list = lambda wildcards: expand("output/mapping/{mapper}/coverage_tables/contigs/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample])
    output:
        bins = directory("output/binning/maxbin2/{mapper}/bin_fastas/{contig_sample}/")
    params:
        basename = "output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/{contig_sample}.maxbin2.bin",
        out_dir = directory("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}"),
        prob = config['params']['maxbin2']['prob_threshold'],  # optional parameters
        min_contig_length = config['params']['maxbin2']['min_contig_length'],
        extra = config['params']['maxbin2']['extra']  # optional parameters
    threads:
        config['threads']['run_maxbin2']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}.log"
    shell:
        """
            mkdir -p {params.out_dir}

            run_MaxBin.pl -thread {threads} -prob_threshold {params.prob} \
            -min_contig_length {params.min_contig_length} {params.extra} \
            -contig {input.contigs} \
            -abund_list {input.abund_list} \
            -out {params.basename} \
            2> {log} 1>&2

            mkdir -p {output.bins}

            mv {params.out_dir}/*.fasta {output.bins}
        """
rule cut_up_fasta:
    """
    Cut up fasta file in non-overlapping or overlapping parts of equal length.
    Optionally creates a BED-file where the cutup contigs are specified in terms
    of the original contigs. This can be used as input to concoct_coverage_table.py.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample)
    output:
        bed="output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.bed",
        contigs_10K="output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.fa"
    conda:
        "../env/concoct_linux.yaml"
    params:
        chunk_size=config['params']['concoct']['chunk_size'],
        overlap_size=config['params']['concoct']['overlap_size']
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/cut_up_fasta/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/cut_up_fasta/{contig_sample}.log"
    shell:
        """
          cut_up_fasta.py {input.contigs} \
          -c {params.chunk_size} \
          -o {params.overlap_size} \
          --merge_last \
          -b {output.bed} > {output.contigs_10K} 2> {log}
        """

rule make_concoct_coverage_table:
    """
    Generates table with per sample coverage depth.
    Assumes the directory "/output/binning/{mapper}/mapped_reads/" contains sorted and indexed bam files where each contig file has has reads mapped against it from the selected prototypes.

    """
    input:
        bed = "output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.bed",
        bam = lambda wildcards: get_bam_list(wildcards.contig_sample, config['mappers'], contig_pairings),
        index = lambda wildcards: get_index_list(wildcards.contig_sample, config['mappers'], contig_pairings)
    output:
        coverage_table = "output/binning/concoct/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}.log"
    shell:
        """
          concoct_coverage_table.py {input.bed} \
          {input.bam} > {output.coverage_table} 2> {log}
        """

rule run_concoct:
    """
    CONCOCT - Clustering cONtigs with COverage and ComposiTion
    CONCOCT does unsupervised binning of metagenomic contigs by using nucleotide composition - kmer frequencies - and coverage data for multiple samples.
    """
    input:
        contigs_10K=rules.cut_up_fasta.output.contigs_10K,
        coverage_table=rules.make_concoct_coverage_table.output.coverage_table
    output:
        clustering = "output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample}_bins_clustering.csv"
    params:
        bins = "output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample,[A-Za-z0-9_]+}_bins",
        min_contig_length=config['params']['concoct']['min_contig_length']
    conda:
        "../env/concoct_linux.yaml"
    threads:
        config['threads']['run_concoct']
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/run_concoct/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/run_concoct/{contig_sample}.log"
    shell:
        """
            concoct --threads {threads} -l {params.min_contig_length} \
            --composition_file {input.contigs_10K} \
            --coverage_file {input.coverage_table} \
            -b {params.bins} \
            2> {log} 1>&2

			# This creates a new concoct clustering.csv file with clusters renamed to reflect their contig_sample name. 
			# Renamed clusters are the basis for concoct bin names in subsequent steps
			awk -F "," '{{OFS=","}} NR==1 {{print}} NR>1{{print $1,"{wildcards.contig_sample}.concoct.bin."$2}}' \
				output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}/{wildcards.contig_sample}_bins_clustering_gt{params.min_contig_length}.csv >\
				output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}/{wildcards.contig_sample}_bins_clustering.csv

			rm output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}/{wildcards.contig_sample}_bins_clustering_gt{params.min_contig_length}.csv
        """

rule merge_cutup_clustering:
    """
    Merges subcontig clustering into original contig clustering.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample}_bins_clustering.csv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        merged = "output/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_clustering_merged.csv"
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}.log"
    shell:
        """
            merge_cutup_clustering.py {input.bins} > {output.merged} 2> {log}
        """

rule extract_fasta_bins:
    """
    Extracts bins as individual FASTA.
    """
    input:
        original_contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                    selected_assembler = selected_assembler,
                    contig_sample = wildcards.contig_sample),
        clustering_merged = rules.merge_cutup_clustering.output.merged
    output:
        fasta_bins = directory("output/binning/concoct/{mapper}/bin_fastas/{contig_sample}/")
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}.log"
    shell:
        """
            mkdir -p {output.fasta_bins}
            extract_fasta_bins.py \
            {input.original_contigs} \
            {input.clustering_merged} \
            --output_path {output.fasta_bins} \
            2> {log}
        """
# TODO: Fix this stupid rule
rule make_vamb_coverage_table:
    input:
        coverages = lambda wildcards: expand("output/mapping/{mapper}/coverage_tables/contigs/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample]),

        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),
    params:
        read_sample = lambda wildcards: contig_pairings[wildcards.contig_sample]
    output:
        abund_list = "output/binning/vamb/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    threads:
        1
    conda:
        "../env/binning.yaml"
    log:
        "output/logs/binning/vamb/{mapper}/make_vamb_coverage_table/{contig_sample}_coverage_table.log"
    script:
        "../scripts/make_vamb_coverage_table.py"

rule run_vamb:
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),

        coverage_table = lambda wildcards: expand("output/binning/vamb/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample),
    params:
        min_contig_length=config['params']['vamb']['min_contig_length'],
        min_bin_length=config['params']['vamb']['min_bin_length'],
        basename = "{contig_sample}.vamb.bin.",
        base_dir = "output/binning/vamb/{mapper}/run_vamb",
        out_dir = "output/binning/vamb/{mapper}/run_vamb/{contig_sample}",
    output:
        bins = directory("output/binning/vamb/{mapper}/bin_fastas/{contig_sample}")
    threads:
        config['threads']['run_vamb']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/vamb/{mapper}/run_vamb/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/vamb/{mapper}/run_vamb/{contig_sample}.log"
    shell:
        """
        mkdir -p {params.base_dir}
        vamb bin default --fasta {input.contigs} --abundance_tsv {input.coverage_table} \
        --outdir {params.out_dir} -m {params.min_contig_length} -p {threads} -o \
        --minfasta {params.min_bin_length} \
        2> {log} 1>&2

        # Vamb bins are renamed to follow the same naming convention as the ouput from the other binners
        mkdir -p {output.bins}
        for bin in {params.out_dir}/bins/*; do mv "${{bin}}" "{output.bins}/{params.basename}${{bin##*/}}"; done
        rmdir {params.out_dir}/bins
        """
rule run_semibin2:
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),
        bams = lambda wildcards: expand("output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample])
    params:
        basename = "{contig_sample}.semibin2.bin.",
        base_dir = "output/binning/semibin2/{mapper}/run_semibin2",
        out_dir = "output/binning/semibin2/{mapper}/run_semibin2/{contig_sample}",
    output:
        bins = directory("output/binning/semibin2/{mapper}/bin_fastas/{contig_sample}")
    threads:
        config['threads']['run_semibin2']
    resources:
        gpu=1,
        cpus_per_gpu=config['threads']['run_semibin2']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/semibin2/{mapper}/run_semibin2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/semibin2/{mapper}/run_semibin2/{contig_sample}.log"
    shell:
        """
        mkdir -p {output.bins}
        SemiBin2 single_easy_bin \
        --tag-output {wildcards.contig_sample} --compression none \
        -t {resources.cpus_per_gpu} \
        -i {input.contigs} \
        -b {input.bams} \
        -o {params.out_dir} \
        --self-supervised \
        2> {log} 1>&2

        for bin in {params.out_dir}/output_bins/*; do mv "${{bin}}" "{output.bins}/{params.basename}${{bin##*{wildcards.contig_sample}_}}"; done
        rmdir {params.out_dir}/output_bins
        """

rule run_comebin:
    input:
        contigs = lambda wildcards: expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                selected_assembler = selected_assembler,
                contig_sample = wildcards.contig_sample),

        bams = lambda wildcards: expand("output/mapping/{mapper}/sorted_bams/contigs/{read_sample}_Mapped_To_{contig_sample}.bam",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample])
    params:
        basename = "{contig_sample}.comebin.bin.",
        base_dir = "output/binning/comebin/{mapper}/run_comebin",
        out_dir = "output/binning/comebin/{mapper}/run_comebin/{contig_sample}",
        tempdir="output/binning/comebin/{mapper}/{contig_sample}_Bins_Temp",
        bamlist=lambda wildcards, input: f'{{{",".join(input.bams)}}}' if len(input.bams) > 1 else input.bams
    output:
        bins = directory("output/binning/comebin/{mapper}/bin_fastas/{contig_sample}")
    threads:
        config['threads']['run_comebin']
    resources:
        gpu=1,
        cpus_per_gpu=config['threads']['run_comebin']
    conda:
        "../env/comebin_binning.yaml"
    benchmark:
        "output/benchmarks/binning/comebin/{mapper}/run_comebin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/comebin/{mapper}/run_comebin/{contig_sample}.log"
    shell:
        """
        mkdir -p {params.tempdir}
        cp {params.bamlist} {params.tempdir}
        mkdir -p {params.out_dir}
        mkdir -p {output.bins}

        CUDA_VISIBLE_DEVICES=0 bash run_comebin.sh -t {resources.cpus_per_gpu} \
        -a {input.contigs} \
        -p {params.tempdir} \
        -o {params.out_dir} \
        2> {log} 1>&2

        for bin in {params.out_dir}/comebin_res/comebin_res_bins/*; do mv "${{bin}}" "{output.bins}/{params.basename}${{bin##*/}}"; done
        rmdir {params.out_dir}/comebin_res/comebin_res_bins
        rm -r {params.tempdir}
        """
rule test_gpu:
    input:
        contigs="output/assemble/megahit/ABX-CJ-IRIS.contigs.fasta",
        bams=["output/mapping/minimap2/sorted_bams/contigs/ABX-CJ-IRIS-14_Mapped_To_ABX-CJ-IRIS.bam", "output/mapping/minimap2/sorted_bams/contigs/ABX-CJ-IRIS-43_Mapped_To_ABX-CJ-IRIS.bam",
              "output/mapping/minimap2/sorted_bams/contigs/ABX-CJ-IRIS-29_Mapped_To_ABX-CJ-IRIS.bam", "output/mapping/minimap2/sorted_bams/contigs/ABX-CJ-IRIS-56_Mapped_To_ABX-CJ-IRIS.bam"]
    params:
        tempdir="iris_bams_temp",
        bamlist=lambda wildcards, input: ','.join(input.bams)
    output:
        "gpu.txt"
    resources:
        gpu=1,
        cpus_per_gpu=1
    conda:
        "../env/comebin_binning.yaml"
    shell:
        """
        python ../scripts/test_cuda_gpu.py
        CUDA_VISIBLE_DEVICES=1 python ../scripts/test_cuda_gpu.py
        CUDA_VISIBLE_DEVICES=0 python ../scripts/test_cuda_gpu.py
        """

