# Overall pipeline idea:
# classify assembly with Whokaryote
# Use classification to select the gene predictor (Prodigal or MetaGeneMark for prokaryotes, and GeneMark-ES for eukaryotes).

# Params CAT uses to call prodigal on contigs
#            path_to_prodigal,
#            "-i", contigs_fasta,
#            "-a", proteins_fasta,
#            "-o", proteins_gff,
#            "-p", "meta",
#            "-g", "11",
#            "-q",
#            "-f", "gff"
#
# General workflow:
# XXX Download cat stuff
# XXX Build cat stuff
# XXX Run cat
# XXX run bat (reuse cat output)
# run rat with mc (maybe r too?)
# cat_pack add names
# cat_pack summarise (one for contigs, one for bins. Maybe can do both in one step?)

# Could do some of these in same rule
# index_chunks = 2 or 4 works for 64Gb of RAM. So mess with block size and
# index chunks (and tmpdir). Make those params for the user. Just have blanket param for any other stuff.
# Diamond memory requirement is roughly 20*b/c, with b=block_size and c=chunks
#TODO: Add the named classification files as output
#TODO: Work on a way to include alternate gene prediction tools to Prodigal (preferably one that's domain-agnostic
#TODO: Get rid of the CAT/BAT/RAT and microbeannotator stuff
rule dereplicate_bins:
    input:
        filt_paths=rules.bin_filter_summary.output.filt_paths,
        stats_csv=rules.bin_filter_summary.output.stats_csv
    params:
        completion=rules.filter_bins.params.completion_cutoff,
        contam=rules.filter_bins.params.contam_cutoff,
        outdir="output/refine_bins/dereplicated_bins"
    output:
        derep_fasta="output/refine_bins/dereplicated_bins/dereplicated_bins.fa",
        derep_bins=directory("output/refine_bins/dereplicated_bins/dereplicated_genomes"),
        clust_reps="output/refine_bins/dereplicated_bins/data_tables/Wdb.csv",
        clusters="output/refine_bins/dereplicated_bins/data_tables/Cdb.csv"
    threads:
        config['threads']['dereplicate_bins']
    conda:
        "../env/annotate_bins.yaml"
    benchmark:
        "output/benchmarks/refine_bins/dereplicate_bins/dereplicate_bins_benchmark.txt"
    log:
        "output/logs/refine_bins/dereplicate_bins/dereplicate_bins.log"
    shell:
        """
        dRep dereplicate {params.outdir} -g {input.filt_paths} \
        --S_algorithm skani -sa 0.95 -nc 0.3 -p {threads} \
        -comp {params.completion} -con {params.contam} --genomeInfo {input.stats_csv} \
        2> {log} 1>&2

        # Representative bins are concatenated together into a single file
        for f in {params.outdir}/dereplicated_genomes/*.fa; do
            filename="${{f##*/}}"
            awk -v name="${{filename%.fa}}" '/^>/ {{$0=">"name"_"substr($0,2)}} 1' "$f" >> {params.outdir}/dereplicated_bins.fa
        done
        """

# In case the environment needs to be recreated but the database is already downloaded, run
# conda env config vars set GTDBTK_DATA_PATH="path/to/gtdb/db" to set the env variable
rule annotate_taxonomy:
    input: 
        db=rules.download_gtdbtk_db.output,
        bins=rules.dereplicate_bins.output.derep_bins
    params:
        outdir=lambda wildcards: f'output/annotate_bins/annotate_bin_taxonomy'
    output:
        "output/annotate_bins/annotate_bin_taxonomy/gtdbtk.bac120.summary.tsv"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_taxonomy/annotate_taxonomy_benchmark.txt"
    log:
        "output/logs/annotate_bins/annotate_bin_taxonomy/annotate_bin_taxonomy.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['gtdbtk']
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.bins} \
        --out_dir {params.outdir} \
        --mash_db {params.outdir}/gtdbtk.msh \
        -x fa --pplacer_cpus {threads} --cpus {threads} \
        2> {log} 1>&2
        """

# TODO: Add a rule before this that selects the best assembler and concatenates the assemblies together.
# That can use the bash code below to accomplish the same thing.
rule predict_genes_prodigal:
    input:
        contigs = lambda wildcards: f"output/assemble/{selected_assembler}/{wildcards.contig_sample}.contigs.fasta"
    output:
        gff="output/refine_bins/predict_genes/{contig_sample}_predicted_genes.gff",
        faa="output/refine_bins/predict_genes/{contig_sample}_predicted_genes.faa",
        fna="output/refine_bins/predict_genes/{contig_sample}_predicted_genes.fna"
    conda:
        "../env/annotate_bins.yaml"
    benchmark:
        "output/benchmarks/refine_bins/predict_genes/{contig_sample}_predict_genes_benchmark.txt"
    log:
        "output/logs/refine_bins/predict_genes/{contig_sample}_predict_genes.log"
    threads:
        1
    shell:
        """
        prodigal -i {input.contigs} -f gff -p meta -a {output.faa} -d {output.fna} -o {output.gff} \
        2> {log} 1>&2
        """

rule derep_genes:
    input:
        contigs = lambda wildcards: expand("output/refine_bins/predict_genes/{contig_sample}_predicted_genes.fna",
                                            contig_sample = contig_pairings.keys()),
    params:
        prefix="output/refine_bins/dereplicated_genes",
        tempdir="output/refine_bins/temp",
        basedir="output/refine_bins/predict_genes"
    output:
        derep_fasta = "output/refine_bins/dereplicated_genes/dereplicated_genes_rep_seq.fasta",
        derep_clust = "output/refine_bins/dereplicated_genes/dereplicated_genes_cluster.tsv",
        concat_fasta = "output/refine_bins/dereplicated_genes/concatenated_genes.fna",
    benchmark:
        "output/benchmarks/refine_bins/dereplicate_genes/dereplicate_genes_benchmark.txt"
    log:
        "output/logs/refine_bins/dereplicate_genes/dereplicate_genes.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['dereplicate_bins']
    shell:
        """
        mkdir -p {params.prefix}

        # The sample name is added to each of the ORF headers and all ORFs are merged into a single file
        for f in {params.basedir}/*.fna; do
            filename="${{f##*/}}"
            awk -v name="${{filename%_predicted_genes.fna}}" '/^>/ {{$0=">"name"_"substr($0,2)}} 1' "$f" >> {output.concat_fasta}
        done

        mmseqs easy-cluster {output.concat_fasta} {params.prefix}/dereplicated_genes {params.tempdir} --min-seq-id 0.95 \
        --cov-mode 1 -c 0.95 --cluster-mode 2 --threads {threads} \
        2> {log} 1>&2
        rm -r {params.tempdir}
        """

rule split_dereplicated_genes:
    input:
        rules.derep_genes.output.derep_fasta
    params:
        contig_sample = contig_groups
    output:
        directory("output/refine_bins/dereplicated_genes/dereplicated_gene_fastas")
    log:
        "output/logs/refine_bins/dereplicate_genes/split_dereplicated_genes.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    script:
        "../scripts/split_dereplicated_genes.py"

#TODO: When joining dereplicated bin cluster files, remove .fa extension. Also remove file prefix from DRAM genes faa output

# This tool has been AWFUL to use, especially during the database setup. Best solution is to just use an
# existing DRAM database if possible. If really needed, can download a pre-formatted database from here:
# https://github.com/WrightonLabCSU/DRAM/tree/dev
# TODO: Add this fix to make this work when installing from Conda: https://github.com/WrightonLabCSU/DRAM/issues/329
rule annotate_genes:
    input:
        db=rules.build_dram_db.output,
        derep_genes=rules.split_dereplicated_genes.output
    params:
        db_path=rules.build_dram_db.params.db_path,
        outdir='output/annotate_bins/annotate_genes'
    output:
        "output/annotate_bins/annotate_genes/annotations.tsv",
    log:
        "output/logs/annotate_bins/annotate_genes/annotate_genes.log"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_genes/annotate_genes_benchmark.txt"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['dram']
    shell:
        """
        rm -rf {params.outdir}
        DRAM.py annotate_genes -i "{input.derep_genes}/*.faa" -o {params.outdir} \
        --config_loc {params.db_path}/dram_configfile.json \
        --threads {threads} --use_vogdb \
        --custom_db_name methyl --custom_fasta_loc {params.db_path}/methyl/methylotrophy.faa \
        --custom_hmm_name camper --custom_hmm_loc {params.db_path}/camper/hmm/camper.hmm \
        --custom_hmm_name canthyd --custom_hmm_loc {params.db_path}/canthyd/hmm/cant_hyd.hmm \
        --custom_hmm_name fegenie --custom_hmm_loc {params.db_path}/fegenie/fegenie.hmm \
        --custom_hmm_name sulfur --custom_hmm_loc {params.db_path}/sulfur/sulfur.hmm \
        --custom_hmm_cutoffs_loc {params.db_path}/camper/hmm/custom_camper_hmm_scores.tsv \
        --custom_hmm_cutoffs_loc {params.db_path}/canthyd/hmm/CANT_HYD_HMM_scores.tsv \
        --custom_hmm_cutoffs_loc {params.db_path}/fegenie/custom_fegenie_iron_cut_offs.txt \
        2> {log} 1>&2
        """

rule annotate_contigs:
    input:
        annotations = rules.annotate_genes.output,
        gene_clust = "output/refine_bins/dereplicated_genes/dereplicated_genes_cluster.tsv",
        concat_fasta = rules.derep_genes.output.concat_fasta,
        contigs = expand("output/assemble/{selected_assembler}/{contig_sample}.contigs.fasta",
                                            selected_assembler = selected_assembler,
                                            contig_sample = contig_pairings.keys()),
        faas = expand("output/refine_bins/predict_genes/{contig_sample}_predicted_genes.faa",
                                            contig_sample = contig_pairings.keys()),
        fnas = expand("output/refine_bins/predict_genes/{contig_sample}_predicted_genes.fna",
                                            contig_sample = contig_pairings.keys()),
        gffs = expand("output/refine_bins/predict_genes/{contig_sample}_predicted_genes.gff",
                                            contig_sample = contig_pairings.keys()),
    params:
        contig_samples = contig_groups,
        outdir = "output/annotate_bins/annotate_contigs",
        tempdir = "output/annotate_bins/annotate_contigs_temp"
    output:
        annotations="output/annotate_bins/annotate_contigs/annotations.tsv",
        annot_contigs="output/annotate_bins/annotate_contigs/scaffolds.fna",
        contig_faa="output/annotate_bins/annotate_contigs/genes.faa",
        contig_fna="output/annotate_bins/annotate_contigs/genes.fna",
        contig_gff="output/annotate_bins/annotate_contigs/genes.gff",
    log:
        "output/logs/annotate_bins/annotate_contigs/annotate_contigs.log"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_contigs/annotate_contigs_benchmark.txt"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    script:
        "../scripts/annotate_contigs.py"

rule annotate_contig_pathways:
    input:
        annotations=rules.annotate_contigs.output.annotations,
    params:
        outdir="output/annotate_bins/annotate_contig_pathways",
        db_path=rules.build_dram_db.params.db_path,
    output:
        genome_stats="output/annotate_bins/annotate_contig_pathways/genome_stats.tsv",
        metab_summary="output/annotate_bins/annotate_contig_pathways/metabolism_summary.xlsx",
        products="output/annotate_bins/annotate_contig_pathways/product.tsv"
    log:
        "output/logs/annotate_bins/annotate_contig_pathways/annotate_pathways.log"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_contig_pathways/annotate_contig_pathways_benchmark.txt"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    shell:
        """
        rm -rf {params.outdir}
        DRAM.py distill -i {input.annotations} -o {params.outdir} \
        --distillate_gene_names --config_loc {params.db_path}/dram_configfile.json \
        --custom_distillate {params.db_path}/distillation/methylotrophy_distillate.tsv \
        --custom_distillate {params.db_path}/distillation/CAMPER_distillate.tsv \
        2> {log} 1>&2
        """

rule annotate_bins:
    input:
        contigs2bins=expand(rules.run_DAS_Tool.output.contigs2bin,
                                               mapper=selected_mapper,
                                               contig_sample=contig_pairings.keys()),
        annot_contigs=rules.annotate_contigs.output.annot_contigs,
        contig_faa=rules.annotate_contigs.output.contig_faa,
        contig_fna=rules.annotate_contigs.output.contig_fna,
        contig_gff=rules.annotate_contigs.output.contig_gff,
        clust_reps=rules.dereplicate_bins.output.clust_reps,
        clusters=rules.dereplicate_bins.output.clusters,
        annot=rules.annotate_contigs.output.annotations,
        bin_stats=rules.bin_filter_summary.output.stats_tsv,
        taxonomy=rules.annotate_taxonomy.output,
    params:
        outdir="output/annotate_bins/annotate_bins"
    output:
        binned_contigs="output/annotate_bins/annotate_bins/scaffolds.fna",
        bin_faa="output/annotate_bins/annotate_bins/genes.faa",
        bin_fna="output/annotate_bins/annotate_bins/genes.fna",
        bin_gff="output/annotate_bins/annotate_bins/genes.gff",
        annotations="output/annotate_bins/annotate_bins/bin_annotations.tsv",
        trnas="output/annotate_bins/annotate_bins/bin_trnas.tsv",
        rrnas="output/annotate_bins/annotate_bins/bin_rrnas.tsv",
    log:
        "output/logs/annotate_bins/annotate_bins/annotate_bins.log"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_bins/annotate_bins_benchmark.txt"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['dereplicate_bins']
    script:
        "../scripts/annotate_bins.py"

rule annotate_bin_pathways:
    input:
        annotations=rules.annotate_bins.output.annotations,
        trnas=rules.annotate_bins.output.trnas,
        rrnas=rules.annotate_bins.output.rrnas,
    params:
        outdir="output/annotate_bins/annotate_bin_pathways",
        db_path=rules.build_dram_db.params.db_path,
    output:
        genome_stats="output/annotate_bins/annotate_bin_pathways/genome_stats.tsv",
        metab_summary="output/annotate_bins/annotate_bin_pathways/metabolism_summary.xlsx",
        products="output/annotate_bins/annotate_bin_pathways/product.tsv"
    log:
        "output/logs/annotate_bins/annotate_bin_pathways/annotate_pathways.log"
    benchmark:
        "output/benchmarks/annotate_bins/annotate_bin_pathways/annotate_bin_pathways_benchmark.txt"

    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    shell:
        """
        rm -rf {params.outdir}
        DRAM.py distill -i {input.annotations} -o {params.outdir} \
        --distillate_gene_names --config_loc {params.db_path}/dram_configfile.json \
        --rrna_path {input.rrnas} --trna_path {input.trnas} \
        --custom_distillate {params.db_path}/distillation/methylotrophy_distillate.tsv \
        --custom_distillate {params.db_path}/distillation/CAMPER_distillate.tsv \
        2> {log} 1>&2
        """

rule annotate_bin_taxonomy:
    input:
        db=rules.download_gtdbtk_db.output,
        bins=lambda wildcards: f"output/refine_bins/filtered_bins/{selected_mapper}/{wildcards.contig_sample}"
    params:
        outdir=lambda wildcards: f'output/annotate_bins/annotate_bin_taxonomy/{wildcards.contig_sample}'
    output:
        "output/annotate_bins/annotate_bin_taxonomy/{contig_sample}/gtdbtk.bac120.summary.tsv"
    log:
        "output/logs/annotate_bins/annotate_bin_taxonomy/{contig_sample}_annotated.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['gtdbtk']
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.bins} \
        --out_dir {params.outdir} \
        --mash_db {params.outdir}/gtdbtk.msh \
        -x fa --pplacer_cpus {threads} --cpus {threads} \
        2> {log} 1>&2
        """

rule annotate_bin_function:
    input:
        db=rules.build_dram_db.output,
        bins=lambda wildcards: f"output/refine_bins/filtered_bins/{selected_mapper}/{wildcards.contig_sample}",
        bin_quality=lambda wildcards: f"output/refine_bins/run_CheckM2/{selected_mapper}/{wildcards.contig_sample}/quality_report.tsv",
        bin_taxonomy=lambda wildcards: f"output/annotate_bins/annotate_bin_taxonomy/{wildcards.contig_sample}/gtdbtk.bac120.summary.tsv"
    params:
        db_path=rules.build_dram_db.params.db_path,
        outdir=lambda wildcards: expand('output/annotate_bins/annotate_bin_function/{contig_sample}',
                                                contig_sample=wildcards),
    output:
        "output/annotate_bins/annotate_bin_function/{contig_sample}/annotations.tsv",
    log:
        "output/logs/annotate_bins/annotate_bin_function/{contig_sample}.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['dram']
    shell:
        """
        rm -rf {params.outdir}
        DRAM.py annotate -i "{input.bins}/*.fa" -o {params.outdir} \
        --checkm_quality {input.bin_quality} \
        --gtdb_taxonomy {input.bin_taxonomy} \
        --config_loc {params.db_path}/dram_configfile.json \
        --prodigal_mode meta \
        --threads {threads} --use_vogdb \
        --custom_db_name methyl --custom_fasta_loc {params.db_path}/methyl/methylotrophy.faa \
        --custom_hmm_name camper --custom_hmm_loc {params.db_path}/camper/hmm/camper.hmm \
        --custom_hmm_name canthyd --custom_hmm_loc {params.db_path}/canthyd/hmm/cant_hyd.hmm \
        --custom_hmm_name fegenie --custom_hmm_loc {params.db_path}/fegenie/fegenie.hmm \
        --custom_hmm_name sulfur --custom_hmm_loc {params.db_path}/sulfur/sulfur.hmm \
        --custom_hmm_cutoffs_loc {params.db_path}/camper/hmm/custom_camper_hmm_scores.tsv \
        --custom_hmm_cutoffs_loc {params.db_path}/canthyd/hmm/CANT_HYD_HMM_scores.tsv \
        --custom_hmm_cutoffs_loc {params.db_path}/fegenie/custom_fegenie_iron_cut_offs.txt \
        2> {log} 1>&2
        """

#rule annotate_bin_pathways:
#    input:
#        #annotations=lambda wildcards: expand("output/annotate_bins/annotate_bin_function/{contig_sample}/annotations.tsv",
#        #        contig_sample=contig_pairings.keys())
#        annotations=lambda wildcards: expand(rules.annotate_bin_function.output,
#                contig_sample=contig_pairings.keys()),
#    params:
#        annot_dir=rules.annotate_bin_function.params.outdir,
#        func_dirs="output/annotate_bins/annotate_bin_function",
#        merged_outdir="output/annotate_bins/annotate_bin_pathways/merged_annotations",
#        pathway_outdir="output/annotate_bins/annotate_bin_pathways/annotated_pathways",
#        db_path=rules.build_dram_db.params.db_path,
#    output:
#        merged_annotations="output/annotate_bins/annotate_bin_pathways/merged_annotations/annotations.tsv",
#        merged_genes="output/annotate_bins/annotate_bin_pathways/merged_annotations/genes.fna",
#        merged_bins="output/annotate_bins/annotate_bin_pathways/merged_annotations/scaffolds.fna",
#        genome_stats="output/annotate_bins/annotate_bin_pathways/annotated_pathways/genome_stats.tsv",
#        metab_summary="output/annotate_bins/annotate_bin_pathways/annotated_pathways/metabolism_summary.xlsx",
#        products="output/annotate_bins/annotate_bin_pathways/annotated_pathways/product.tsv"
#    log:
#        "output/logs/annotate_bins/annotate_bin_pathways/annotate_pathways.log"
#    conda:
#        "../env/annotate_bins.yaml"
#    threads:
#        1
#    shell:
#        """
#        rm -rf {{{params.merged_outdir},{params.pathway_outdir}}}
#        DRAM.py merge_annotations -i "{params.func_dirs}/*" -o {params.merged_outdir} \
#        2> {log} 1>&2
#        DRAM.py distill -i {params.merged_outdir}/annotations.tsv -o {params.pathway_outdir} \
#        --distillate_gene_names --config_loc {params.db_path}/dram_configfile.json \
#        --rrna_path {params.merged_outdir}/rrnas.tsv --trna_path {params.merged_outdir}/trnas.tsv \
#        --custom_distillate {params.db_path}/distillation/methylotrophy_distillate.tsv \
#        --custom_distillate {params.db_path}/distillation/CAMPER_distillate.tsv \
#        2>> {log} 1>&2
#        """

#rule dereplicate_genes:
#    input:
#        rules.annotate_bin_pathways.output.merged_genes
#    params:
#        prefix="output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes",
#        tempdir="output/annotate_bins/temp"
#    output:
#        derep_fasta = "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_rep_seq.fasta",
#        derep_clust = "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_cluster.tsv"
#    log:
#        "output/logs/annotate_bins/annotate_bin_pathways/dereplicate_genes.log"
#    conda:
#        "../env/annotate_bins.yaml"
#    threads:
#        config['threads']['dram']
#    shell:
#        """
#        mmseqs easy-cluster {input} {params.prefix} {params.tempdir} --min-seq-id 0.95 \
#        --cov-mode 1 -c 0.95 --cluster-mode 2 --threads {threads} \
#        2> {log} 1>&2
#        rm -r {params.tempdir}
#        """

#TODO: Make it so the pipeline doesn't fail if rrnas.tsv or trnas.tsv don't exist
#rule annotate_bin_pathways:
#    input:
#        annotations=rules.annotate_bin_function.output.annotations,
#    params:
#        annot_dir=rules.annotate_bin_function.params.outdir,
#        outdir=lambda wildcards: f'output/annotate_bins/annotate_bin_pathways/{wildcards.contig_sample}',
#        db_path=rules.build_dram_db.params.db_path,
#    output:
#        genome_stats="output/annotate_bins/annotate_bin_pathways/{contig_sample}/genome_stats.tsv",
#        metab_summary="output/annotate_bins/annotate_bin_pathways/{contig_sample}/metabolism_summary.xlsx",
#        products="output/annotate_bins/annotate_bin_pathways/{contig_sample}/product.tsv"
#    log:
#        "output/logs/annotate_bins/annotate_bin_pathways/{contig_sample}.log"
#    conda:
#        "../env/annotate_bins.yaml"
#    threads:
#        config['threads']['dram']
#    shell:
#        """
#        rm -rf {params.outdir}
#        DRAM.py distill -i {input.annotations} -o {params.outdir} --distillate_gene_names \
#        --config_loc {params.db_path}/dram_configfile.json \
#        --rrna_path {params.annot_dir}/rrnas.tsv --trna_path {params.annot_dir}/trnas.tsv \
#        --custom_distillate {params.db_path}/distillation/methylotrophy_distillate.tsv \
#        --custom_distillate {params.db_path}/distillation/CAMPER_distillate.tsv \
#        2> {log} 1>&2
#        """

#microbeannotator [-h] [-i INPUT_LIST [INPUT_LIST ...]] [-l FILE_LIST]
#                        -o OUTDIR -d DATABASE -m METHOD
#                        [--method_bin METHOD_BIN] [--id_perc ID_PERC]
#                        [--bitscore BITSCORE] [--evalue EVALUE]
#                        [--aln_percent ALN_PERCENT] [--cluster CLUSTER]
#                        [--filename PLOT_FILENAME] [-t THREADS] [-p PROCESSES]
#                        [--light] [--eukaryote] [--full] [--continue_run]
#                        [--refine] [--version]

#
#    shell:
#        """
#        {params.cat_prefix}/CAT_pack reads \
#        --mode mcr \
#        -1 {input.fastq1} \
#        -2 {input.fastq2} \
#        --bam1 {input.bam} \
#        -c {input.contigs} \
#        -b {input.bins} \
#        -d {params.db_path}/db \
#        -t {params.db_path}/tax \
#        -o {params.outdir}/{wildcards.read_sample} \
#        -s fa \
#        --c2c {input.c2c} \
#        --b2c {input.b2c} \
#        -p {input.pred_prots} \
#        -a {input.diamond_align} \
#        --verbose \
#        --force \
#        --block_size {params.block_size} \
#        --index_chunks {params.index_chunks} \
#        -n {threads}
#        mv {params.outdir}/{wildcards.contig_sample}.log {log}
#        """
#

#TODO: Important: Make sure to only annotate high-quality MAGS with RAT. Even the authors recommend that.
#
#  -c , --contigs_fasta
#                        Path to contigs fasta file.
#  -t , --taxonomy_folder
#                        Path to directory that contains taxonomy files.
#  -m , --mode           classification mode. "mcr": integrate annotations from MAGs, contigs, and reads; "cr": integrate annotations from contigs and reads; "mr": integrate annotations from MAGs and reads.
#Optional arguments:
#  -o , --out_prefix     Prefix for output files (default: ./out.RAT).
#  -1 , --read_file1     Path to (forward) read file. Please note that RAT does not currently support interlaced read files. Please supply a single read file or two files for paired-end reads.
#  -2 , --read_file2     Path to reverse read file.
#  --bam1                Path to sorted mapping file.
#  --bam2                Path to second sorted mapping file (not recommended).
#  --alignment_unmapped
#                        Path to alignment file of reads and contigs that couldnot be classified by CAT/BAT.
#  -b , --bin_fasta , --bin_folder
#                        Path to bin fasta file or to directory containing bins.
#  -s , --bin_suffix     Suffix of bins in bin directory (default: .fna).
#  --c2c                 Path to contig2classification file.
#  --b2c                 Path to bin2classification file.
#  --read2classification
#                        Includes read classification step.
#  --u2c                 Path to bin2classification file.
#  --mapping_quality     Minimum mapping quality phred score (default: 2)
#  --path_to_bwa         Path to bwa binaries. Supply if RAT cannot find bwa.
#  --path_to_samtools    Path to samtools binaries. Supply if RAT cannot find samtools.
#  --force               Force overwrite existing files.
#  -q, --quiet           Suppress verbosity.
#  --verbose             Increase verbosity.
#  --no_log              Suppress log file.
#  -h, --help            Show this help message and exit.
#  -p , --proteins_fasta
#                        Path to predicted proteins fasta file. If supplied, the protein prediction step is skipped.
#  -a , --diamond_alignment
#                        Path to alignment table. If supplied, the alignment step is skipped and classification is carried out directly. A predicted proteins fasta file should also be supplied with argument
#                        [-p / --proteins].
#
#CAT/BAT-specific arguments:
#  -d , --database_folder
#                        Path to directory that contains database files.
#  -r , --range          r parameter [0-100] (default: 10).
#  -f , --fraction       f parameter [0-0.99] (default: 0.50).
#  --path_to_prodigal    Path to Prodigal binaries. Supply if CAT/BAT/RAT cannot find Prodigal
#  --path_to_diamond     Path to DIAMOND binaries. Supply if CAT/BAT/RAT cannot find DIAMOND.
#  --no_stars            Suppress marking of suggestive taxonomic assignments.
#  --I_know_what_Im_doing
#                        Flag for experimental features.
#
#DIAMOND specific optional arguments:
#  -n , --nproc          Number of cores to deploy by DIAMOND (default: maximum).
#  --sensitive           Run DIAMOND in sensitive mode (default: not enabled).
#  --no_self_hits        Do not report identical self hits by DIAMOND (default: not enabled).
#  --block_size          DIAMOND block-size parameter (default: 12.0). Lower numbers will decrease memory and temporary disk space usage.
#  --index_chunks        DIAMOND index-chunks parameter (default: 1). Set to 4 on low memory machines. The parameter has no effect on temporary disk space usage.
#  --tmpdir              Directory for temporary DIAMOND files (default: directory to which output files are written).
#  --compress            Compress DIAMOND alignment file (default: not enabled).
#  --top                 DIAMOND top parameter [0-100] (default: 11). Governs hits within range of best hit that are written to the alignment file. This is not the [-r / --range] parameter! See README.md
#
# For annotate reads, will need sequence files, BAM files, contigs, and bins. Will have to do something
# like in the read mapping step, where the read_sample and contig_sample are concurrent wildcards in the
# same rule.

rule summarize_taxonomy:
    input:
        lambda wildcards: expand(rules.annotate_taxonomy_reads.output,
                    contig_sample=wildcards.contig_sample,
                    read_sample=contig_pairings[wildcards.contig_sample])
    output:
        "output/annotate_bins/{contig_sample}_taxonomy.txt"
    shell:
        """
        touch {output}
        """

#rule compile_classified_bins:
#    input:
#        rules.compile_quantified_bins.output,

        #expand(rules.summarize_taxonomy.output,
        #        contig_sample=contig_pairings.keys()),
        #expand(rules.annotate_taxonomy_bins.output,
        #        contig_sample=contig_pairings.keys()),

        #expand(rules.annotate_function.output,
        #        contig_sample=contig_pairings.keys())
        #rules.profile_bins.output,
        #lambda wildcards: get_pairing_list(sample, wildcards.mapper, contig_pairings)
#        coverage_tables=lambda wildcards: expand("output/mapping/{selected_mapper}/coverage_tables/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
            #selected_mapper=selected_mapper,
            #read_sample=wildcards.read_sample,
            #contig_sample=wildcards.contig_sample)
        #contig_lengths=expand("output/mapping/{mapper}/contig_lengths/{read_sample}_Mapped_To_{contig_sample}_contig_lengths.txt",
#    output:
#        "output/classify_bins/classified_bins.txt"
#    shell:
#        """
#            touch {output}
#        """
# Modified diamond search format by microbeannotator shouldn't affect cat diamond parsing.
# Download nr (CAT) > build diamond db (CAT) > diamond search (modified output format from microbeannotator)
