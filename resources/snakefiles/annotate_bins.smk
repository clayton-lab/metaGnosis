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
rule annotate_taxonomy_contigs:
    input: 
        built_cat_db=rules.build_cat_pack_db.output,
        contigs="output/assemble/megahit/{contig_sample}.contigs.fasta"
    params:
        cat_prefix=rules.download_cat_pack_db.params.cat_prefix,
        db_path=rules.download_cat_pack_db.params.db_path,
        outdir=lambda wildcards: f'output/annotate_bins/annotated_contigs/{wildcards.contig_sample}',
        r_param=config['params']['cat_pack']['contigs_r_param'],
        f_param=config['params']['cat_pack']['contigs_f_param'],
        block_size=config['params']['cat_pack']['block_size'],
        index_chunks=config['params']['cat_pack']['index_chunks']
    output:
        pred_prots="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.predicted_proteins.faa",
        diamond_align="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.alignment.diamond",
        orf2lca="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.ORF2LCA.txt",
        c2c="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.contig2classification.txt",
        named_c2c="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.contig2classification.names.txt"
    log:
        "output/logs/annotate_bins/annotate_contigs/{contig_sample}_annotated.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['cat_pack']

    shell:
        """
        {params.cat_prefix}/CAT_pack contigs \
        -c {input.contigs} \
        -d {params.db_path}/db \
        -t {params.db_path}/tax \
        -o {params.outdir}/{wildcards.contig_sample} \
        --range {params.r_param} \
        --fraction {params.f_param} \
        --verbose \
        --force \
        --block_size {params.block_size} \
        --index_chunks {params.index_chunks} \
        -n {threads} \
        2> {log} 1>&2

        {params.cat_prefix}/CAT_pack add_names \
        -i {output.c2c} \
        -t {params.db_path}/tax \
        -o {output.named_c2c} \
        --only_official \
        2> {log} 1>&2

        mv {params.outdir}/{wildcards.contig_sample}.log {log}

        """

rule annotate_taxonomy_bins:
    input: 
        pred_prots=rules.annotate_taxonomy_contigs.output.pred_prots,
        diamond_align=rules.annotate_taxonomy_contigs.output.diamond_align,
        bins="output/refine_bins/filtered_bins/minimap2/{contig_sample}"
    params:
        cat_prefix=rules.download_cat_pack_db.params.cat_prefix,
        db_path=rules.download_cat_pack_db.params.db_path,
        outdir=lambda wildcards: f'output/annotate_bins/annotated_bins/{wildcards.contig_sample}',
        r_param=config['params']['cat_pack']['bins_r_param'],
        f_param=config['params']['cat_pack']['bins_f_param'],
        block_size=config['params']['cat_pack']['block_size'],
        index_chunks=config['params']['cat_pack']['index_chunks']

    output:
        orf2lca="output/annotate_bins/annotated_bins/{contig_sample}/{contig_sample}.ORF2LCA.txt",
        b2c="output/annotate_bins/annotated_bins/{contig_sample}/{contig_sample}.bin2classification.txt",
        named_b2c="output/annotate_bins/annotated_contigs/{contig_sample}/{contig_sample}.bin2classification.names.txt"
    log:
        "output/logs/annotate_bins/annotate_bins/{contig_sample}_annotated.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1

    shell:
        """
        {params.cat_prefix}/CAT_pack bins \
        -b {input.bins} \
        -d {params.db_path}/db \
        -t {params.db_path}/tax \
        -o {params.outdir}/{wildcards.contig_sample} \
        -p {input.pred_prots} \
        -a {input.diamond_align} \
        -s fa \
        --range {params.r_param} \
        --fraction {params.f_param} \
        --verbose \
        --force \
        --block_size {params.block_size} \
        --index_chunks {params.index_chunks} \
        -n {threads} \
        2> {log} 1>&2

        {params.cat_prefix}/CAT_pack add_names \
        -i {output.b2c} \
        -t {params.db_path}/tax \
        -o {output.named_b2c} \
        --only_official \
        2> {log} 1>&2

        mv {params.outdir}/{wildcards.contig_sample}.log {log}
        """

rule annotate_taxonomy_reads:
    input: 
        pred_prots=rules.annotate_taxonomy_contigs.output.pred_prots,
        diamond_align=rules.annotate_taxonomy_contigs.output.diamond_align,
        contigs=rules.annotate_taxonomy_contigs.input.contigs,
        c2c=rules.annotate_taxonomy_contigs.output.c2c,
        bins=rules.annotate_taxonomy_bins.input.bins,
        b2c=rules.annotate_taxonomy_bins.output.b2c,
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
        bam="output/mapping/minimap2/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.bam"
    params:
        cat_prefix=rules.download_cat_pack_db.params.cat_prefix,
        db_path=rules.download_cat_pack_db.params.db_path,
        outdir=lambda wildcards: f'output/annotate_bins/annotated_reads/{wildcards.contig_sample}',
        block_size=config['params']['cat_pack']['block_size'],
        index_chunks=config['params']['cat_pack']['index_chunks']

    output:
        "output/annotate_bins/annotated_reads/{contig_sample}/{read_sample}.unmapped2classification.txt"
    log:
        "output/logs/annotate_bins/annotate_reads/{contig_sample}/{read_sample}_annotated.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['cat_pack']
    shell:
        """
        {params.cat_prefix}/CAT_pack reads \
        --mode m \
        -1 {input.fastq1} \
        -2 {input.fastq2} \
        --bam1 {input.bam} \
        -b {input.bins} \
        -d {params.db_path}/db \
        -t {params.db_path}/tax \
        -o {params.outdir}/{wildcards.read_sample} \
        -s fa \
        --b2c {input.b2c} \
        --verbose \
        --force \
        --block_size {params.block_size} \
        --index_chunks {params.index_chunks} \
        -n {threads}
        mv {params.outdir}/{wildcards.contig_sample}.log {log}
        """
rule annotate_function:
    input: 
        built_microbe_db=rules.build_microbeannotator_db.output,
        pred_prots=rules.predict_genes_prodigal.output.faa
    params:
        db_path=rules.build_microbeannotator_db.params.db_path,
        outdir=lambda wildcards: f'output/annotate_bins/annotated_function/{wildcards.contig_sample}',
        method=config['params']['microbeannotator']['method']
    output:
        "output/annotate_bins/annotated_function/{contig_sample}/metabolic_summary__heatmap.pdf"
    log:
        "output/logs/annotate_bins/annotate_function/{contig_sample}_annotated.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['microbeannotator']

    shell:
        """
        microbeannotator -i {input.pred_prots} -o {params.outdir} -d {params.db_path} -m {params.method} -t {threads} \
        2> {log}
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

# This tool has been AWFUL to use, especially during the database setup. Best solution is to just use an
# existing DRAM database if possible. If really needed, can download a pre-formatted database from here:
# https://github.com/WrightonLabCSU/DRAM/tree/dev
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
        --prodigal_mode single \
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

rule annotate_bin_pathways:
    input:
        #annotations=lambda wildcards: expand("output/annotate_bins/annotate_bin_function/{contig_sample}/annotations.tsv",
        #        contig_sample=contig_pairings.keys())
        annotations=lambda wildcards: expand(rules.annotate_bin_function.output,
                contig_sample=contig_pairings.keys()),
    params:
        annot_dir=rules.annotate_bin_function.params.outdir,
        func_dirs="output/annotate_bins/annotate_bin_function",
        merged_outdir="output/annotate_bins/annotate_bin_pathways/merged_annotations",
        pathway_outdir="output/annotate_bins/annotate_bin_pathways/annotated_pathways",
        db_path=rules.build_dram_db.params.db_path,
    output:
        merged_annotations="output/annotate_bins/annotate_bin_pathways/merged_annotations/annotations.tsv",
        merged_genes="output/annotate_bins/annotate_bin_pathways/merged_annotations/genes.fna",
        merged_bins="output/annotate_bins/annotate_bin_pathways/merged_annotations/scaffolds.fna",
        genome_stats="output/annotate_bins/annotate_bin_pathways/annotated_pathways/genome_stats.tsv",
        metab_summary="output/annotate_bins/annotate_bin_pathways/annotated_pathways/metabolism_summary.xlsx",
        products="output/annotate_bins/annotate_bin_pathways/annotated_pathways/product.tsv"
    log:
        "output/logs/annotate_bins/annotate_bin_pathways/annotate_pathways.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    shell:
        """
        rm -rf {{{params.merged_outdir},{params.pathway_outdir}}}
        DRAM.py merge_annotations -i "{params.func_dirs}/*" -o {params.merged_outdir} \
        2> {log} 1>&2
        DRAM.py distill -i {params.merged_outdir}/annotations.tsv -o {params.pathway_outdir} \
        --distillate_gene_names --config_loc {params.db_path}/dram_configfile.json \
        --rrna_path {params.merged_outdir}/rrnas.tsv --trna_path {params.merged_outdir}/trnas.tsv \
        --custom_distillate {params.db_path}/distillation/methylotrophy_distillate.tsv \
        --custom_distillate {params.db_path}/distillation/CAMPER_distillate.tsv \
        2>> {log} 1>&2
        """

rule dereplicate_genes:
    input:
        rules.annotate_bin_pathways.output.merged_genes
    params:
        prefix="output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes",
        tempdir="output/annotate_bins/temp"
    output:
        derep_fasta = "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_rep_seq.fasta",
        derep_clust = "output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_cluster.tsv"
    log:
        "output/logs/annotate_bins/annotate_bin_pathways/dereplicate_genes.log"
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['dram']
    shell:
        """
        mmseqs easy-cluster {input} {params.prefix} {params.tempdir} --min-seq-id 0.95 \
        --cov-mode 1 -c 0.95 --cluster-mode 2 --threads {threads} \
        2> {log} 1>&2
        """
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
