rule make_ktaxonomy:
    input:
        seqid2taxid=rules.kraken2_build_db.output.seqid2taxid,
        taxnames=rules.kraken2_download_taxonomy.output.taxnames,
        taxnodes=rules.kraken2_download_taxonomy.output.taxnodes
    output:
        "output/profile/kraken2/ktaxonomy.txt"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/kraken2/make_ktaxonomy.log"
    benchmark:
        "output/benchmarks/profile/kraken2/make_ktaxonomy_benchmark.txt"
    shell:
        """
        make_ktaxonomy.py \
            --nodes {input.taxnodes} \
            --names {input.taxnames} \
            --seqid2taxid {input.seqid2taxid} \
            --output {output} \
            2> {log} 1>&2
        """

# If using --classified-out or --unclassified-out, the file name needs a #, which will be replaced by _1/_2 by kraken.
# --memory-mapping can be used to prevent kraken from loading the entire database into RAM (which usually crashes the
# program unless only a single process with 120gb memory is used). However, it runs faster without memory-mapping, so
# there's a tradeoff between speed and bugginess.
rule taxonomy_kraken2:
    """
    Runs Kraken2 to construct taxonomic profiles for sequence reads.
    """
    input:
        rules.kraken2_build_db.output,
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
    output:
        report = "output/profile/kraken2/{read_sample}/{read_sample}.report.txt",
        cread1="output/profile/kraken2/{read_sample}/{read_sample}.classified_1.fastq",
        cread2="output/profile/kraken2/{read_sample}/{read_sample}.classified_2.fastq",
        uread1="output/profile/kraken2/{read_sample}/{read_sample}.unclassified_1.fastq",
        uread2="output/profile/kraken2/{read_sample}/{read_sample}.unclassified_2.fastq",
        outfile="output/profile/kraken2/{read_sample}/{read_sample}.output.txt",
    params:
        db_path=rules.kraken2_build_db.params.db_path,
        cread=lambda wildcards, output: output.cread1.replace('_1', '#'),
        uread=lambda wildcards, output: output.uread1.replace('_1', '#')
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['kraken2']
    log:
        "output/logs/profile/kraken2/taxonomy_kraken2/{read_sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/taxonomy_kraken2/{read_sample}_benchmark.txt"
    shell:
        """
        # run Kraken to align reads against reference genomes
        kraken2 {input.fastq1} {input.fastq2} \
            --db {params.db_path} \
            --paired \
            --gzip-compressed \
            --threads {threads} \
            --report {output.report} \
            --output {output.outfile} \
            --classified-out {params.cread} \
            --unclassified-out {params.uread} \
            2> {log}
       """

rule bracken_abundance:
    input:
        rules.bracken_build.output,
        report=rules.taxonomy_kraken2.output.report,
    output:
        "output/profile/bracken/{read_sample}.bracken.txt",
    params:
        levels = config['params']['kraken2']['levels'],
        db_path = rules.kraken2_build_db.params.db_path
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/bracken/bracken_abundance/{read_sample}.log"
    benchmark:
        "output/benchmarks/profile/bracken/bracken_abundance/{read_sample}_benchmark.txt"
    shell:
        """
        bracken \
            -d {params.db_path} \
            -i {input.report} \
            -t 10 \
            -l 'S' \
            -o {output} \
            2>> {log} 1>&2
        """
# This is part of the bracken script that is untested
## get stem file path
#        #stem={wildcards.sample}
#
#        # run Bracken to re-estimate abundance at given rank
#        #if [[ ! -z {params.levels} ]]
#        #then
#            #IFS=',' read -r -a levels <<< "{params.levels}"
#            #for level in "${{levels[@]}}"
#            #do
#                #bracken \
#                        #-d {params.db_path} \
#                        #-i {input.report} \
#                        #-t 10 \
#                        #-l $(echo $level | head -c 1 | tr a-z A-Z) \
#                        #-o $stem.redist.$level.txt \
#                        #2>> {log} 1>&2
#            #done
#        #fi
#        #mv ${{stem}} {output}
#        """
#
## For krona rule:  make_ktaxonomy.py > make_kreport.py > kreport2krona.py
## The current rule seems to use a custom script to convert kraken2 output to krona-compatible input
rule krona:
    input:
        rules.taxonomy_kraken2.output.report
    output:
        "output/profile/krona/{read_sample}.report.html"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/krona/{read_sample}.log"
    benchmark:
        "output/benchmarks/profile/krona/{read_sample}_benchmark.txt"
    shell:
        """
        kreport2krona.py -r {input} -o {input}.krona.temp \
        2> {log} 1>&2
        ktImportText {input}.krona.temp -o {output} 
        2>> {log} 1>&2
        rm {input}.krona.temp \
        """
## For final rule of combining kraken (and maybe metaphlan) into single output file:
## make_ktaxonomy > make_kreport.py > kreport2mpa.py > combine_mpa.py
#rule make_kreport:
#    input:
#        kraken_output=rules.taxonomy_kraken2.output.outfile,
#        ktaxonomy=rules.make_ktaxonomy.output
#    output:
#        "output/profile/kraken2/{sample}.kreport.txt"
#    conda:
#        "../env/profile.yaml"
#    threads:
#        1
#    log:
#        "output/logs/profile/kraken2/make_kreport/{sample}.log"
#    benchmark:
#        "output/benchmarks/profile/kraken2/make_kreport/{sample}_benchmark.txt"
#    shell:
#        """
#        make_kreport.py \
#            --input {input.kraken_output} \
#            --taxonomy {input.ktaxonomy} \
#            --output {output} \
#            2> {log} 1>&2
#        """

rule kreport2mpa:
    input:
        rules.taxonomy_kraken2.output.report
    output:
        temp("output/profile/kraken2/{read_sample}/{read_sample}.kreport2mpa.txt")
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/kraken2/make_kreport/{read_sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/make_kreport/{read_sample}_benchmark.txt"
    shell:
        """
        kreport2mpa.py \
            --report {input} \
            --output {output} \
            --display-header \
            2> {log} 1>&2
        """

rule combine_kreport2mpa_tables:
    input:
        kraken=expand(rules.kreport2mpa.output,
               read_sample=read_groups),
        bracken=expand(rules.bracken_abundance.output,
                read_sample=read_groups),
        krona=expand(rules.krona.output,
                        read_sample=read_groups)

    output:
        "output/profile/kraken2/merged_kreport2mpa_table.txt"
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/profile/kraken2/merge_kreport2mpa_tables/merged_kreport2mpa_table.log"
    benchmark:
        "output/benchmarks/profile/metaphlan/merge_kreport2mpa_tables/merged_kreport2mpa_table_benchmark.txt"
    shell:
        """
        combine_mpa.py \
        --input {input.kraken} \
        --output {output} \
        2> {log} 1>&2
        """

rule metaphlan:
    """

    Performs taxonomic profiling using MetaPhlAn3.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
        downloaded_db=rules.download_metaphlan_db.output
    params:
        db_path = rules.download_metaphlan_db.params.db_path,
        db_index = config['params']['metaphlan']['db_index'] if config['params']['metaphlan']['db_index'] != 'latest' else '',
        extra=config['params']['metaphlan']['extra']
    output:
        bt2="output/profile/metaphlan/bowtie2s/{read_sample}.bowtie2.bz2",
        sam="output/profile/metaphlan/sams/{read_sample}.sam.bz2",
        profile="output/profile/metaphlan/profiles/{read_sample}.tsv"
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['metaphlan']
    benchmark:
        "output/benchmarks/profile/metaphlan/{read_sample}_benchmark.txt"
    log:
        "output/logs/profile/metaphlan/{read_sample}.log"
    shell:
        """
        metaphlan {input.fastq1},{input.fastq2} \
        --input_type fastq \
        --sample_id {wildcards.read_sample} \
        --nproc {threads} {params.extra} \
        --bowtie2db {params.db_path} \
        --index {params.db_index} \
        --bowtie2out {output.bt2} \
        -s {output.sam} \
        -o {output.profile} \
        --offline \
        2> {log} 1>&2
        """

rule merge_metaphlan_tables:
    """

    Merges MetaPhlAn3 profiles into a single table.

    """
    input:
        expand(rules.metaphlan.output.profile,
               read_sample=read_groups)
    output:
        "output/profile/metaphlan/merged_abundance_table.tsv"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table.log"
    benchmark:
        "output/benchmarks/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table_benchmark.txt"
    shell:
        """
        merge_metaphlan_tables.py {input} \
        -o {output} \
        2> {log} 1>&2
        """
rule humann:
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
        profile=rules.metaphlan.output.profile,
    params:
        outdir = "output/profile/humann/sample_tables/{read_sample}",
        temp_fastq = "output/profile/humann/{read_sample}.fastq.gz",
        prot_db_path = rules.download_humann_db.params.prot_db_path,
        nuc_db_path=rules.download_humann_db.params.nuc_db_path,
        extra=config['params']['humann']['extra']
    output:
        gene_fams = "output/profile/humann/sample_tables/{read_sample}/{read_sample}_genefamilies.tsv",
        path_abund = "output/profile/humann/sample_tables/{read_sample}/{read_sample}_pathabundance.tsv",
        path_cov = "output/profile/humann/sample_tables/{read_sample}/{read_sample}_pathcoverage.tsv"
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['humann']
    log:
        "output/logs/profile/humann/humann/{read_sample}_humann.log"
    benchmark:
        "output/benchmarks/profile/humann/humann/{read_sample}_humann.txt"
    shell:
        """
        cat {input.fastq1} {input.fastq2} > {params.temp_fastq}
        humann -i {params.temp_fastq} \
               -o {params.outdir} \
               --nucleotide-database {params.nuc_db_path} \
               --protein-database {params.prot_db_path} \
               --taxonomic-profile {input.profile} \
               --memory-use maximum \
               --remove-temp-output \
               --threads {threads} {params.extra} \
               --o-log {log} \
               2>> {log} 1>&2
        rm {params.temp_fastq}
        """
rule join_humann_tables:
    input:
        gene_fams = expand(rules.humann.output.gene_fams, 
                           read_sample=read_groups),
        path_abund = expand(rules.humann.output.path_abund, 
                           read_sample=read_groups),
        path_cov = expand(rules.humann.output.path_cov, 
                           read_sample=read_groups),
    params:
        searchdir = "output/profile/humann/sample_tables",
    output:
        gene_fams_merged = "output/profile/humann/merged_tables/merged_genefamilies.tsv",
        path_abund_merged = "output/profile/humann/merged_tables/merged_pathabundance.tsv",
        path_abund_cpm = "output/profile/humann/merged_tables/normed_pathabundance.tsv",
        path_cov_merged = "output/profile/humann/merged_tables/merged_pathcoverage.tsv"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/humann/join_humann_tables/join_humann_tables.log"
    benchmark:
        "output/benchmarks/profile/humann/join_humann_tables/join_humann_tables.txt"
    shell:
        """
        humann_join_tables -i {params.searchdir} \
        -o {output.gene_fams_merged} -s \
        --file_name genefamilies.tsv \
        2> {log} 1>&2

        humann_join_tables -i {params.searchdir} \
        -o {output.path_abund_merged} \
        -s --file_name pathabundance.tsv \
        2>> {log} 1>&2

        humann_renorm_table \
        -i {output.path_abund_merged} \
        -o {output.path_abund_cpm} \
        2>> {log} 1>&2

        humann_join_tables -i {params.searchdir} \
        -o {output.path_cov_merged} -s \
        --file_name pathcoverage.tsv \
        2>> {log} 1>&2
        """
