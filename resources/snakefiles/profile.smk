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
        report = "output/profile/kraken2/{sample}.report.txt",
        outfile="output/profile/kraken2/{sample}.output.txt",
    params:
        db_path=rules.kraken2_build_db.params.db_path,
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['kraken2']
    log:
        "output/logs/profile/kraken2/taxonomy_kraken2/{sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/taxonomy_kraken2/{sample}_benchmark.txt"
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
            --only-classified-output \
            2> {log}
       """

rule bracken_abundance:
    input:
        rules.bracken_build.output,
        report=rules.taxonomy_kraken2.output.report,
    output:
        "output/profile/bracken/{sample}.bracken.txt",
    params:
        levels = config['params']['kraken2']['levels'],
        db_path = rules.kraken2_build_db.params.db_path
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/bracken/bracken_abundance/{sample}.log"
    benchmark:
        "output/benchmarks/profile/bracken/bracken_abundance/{sample}_benchmark.txt"
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
        "output/profile/krona/{sample}.report.html"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/krona/{sample}.log"
    benchmark:
        "output/benchmarks/profile/krona/{sample}_benchmark.txt"
    shell:
        """
        perl resources/scripts/kraken2-translate.pl {input} > {input}.temp
        ktImportText -o {output} {input}.temp
        rm {input}.temp
        """
## For final rule of combining kraken (and maybe metaphlan) into single output file:
## make_ktaxonomy > make_kreport.py > kreport2mpa.py > combine_mpa.py
rule make_kreport:
    input:
        kraken_output=rules.taxonomy_kraken2.output.outfile,
        ktaxonomy=rules.make_ktaxonomy.output
    output:
        "output/profile/kraken2/{sample}.kreport.txt"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/kraken2/make_kreport/{sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/make_kreport/{sample}_benchmark.txt"
    shell:
        """
        make_kreport.py \
            --input {input.kraken_output} \
            --taxonomy {input.ktaxonomy} \
            --output {output} \
            2> {log} 1>&2
        """

rule kreport2mpa:
    input:
        rules.make_kreport.output
    output:
        "output/profile/kraken2/{sample}.kreport2mpa.txt"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/kraken2/make_kreport/{sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/make_kreport/{sample}_benchmark.txt"
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
               sample=samples),
        bracken=expand(rules.bracken_abundance.output,
                sample=samples),
        krona=expand(rules.krona.output,
                        sample=samples)

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
## Okay. This is probably the main trigger for the rest of the rules in the snakefile
#rule kraken:
#    input:
#        "output/profile/kraken2/merged_kreport2mpa_table.txt",
#        expand("output/profile/krona/{sample}.report.html",
#                sample=samples),
#        expand("output/profile/bracken/{sample}.bracken.txt",
#                sample=samples)
#
rule metaphlan:
    """

    Performs taxonomic profiling using MetaPhlAn3.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
        db_path=rules.download_metaphlan_db.output
    output:
        bt2="output/profile/metaphlan/bowtie2s/{sample}.bowtie2.bz2",
        sam="output/profile/metaphlan/sams/{sample}.sam.bz2",
        profile="output/profile/metaphlan/profiles/{sample}.txt"
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['metaphlan']
    params:
        other=config['params']['metaphlan']['other']
    benchmark:
        "output/benchmarks/profile/metaphlan/{sample}_benchmark.txt"
    log:
        "output/logs/profile/metaphlan/{sample}.log"
    shell:
        """
        metaphlan {input.fastq1},{input.fastq2} \
        --input_type fastq \
        --nproc {threads} {params.other} \
        --bowtie2db {input.db_path} \
        --bowtie2out {output.bt2}  \
        -s {output.sam}  \
        -o {output.profile} \
        2> {log} 1>&2
        """

rule merge_metaphlan_tables:
    """

    Merges MetaPhlAn3 profiles into a single table.

    """
    input:
        expand(rules.metaphlan.output.profile,
               sample=samples)
    output:
        merged_abundance_table="output/profile/metaphlan/merged_abundance_table.txt"
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table.log"
    benchmark:
        "output/benchmarks/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table_benchmark.txt"
    shell:
        """
        merge_metaphlan_tables.py {input} \
        -o {output.merged_abundance_table} \
        2> {log} 1>&2
        """
