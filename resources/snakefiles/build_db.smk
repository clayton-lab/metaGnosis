## To make Kraken2 work, it needs at least 100 Gb of memory to build the database.
# TODO: Include an option to remove the extra files (reference_sequences and taxonomy mapping info) for users to reduce space if they want.
# Make it so the temporary ncbi_taxonomy_db is removed, otherwise it probably won't be
rule get_taxid_list:
    output:
        temp("output/build_db/build_genome_db/taxid_list.txt")
    params:
        taxids=config['add_genomes']['genome_taxids'],
        temp_taxdir=temp(directory("output/build_db/build_genome_db/ncbi_taxonomy_db"))
    conda:
        "../env/qc.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_genome_db/get_taxid_list.log"
    shell:
        """
        gimme_taxa.py -v --just-taxids \
        -d {params.temp_taxdir} \
        -o {output} \
        {params.taxids} \
        2> {log} 1>&2
        """

# TODO: Add option to download protein sequences instead. Will also need to edit kraken2_addto_libraries
# to enable the addition of protein sequences
rule build_genome_db:
    input:
        taxid_list=rules.get_taxid_list.output if rules.get_taxid_list.params.taxids else [],
    output:
        join(config['user_paths']['genome_db_path'], "build_genome_db_metadata.tsv"),
    params:
        db_path=config['user_paths']['genome_db_path'],
        taxids=config['add_genomes']['genome_taxids'],
        # If genome_accessions are specified and host_filter_accns isn't part of it, the host_filter_accns is automatically
        # added. This ensures that the specified host_filter genome will be fetched for filtering in qc.smk
        accns=config['add_genomes']['genome_accessions'] + \
              config['host_filter_accn'] if (config['host_filter_accn'] and \
                                            config['host_filter_accn'] not in config['add_genomes']['genome_accessions']) \
                                         else config['add_genomes']['genome_accessions']
    conda:
        "../env/qc.yaml"
    threads:
        config['threads']['build_db']
    log:
        "output/logs/build_db/build_genome_db/build_genome_db.log"
    shell:
        """
        main() {{
            # -p for parallel downloads
            metadata=()
            pwd
            echo "{input}"
            echo "Checking for specified genome files" 2> {log} 1>&2
            if [[ ! -z {params.accns} ]]; then
                echo "Starting with accession numbers" 2>> {log} 1>&2
                if [[ -s {output} ]]; then
                    touch {output}
                    echo "Genomes for specified accession numbers already present at {params.db_path}" \
                    2>> {log} 1>&2
                else
                    echo "Downloading specified accession numbers" 2>> {log} 1>&2
                    ncbi-genome-download -s 'genbank' \
                    -N -F 'fasta' -R 'reference,representative' 'all' -r 3 --flat-output \
                    -p {threads} \
                    -A {params.accns} \
                    -o {params.db_path} \
                    -m {params.db_path}/accn_genome_metadata.tsv \
                    2>> {log} 1>&2
                fi
                metadata+=({params.db_path}/accn_genome_metadata.tsv)
                echo "${{metadata}}"
            fi
            if [[ ! -z {params.taxids} ]]; then
                echo "Moving on to taxids" 2>> {log} 1>&2
                if [[ -s {output} ]]; then
                    touch {output}
                    echo "Genomes for specified taxids already present at {params.db_path}" \
                    2>> {log} 1>&2
                else
                    echo "Downloading specified taxids" 2>> {log} 1>&2
                    ncbi-genome-download -s 'genbank' \
                    -N -F 'fasta' -R 'reference,representative' 'all' -r 3 --flat-output \
                    -p {threads} \
                    -t {input.taxid_list} \
                    -o {params.db_path} \
                    -m {params.db_path}/taxid_genome_metadata.tsv \
                    2>> {log} 1>&2
                fi
                metadata+=({params.db_path}/taxid_genome_metadata.tsv)
            fi
            # Reads the header to the output file for all specified accessions/taxids,
            # then combines all metadata into a single file for easier downstream processing.
            # Accession/taxid metadata files are subsequently removed
            head -n 1 "${{metadata[0]}}" > {output} 
            for metafile in "${{metadata[@]}}"
            do
                cat "${{metafile}}"
                tail -n +2 "${{metafile}}" >> {output}
                rm "${{metafile}}"
            done
        }}
        main 2> {log} 1>&2

        """

rule kraken2_reformat_genome_fastas:
    input:
        rules.build_genome_db.output,
    output:
        temp("output/build_db/build_kraken2_db/kraken2_reformat_genome_fastas/kraken2_reformat_genome_fastas.done")
    params:
        temp(directory("output/build_db/build_kraken2_db/kraken2_reformat_genome_fastas"))
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_reformat_genome_fastas.log"
    benchmark:
        "output/benchmarks/build_db/build_kraken2_db/kraken2_reformat_genome_fastas.txt"
    shell:
        """
        python resources/scripts/reformat_genome_fastas_kraken2.py \
        --genome_metadata {input} \
        --output_path {params} \
        2> {log} 1>&2
        touch {output}
        """
 
rule kraken2_download_taxonomy:
    output:
        taxnames=join(config['user_paths']['kraken2_db_path'], "taxonomy/names.dmp"),
        taxnodes=join(config['user_paths']['kraken2_db_path'], "taxonomy/nodes.dmp"),
    params:
        db_path=directory(config['user_paths']['kraken2_db_path']),
        db_type='' if config['params']['kraken2']['db_type'] == 'nucleotide' else '--protein'
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_download_taxonomy.log"
    benchmark:
        "output/benchmarks/build_db/build_kraken2_db/kraken2_download_taxonomy_benchmark.txt"
    shell:
        """
        kraken2-build --download-taxonomy {params.db_type} \
        --threads {threads} \
        --db {params.db_path} \
        2> {log} 1>&2
        """

rule kraken2_addto_libraries:
    input:
        rules.kraken2_download_taxonomy.output,
        rules.kraken2_reformat_genome_fastas.output
    output:
        temp("output/build_db/build_kraken2_db/kraken2_addto_libraries/kraken2_addto_libraries.done")
    params:
        fasta_dir=rules.kraken2_reformat_genome_fastas.params,
        db_path=rules.kraken2_download_taxonomy.params.db_path
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['build_db'] * 0.5
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_addto_libraries.log"
    benchmark:
        "output/benchmarks/build_db/build_kraken2_db/kraken2_addto_libraries_benchmark.txt"
    shell:
        """
        find {params.fasta_dir}/ -name '*.fna' -print0 | \
             xargs -0 -I{{}} -n1 -P {threads} \
             kraken2-build --add-to-library {{}} \
             --db {params.db_path} \
             2> {log} 1>&2
        rm -r {params.fasta_dir}
        touch {output}
        """

rule kraken2_download_libraries:
    input:
        rules.kraken2_download_taxonomy.output
    output:
        join(rules.kraken2_download_taxonomy.params.db_path, "library/{lib}/prelim_map.txt")
    params:
        db_path=rules.kraken2_download_taxonomy.params.db_path,
        db_type=rules.kraken2_download_taxonomy.params.db_type,
        library=lambda wildcards: wildcards.lib
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_download_libraries_{lib}.log"
    benchmark:
        "output/benchmarks/build_db/build_kraken2_db/kraken2_download_libraries_{lib}_benchmark.txt"
    shell:
        """
        kraken2-build --download-library {params.library} {params.db_type} \
        --threads {threads} \
        --db {params.db_path} \
        2> {log} 1>&2
        """

rule kraken2_build_db:
    # Confusing syntax, but this makes it so this rule is called on its own if the user specifies "standard" in config.yaml
    # (i.e., just kraken2-build standard is called; no downloading taxonomy or downloading libraries)
    input:
        expand(join(rules.kraken2_download_taxonomy.params.db_path, "library/{lib}/prelim_map.txt"),
                    lib=config['params']['kraken2']['ref_libs']
              ) if 'standard' not in config['params']['kraken2']['ref_libs'] else [],
        rules.kraken2_addto_libraries.output if 'standard' not in config['params']['kraken2']['ref_libs'] else []
    output:
        hashfile=join(rules.kraken2_download_taxonomy.params.db_path, "hash.k2d"),
        optfile=join(rules.kraken2_download_taxonomy.params.db_path, "opts.k2d"),
        taxofile=join(rules.kraken2_download_taxonomy.params.db_path, "taxo.k2d"),
        seqid2taxid=join(rules.kraken2_download_taxonomy.params.db_path, "seqid2taxid.map")
    params:
        # Going with conditional input, this makes it so a custom library is built unless "standard" is specified in config.yaml
        build_type=lambda wildcards, input: "--build" if input else "--standard",
        db_path=rules.kraken2_download_taxonomy.params.db_path,
        db_type=rules.kraken2_download_taxonomy.params.db_type
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['build_db']
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_build_db.log"
    benchmark:
        "output/benchmarks/build_db/build_kraken2_db/kraken2_build_db_benchmark.txt"
    shell:
        """
        kraken2-build {params.build_type} {params.db_type} \
        --threads {threads} \
        --db {params.db_path} \
        2> {log} 1>&2
        """

# ${KRAKEN_DB}/database${READ_LEN}mers.kmer_distrib
rule bracken_build:
    input:
        rules.kraken2_build_db.output
    params:
        kmer_len=config['params']['bracken']['kmer_len'],
        read_len=config['params']['bracken']['read_len'],
        db_path=rules.kraken2_build_db.params.db_path
    output:
        join(rules.kraken2_build_db.params.db_path, f"database{config['params']['bracken']['read_len']}mers.kmer_distrib")
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['build_db']
    log:
        "output/logs/build_db/bracken_build/bracken_build.log"
    benchmark:
        "output/benchmarks/build_db/bracken_build/bracken_build.txt"
    shell:
        """
        bracken-build \
            -d {params.db_path} \
            -t {threads} \
            -k {params.kmer_len} \
            -l {params.read_len} \
            2> {log} 1>&2
        """
 
rule download_metaphlan_db:
    output:
        directory(config['params']['metaphlan']['db_path'])
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['build_db']
    log:
        "output/logs/build_db/download_metaphlan_db/download_metaphlan_db.log"
    benchmark:
        "output/benchmarks/build_db/download_metaphlan_db/download_metaphlan_db_benchmark.txt"
    shell:
        """
        if test -f "{output}/mpa_latest"; then
            touch {output}
            echo "DB already installed at {output}"
        else
            metaphlan --install \
            --bowtie2db {output} \
            --nproc {threads} \
            2> {log} 1>&2
        fi
        """
