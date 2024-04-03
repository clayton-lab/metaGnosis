## To make Kraken2 work, it needs at least 100 Gb of memory to build the database.
# TODO: Include an option to remove the extra files (reference_sequences and taxonomy mapping info) for users to reduce space if they want.
rule kraken2_download_taxonomy:
    output:
        taxnames=join(config['user_paths']['kraken2_db_path'], "taxonomy/names.dmp"),
        taxnodes=join(config['user_paths']['kraken2_db_path'], "taxonomy/nodes.dmp"),
    params:
        db_path=config['user_paths']['kraken2_db_path'],
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

def create_accn_list(host_accn, genome_accns):
    if host_accn:
        if genome_accns:
            if host_accn not in genome_accns:
                return host_accn + ',' + genome_accns
            else:
                return genome_accns
        else:
            return host_accn
    else:
        if genome_accns:
            return genome_accns
        else:
            return []
   
checkpoint create_genome_metadata:
    output:
        join(config['user_paths']['genome_db_path'], "build_genome_db_metadata.tsv")
    shadow: "shallow"
    params:
        taxids=config['add_genomes']['genome_taxids'],
        accns=create_accn_list(config['host_filter_accn'], config['add_genomes']['genome_accessions']),
        genome_db=config['user_paths']['genome_db_path'],
        file_formats='fasta' # TODO: Make this a parameter later, which would allow the user to download protein fastas if they want
    conda:
        "../env/qc.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_genome_db/create_genome_meta.log"
    script:
        "../scripts/create_genome_metadata.py"

rule download_genome_fasta:
    input:
        join(config['user_paths']['genome_db_path'], "build_genome_db_metadata.tsv")
    output:
        join(config['user_paths']['genome_db_path'], "{accn}_{asm_name}_genomic.fna.gz")
    # Keep tweaking this, put in in different places, etc.
    params:
        genome_db=config['user_paths']['genome_db_path'],
        file_formats='fasta' # TODO: Make this a parameter later, which would allow the user to download protein fastas if they want
    conda:
        "../env/qc.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_genome_db/download_genome_fasta.{accn}_{asm_name}.log"
    retries: 3
    script:
        "../scripts/download_genome_fasta.py"

rule kraken2_reformat_genome_fasta:
    input:
        genome_meta=join(config['user_paths']['genome_db_path'], "build_genome_db_metadata.tsv"),
        seqfile=rules.download_genome_fasta.output
    output:
        temp("output/build_db/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna")
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_kraken2_db/reformat_genome_fasta/reformat_genome_fasta.{accn}_{asm_name}.log"
    script:
        "../scripts/reformat_genome_fasta_kraken2.py"

rule kraken2_addto_libraries:
    input:
        taxfile=rules.kraken2_download_taxonomy.output,
        addfile="output/build_db/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna"
    output:
        temp("output/build_db/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}.added")
    params:
        fasta_dir=rules.download_genome_fasta.params.genome_db,
        db_path=rules.kraken2_download_taxonomy.params.db_path
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_kraken2_db/kraken2_addto_libraries/kraken2_addto_libraries.{accn}_{asm_name}.log"
    shell:
        """
        kraken2-build --add-to-library {input.addfile} --db {params.db_path} \
        2> {log} 1>&2
        tail -n 1 {log} > {output}
        """

# This could be even further optimized with snakemake's "lookup" function, which handles queries to pandas dataframes
# or python data structures to get associated wildcard info. It works in combination with checkpoints, but only if
# used via an input function like this one. In theory though, the workflow could start by iterating over the
# 'config['add_genomes'] options and retrieve wildcard info associated with each user input (taxid, accession number, etc.)
# needed for each step by passing those as wildcards to this function, which could retrieve specific info depending on
# which wildcard it gets.  Like it could take a family-level taxid as a wildcard and return a list of species-level
# taxids and accession numbers for use in the rules above.
def read_accns(wildcards):
    with checkpoints.create_genome_metadata.get(**wildcards).output[0].open() as f:
        genome_meta = pd.read_csv(f, sep='\t', header=0, na_filter=False, usecols=['assembly_accession', 'asm_name'])
        accns=genome_meta['assembly_accession'].values
        asm_names=genome_meta['asm_name'].replace('_{2,}', '_', regex=True).values
        wildcard_constraints:
            accn="\w{2,3}_\d+\.\d+"
        return expand(["output/build_db/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}.added"], zip, accn=accns, asm_name=asm_names)

rule kraken2_summarize_added_seqs:
    input:
        read_accns
    output:
        "output/build_db/build_kraken2_db/kraken2_addto_libraries_summary.txt"
    shell:
        """
        cat {input} > {output}
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
        rules.kraken2_summarize_added_seqs.output if 'standard' not in config['params']['kraken2']['ref_libs'] else []
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
