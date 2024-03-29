## To make Kraken2 work, it needs at least 100 Gb of memory to build the database.
# TODO: Include an option to remove the extra files (reference_sequences and taxonomy mapping info) for users to reduce space if they want.
# Make it so the temporary ncbi_taxonomy_db is removed, otherwise it probably won't be

#| cut -f 1 | tail -n +1

# Gameplan: Use NGD python calls. Download summary files and create metadata. Use metadata to create checkpoint.
# Checksums for filenames of of directory (maybe stored in dict) provide way to ensure all downloads are complete. Python calls might(?) allow
# parallel downloading. Could allow errors so the pipeline doesn't break if a single file isn't found.

# Oooh, could also check metadata so the right subdirectory (mammal, bacteria, etc) is chosen instead of 'all', speeding that up. This could be a
# dict to avoid the overhead of pandas lookups. Best way to do this is to just include an extra column in the created metadata file from
# the create_metadata checkpoint

# So far this works when there are no taxids in config.yaml, but not with there are no accession numbers.
# It could be further optimized by creating a checkpoint from the output of gimme_taxa.py & the ncbi-genome-download (NGD)
# dry run, then actually downloading those files in parallel in a subsequent rule. However, NGD is restricting because 
# it doesn't allow the option to change its cache settings, and there's no easy way to access the cache from here,
# meaning that the only workaround is to disable the cache altogether. Parellel downloads would then be slower than
# doing everything in a single job because the cache (assembly summary files) would have to be re-downloaded each time.
# Downloaded genome files follow this naming convention: assembly_accession + '_' + asm_name + '_genomic.fna.gz'

# Without knowing the asm_name (which can be created through NGD's cache or its metadata file that it creates after 
# downloading something), it's impossible to create wildcards for the genome files without downloading them explicitly.
# An alternative solution is to download files via the Snakemake efetch wrapper. This has potential, but
# doesn't work with RefSeq or GenBank accession numbers (or taxonomy IDs) to fetch reference genomes assemblies.
# So wrappers to the other E-Utilities need to be created.

# Do rules that use the input generated from the checkpoint work if their output is specified for a target?
#rule post_process_accns:
#    input: 
#        expand(join(config['user_paths']['genome_db_path'], "{accn}_{asm_name}_genomic.fna.gz"),
#                accn=config['add_genomes']['genome_accessions'],
#                asm_name=lookup(

#rule process_accns2:
#    input:
#        read_accns
#    output:
#        "output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna.gz"
#    shell:
#        """
#        echo {input}
#        touch {output}
#        """
#rule summarize_accns:
#    input:
#        expand("output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna.gz",
#                accn="bob_accn",
#                asm_name="bob_asm_name"),
#    output:
#        "output/bob2.txt"
#    shell:
#        """
#        echo {input}
#        touch {output}
#        """
#def lookup_wildcards(wildcards):
#    with checkpoints.create_genome_metadata.get(**wildcards).output[0].open() as f:
#        print('function wildcards: ', wildcards)
#        genome_meta = pd.read_csv(f, sep='\t', header=0, na_filter=False)
#        q = lookup(
#                     query=f"assembly_accession == '{config['host_filter_accn']}'",
#                     cols='assembly_accession',
#                     within=genome_meta,
#                 )
#        print(q)
#        accns=genome_meta['assembly_accession'].values
#        asm_names=genome_meta['asm_name'].values
#        print('Funtion Wildcards: ', wildcards)
#        return q
  
#rule temp_rule:
#    input:
#        lookup_wildcards
        #lambda wildcards: expand("output/{genome_idtype}",
        #            genome_idtype=wildcards.genome_idtype,
        #            genome_id=config['add_genomes'][wildcards.genome_idtype].split(','),
        #            accn=lookup_wildcards)
#    output:
#        "output/build_db/{genome_idtype}/{genome_id}.txt"

# Can use the lookup rule for input functions, then pass the wildcards
# to the lookup rule to get the desired attributes. Also works with
# checkpoint-dependent input functions
#rule build_db_summary:
#    input:
#        expand("output/build_db/{genome_idtype}/{genome_id}.txt",
#                genome_idtype=[key for key,value in config['add_genomes'].items() if value],
#                genome_id=config['add_genomes']['genome_accessions'])
#    output:
#        "output/build_db/build_db_summary.txt"

# This rule and build_genome_db can be further optimized, but no solid ways to
# optimize them exist yet. The Snakemake efetch wrapper has potential, but
# doesn't work with RefSeq or GenBank accession numbers to fetch reference
# genomes. So wrappers to the other E-Utilities need to be created. It's also
# possible to create a checkpoint by combining get_taxid_list and a dry run of
# ncbi-genome-download (NGD) to get a list of accession numbers for further 
# processing. However, NGD is restricting because it doesn't allow the option
# to change its cache settings.
# Downloaded genome files follow this naming convention:
# assembly_accession + '_' + asm_name + '_genomic.fna.gz'

# Without knowing the asm_name (which can be created through NGD's cache or its
# metadata file that it creates after downloading something), it's impossible 
# to create wildcards for the genome files without downloading them explicitly.
rule get_taxid_list:
    output:
        temp("output/build_db/build_genome_db/taxid_list.txt")
    shadow: "shallow"
    params:
        taxids=config['add_genomes']['genome_taxids'],
        #temp_taxdir=temp(directory("output/build_db/build_genome_db/ncbi_taxonomy_db"))
    conda:
        "../env/qc.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_genome_db/get_taxid_list.log"
    shell:
        """
        gimme_taxa.py -v --just-taxids \
        -d ncbi_taxonomy_db \
        -o {output} \
        {params.taxids} \
        2> {log} 1>&2
        """

#def read_tissues_output(wildcards):
#    with open('tissuesused.txt') as f:
#        samples = [sample for sample in f.read().split('\n') if len(sample) > 0]  # we dont want empty lines
#        return expand("tissue_{sample}.txt", sample=samples)

# Basically, create a param that opens the input file and downloads accessions. Then
# an input function can create wildcards from those, allowing them to be processed in parallel

# TODO: Add option to download protein sequences instead. Will also need to edit kraken2_addto_libraries
# to enable the addition of protein sequences
rule build_genome_db:
    input:
        taxid_list=rules.get_taxid_list.output if rules.get_taxid_list.params.taxids else [],
    output:
        join("output/old_genome_db", "build_genome_db_metadata.tsv"),
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
        """

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
        temp("output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna")
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_kraken2_db/reformat_genome_fasta/reformat_genome_fasta.{accn}_{asm_name}.log"
    script:
        "../scripts/reformat_genome_fasta_kraken2.py"


# rule kraken2_reformat_genome_fastas:
#     input:
#         rules.build_genome_db.output,
#     output:
#         temp("output/build_db/build_kraken2_db/kraken2_reformat_genome_fastas/kraken2_reformat_genome_fastas.done")
#     params:
#         temp(directory("output/build_db/build_kraken2_db/kraken2_reformat_genome_fastas"))
#     conda:
#         "../env/profile.yaml"
#     threads:
#         1
#     log:
#         "output/logs/build_db/build_kraken2_db/kraken2_reformat_genome_fastas.log"
#     benchmark:
#         "output/benchmarks/build_db/build_kraken2_db/kraken2_reformat_genome_fastas.txt"
#     shell:
#         """
#         python resources/scripts/reformat_genome_fastas_kraken2.py \
#         --genome_metadata {input} \
#         --output_path {params} \
#         2> {log} 1>&2
#         touch {output}
#         """
 

rule kraken2_addto_libraries:
    input:
        taxfile=rules.kraken2_download_taxonomy.output,
        addfile="output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}_reformatted.fna"
    output:
        temp("output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}.added")
    params:
        fasta_dir=rules.download_genome_fasta.params.genome_db,
        db_path=rules.kraken2_download_taxonomy.params.db_path
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_kraken2_db/kraken2_addto_libraries.{accn}_{asm_name}.log"
    shell:
        """
        kraken2-build --add-to-library {input.addfile} --db {params.db_path} \
        2> {log} 1>&2
        tail -n 1 {log} > {output}
        """

# rule kraken2_addto_libraries:
#     input:
#         rules.kraken2_download_taxonomy.output,
#         rules.kraken2_reformat_genome_fastas.output
#     output:
#         temp("output/build_db/build_kraken2_db/kraken2_addto_libraries/kraken2_addto_libraries.done")
#     params:
#         fasta_dir=rules.kraken2_reformat_genome_fastas.params,
#         db_path=rules.kraken2_download_taxonomy.params.db_path
#     conda:
#         "../env/profile.yaml"
#     threads:
#         config['threads']['build_db'] * 0.5
#     log:
#         "output/logs/build_db/build_kraken2_db/kraken2_addto_libraries.log"
#     benchmark:
#         "output/benchmarks/build_db/build_kraken2_db/kraken2_addto_libraries_benchmark.txt"
#     shell:
#         """
#         find {params.fasta_dir}/ -name '*.fna' -print0 | \
#              xargs -0 -I{{}} -n1 -P {threads} \
#              kraken2-build --add-to-library {{}} \
#              --db {params.db_path} \
#              2> {log} 1>&2
#         rm -r {params.fasta_dir}
#         touch {output}
#         """
# 
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
        print(accns)
        asm_names=genome_meta['asm_name'].values
        print(asm_names)
        print('Funtion Wildcards: ', wildcards)
        wildcard_constraints:
            accn="\w{2,3}_\d+\.\d+"
        return expand(["output/build_kraken2_db/reformat_genome_fasta/{accn}_{asm_name}.added"], zip, accn=accns, asm_name=asm_names)
        #accns = [accn for accn in f.read().split('\n') if len(accn) > 0]
        #return expand("{accn}.fna", accn=accns)

# Now just rebuild this with new versions of the kraken2 stuff below. That's easier then trying to edit the existing rules.   
# The only way to do it for now is to have a summary rule aggragating everything that was added to kraken.
# As long as that output file is there (maybe only a log instead of output?), the rest shouldn't be retriggered.
# The summary rule can be used as input for kraken2 build.
# Wait... could kraken2 build be the thing that aggregates them? I guess I could make temp seq.added.done for everything as the
# input for kraken2 build, which would guarantee that it gets added before building. Or I could just summarize that in a rule
# directly upstream, so kraken2 build only takes a single input file. That's a better idea.

# Is is possible to find out what the new sequence names are after adding to kraken?
# TODO: For some reason the wildcards don't get parsed correctly. It's taking more than just
# the accession for "{accn}", but no clear cause yet. "read_accns" output looks good.
rule kraken2_summarize_added_seqs:
    input:
        read_accns
    output:
        "output/build_kraken2_db/kraken2_addto_libraries_summary.txt"
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
        #rules.kraken2_addto_libraries.output if 'standard' not in config['params']['kraken2']['ref_libs'] else []
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
