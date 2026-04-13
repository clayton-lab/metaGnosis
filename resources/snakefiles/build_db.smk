## To make Kraken2 work, it needs at least 100 Gb of memory to build the database.
# TODO: Include an option to remove the extra files (reference_sequences and taxonomy mapping info) for users to reduce space if they want.
# TODO: Consider making the database construction a separate run, since so much extra memory is required compared to post-construction. Then
# all databases could be constructed simultaneously.
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
        taxids=config['add_genomes']['genome_taxids'] if config.get('add_genomes') and config['add_genomes'].get('genome_taxids') else [],
        accns=create_accn_list(config['host_filter_accn'], config['add_genomes']['genome_accessions']) if config.get('host_filter_accn') else [],
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
    params:
        db_path=config['user_paths']['metaphlan_db_path'],
        db_index=config['params']['metaphlan']['db_index']
    output:
        join(config['user_paths']['metaphlan_db_path'], 'mpa_latest')
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
        metaphlan --install \
        --bowtie2db {params.db_path} \
        --index {params.db_index} \
        --nproc {threads} \
        2> {log} 1>&2

        wget -P {params.db_path} http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest
        """

# Download humann_db rule:
rule download_humann_db:
    params:
        db_path=config['user_paths']['humann_db_path'],
        nuc_db_path=join(config['user_paths']['humann_db_path'], 'chocophlan'),
        prot_db_path=join(config['user_paths']['humann_db_path'], 'uniref'),
        map_db_path=join(config['user_paths']['humann_db_path'], 'utility_mapping')
    output:
        join(config['user_paths']['humann_db_path'], 'utility_mapping', 'map_ec_name.txt.gz')
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/build_db/download_humann_db/download_humann_db.log"
    benchmark:
        "output/benchmarks/build_db/download_humann_db/download_humann_db_benchmark.txt"
    shell:
        """
        humann_databases --download chocophlan full {params.db_path} --update-config no \
        2> {log} 1>&2

        humann_databases --download uniref uniref90_diamond {params.db_path} --update-config no \
        2>> {log} 1>&2

        humann_databases --download utility_mapping full {params.db_path} --update-config no \
        2>> {log} 1>&2
        """


rule download_checkm_db:
    output:
        directory(config['user_paths']['checkm_db_path'])
    params:
        checkm_db_file="checkm_data_2015_01_16.tar.gz"
    conda:
        "../env/refine_bins.yaml"
    threads:
        1
    log:
        "output/logs/build_db/download_checkm_db/download_checkm_db.log"
    shell:
        """
        wget -P {output} https://data.ace.uq.edu.au/public/CheckM_databases/{params.checkm_db_file}
        tar -xvf {output}/{params.checkm_db_file} --one-top-level={output}
        rm {output}/{params.checkm_db_file} 
        checkm data setRoot {output}
        2> {log} 1>&2
        """

rule download_checkm2_db:
    output:
        join(config['user_paths']['checkm2_db_path'], 'CheckM2_database/uniref100.KO.1.dmnd')
    conda:
        "../env/refine_bins.yaml"
    threads:
        1
    log:
        "output/logs/build_db/download_checkm_db/download_checkm2_db.log"
    shell:
        """
        checkm2 database --download --path {output} \
        2> {log} 1>&2
        """

rule download_gtdbtk_db:
    output:
        join(config['user_paths']['gtdbtk_db_path'], 'taxonomy/gtdb_taxonomy.tsv')
    params:
        config['user_paths']['gtdbtk_db_path']
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    log:
        "output/logs/build_db/download_gtdbtk_db/download_gtdbtk_db.log"
    shell:
        """
        download-db.sh {params} -p \
        2> {log} 1>&2
        """

#$WORK/projects/nih_k01_pilot_aim1/full_shotgun_analysis/db/dram/database_files/
#dbCAN-HMMdb-V11.txt   viral.1.protein.faa.gz 

#   --kofam_hmm_loc /my/path/database_files/kofam_profiles.tar.gz                          //# hmm file for KOfam (profiles.tar.gz) (default: None)
#   --kofam_ko_list_loc /my/path/database_files/kofam_ko_list.tsv.gz                       //# KOfam ko list file (ko_list.gz) (default: None)
#   --uniref_loc /my/path/database_files/uniref90.fasta.gz                                 //# File path to uniref, if already downloaded (uniref90.fasta.gz) (default: None)
#   --pfam_loc /my/path/database_files/Pfam-A.full.gz                                      //# File path to pfam-A full file, if already downloaded (Pfam-A.full.gz) (default: None)
#   --pfam_hmm_dat /my/path/Pfam-A.hmm.dat.gz                                              //# pfam hmm .dat file to get PF descriptions, if already downloaded (Pfam-A.hmm.dat.gz) (default: None)
#   --dbcan_loc /my/path/database_files/CAMPER_v1.0.0-beta.1.tar.gz                        //# File path to dbCAN, if already downloaded (dbCAN-HMMdb-V9.txt) (default: None)
#   --dbcan_fam_activities /my/path/CAZyDB.07292021.fam-activities.txt                     //# CAZY family activities file, if already downloaded (CAZyDB.07302020.fam-activities.txt) (default: None)
#   --dbcan_sub_fam_activities /my/path/CAZyDB.07292021.fam.subfam.ec.txt                  //# CAZY subfamily activities file, if already downloaded (CAZyDB.07292021.fam.subfam.ec.txt) (default: None)
#   --vogdb_loc /my/path/database_files/vog.hmm.tar.gz                                     //# hmm file for vogdb, if already downloaded (vog.hmm.tar.gz) (default: None)
#   --vog_annotations /my/path/vog_annotations_latest.tsv.gz                               //# vogdb annotations file, if already downloaded (vog.annotations.tsv.gz) (default: None)
#   --camper_tar_gz_loc /my/path/database_files/CAMPER_v1.0.0-beta.1.tar.gz                //#
#   --viral_loc /my/path/database_files/viral.merged.protein.faa.gz                        //# File path to merged viral protein faa, if already downloaded (viral.x.protein.faa.gz) (default: None)
#   --peptidase_loc /my/path/database_files/merops_peptidases_nr.faa                       //# File path to MEROPS peptidase fasta, if already downloaded (pepunit.lib) (default: None)
#   --genome_summary_form_loc /my/path/database_files/genome_summary_form.20220504.tsv     //# File path to genome summary form,if already downloaded (default: None)
#   --module_step_form_loc /my/path/database_files/module_step_form.20220504.tsv           //# File path to module step form, ifalready downloaded (default: None)
#   --etc_module_database_loc /my/path/database_files/etc_mdoule_database.20220504.tsv     //# File path to etc module database, if already downloaded (default: None)
#   --function_heatmap_form_loc /my/path/database_files/function_heatmap_form.20220504.tsv //# File path to function heatmap form, if already downloaded (default: None)
#   --amg_database_loc /my/path/database_files/amg_database.20220504.tsv                     # File path to amg database, if already downloaded (default: None)
rule modify_kegg_fasta:
    output:
        join(config['user_paths']['dram_db_path'], "kegg_modified.pep")
    params:
        kegg_db=config['user_paths']['kegg_genes_db_path'] if config['user_paths'].get('kegg_genes_db_path') else "",
        kegg_link=config['user_paths']['kegg_genes_link_path'] if config['user_paths'].get('kegg_genes_link_path') else ""
    conda:
        "../env/annotate_bins.yaml"
    threads:
        1
    log:
        "output/logs/build_db/build_dram_db/modify_kegg_fasta.log"
    script:
        "../scripts/modify_kegg_fasta.py"

#TODO:
# Document that I had to fix the source code with this solution: https://github.com/WrightonLabCSU/DRAM/issues/340
rule build_dram_db:
    input:
        lambda wildcards: rules.modify_kegg_fasta.output if rules.modify_kegg_fasta.params.kegg_link else []
    output:
        join(config['user_paths']['dram_db_path'], 'db_descriptions/description_db.sqlite')
    params:
        kegg_db=lambda wildcards, input: f"--kegg_loc {input}" if input else "",
        db_path=config['user_paths']['dram_db_path']
    conda:
        "../env/annotate_bins.yaml"
    threads:
        config['threads']['build_db']
    log:
        "output/logs/build_db/build_dram_db/build_dram_db.log"

        # Extra thing you might need
        #--dbcan_sub_fam_activities {params.db_path}/database_files/CAZyDB.08062022.fam.subfam.ec.txt \
    shell:
        """
        DRAM-setup.py prepare_databases --threads {threads} \
        {params.kegg_db} --output_dir {params.db_path} \
        2> {log} 1>&2
        mkdir {params.db_path}/db_descriptions
        mv {params.db_path}/description_db.sqlite {params.db_path}/db_descriptions
        DRAM-setup.py export_config --output_file {params.db_path}/dram_configfile.json
        echo "Database setup complete!" >> {log}
        """
