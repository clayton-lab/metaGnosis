from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob
import pathlib

# So DASTool and CAT's prodigal runs are (probably) identical.

#TODO: Include BUSCO, which is comparable to CheckM and has a “–auto-lineage” function that can be used for non-bacteria
# Snakemake workflow here! (https://gitlab.com/ezlab/plugins_buscov5/-/tree/master/BUSCO_batch_analysis_with_snakemake)
#TODO: Implement ACR, which works for prokaryote and eukaryotes

######
# ACR Config requirements (in YAML format). 
#  - eukcc
#  - eukrep
#  - scikit-learn
#  - pandas
#  - perl
#  - perl-app-cpanminus
#  - pip
#  - pip:
#      - kmeans1d

# In additon to these, ACR needs to use cpanm to install perl packages with this command: 
# cpanm YAML Hash::Merge Parallel::ForkManager MCE::Mutex Thread::Queue Math::Utils threads

# And it needs an installation of GeneMark-ES, plus a copy of .gm_key in the home directory (the location of which can be
# changed with change_path_in_perl_scripts.pl, though this is annoying).
######

#TODO: Stop the binner wildcard after dastool.
# So, bin individually, refine bins with DAS_Tool/ACR, then Busco (or whatever ACR uses) to filter bins.
# Then Whokaryote to decide which gene predictor to use per contig. Predict the genes per contig, and cat predictions to
# a single predicted_proteins fasta and gff for each assembly. Then use CAT/BAT/RAT and MicrobeAnnotator.

#Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
#                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
#                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

rule Fasta_to_Contig2Bin:
    """
    Uses Fasta_to_Contig2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins="output/binning/{binner}/{mapper}/bin_fastas/{contig_sample}"
    output:
        scaffolds="output/refine_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}_scaffolds2bin.tsv"
    params:
        extension=lambda wildcards, input: list(pathlib.Path(str(input)).iterdir())[0].suffix[1:]
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/refine_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}.log"
    threads:
        1
    shell:
       """
           Fasta_to_Scaffolds2Bin.sh \
           -i {input.bins} \
           -e {params.extension} > {output.scaffolds}
       """

rule predict_genes_prodigal:
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                                            assembler = config['assemblers'],
                                            contig_sample = wildcards.contig_sample),
    output:
        gff="output/refine_bins/predict_genes/{contig_sample}/{contig_sample}_predicted_proteins.gff",
        faa="output/refine_bins/predict_genes/{contig_sample}/{contig_sample}_predicted_proteins.faa",
    conda:
        "../env/refine_bins.yaml"
    benchmark:
        "output/benchmarks/refine_bins/predict_genes/{contig_sample}/{contig_sample}_benchmark.txt"
    log:
        "output/logs/refine_bins/predict_genes/{contig_sample}/{contig_sample}_predicted_genes.log"
    threads:
        1
    shell:
        """
        prodigal -i {input} -f gff -p meta -a {output.faa} -o {output.gff} \
        2> {log} 1>&2
        """

#Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
#                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
#                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                                            assembler = config['assemblers'],
                                            contig_sample = wildcards.contig_sample),
        bins = lambda wildcards: expand("output/refine_bins/scaffolds2bin/{binner}/{mapper}/{contig_sample}_scaffolds2bin.tsv",
                                         binner = config['binners'],
                                         mapper = config['mappers'],
                                         contig_sample = wildcards.contig_sample)
    output:
        summary="output/refine_bins/DAS_Tool/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}_DASTool_summary.txt",
        contigs2bin="output/refine_bins/DAS_Tool/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}_DASTool_scaffolds2bin.txt"
    params:
        basename = "output/refine_bins/DAS_Tool/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}",
        fasta_dir = "output/refine_bins/DAS_Tool/{mapper}/bin_fastas/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine'],
        threshold = config['params']['das_tool']['score_threshold'],
        binner_paths=lambda wildcards, input: ",".join(input.bins),
        binners=",".join(config['binners']),
    conda:
        "../env/selected_bins.yaml"
    threads:
        config['threads']['run_DAS_Tool']
    benchmark:
        "output/benchmarks/refine_bins/DAS_Tool/{mapper}/{contig_sample}_benchmark.txt"
    log:
        "output/logs/refine_bins/DAS_Tool/{mapper}/{contig_sample}.log"
    shell:
        """ 
           DAS_Tool \
           --bins {params.binner_paths} \
           --contigs {input.contigs} \
           --outputbasename {params.basename} \
           --labels {params.binners} \
           --write_bins 1 \
           --write_bin_evals 1 \
           --threads {threads} \
           --search_engine {params.search_engine} \
           --score_threshold {params.threshold} \
           2> {log} 1>&2

           mkdir -p {params.fasta_dir}
           mv {params.basename}_DASTool_bins/* {params.fasta_dir}
           rmdir {params.basename}_DASTool_bins
        """

rule download_acr_db:
    params:
        db_path=config['user_paths']['acr_db_path']
    output:
        directory(config['user_paths']['acr_clone_path'])
    log:
        "output/logs/refine_bins/clone_acr.log"
    conda:
        "../env/refine_bins.yaml"
    shell:
        """
        git clone https://github.com/hoonjeseong/acr.git {output}
        wget -O {params.db_path}/data.tar.gz https://figshare.com/ndownloader/files/41282157
        tar -zxvf {params.db_path}/data.tar.gz -C {params.db_path} \
        2> {log} 1>&2
        """
 
rule run_acr:
    input:
        built_db=rules.download_acr_db.output
    params:
        acr_prefix=rules.download_acr_db.params.db_path,
        db_path=rules.download_acr_db.params.db_path
    output:
        directory('output/refine_bins/acr')
    threads:
        16
    conda:
        "../env/refine_bins.yaml"
    log:
        "output/logs/refine_bins/run_acr.log"
    shell:
        """
        python {params.acr_prefix}/acr.py -g output/binning/concoct/minimap2/bin_fastas/ABX-CJ-IRIS \
        -c output/binning/metabat2/minimap2/coverage_tables/ABX-CJ-IRIS_coverage_table.txt \
        -e fa -o {output} -p bobert --comp=50 --cont=10 -t {threads} \
        2> {log} 1>&2
        """

#python acr.py -g ../metaGnosis/output/binning/concoct/minimap2/bin_fastas/ABX-CJ-IRIS -c ../metaGnosis/output/binning/concoct/minimap2/coverage_tables/ABX-CJ
#-IRIS_coverage_table.txt -e fa -o ../bob -p bobert --comp=50 --cont=10


# # tree, analyze, and qa use threads. Additionally, tree uses pplacer threads
# [-h] [-r] [--ali] [--nt] [-g] [-x EXTENSION] [-t THREADS] [--pplacer_threads PPLACER_THREADS] [-q] [--tmpdir TMPDIR] bin_input output_dir
# If running multiple instances of this rule, it's helpful to have at least 90Gb of RAM available
#TODO: Make it so this rule can read the amount of memory available and adapt accordingly
rule run_checkm:
    input:
        db = rules.download_checkm_db.output,
        bins = rules.run_DAS_Tool.output.summary,
        pred_genes = rules.predict_genes_prodigal.output.faa,
    params:
        bin_dir=rules.run_DAS_Tool.params.fasta_dir,
        out_dir="output/refine_bins/run_CheckM/{mapper}/{contig_sample}"
    output:
        "output/refine_bins/run_CheckM/{mapper}/{contig_sample}/{contig_sample}_output.tsv"
    threads:
        config['threads']['run_checkm']
    conda:
        "../env/refine_bins.yaml"
    log:
        "output/logs/refine_bins/run_CheckM/{mapper}/{contig_sample}.log"
    shell:
        """
        checkm lineage_wf --ali --nt --tab_table -x fa \
        -t {threads} --pplacer_threads {threads} \
        -a {params.out_dir}/{wildcards.contig_sample}_alignment.tsv \
        -f {params.out_dir}/{wildcards.contig_sample}_output.tsv \
        {params.bin_dir} {params.out_dir} \
        2> {log} 1>&2
        """

# Need the summary file and the prodigal predicted file as output
#TODO: Make it so the pipeline doesn't use the pred_prots (i.e., remove the CAT/BAT rules)
rule run_checkm2:
    input:
        db = rules.download_checkm2_db.output,
        bins = rules.run_DAS_Tool.output.summary,
    params:
        bin_dir=rules.run_DAS_Tool.params.fasta_dir,
        out_dir="output/refine_bins/run_CheckM2/{mapper}/{contig_sample}"
    output:
        summary_file="output/refine_bins/run_CheckM2/{mapper}/{contig_sample}/quality_report.tsv",
        pred_prots=directory("output/refine_bins/run_CheckM2/{mapper}/{contig_sample}/protein_files")
    threads:
        config['threads']['run_checkm2']
    conda:
        "../env/refine_bins.yaml"
    log:
        "output/logs/refine_bins/run_CheckM2/{mapper}/{contig_sample}.log"
    shell:
        """
        checkm2 predict -x fa --database_path {input.db} \
        --input {params.bin_dir} \
        --output-directory {params.out_dir} \
        -t {threads} \
        2> {log} 1>&2
        """


rule filter_bins:
    input:
        #rules.run_checkm.output
        rules.run_checkm2.output.summary_file
    params:
        bin_dir=rules.run_DAS_Tool.params.fasta_dir,
        completion_cutoff=config['params']['checkm']['percent_completion_filter'],
        contam_cutoff=config['params']['checkm']['percent_contam_filter']
    output:
        filtered_dir=directory("output/refine_bins/filtered_bins/{mapper}/{contig_sample}"),
        temp_summary=temp("output/refine_bins/summarize_bins/{mapper}/{contig_sample}_bin_summary.tsv"),
        bin_paths=temp("output/refine_bins/summarize_bins/{mapper}/{contig_sample}_filtered_bin_paths.txt")
    log:
        "output/logs/refine_bins/filter_bins/{mapper}/{contig_sample}_filter_bins.log"
    script:
        "../scripts/filter_bins.py"

#TODO: Add a thing here that chooses the best mapper (would output a file named best_mapper.txt that could be used by later rules)
rule bin_filter_summary:
    input:
        sum_files = expand("output/refine_bins/summarize_bins/{mapper}/{contig_sample}_bin_summary.tsv",
                    mapper = config['mappers'],
                    contig_sample = contig_pairings.keys()),
        path_files = expand("output/refine_bins/summarize_bins/{mapper}/{contig_sample}_filtered_bin_paths.txt",
                     mapper = config['mappers'],
                     contig_sample = contig_pairings.keys())

    output:
        filt_paths=temp("output/refine_bins/summarize_bins/filtered_bin_paths.txt"),
        summary="output/refine_bins/summarize_bins/bin_filter_summary.txt",
        stats_tsv="output/refine_bins/summarize_bins/compiled_bin_statistics.tsv",
        stats_csv="output/refine_bins/summarize_bins/compiled_bin_statistics.csv",
    log:
       "output/logs/refine_bins/bin_filter_summary.log"
    script:
        "../scripts/summarize_bins.py"

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
        clust_reps="output/refine_bins/data_tables/Wdb.csv",
        clusters="output/refine_bins/data_tables/Cdb.csv"
    threads:
        config['threads']['dereplicate_bins']
    conda:
        "../env/annotate_bins.yaml"
    log:
        "output/logs/refine_bins/dereplicate_bins/dereplicate_bins.log"
    shell:
        """
        dRep dereplicate --S_algorithm skani -sa 0.95 -nc 0.3 -p {threads} \
        -comp {params.completion} -con {params.contam} --genomeInfo {input.stats_csv} \
        -g {input.filt_paths} {params.outdir} \
        2> {log} 1>&2

        # Representative bins are concatenated together into a single file
        for f in {params.outdir}/dereplicated_genomes/*.fa; do
            filename="${{f##*/}}"
            awk -v name="${filename%.fa}" '/^>/ {$0=">"name"_"substr($0,2)} 1' "$f" >> {params.outdir}/dereplicated_bins.fa
        done
        """

# Checkm quality file needs to be in csv format and bins need the extension (.fa) in the csv. Header needs to be 'genome', 'completeness', 'contamination'
# Dereplicate bins: dRep dereplicate -g output/refine_bins/filtered_bins/minimap2/*/* --S_algorithm skani --genomeInfo output/refine_bins/summarize_bins/compiled_bin_statistics.csv -sa 0.95 -nc 0.3 -p 1 -comp 90 -con 5 all_derep2
# After derepelicated genomes are created, this cats the output to a new file with all bin_ids at the start of each header
# for f in copy_derep_genomes/*.fa; do filename="${f##*/}"; awk -v name="${filename%.fa}" '/^>/ {$0=">"name"_"substr($0,2)} 1' "$f" >> mod.fa; done
