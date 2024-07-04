from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob

# This doesn't seem to be ran in the pipeline, and it seems to build from the input of the previous rule
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
# # tree, analyze, and qa use threads. Additionally, tree uses pplacer threads
# [-h] [-r] [--ali] [--nt] [-g] [-x EXTENSION] [-t THREADS] [--pplacer_threads PPLACER_THREADS] [-q] [--tmpdir TMPDIR] bin_input output_dir
rule run_checkm:
    input:
        rules.download_checkm_db.output,
#        #lambda wildcards: "output/selected_bins/{wildcards.mapper}/DAS_Tool_Fastas/{wildcards.contig_sample}.done",
        rules.consolidate_DAS_Tool_bins.output,
    params:
#        #bin_dir=lambda wildcards: "output/selected_bins/{wildcards.mapper}/Run_DAS_Tool/{wildcards.contig_sample}_DASTool_bins"
        bin_dir="output/selected_bins/{mapper}/Run_DAS_Tool/{contig_sample}_DASTool_bins"
    output:
        directory("output/refine_bins/{mapper}/Run_CheckM/run_checkm/{contig_sample}")
    threads:
        config['threads']['run_checkm']
    conda:
        "../env/refine_bins.yaml"
    log:
        "output/logs/refine_bins/{mapper}/checkm_tree/{contig_sample}.log"
    shell:
        """
        checkm lineage_wf --ali --nt --tab_table -x fa \
        -t {threads} --pplacer_threads {threads} \
        -a {wildcards.contig_sample}_alignment \
        -f {wildcards.contig_sample}_output \
        {params.bin_dir} {output}
        2> {log} 1>&2
        """
#
#rule checkm_tree:
#    input:
#        rules.download_checkm_db.output,
#        #lambda wildcards: "output/selected_bins/{wildcards.mapper}/DAS_Tool_Fastas/{wildcards.contig_sample}.done",
#        rules.consolidate_DAS_Tool_bins.output,
#    params:
#        #bin_dir=lambda wildcards: "output/selected_bins/{wildcards.mapper}/Run_DAS_Tool/{wildcards.contig_sample}_DASTool_bins"
#        bin_dir="output/selected_bins/{mapper}/Run_DAS_Tool/{contig_sample}_DASTool_bins"
#    output:
#        "output/refine_bins/{mapper}/Run_CheckM/checkm_tree/{contig_sample}"
#    threads:
#        config['threads']['run_checkm']
#    conda:
#        "../env/refine_bins.yaml"
#    log:
#        "output/logs/refine_bins/{mapper}/checkm_tree/{contig_sample}.log"
#    shell:
#        """
#        checkm tree -x fa -t {threads} --pplacer_threads {threads} \
#        --ali --nt {params.bin_dir} {output}
#        2> {log} 1>&2
#        """
##
#    input:
#        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
#                                 mapper=config['mappers'],
#                                 contig_sample=contig_pairings.keys())

