selected_mapper='minimap2'

#samtools view -c {read_sample}.bam to count host reads
#grep 'Total Sequences' to find nonhost reads
# Best way so far to compile nonhost reads: echo -e "Nonhost Reads:\t$(grep -m 1 -e "Total\WSequences" bob/file1.txt bob/file2.txt | cut -d : -f 3 | paste - - -sd ,)"
# Or grep -hm 1  -e "Total\WSequences" bob/file1.txt bob/file2.txt | paste - - -d , | tr -d -c [:digit:],

# echo "$(gunzip -c qc/host_filter/nonhost/ABX-CJ-CLEM-14.R1.fastq.gz qc/host_filter/nonhost/ABX-CJ-CLEM-14.R1.fastq.gz | wc -l )/(4 * 2)" | bc

# This can be greatly sped up by reading from the output of fastqc, but that would require more complicated script logic
# since the only way to access that info is by opening a zipfile. And unzip only works on one file at a time, meaning that
# it wouldn't work with paired end data (which comes in 2 separate zipfiles). And it would force the user to have used
# fastqc beforehand, which restricts freedom unneccessarily.
#echo "Some text: $(($(gunzip -c qc/host_filter/nonhost/ABX-CJ-CLEM-14.R1.fastq.gz qc/host_filter/nonhost/ABX-CJ-CLEM-14.R2.fastq.gz | wc -l ) / (4 * 2)))"
#TODO: Make the output temp when finished
rule count_sample_reads:
    input:
        nonhost_reads=lambda wildcards: expand([rules.host_filter.output.nonhost_R1, rules.host_filter.output.nonhost_R2],
                             read_sample=wildcards.read_sample),
        host_reads=rules.host_filter.output.host
    params:
        # This will allow easy switching later to allow single-end reads to be used
        readlength=lambda wildcards, input: len(input.nonhost_reads)
    output:
        "output/quant_bins/sample_read_counts/{read_sample}_read_counts.csv"
    log:
        "output/logs/quant_bins/sample_read_counts/{read_sample}_count_sample_reads.log"
    conda:
            "../env/qc.yaml"
    shell:
        """
        echo "{wildcards.read_sample}" > {output}
        echo "Nonhost readcount,$(($(gunzip -c {input.nonhost_reads} | wc -l ) / (4 * {params.readlength})))" >> {output}
        echo "Host readcount,$(samtools view -c {input.host_reads})" >> {output}
        """

# Take the coverage and length and combine them into a contig_sample-wise dataframe
# This allows wildcards to be collapsed so the read_sample can be ignored after this
# Also, grab the host and non-host reads
rule calculate_bin_abundance:
    input:
        coverages=expand("output/mapping/minimap2/coverage_tables/bins/{read_sample}_bin_coverage.txt",
                                          read_sample = read_groups),

        contigs2bins=expand(rules.run_DAS_Tool.output.contigs2bin,
                                               mapper=selected_mapper,
                                               contig_sample=contig_groups),
        read_counts=expand(rules.count_sample_reads.output,
                                             read_sample=read_groups),

        genome_stats="output/annotate_bins/annotate_bin_pathways/annotated_pathways/genome_stats.tsv",

    params:
        read_sample=read_groups,
        contig_sample=contig_groups
    output:
        quant_bins = "output/quant_bins/quantified_bins.tsv",
        bin_info = "output/quant_bins/quantified_bin_info.tsv"
    log:
        "output/logs/quant_bins/quant_bins.log"
    conda:
        "../env/qc.yaml"
#    shell:
#        """
#            touch {output}
#        """

    script:
        "../scripts/quantify_bins.py"

#module profile_bins:
#    github("MGXlab/CAT_pack", path="CAT_pack/CAT_pack", tag="v6.0.1")
#    output:
#       "output/testfile.txt"
#profile_bins.CAT_pack()
   
# Then actually quantify bins for the {contig_sample}. Grab the scaffold2bin info in this step as well.
rule calculate_gene_abundance:
     input:
         coverages=expand("output/mapping/minimap2/coverage_tables/genes/{read_sample}_gene_coverage.txt",
                 read_sample=read_groups),
    
         lengths=expand("output/mapping/minimap2/lengths/genes/{read_sample}_gene_lengths.txt",
                 read_sample=read_groups),

         derep_clust="output/annotate_bins/annotate_bin_pathways/merged_annotations/dereplicated_genes_cluster.tsv",

         function="output/annotate_bins/annotate_bin_pathways/merged_annotations/annotations.tsv",

         pathways="output/annotate_bins/annotate_bin_pathways/annotated_pathways/metabolism_summary.xlsx",

     params:
         read_sample=read_groups
     output:
        quant_genes = "output/quant_bins/quantified_genes.tsv",
        gene_info = "output/quant_bins/quantified_gene_info.tsv"
     log:
         "output/logs/quant_bins/quant_genes.log"
     conda:
         "../env/qc.yaml"
     #shell:
         #"""
         #    touch {output}
         #"""
     script:
         "../scripts/quantify_genes.py"

#rule compile_quantified_bins:
#    input:
#        #rules.profile_bins.output,
#        "output/quant_bins/quantified_bins.txt",
#        "output/quant_bins/quantified_genes.txt",
#
#        #lambda wildcards: get_pairing_list(sample, wildcards.mapper, contig_pairings)
##        coverage_tables=lambda wildcards: expand("output/mapping/{selected_mapper}/coverage_tables/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
#            #selected_mapper=selected_mapper,
#            #read_sample=wildcards.read_sample,
#            #contig_sample=wildcards.contig_sample)
#        #contig_lengths=expand("output/mapping/{mapper}/contig_lengths/{read_sample}_Mapped_To_{contig_sample}_contig_lengths.txt",
#    output:
#        "output/quant_bins/quantified_mags.txt"
#    shell:
#        """
#            touch {output}
#        """
