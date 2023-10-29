# BugSeq-er2 Pipeline Notes
## Purpose
These notes serve as a general summary of the current state of BugSeq-er2 and future goals. In other words, they function to say a lot in a short amount of space instead of blowing the repo up with 1000 issues. The notes consist of 3 main sections:
1. Current Progress and Immediate Next Steps
  - Current state of BugSeq-er2. This section will become outdated and might be removed from the notes entirely. For now, it functions to provide background info that is useful to keep in mind when making changes.

2. Past Changes
  - A separate fork of this repo was briefly used to attempt updating package info. The changes made to that fork will be documented here. This section will also become outdated eventually and subsequently removed.

3. Future Fixes and Additions
  - Heading should be straighforward. In this section, requested fixes and additions will have a subheading, followed by an optional description in plain text if further clarification is needed. This section will be used as a source for creating and assigning issues.

# Current Progress and Next Steps
## Current Progress
This repo was cloned from an existing one (source [here](https://github.com/CUMoellerLab/sn-mg-pipeline)). As of now, the repo works with real data (the samples added by me) in both the Snakefile and Snakefile-bin workflows. The provided test data (the samples that initially came with the cloned repo) work with the Snakefile workflow, but not the Snakefile-bin workflow. The real data is therefore preferable to use, even if it takes longer to run. The test data were kept to make it easier for everyone to get a working version of the pipeline, but will be removed after that.

## Immediate Next Steps
This is mostly a copy of the current issues on GitHub. This section will be removed when these steps are completed.
- Finish the incomplete binning pipeline
- Analyze real data to optimize tool parameters
- Containerize the existing packages
- Merge the assembly and binning pipelines
- Rewrite existing code to match the Snakemake-workflows template

# Past Changes
Again, these changes will need to be made when updating package versions, but for now the pipeline works with minimal changes to actual code.
1. Changed "--unknown_estimation" to "--unclassified_estimation" in config.yaml. This change occured in MetaPhlAn from version 3.1.0 and onward.
2. Changed "Fasta_to_Scaffolds2Bin.sh" to "Fasta_to_Contig2Bin.sh" in selected_bins.smk. Removed 1's in DAS_Tool's parameters (also in selected_bins.smk). Both of these changes are contingent on DAS_Tool's version, and must be made past version 1.1.7. Not sure about anything between 1.1.3 (the current pipeline version) and 1.1.7.
3. Changed the iteration method in Snakefile-bin. Changed ".iteritems()" (deprecated in pandas) to ".items()".
4. Messed with conda's channel_preference settings (as per snakemake's recommendations). So far, "disabled" gets the best results, "strict" crashes before solving any environments and "flexible" might(?) work. Haven't compared it to "disabled".

In addition to these changes, a brief note about Concoct: It needs libgsl.so.25 (a C dependency), which isn't included in Conda's Concoct packages. Best workaround so far is to add gsl<2.6 and scikit-learn<1.2 to concoct_linux.yaml. Gsl includes libgsl.so.25 unless it's >= 2.6, when it switches to libgsl.so.27. And scikit-learn > 1.1 throws an error because something Concoct uses is deprecated. The oldest version of Concoct doesn't install correctly and 1.0 migh or might not have worked. Not sure. Regardless, this package is a pain to maintain and future updates will likely be challenging.

# Future Fixes and Additions
From here on, each proposed change will have a header, followed by an optional description in plain text if further clarification is needed.

## Fix recursive genome location in db
When bowtie2 constructs a host genome, it creates a duplicate directory in resources/db/bt2, so that the constructed host genome is in resources/db/bt2/resources/db/bt2. Changing the settings config.yaml (which currently has db_dir and genome parameters) breaks the pipeline. This can be fixed by either making it so constructed host genomes go in a different directory, or putting everything in the bt2 directory.

## Add tool to automatically download NCBI GenBank Accession number for the host genome assembly
The comments of config.yaml look like this was an intended feature that was never added. An optional 'accn' parameter could be added for the user to specify which accession number they wanted to download. This could be accomplished with existing NCBI command-line tools.

## Create new container with updated packages
Previous containers can be kept as well, with a description of the packages and their versions available in each container.

## Add capability to handle long reads
The future of shotgun_sequencing will likely shift toward long read data. Including tools to handle that will be a useful addition.

## Add additional tools based on KBase Pipeline
A KBase metagenomics pipeline was published with assembly-based metagenomics tools including Kaiju, CheckM, and GTDBTk ([article](https://www.nature.com/articles/s41596-022-00747-x)). These will be useful to add to the current pipeline.

## Implement unused tools that are currently specified in config.yaml
Tools like kracken2 and bracken have parameters in config.yaml, but are unused in the current pipeline. 

## Add analysis of non-bacteria
A tool called [EukDetect](https://github.com/allind/EukDetect) allows microbial eukaryotes to be identified. MetaPhlAn includes an option to detect viruses (included in the current pipeline). More tools of this nature will allow all the bugs to be seq-ed.

## Add processing of host reads
There are many things that can be done with host reads instead of simply throwing them out. Some studies ([example_1](https://www.nature.com/articles/s41522-020-00171-7), [example_2](https://www.nature.com/articles/s41396-020-0665-8)) have used the ratio of microbial and host reads (microbial load) to approximate absolute microbial abundance, which would be useful going forward. Another idea is to use the host reads to create full-on host-specific assemblies, paving the way for combined genomic and metagenomic analyses.

## Expand on current parsing of samples.txt, units.txt, and binning.txt
Ideally, experimental metadata variables like subject identifiers, treatment conditions, and time points could be parsed by BugSeq-er2. This would allow many possibilities such as comparing individual and co-assemblies based on the aforementioned metadata variables to find which combination produces bins with the highest quality. Combined with the addition of proposed tools above, statistical analysis could also be performed to find differentially abundant taxa, genes, etc. between groups.

## Add a development tool to summarize Conda environments in .snakefile directory
A major problem is figuring out which packages and versions are installed in which of the multiple Conda environments that are created to run BugSeq-er2. This functionality might currently be present in snakemake, but I have yet to find it. It could be as simple as a function that goes through every installed conda environment, lists their packages (via conda list), and prints which of the .yaml files in the env directory was used to create said environment. It might be possible to use Git to find which version of the .yaml file produced the environment, though this would be messy. Regardless, adding this tool would make it easier to detect problems related to package versioning and update packages into stable containers. Snakemake has an option to clear unused environments. That will be essential as well, since it's possible to have old environments from previous versions of .yaml files and a new environment is created each time a .yaml file is updated. When an environment is created, the source snakefile, the env/.yaml file, and the hashed .snakemake/conda directory are printed to the screen, which could also be used. I believe snakemake allows users to create named environments, but this isn't ideal and would get complicated when changing env/.yaml files.

## Consolidate logs and reports
Currently, there multiple log directories created after running the pipeline (e.g., .snakemake/logs, output/log, output/benchmarks). These could be consolidated to a single output directory. Snakemake has a functionality for something called 'reports' that could be useful as well.

## Add color functionality to output
While snakemake produces colors when printing to std_out, this functionality is lost in snakemake logs and in slurm files. It might be possible either within snakemake or with an external function to add colors to these as well. 

## Render DAG of jobs
Snakemake produces a directed acyclic graph (DAG) of job output that can be rendered. Details can be found in the snakemake documentation. DAGs would be useful to describe the pipeline in presentations and publications.

## Figure out how to deal with output between runs
Ideally, an output directory could be specified for snakemake to write everything to. This could be added to config.yaml. This functionality will be important to compare output between multiple runs of the pipeline with different tool parameters, so it's a high priority thing to add.

## Add Slurm integration
Snakemake can be integrated directly with slurm, removing the need to create slurm files. Can be used to create multiple jobs and job arrays, which would be challenging to coordinate otherwise. This would be a great addition, and could pave the way for future performance optimization (see below).

## Optimize parallelization and pipeline performance
This will require thorough research into parallelization methods. This is a good article on the topic ([source](https://arxiv.org/pdf/1303.7195.pdf)). Especially powerful with slurm integration, as that would allow multiple nodes to be utilized simultaneously, exponentially speeding up the pipeline. It's already to the point where a single node takes > 3 days to analyze 64 samples with this pipeline, even with max resources allocated (64 GB memory and 16 threads). So this is a high-priority addition. GPUs could provide a further performance boost, and are available on HCC, though this would be further down the line.
