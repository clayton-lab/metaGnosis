# metaGnosis Pipeline Notes
## Purpose
These notes serve as a general summary of the current state of metaGnosis and future goals. In other words, they function to say a lot in a short amount of space instead of blowing the repo up with 1000 issues. The notes consist of 3 main sections:
1. Current Progress and Immediate Next Steps
  - Current state of metaGnosis. This section will become outdated and might be removed from the notes entirely. For now, it functions to provide background info that is useful to keep in mind when making changes.

2. Past Changes
  - A separate fork of this repo was briefly used to attempt updating package info. The changes made to that fork will be documented here. This section will also become outdated eventually and subsequently removed.

3. Future Fixes and Additions
  - Heading should be straighforward. In this section, requested fixes and additions will have a subheading, followed by an optional description in plain text if further clarification is needed. This section will be used as a source for creating and assigning issues.

# Current Progress and Next Steps
## Current Progress
This repo was cloned from an existing one (source [here](https://github.com/CUMoellerLab/sn-mg-pipeline)). As of now, the repo works with real data (the samples added by me) in the Snakefile and workflow. The provided test data (the samples that initially came with the cloned repo) work with the Snakefile workflow, but not the Snakefile-bin workflow. The real data is therefore preferable to use, even if it takes longer to run. The test data were kept to make it easier for everyone to get a working version of the pipeline, but will be removed after that.

## Immediate Next Steps
This is mostly a copy of the current issues on GitHub. This section will be removed when these steps are completed.
- Containerize the existing packages
- Rewrite existing code to match the Snakemake-workflows template

# Past Changes
Again, these changes will need to be made when updating package versions, but for now the pipeline works with minimal changes to actual code.
1. Changed "--unknown\_estimation" to "--unclassified\_estimation" in config.yaml. This change occured in MetaPhlAn from version 3.1.0 and onward.
2. Changed "Fasta\_to\_Scaffolds2Bin.sh" to "Fasta\_to\_Contig2Bin.sh" in selected\_bins.smk. Removed 1's in DAS\_Tool's parameters (also in selected\_bins.smk). Both of these changes are contingent on DAS\_Tool's version, and must be made past version 1.1.7. Not sure about anything between 1.1.3 (the current pipeline version) and 1.1.7.
3. Changed the iteration method in Snakefile-bin. Changed ".iteritems()" (deprecated in pandas) to ".items()".
4. Messed with conda's channel\_preference settings (as per snakemake's recommendations). So far, "disabled" gets the best results, "strict" crashes before solving any environments and "flexible" might(?) work. Haven't compared it to "disabled".
5. Researched kraken2-build, which occasionally fails (examples [here](https://github.com/DerrickWood/kraken2/issues/571) and [here](https://github.com/DerrickWood/kraken2/issues/465)). This is usually due to things on NCBI's end (files that don't exist anymore, etc) and is outside the control of Kraken2 or us. If this happends, script files can be manually modified to ignore missing files, though the issue usually resolves itself in a few days tops.

In addition to these changes, a brief note about Concoct: It needs libgsl.so.25 (a C dependency), which isn't included in Conda's Concoct packages. Best workaround so far is to add gsl<2.6 and scikit-learn<1.2 to concoct\_linux.yaml. Gsl includes libgsl.so.25 unless it's >= 2.6, when it switches to libgsl.so.27. And scikit-learn > 1.1 throws an error because something Concoct uses is deprecated. The oldest version of Concoct doesn't install correctly and 1.0 migh or might not have worked. Not sure. Regardless, this package is a pain to maintain and future updates will likely be challenging.

# Future Fixes and Additions
From here on, each proposed change will have a header, followed by an optional description in plain text if further clarification is needed.

## Provide more user customization and pipeline flexibility
The pipeline currently only works with paired-end reads, but could be tweaked to work with single-end reads as well. Other steps like host filtering are currently mandatory but should also be made optional in the future, especially for handling environmental samples. These are just a few examples of how the pipeline can be improved in this area.

## Allow Kraken2 workflow to be run for both reads and assemblies
This would involve restructuring the Kraken2 workflow from its current state, where it defaults to only using nonhost reads.

## Create new container with updated packages
Previous containers could be kept as well, with a description of the packages and their versions available in each container. These could be hosted somewhere for users to download

## Add capability to handle long reads
The future of shotgun sequencing will likely shift toward long read data. Although the pipeline currently incorporates Minimap2, which is optimized for long reads, later steps may not work with long reads yet. More testing is needed in this area.

## Allow Bowtie2 and Minimap2 to be used simultaneously
Both mappers work on their own, but the pipeline breaks when both are simultaneously used. A possible fix is to use CheckM/CheckM2 to compare overall bin quality between mappers and select the mapper that produces the highest-quality bins.

## Fix multi-assembler error
Like the mapper problem above, MEGAHIT and Metaspades work fine if only one is specified, but using both simultaneously breaks the pipeline. This could be resolved by adding a step that aggregates assembler statistics (a.k.a., Quast output) across all samples, compares them, and selects the best-performing assembler for all subsequent steps (binning and beyond).

## Add additional binners and bin refiners
A strength of DAS\_Tool is that it can potentially create bins of even higher quality with the addition of more binners. DAS\_Tool isn't the only tool of its kind either, so similar tools like MetaWRAP would make nice additions to the pipeline. Useful review article [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03667-3).

## Create a domain-agnostic pipeline
Several steps in the conventional shotgun sequencing pipeline (Prodigal gene calling, DASTool for bin refinement, CheckM for bin qc, GTDB for taxonomy classification) have an inherent bias toward prokaryotes. A domain-agnostic workflow could be created to process eukaryote, prokaryote, and archaeal reads simultaneously. Gene calling is a nightmare for eukaryotes, so one way to accomplish this is to use [ACR](https://academic.oup.com/bib/article/24/6/bbad381/7330501) to refine bins (possibly for whatever bins weren't used by DASTool) and [Whokaryote](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000823) to classify the bins/contigs as bacterial or eukaryotic. Then specialized eukaryotic gene callers like GeneMark-ES and Augustus to call eukaryotic genes while Prodical and MetaGeneMark can be used for prokaryotic/archael genes. Taxonomy and functional classification could then be done by aligning them to NCBI's NR database (via DIAMOND/BAT). Unbinned contigs and unassembled reads could also be processed with [RAT](https://www.nature.com/articles/s41467-024-47155-1). A serious attempt at integrating said pipeline into this one was made, but progress stalled because several of the tools in this paragraph are poorly documented and not well maintained. Although viral reads weren't mentioned yet, additional tools can be used to process these. Another approach could be incorporating [BUSCO](https://academic.oup.com/mbe/article/38/10/4647/6329644) into the pipeline, as they seem to have a similar idea.

## Add processing of host reads
There are many things that can be done with host reads instead of simply throwing them out. Some studies ([example\_1](https://www.nature.com/articles/s41522-020-00171-7), [example\_2](https://www.nature.com/articles/s41396-020-0665-8)) have used the ratio of microbial and host reads (microbial load) to approximate absolute microbial abundance, which the current iteration of the pipeline also does. The next logical step from here is to use the host reads to create full-on host-specific assemblies, paving the way for combined genomic and metagenomic analyses.

## Expand on current parsing of samples.txt, units.txt, and binning.txt
Streamlining the samples.txt file so that duplicate samples are allowed (as in with units.txt) would be an important change, as this will allow samples to be part of multiple groups. This in turn will give the user much more flexibility in terms of things they could do (e.g., co-assembly + individual assembly, specifying which samples are mapped to respective bins with the mapping\_group variable). Ideally, experimental metadata variables like subject identifiers, treatment conditions, and time points could be parsed by metaGnosis. This would allow many possibilities such as comparing individual and co-assemblies based on the aforementioned metadata variables to find which combination produces bins with the highest quality. Combined with the addition of proposed tools above, statistical analysis could also be performed to find differentially abundant taxa, genes, etc. between groups.

## Fix or move away from DRAM and KEGG for functional annotation
This is one of the most important items on this list because DRAM has been an absolute nightmare to use. DRAM1 isn't supported anymore because the developers are trying to implement DRAM2 in Nextflow, but DRAM2 is still in beta. So it's like choosing between a crumbling sand castle and a house of cards. DRAM1 database setup is impossible, requiring 500+ GB of RAM and 16+ threads for more than 3 days straight. Even then, it still didn't work correctly. The only way we got it to work was by downloading pre-formatted databases [here](https://app.globus.org/file-manager?origin_id=97ed64b9-dea0-4eb5-a7e0-0b50ac94e889) with Globus. However, users won't want to do this if they need to create an account. KEGG, while not required, improves DRAM annotation accuracy, but is an even bigger hurdle for users because they charge a large yearly fee to use their database. An attempt was made to use DRAM2 instead, but this ultimately proved fruitless. In summary, the entire taxonomic and functional annotation section needs a massive overhaul, preferrably to a more domain-agnostic approach.

## Consolidate logs and reports
Currently, there multiple log directories created after running the pipeline (e.g., .snakemake/logs, output/log, output/benchmarks). These could be consolidated to a single output directory. Snakemake has a functionality for something called 'reports' that could be useful as well. There is also a --summarize option that can be passed, which might be what I'm looking for.

## Add color functionality to output
While snakemake produces colors when printing to std\_out, this functionality is lost in snakemake logs and in slurm files. It might be possible either within snakemake or with an external function to add colors to these as well. 

## Add Slurm integration
Snakemake can be integrated directly with slurm, removing the need to create slurm files. Can be used to create multiple jobs and job arrays, which would be challenging to coordinate otherwise. This would be a great addition, and could pave the way for future performance optimization (see below).

## Optimize parallelization and pipeline performance
This will require thorough research into parallelization methods. This is a good article on the topic ([source](https://arxiv.org/pdf/1303.7195.pdf)). Especially powerful with slurm integration, as that would allow multiple nodes to be utilized simultaneously, exponentially speeding up the pipeline. It's already to the point where a single node takes > 3 days to analyze 64 samples with this pipeline, even with max resources allocated (64 GB memory and 16 threads). So this is a high-priority addition. GPUs could provide a further performance boost, and are available on HCC, though this would be further down the line.

## Implement memory and resource optimization
Snakemake has options to not run certain rules unless a minimum amount of memory is available, or to partition memory between rules. Along with other resource management capabilities, this pipeline can be further optimized run efficiently on outdated laptops, HPC clusters, and cloud computing environments. Together with containerization, it could potentially run on any operating system or software environment imaginable.
