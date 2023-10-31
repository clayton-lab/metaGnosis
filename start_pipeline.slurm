#!/bin/bash


#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=12:00:00
#SBATCH --job-name=aa117
#SBATCH --open-mode=truncate    # Can change to 'append' to append to existing file
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --licenses=common		# Required for /common to be accessible by SLURM


set -Eeuo pipefail

### Main script logic starts here ###
echo "Activating environment"
module purge
module load mamba
conda activate snakemakeenv

echo "Executing pipeline"

snakemake --verbose --cores 16 --use-conda --rerun-incomplete
snakemake --verbose -s Snakefile-bin --cores 16 --use-conda

#snakemake --cores 16 --use-conda --until multiqc --rerun-incomplete
#snakemake --cores 16 --use-conda
#snakemake -s Snakefile-bin --cores 16 --use-conda
#snakemake -n -f -s Snakefile-bin --config samples=resources/config/samples2.txt units=resources/config/units2.txt binning=resources/config/binning2.txt
echo "Pipeline finished!"
