#!/bin/sh
#SBATCH --job-name=SIM_set1_centrifuge_example_run
#SBATCH --account=co_moorjani
#SBATCH --partition=savio4_htc
#SBATCH --time=48:00:00
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=32

# Load necessary modules
module load samtools/1.14 # samtools need to be higher than a certain version, otherwise -N flag cannot be used to select reads
module load bedtools
module load gsl
module load bwa/0.7.17
module load r/4.4.0-gcc-11.4.0
module load python/3.11.6-gcc-11.4.0

# Unlock the Snakemake workflow if locked (sometimes if the pipeline is interrupted, it will be locked. This line of code will unlock it, and it never hurts to include it.)
snakemake -s ../rules/snakefile_sedimix --unlock

snakemake -s ../rules/snakefile_sedimix \
    --cores 32 --resources mem_mb=200000 --jobs 1 --rerun-incomplete 