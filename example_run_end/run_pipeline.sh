#!/bin/sh
#SBATCH --job-name=example_run
#SBATCH --account=xxx
#SBATCH --partition=xxx
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32

# Unlock the Snakemake workflow if locked (sometimes if the pipeline is interrupted, it will be locked. Then execute the workflow. 
singularity exec --cleanenv --no-home \
   --bind "$(pwd)/..":/workdir \
   ../sedimix-v1.0_latest.sif \
   bash -c '\
      snakemake -s /workdir/rules/snakefile_sedimix --unlock && \
      snakemake -s /workdir/rules/snakefile_sedimix \
         --cores 16 \
         --resources mem_mb=200000 \
         --jobs 1 \
         --rerun-incomplete
   '