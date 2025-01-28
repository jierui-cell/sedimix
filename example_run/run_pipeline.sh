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

# first build alternative reference genome - DONE
# python ../scripts/create_modified_hg19.py ../vernot_snp_panels/big_steffi_AA211_def.txt /global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa /global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human_alternative_allele_vernot_ssAA211/modified_hg19_vernot_ssAA211.fa
# cd /global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human_alternative_allele_vernot_ssAA211/
# bwa index modified_hg19_vernot_ssAA211.fa

# Unlock the Snakemake workflow if locked
snakemake -s ../rules/snakefile_sedimix --unlock

# mapdamage seems to be suddenly not working, remove it from the rules for now 
# the rule all is also changed for ssAA211 in the Snakefile in rules 

# Loop through each sample and run Snakemake for each
snakemake -s ../rules/snakefile_sedimix \
    --cores 32 --resources mem_mb=200000 --jobs 1 --rerun-incomplete 

# # Run the Snakemake pipeline
# snakemake -s ../rules/Snakefile_centrifuge_vernot_third_allele --cores 16 --jobs 1 --rerun-incomplete

# Combine final report TSV files if needed
# (head -n 1 final_report/A15919_humanNucCapture_ssAA197-200_rawReads.tsv && tail -n +2 -q final_report/*humanNucCapture*.tsv) > final_report/combined_final_report.tsv && rm final_report/final_report_*.tsv
# (head -n 1 final_report/A15919_humanNucCapture_ssAA211_rawReads.tsv && tail -n +2 -q final_report/*humanNucCapture*.tsv) > final_report/combined_final_report.tsv
