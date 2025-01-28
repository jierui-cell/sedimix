#!/bin/sh
#SBATCH --job-name=create_alternate_reference
#SBATCH --account=co_moorjani
#SBATCH --partition=savio4_htc 
#SBATCH --time=72:00:00 --mem 50gb 
#SBATCH --cpus-per-task=4

module load python/3.11.4 
python create_modified_hg19.py 
seqtk seq -l 50 modified_hg19.fa > modified_hg19_l50.fa
rm modified_hg19.fa 
mv modified_hg19_l50.fa modified_hg19.fa