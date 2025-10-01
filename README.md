# *sedimix*

*sedimix*: A workflow for the analysis of hominin nuclear DNA sequences from sediments

## Overview
Here we present an open-source snakemake workflow that identifies human sequences from sequencing data and provides relevant summary statistics. The final tool prioritizes the retention of human DNA while minimizing detection errors, offering a robust and accessible solution to support the growing needs of human evolutionary research. See [paper](https://www.biorxiv.org/content/10.1101/2025.02.28.640818v1) for details.

**Key Features:**
- Utilizes Snakemake, a workflow management system for Python
- Identifies ancient hominin reads in the samples through taxonomic classification followed by mapping 
- Generates a report file with summary statistics 

## Requirements

### 0. Clone the repository
```bash
git clone https://github.com/jierui-cell/sedimix.git
cd sedimix
```

### 1. Using container 
We provide a fully-baked container image you can pull & run in one step—no installs required beyond Singularity or Apptainer:

#### ▶️ Option A: Using **Singularity**
```bash
singularity pull library://jieruixu/sedimix/sedimix-v1.0:latest
```

#### ▶️ Option B: Using **Apptainer** 
Apptainer is the community-maintained successor to Singularity. To pull the same container image:
```bash
apptainer pull sedimix-v1.0_latest.sif library://jieruixu/sedimix/sedimix-v1.0:latest
```

This should download the container as `sedimix-v1.0_latest.sif` in your current working directory.        

If you prefer to install every dependency yourself and use conda environment, see [manual install](./MANUAL_INSTALL.md).

### 2. Index Files
You only need to download one of the following files below to run *sedimix*. We recommend Centrifuge as the default classification software.            
             
Download index files for Centrifuge and Kraken2 from the following:
- **Recommended**: [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)  
  We recommend NCBI: nucleotide non-redundant sequences (64GB) for Centrifuge.  
  ```bash
  wget https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
  mkdir centrifuge
  tar -xvzf nt_2018_3_3.tar.gz -C centrifuge
  rm nt_2018_3_3.tar.gz
  ```

- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)  
  The paper was tested out using nt Database version 11/29/2023 (710GB) for Kraken2. A later version should work as well though not tested. Change the code accordingly for the latest version.
  ```bash
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz
  mkdir kraken2
  tar -xvzf k2_nt_20231129.tar.gz -C kraken2
  rm k2_nt_20231129.tar.gz
  ```

Alternatively, you can build Centrifuge and Kraken2 indexes yourself by following the instructions provided on their respective GitHub repositories.

### 3. Human Reference Genome 
Download the human reference genome (e.g. [hg19.fa.gz from UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz)), unzip it, and build the BWA index.

```bash
mkdir human_ref && cd human_ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
gunzip hg19.fa.gz
bwa index hg19.fa
samtools faidx hg19.fa
```

Your folder should look like this:
```human_ref/
├── hg19.fa
├── hg19.fa.amb
├── hg19.fa.ann
├── hg19.fa.bwt
├── hg19.fa.pac
├── hg19.fa.sa
├── hg19.fa.fai
```

## Usage Instructions
Once you have:

- Pulled the Singularity image  
- Downloaded any required pre-built indices for Centrifuge or Kraken2 
- Built the human reference genome index       

Your folder should look like this:
```sedimix/
├── sedimix-v1.0_latest.sif
├── centrifuge (or kraken2)
├── human_ref
├── rules
├── scripts
├── hominin_informative_sites.txt
├── my_run (folder you create to run sedimix)
│   ├── config.yaml
│   └── 0_data
│       ├── sample_1.fq
│       └── sample_2.fq
```
 
 you can run *sedimix* in four simple steps:     
 
1. Create and enter a run directory (same level as `scripts` and `rules`) 
2. Create a folder named `0_data` within the folder you just created, then place your input FASTQ files in it. **Input files must be in format that is either .fq or .fq.gz.** 
3. Update the `config.yaml` file to set parameters (see [details](./CONFIG_PARAMETERS.md)).     
4. Execute the workflow within the run directory with the following command:

   ```bash
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
   ```
- `--cores` and `--resources` should equal to number of `threads` and `memory_mb` specified in `config.yaml`        
- Append `-n` (dry run) to the second snakemake call to preview execution without running.     
- For reproducibility, it is recommended that you always define all necessary parameters in `config.yaml`. 

An example folder can be found in `example_run/`. The start folder is only `0_data`, and start file is `config.yaml` and a bash script `run_pipeline.sh` to submit the snakemake command. All the other folders and files are generated by *sedimix*. 

## Retrieve Your Results
- **Classified hominin reads**: Located in the `3_final_reads` folder, ending with `{sample_name}_final.bam`. 
- **Classified hominin reads that have deamination**: Located in the `3_final_reads` folder, ending with `{sample_name}_final_deaminated.bam`.
- **Classified hominin reads that do not have deamination**: Located in the `3_final_reads` folder, ending with `{sample_name}_final_non_deaminated.bam`.
- **Classified non-hominin reads (if specified in config.yaml)**: Located in the `3_final_reads` folder, ending with `{sample_name}_non_hominin.fq`. 
- **Data summary report**: Located in the `4_final_report` folder. `combined_final_report.tsv` contains results for all samples.
- **Deamination profile**: Located in the `4_mapdamage_results` folder.  

## License
This project is licensed under the [MIT License](LICENSE).
