# Sedimix

**Sedimix**: A workflow for the analysis of human nuclear DNA from sediments

## Overview
Here we present an open-source snakemake workflow that identifies human sequences from sequencing data and provides relevant summary statistics. The final tool prioritizes the retention of human DNA while minimizing detection errors, offering a robust, accessible, and adaptable solution to support the growing needs of human evolutionary research. See paper for details. 

**Key Features:**
- Utilizes Snakemake, a workflow management system for Python
- Identifies ancient hominin reads in the samples through mapping and taxonomic classification
- Generates a report file with summary statistics 

## Requirements
Ensure the following tools are installed and configured before running the pipeline:

#### [Centrifuge (default)](https://ccb.jhu.edu/software/centrifuge/manual.shtml)
```bash
git clone https://github.com/DaehwanKimLab/centrifuge
make -C centrifuge
echo 'export PATH=$PATH:$(pwd)/centrifuge' >> ~/.bashrc
source ~/.bashrc
```

#### [Kraken2 (optional)](https://github.com/DerrickWood/kraken2/wiki/Manual)
```bash
git clone https://github.com/DerrickWood/kraken2.git
./kraken2/install_kraken2.sh kraken2
echo 'export PATH=$PATH:$(pwd)/kraken2' >> ~/.bashrc
source ~/.bashrc
```

#### Seqtk 
```bash
git clone https://github.com/lh3/seqtk.git
make -C seqtk
echo 'export PATH=$PATH:$(pwd)/seqtk' >> ~/.bashrc
source ~/.bashrc
```

#### BWA 
```bash
git clone https://github.com/lh3/bwa.git
make -C bwa
echo 'export PATH=$PATH:$(pwd)/bwa' >> ~/.bashrc
source ~/.bashrc
```

#### Samtools
```bash
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xvjf samtools-1.20.tar.bz2
rm samtools-1.20.tar.bz2
./samtools-1.20/configure
make -C samtools-1.20
echo 'export PATH=$PATH:$(pwd)/samtools-1.20' >> ~/.bashrc
source ~/.bashrc
```

#### Bedtools 
```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```

#### [MapDamage2](https://ginolhac.github.io/mapDamage/)
The following requirements are for **mapDamage2** to successfully run:

1. Python (version >= 3.5)
2. Git
3. R (version >= 3.1) must be present in your `$PATH`.
4. GNU Scientific Library (GSL)

#### R Libraries:
- `inline`
- `gam`
- `Rcpp`
- `RcppGSL`
- `ggplot2` (>= 0.9.2)

If you're using `zsh` as your shell, replace `~/.bashrc` with `~/.zshrc` in the commands above.

#### Python Packages

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (version 7.32.4)
2. `pysam` (version >= 0.6): Python package, an interface for reading/writing SAM/BAM files
3. `Cython`: Only necessary if `pysam` needs to be built
4. `numpy` (version >=1.24.4)
5. `tqdm` (version 4.66.1)
6. `pandas` (version >=2.1.2)
7. `pyfaidx` (version >=0.8.1.1)
8. `biopython` (version >=1.84)
9. `scipy` (version >=1.11.3)

To ensure all necessary dependencies are installed, you can create a conda environment using the provided `environment.yaml` file.

#### Conda Environment

Create a conda environment with the following command:

```bash
conda env create -f environment.yaml
```

Activate the environment with:

```bash
conda activate sedimix
```

### Dependency Check Script

To ensure all required tools, Python packages, and R libraries are installed and correctly configured, we provide a script named `check_dependencies.py`.

#### Usage

1. Save the script `check_dependencies.py` in the root directory of your project.
2. Run the script using Python:
   ```bash
   python check_dependencies.py

#### Index Files
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)  
  We recommend NCBI: nucleotide non-redundant sequences (64GB) for Centrifuge.  
  ```bash
  wget https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz
  tar -xvzf p%2Bh%2Bv.tar.gz -C centrifuge
  rm p%2Bh%2Bv.tar.gz
  ```

- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)  
  The paper was tested out using nt Database version 11/29/2023 (710GB) for Kraken2. A later version should work as well though not tested. Change the code accordingly for the latest version.
  ```bash
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz
  tar -xvzf k2_standard_20240605.tar.gz -C kraken2
  rm k2_standard_20240605.tar.gz
  ```

Alternatively, you can build Centrifuge and Kraken2 indexes yourself by following the instructions provided on their respective GitHub repositories.

#### Human Reference Genome 
Download the human reference genome hg19.fq.gz and build the BWA index. 

## Pipeline Functionality
1. Takes raw FASTQ file(s) as input.
2. Classifies Homo sapiens or Primate reads using Centrifuge or Kraken2.
2. Processes through BWA to further identify reads of interest.
4. Outputs a comprehensive taxonomic report and a folder containing classified hominin reads.

## Usage Instructions
1. Create a new folder for your current run (same level as `scripts` and `rules`)
2. Create a folder named `0_data` within the folder you just created, then place your input FASTQ files in it. 
3. Update the `config.yaml` file to select parameters (see details below) 
2. Run the pipeline with the following command:

   ```bash
   snakemake -s ../rules/snakefile_sedimix --cores {n_cores} --resources mem_mb={max_memory} --jobs 1 --rerun-incomplete 
   ```
An example folder can be found in `example_run/`.

# Path to the reference genome FASTA file.
ref_genome: "/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa"

# Number of threads to use for parallel processing.
threads: 32

# Minimum read length to retain after filtering.
min_length: 35

# Minimum base quality score for reads. Reads below this threshold will be filtered out.
min_quality: 25

# Boolean (True or False) to indicate whether to use a SNP panel for read filtering.
use_snp_panel: False 

# (Optional) Path to an alternative reference genome with additional or modified alleles.
# alt_ref_genome: "/global/scratch/users/jieruixu/sediment_dna/sedimix/reference_data/human_third_allele/modified_hg19.fa"

# (Optional) Path to a BED file defining regions of interest based on the SNP panel.
# snp_panel_bed: "/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/final_pipeline/pendant_2023/probes_reich_n3_b8_CONTROL_BV09-BV09.bed"

# Boolean (True or False) to indicate whether to mask SNPs from the targeted panel regions during processing.
mask_on_target_snps: False  

# Classification tool to use for taxonomic assignment. Options: "centrifuge" or "kraken2".
classification_software: "centrifuge" 

# Path to the directory containing the classification software. The path should not end with a '/'.
classification_software_path: "/global/home/users/jieruixu/jieruixu/sediment_dna/peerj_replication/centrifuge"

# Name of the classification index to use.
classification_index: "nt"

# Maximum memory (in megabytes) allocated to the pipeline.
memory_mb: 200000

# Path to a CSV file containing taxonomic IDs of interest for classification.
taxID: "/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/final_pipeline/temp/primates_taxids.csv"

# Boolean (True or False) to indicate whether to perform additional analysis using mapDamage results.
calculate_from_mapdamage: True

# Path to a file containing sites of interest for lineage-specific analysis.
lineage_sites: "/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/final_pipeline/pendant_2023/probes_reich_n3_b8_CONTROL_BV09-BV09.transformed.txt"

# Analysis type or category for output classification.
types: "hominin_informative"

## Retrieve Your Results
- **Classified Reads**: Located in the `3_final_reads` folder
- **Nonclassified Reads (if specified in congig.yaml)**: Located in the `3_non_classified_reads` folder
- **Data Summary Report**: Located in the `4_final_report` folder, `combined_final_report.tsv` 
- **Deamination Profile**: Located in the `4_mapdamage_results` folder 