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

### Centrifuge (default)
```bash
git clone https://github.com/DaehwanKimLab/centrifuge
make -C centrifuge
echo 'export PATH=$PATH:$(pwd)/centrifuge' >> ~/.bashrc
source ~/.bashrc
```

### Kraken2 (optional)
```bash
git clone https://github.com/DerrickWood/kraken2.git
./kraken2/install_kraken2.sh kraken2
echo 'export PATH=$PATH:$(pwd)/kraken2' >> ~/.bashrc
source ~/.bashrc
```

### Seqtk
```bash
git clone https://github.com/lh3/seqtk.git
make -C seqtk
echo 'export PATH=$PATH:$(pwd)/seqtk' >> ~/.bashrc
source ~/.bashrc
```

### BWA 
```bash
git clone https://github.com/lh3/bwa.git
make -C bwa
echo 'export PATH=$PATH:$(pwd)/bwa' >> ~/.bashrc
source ~/.bashrc
```

### Samtools
```bash
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xvjf samtools-1.20.tar.bz2
rm samtools-1.20.tar.bz2
./samtools-1.20/configure
make -C samtools-1.20
echo 'export PATH=$PATH:$(pwd)/samtools-1.20' >> ~/.bashrc
source ~/.bashrc
```

### Bedtools 
```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```

### MapDamage2
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

### Python Packages

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. `pysam` (version >= 0.6): Python package, an interface for reading/writing SAM/BAM files
3. `Cython`: Only necessary if `pysam` needs to be built
4. `numpy` (version >= 1.24.4)
5. `tqdm`
6. `pandas` (version 2.1.2)
7. `pyfaidx` (version 0.8.1.1)
8. `biopython` (version 1.84)
9. `scipy` (version 1.11.3)

**Index Files**  
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)  
  We recommend Refseq: bacteria, archaea, viral, human (7.9GB) for Centrifuge.  
  ```bash
  wget https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz
  tar -xvzf p%2Bh%2Bv.tar.gz -C centrifuge
  rm p%2Bh%2Bv.tar.gz
  ```

- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)  
  We recommend Refseq: archaea, bacteria, viral, plasmid, human1, UniVec_Core (60GB) for Kraken2.  
  ```bash
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
  tar -xvzf k2_standard_20240605.tar.gz -C kraken2
  rm k2_standard_20240605.tar.gz
  ```

Alternatively, you can build Centrifuge and Kraken2 indexes yourself by following the instructions provided on their respective GitHub repositories.

**Human Reference Genome**  
Download the human reference genome hg19.fq.gz and build the BWA index:  
```bash
mkdir hg19
cd hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz
gunzip hg19.fa.gz
../bwa/bwa index hg19.fa
cd ..
```

## Python and Other Dependencies
To ensure all necessary dependencies are installed, create a conda environment using the provided `environment.yaml` file.

### Conda Environment

Create a conda environment with the following command:

```bash
conda env create -f environment.yaml
```

Alternatively, you can use mamba for faster environment creation:

```bash
mamba env create -f environment.yaml
```

Activate the environment with:

```bash
conda activate sedimix
```

## Pipeline Functionality
1. Takes raw FASTQ files as input.
2. Processes these through BWA to identify reads of interest.
3. Classifies Homo sapiens reads using Centrifuge.
4. Outputs a comprehensive taxonomic report and a folder containing classified hominin reads.

## Usage Instructions
1. Place your input FASTQ files in the `data` folder.
2. Run the pipeline with the following command:

   ```bash
   snakemake -s ../rules/centrifuge.smk --cores {n_cores}
   ```

## Retrieve Your Results
- **Classified Reads**: Located in the `final_reads` folder
- **Data Summary Report**: Located in the `final_report` folder
- **Example Folder**: An example folder can be found in `SIM_Set3_centrifuge`.