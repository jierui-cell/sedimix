# *sedimix*

*sedimix*: A workflow for the analysis of hominin nuclear DNA sequences from sediments

## Overview
Here we present an open-source snakemake workflow that identifies human sequences from sequencing data and provides relevant summary statistics. The final tool prioritizes the retention of human DNA while minimizing detection errors, offering a robust and accessible solution to support the growing needs of human evolutionary research. See [paper](https://www.biorxiv.org/content/10.1101/2025.02.28.640818v1) for details.

**Key Features:**
- Utilizes Snakemake, a workflow management system for Python
- Identifies ancient hominin reads in the samples through taxonomic classification followed by mapping 
- Generates a report file with summary statistics 

## Requirements
Ensure the following tools are installed and configured before running the pipeline:

If you're using `zsh` as your shell, replace `~/.bashrc` with `~/.zshrc` in the commands below.

### [Centrifuge (default)](https://ccb.jhu.edu/software/centrifuge/manual.shtml)
```bash
git clone https://github.com/DaehwanKimLab/centrifuge
make -C centrifuge
echo 'export PATH=$PATH:$(pwd)/centrifuge' >> ~/.bashrc
source ~/.bashrc
```

### [Kraken2 (optional)](https://github.com/DerrickWood/kraken2/wiki/Manual)
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

### [MapDamage2](https://ginolhac.github.io/mapDamage/)
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

### Python Packages

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

### Conda Environment

Create a conda environment with the following command:

```bash
conda env create -f environment.yaml
```

Activate the environment with:

```bash
conda activate sedimix
```

Note that r-gam might show errors for the latest mac users that are on arm64 architecture. If you encounter this issue, you can remove the line for r-gam from the environment.yaml file and install it manually in your R session. 

### Dependency Check Script

To ensure all required tools, Python packages, and R libraries are installed and correctly configured, we provide a script named `check_dependencies.py`.

1. Save the script `check_dependencies.py` in the root directory of your project.
2. Run the script using Python:
   ```bash
   python check_dependencies.py
   ```

### Index Files
Download index files for Centrifuge and Kraken2 from the following:
- [AWS Indexes for Centrifuge](https://benlangmead.github.io/aws-indexes/centrifuge)  
  We recommend NCBI: nucleotide non-redundant sequences (64GB) for Centrifuge.  
  ```bash
  wget https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
  tar -xvzf nt_2018_3_3.tar.gz -C centrifuge
  rm nt_2018_3_3.tar.gz
  ```

- [AWS Indexes for Kraken2](https://benlangmead.github.io/aws-indexes/k2)  
  The paper was tested out using nt Database version 11/29/2023 (710GB) for Kraken2. A later version should work as well though not tested. Change the code accordingly for the latest version.
  ```bash
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_nt_20231129.tar.gz
  tar -xvzf k2_nt_20231129.tar.gz -C kraken2
  rm k2_nt_20231129.tar.gz
  ```

Alternatively, you can build Centrifuge and Kraken2 indexes yourself by following the instructions provided on their respective GitHub repositories.

### Human Reference Genome 
Download the human reference genome (e.g. hg19.fq.gz) and build the BWA index. 

If you have a specified SNP panel, you can generate an alternative reference genome to minimize reference bias during sequence mapping.  
Use the script `scripts/generate_alternative_ref.py` with three arguments:  

1. The SNP panel/probe file (e.g. containing positions and alleles to modify)  
2. Your original reference genome (FASTA format)  
3. A name for your output modified reference genome file  

**SNP Panel/Probe File Format:**  
- Tab-delimited text file with **no header**.  
- Must contain **five columns** in the order: `chrom`, `pos`, `ref`, `a1`, `a2`, for example:  
   ```
   1 10000 A G G
   1 20000 T C C
   2 30000 C T T
   ```
- `chrom`: Chromosome number or letter (e.g., `1`, `2`, ..., `X`, `Y`). The script automatically adds `chr` to this value.  
- `pos`: 1-based genomic position.  
- `ref`: The reference allele at this position.  
- `a1` and `a2`: Observed alternate alleles at this position.  

The script replaces bases in the reference genome at specified SNP positions with a "third allele," ensuring it differs from both the original reference and provided SNP alleles. This helps reduce reference bias when mapping ancient or metagenomic reads.

**Example usage:**
```bash
python scripts/generate_alternative_ref.py \
    <snp_panel_file> \
    <hg19.fasta> \
    <modified_hg19.fasta>

# Build a BWA index on the newly created reference
bwa index <modified_hg19.fasta>
```

## Usage Instructions
After you have run `check_dependencies.py` and it shows no error, downloaded the pre-built index, and built the alternative humen reference (optional, if a SNP panel is specified), follow the steps below to run Sedimix. 
 
1. Create a new folder for your current run (same level as `scripts` and `rules`)
2. Create a folder named `0_data` within the folder you just created, then place your input FASTQ files in it. **Input files must be in format that is either .fq or .fq.gz.** 
3. Update the `config.yaml` file to select parameters (see details below) 
2. Run the pipeline with the following command:

   ```bash
   snakemake -s ../rules/snakefile_sedimix --cores {n_cores} --resources mem_mb={max_memory} --jobs 1 --rerun-incomplete 
   ```
An example folder can be found in `example_run/`. The start folder is only 0_data, and start file is config.yaml and a bash script run_pipeline.sh to submit the snakemake command. All the other folders and files are generated by *sedimix*. 

### Explanation of `config.yaml` Parameters:
 
1. **memory_mb**:
   - Maximum memory (in megabytes) allocated to the pipeline. This should align with the max_memory parameter in the snakemake command. 
   - Example: 200000 (around 200 GB)

2. **threads**:
   - Number of threads to use for parallel processing. This should align with the n_cores parameter in the snakemake command.
   - Example: 32

3. **min_length**:
   - Minimum read length to retain after filtering.
   - Default: 30

4. **min_quality**:
   - Minimum base quality score after mapping for reads. Reads below this threshold will be filtered out.
   - Default: 25

5. **classification_software**:
   - Classification tool to use for taxonomic assignment. 
   - Options: "centrifuge" or "kraken2".
   - Default: "centrifuge"

6. **classification_index**:
   - Folder path and name of the classification index to use.
   - Example: "/path/to/index/nt" (after untaring for Centrifuge, nt is the index filename prefix (minus trailing .X.cf); for Kraken2, there is no need for this extra reference, and the path to the untar index folder is suffix)

7. **taxID**:
   - Path to a CSV file containing taxonomic IDs of interest for classification.
   - Example: "/path/to/primates_taxids.csv"

8. **use_snp_panel**:
   - Boolean (True or False) to indicate whether to use a SNP panel for read filtering.
   - Default: False

9. **ref_genome**:
   - Path to the reference genome FASTA file (or to the pre-built alternative reference genome if use_snp_panel is set to True.). 
   - Example: "/path/to/reference.fa"

10. **snp_panel_bed** (Optional):
   - Path to a BED file defining regions of interest based on the SNP panel.
   - Only required if use_snp_panel is set to True. Commented out by default. 

11. **calculate_from_mapdamage**:
   - Boolean (True or False) to indicate whether to perform deamination analysis using mapDamage2's results.
   - Default: True

12. **lineage_sites**:
   - Path to a file containing sites of interest for lineage-specific analysis.
   - This file should be a tab-delimited text file with the following columns with headers:
     ```
     Chromosome   Start   End   Reference   Alternate   Type
     ```
   - Example:
     ```
     1       949200  949200  C       G       hominin_informative
     1       1500380 1500380 G       C       hominin_informative
     1       1500941 1500941 G       A       neanderthal
     ```
   - The `Type` column defines the classification of each site. 

13. **types**:
   - One or more site types to be analyzed, corresponding to the `Type` column in the `lineage_sites` file.
   - Multiple types can be specified as a space-separated string.
   - Example: "hominin_informative", OR "hominin_informative neanderthal denisova"

14. **to_clean**:
   - Boolean (True or False) indicating whether to remove all intermediate files once the final output is generated.
   - Default: False

15. **keep_non_hominin_reads**:
   - Boolean (True or False) determining whether to save reads classified as non-hominin into a separate FASTQ file.
   - Default: False

## Retrieve Your Results
- **Classified hominin reads**: Located in the `3_final_reads` folder, ending with `{sample_name}_final.bam`. 
- **Classified hominin reads that have deamination**: Located in the `3_final_reads` folder, ending with `{sample_name}_final_deaminated.bam`.
- **Classified hominin reads that do not have deamination**: Located in the `3_final_reads` folder, ending with `{sample_name}_final_non_deaminated.bam`.
- **Classified Non-hominin Reads (if specified in config.yaml)**: Located in the `3_final_reads` folder, ending with `{sample_name}_non_hominin.fq`. 
- **Data Summary Report**: Located in the `4_final_report` folder. `combined_final_report.tsv` contains results for all samples.
- **Deamination Profile**: Located in the `4_mapdamage_results` folder. 