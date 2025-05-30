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