## Explanation of `config.yaml` Parameters:
 
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
   - Path to a CSV file containing taxonomic IDs of interest for classification. A curated list can be found in the example_run folder.
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
   - This file should be a tab-delimited text file with the following columns **but no header**:
     ```
     Chromosome   Start   End   Reference   Alternate   Type
     ```
   - File example:
     ```
     1       949200  949200  C       G       hominin_informative
     1       1500380 1500380 G       C       hominin_informative
     1       1500941 1500941 G       A       neanderthal
     ```
   - The `Type` column defines the classification of each site. 
   - The file used in this paper is [hominin_informative_sites.txt](../hominin_informative_sites.txt). Original file can be downloaded from https://datadryad.org/dataset/doi:10.5061/dryad.41ns1rnj1, subsetting to hominin_informative only. 

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