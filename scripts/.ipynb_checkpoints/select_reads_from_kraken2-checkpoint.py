import pandas as pd 
import sys 

_, taxID, kraken2_result, sample = sys.argv 
taxID = int(taxID) 
output_file = f"{sample}_reads.lst"

"""
Kraken2 Output format: 
"C"/"U": a one letter code indicating that the sequence was either classified or unclassified.

The sequence ID, obtained from the FASTA/FASTQ header.

The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.

The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".

A space-delimited list indicating the LCA mapping of each k-mer in the sequence(s). 
""" 

def extract_homo_sapiens_reads(taxID, kraken2_file, output_file):
    # Define column names based on Kraken2 output structure
    col_names = ['classification', 'readID', 'taxID', 'seqLength', 'kmerMappings']

    # Read kraken2 output into pandas DataFrame, specify only the first four columns to avoid parsing issues
    kraken2_df = pd.read_csv(kraken2_file, sep='\t', names=col_names, usecols=range(4))
    
    # Convert the 'taxID' column to int for accurate filtering
    kraken2_df['taxID'] = kraken2_df['taxID'].astype(int)

    # Filter reads classified as Homo sapiens
    homo_sapiens_reads = kraken2_df[kraken2_df['taxID'] == taxID]['readID'].tolist()
    
    print(f'Number of total reads: {len(kraken2_df)}')
    print(f'Number of homo sapien reads classified: {len(homo_sapiens_reads)}')
    
    # Write homo_sapiens_reads to output file, each ID on a new line
    with open(output_file, 'w') as file:
        for read_id in homo_sapiens_reads:
            file.write(f'{read_id}\n')
            
extract_homo_sapiens_reads(taxID, kraken2_result, output_file) 
