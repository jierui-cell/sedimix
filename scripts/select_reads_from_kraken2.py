import pandas as pd
import sys
import os

_, taxID_input, kraken2_result, sample = sys.argv
output_file = f"{sample}_reads.lst"

"""
Kraken2 Output format: 
"C"/"U": a one letter code indicating that the sequence was either classified or unclassified.

The sequence ID, obtained from the FASTA/FASTQ header.

The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.

The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".

A space-delimited list indicating the LCA mapping of each k-mer in the sequence(s). 
"""

def extract_reads(taxID_input, kraken2_file, output_file):
    # Determine if taxID_input is a file or a single number
    if os.path.isfile(taxID_input):
        # Read the taxID database
        taxid_df = pd.read_csv(taxID_input)
        taxid_list = taxid_df['TaxID'].tolist()
    else:
        # Single taxID provided
        taxid_list = [int(taxID_input)]
    
    # Define column names based on Kraken2 output structure
    col_names = ['classification', 'readID', 'taxID', 'seqLength', 'kmerMappings']

    # Read kraken2 output into pandas DataFrame, specify only the first four columns to avoid parsing issues
    kraken2_df = pd.read_csv(kraken2_file, sep='\t', names=col_names, usecols=range(4))

    # Convert the 'taxID' column to int for accurate filtering
    kraken2_df['taxID'] = kraken2_df['taxID'].astype(int)

    # Filter reads classified as any of the taxIDs in the list
    filtered_reads = kraken2_df[kraken2_df['taxID'].isin(taxid_list)]['readID'].tolist()
    
    print(f'Number of total reads: {len(kraken2_df)}')
    print(f'Number of reads classified in taxID database: {len(filtered_reads)}')
    
    # Write filtered_reads to output file, each ID on a new line
    with open(output_file, 'w') as file:
        for read_id in filtered_reads:
            file.write(f'{read_id}\n')
            
extract_reads(taxID_input, kraken2_result, output_file)