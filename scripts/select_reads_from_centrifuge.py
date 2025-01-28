import pandas as pd
import sys
import os

_, taxID_input, centrifuge_result, sample = sys.argv
output_file = f"{sample}_reads.lst"

def extract_reads(taxID_input, centrifuge_file, output_file):
    # Determine if taxID_input is a file or a single number
    if os.path.isfile(taxID_input):
        # Read the taxID database
        taxid_df = pd.read_csv(taxID_input)
        taxid_list = taxid_df['TaxID'].tolist()
    else:
        # Single taxID provided
        taxid_list = [int(taxID_input)]
    
    # Read centrifuge output into pandas DataFrame
    centrifuge_df = pd.read_csv(centrifuge_file, sep='\t', header=0)

    # Filter reads classified as any of the taxIDs in the list
    filtered_reads = set(centrifuge_df[centrifuge_df['taxID'].isin(taxid_list)]['readID'].tolist())
    
    print(f'Number of total reads: {len(centrifuge_df)}')
    print(f'Number of reads classified in taxID database: {len(filtered_reads)}')
    
    # Prepend \n to output 
    filtered_reads_output = [f'{read_id}\n' for read_id in filtered_reads]

    # Write filtered_reads to output file
    with open(output_file, 'w') as file:
        file.writelines(filtered_reads_output)

extract_reads(taxID_input, centrifuge_result, output_file)