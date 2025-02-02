import argparse
import pysam
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import random

def parse_args():
    parser = argparse.ArgumentParser(description="Check for alternate alleles in BAM file at specific positions.")
    parser.add_argument("-b", "--bamfile", required=True, help="Input BAM file")
    parser.add_argument("-v", "--variants", required=True, help="Variant positions file")
    parser.add_argument("-t", "--type", required=True, help="Variant type to consider")
    parser.add_argument("-o", "--output", help="Output BAM file for reads not matching the alternate allele")
    parser.add_argument("-p", "--print", action="store_true", help="Activate print mode to display output messages") 
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load the variant positions file into a DataFrame with specified data types
    variants = pd.read_csv(args.variants, sep='\t', header=None, names=['chrom', 'start', 'end', 'ref', 'alt', 'type'], dtype={'chrom': str, 'start': int, 'end': int, 'ref': str, 'alt': str, 'type': str})

    # Filter to include only variants with the specified type
    variants = variants[variants['type'] == args.type]

    # Create a dictionary to store the variant positions and their alternate alleles with progress tracking
    variant_dict = variants.groupby('chrom').apply(lambda df: df.set_index('start')['alt'].to_dict()).to_dict()
    ref_dict = variants.groupby('chrom').apply(lambda df: df.set_index('start')['ref'].to_dict()).to_dict()
    if args.print: 
        print("Finish creating variant dictionary.")
    
    # Open the BAM file
    bamfile = pysam.AlignmentFile(args.bamfile, 'rb')
    
    # check if we need to output reads with non-derived allele (the wrong ones)
    if args.output:
        output_bamfile = pysam.AlignmentFile(args.output, 'wb', header=bamfile.header)
    else:
        output_bamfile = None

    # Initialize counters and storage for reads at each position
    matched_alternate_allele = 0
    reads_not_matching = 0
    total_allele_positions = 0
    reads_matching_alt = 0
    position_reads = defaultdict(list)

    # Iterate through each read in the BAM file with a progress bar
    for read in tqdm(bamfile.fetch(), desc="Processing reads", unit="read"):
        if read.is_unmapped:
            continue

        read_seq = read.query_sequence
        read_pos = read.get_reference_positions()  # This is 0-based
        chrom = bamfile.get_reference_name(read.reference_id).lstrip('chr')  # Remove 'chr' prefix

        # Store read sequences at each position
        for pos in read_pos:
            ref_pos = pos + 1  # Convert 0-based to 1-based
            if ref_pos in variant_dict.get(chrom, {}):
                ref_allele = ref_dict[chrom][ref_pos]
                alt_allele = variant_dict[chrom][ref_pos]
                read_index = read_pos.index(pos)

                # Check if the base at this position matches the reference or alternate allele
                total_allele_positions += 1
                if read_seq[read_index] == alt_allele:
                    matched_alternate_allele += 1
                elif output_bamfile:
                    output_bamfile.write(read)
                    reads_not_matching += 1

    # Calculate the percentage of sites with the alternate allele
    if total_allele_positions > 0:
        percentage_sites = (matched_alternate_allele / total_allele_positions) * 100
        result = f"{percentage_sites:.2f}% ({matched_alternate_allele}/{total_allele_positions})"
        if args.print:
            print(f"Percentage of sites with the alternate allele: {result}")
    else:
        percentage_sites = 0
        result = "0.00% (0/0)"
        if args.print:
            print("No matching positions found.")

    if args.print:
        print(f"Total reads not matching the alternate allele: {reads_not_matching}")
    
    # Close the BAM files
    bamfile.close()
    if output_bamfile:
        output_bamfile.close()
    
    return result

if __name__ == "__main__":
    percentage_result = main()
    print(percentage_result, end='')