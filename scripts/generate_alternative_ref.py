import argparse
import random
import pandas as pd
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description="Generate modified hg19 genome with third alleles")
    parser.add_argument("input_file", help="Path to the input probes file")
    parser.add_argument("hg19_file", help="Path to the hg19 reference genome file")
    parser.add_argument("output_file", help="Path to the output modified hg19 file")
    return parser.parse_args()

def third_allele(ref, alt1, alt2):
    """Generate a third allele that is neither ancestral nor derived."""
    bases = {'A', 'C', 'G', 'T'}
    ref_alt_bases = {ref, alt1, alt2}
    remaining_bases = list(bases - ref_alt_bases)
    return random.choice(remaining_bases) if remaining_bases else 'N'

def main():
    args = parse_args()

    # Specify data types for the input file
    dtypes_probes = {
        'chrom': str,
        'pos': int,
        'ref': str,
        'a1': str,
        'a2': str
    }

    # Load the targeted sites from the input file
    targeted_sites = pd.read_csv(args.input_file, sep='\t', dtype=dtypes_probes)

    # Normalize chromosome names to match hg19 (e.g., '1' to 'chr1')
    targeted_sites['chrom'] = 'chr' + targeted_sites['chrom'].astype(str)

    # Create a dictionary to store the modifications
    modifications = {}

    # Iterate over the targeted sites and generate the third allele
    for _, row in targeted_sites.iterrows():
        chrom = row['chrom']
        pos = row['pos']  # 1-based
        ref = row['ref']
        a1 = row['a1']
        a2 = row['a2']

        new_allele = third_allele(ref, a1, a2)
        
        if chrom not in modifications:
            modifications[chrom] = []

        modifications[chrom].append((pos, new_allele))

    print("Finish generating dictionary with third allele")

    # Load the hg19 reference genome
    hg19 = Fasta(args.hg19_file)

    # Open the output file for writing
    with open(args.output_file, 'w') as out_fasta:
        # Process each chromosome one by one
        for chrom in hg19.keys():
            seq = str(hg19[chrom])
            mod_seq = list(seq)
            if chrom in modifications:
                for pos, new_allele in modifications[chrom]:
                    mod_seq[pos - 1] = new_allele  # Convert 1-based to 0-based
            mod_seq = ''.join(mod_seq)
            mod_record = SeqRecord(Seq(mod_seq), id=chrom, description='')
            SeqIO.write(mod_record, out_fasta, 'fasta')

    print(f"Modified hg19 genome saved to {args.output_file}")

if __name__ == "__main__":
    main()