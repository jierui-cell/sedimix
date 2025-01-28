import random
import pandas as pd
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Specify the input and output files
input_file = '../probes_reich_n3_b8_CONTROL_BV09-BV09.txt.gz.bam_sorted.txt'
tab_file = '../lineage_assignment_sites.Mbuti_Denisova_Altai_Chagyrskaya_Vindija.tab'
hg19_file = '/global/home/users/jieruixu/jieruixu/sediment_dna/sedimix/reference_data/human/hg19.fa'  
output_file = 'modified_hg19.fa'

# Specify data types for the input file
dtypes_probes = {
    'chrom': str,
    'pos0': int,
    'pos1': int,
    'hg19_ref': str,
    'allele1': str,
    'allele2': str,
    'capture_region': str,
    'snp_category_large': str,
    'snp_category': str
}

dtypes_tab = {
    'chrom': str,
    'pos0': int,
    'pos1': int,
    'ref': str,
    'derived': str,
    'category': str
}

# Load the targeted sites from the input file
targeted_sites = pd.read_csv(input_file, sep='\t', dtype=dtypes_probes)

# Load the allele information from the .tab file
allele_info = pd.read_csv(tab_file, sep='\t', header=None, names=['chrom', 'pos0', 'pos1', 'ref', 'derived', 'category'], dtype=dtypes_tab)

# Normalize chromosome names to match hg19 (e.g., '1' to 'chr1')
targeted_sites['chrom'] = 'chr' + targeted_sites['chrom']
allele_info['chrom'] = 'chr' + allele_info['chrom']

# Create a dictionary to map positions to alleles
allele_dict = {(row['chrom'], row['pos1']): (row['ref'], row['derived']) for _, row in allele_info.iterrows()}

def third_allele(ref, alt1, alt2):
    """Generate a third allele that is neither ancestral nor derived."""
    bases = {'A', 'C', 'G', 'T'}
    ref_alt_bases = {ref, alt1, alt2}
    remaining_bases = list(bases - ref_alt_bases)
    return random.choice(remaining_bases) if remaining_bases else 'N'

# Create a dictionary to store the modifications
modifications = {}

# Iterate over the targeted sites and generate the third allele
for _, row in targeted_sites.iterrows():
    chrom = row['chrom']
    pos = row['pos1']  # 1-based

    # Get the reference and derived alleles from the .tab file
    if (chrom, pos) in allele_dict:
        ref, derived = allele_dict[(chrom, pos)]
    else:
        continue  # Skip if the position is not found in the .tab file
    
    new_allele = third_allele(ref, derived, derived)  # Use derived twice because we're considering ref and derived
    
    if chrom not in modifications:
        modifications[chrom] = []
    
    modifications[chrom].append((pos, new_allele))

print("Finish generating dictionary with third allele")

# Load the hg19 reference genome
hg19 = Fasta(hg19_file)

# Open the output file for writing
with open(output_file, 'w') as out_fasta:
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

print(f"Modified hg19 genome saved to {output_file}")
