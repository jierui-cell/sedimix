import argparse
import pysam

def load_snp_positions(snp_file):
    snp_positions = {}
    with open(snp_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split()
            chrom = 'chr' + parts[0]
            pos = int(parts[1])  # 1-based position
            alleles = (parts[3], parts[4])  # expected allele states
            snp_positions.setdefault(chrom, {})[pos] = alleles
    return snp_positions

def filter_reads(bam_file, snp_positions, output_bam, filter_allele_mismatch):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    out_bam = pysam.AlignmentFile(output_bam, 'wb', header=bam.header)
    
    for read in bam.fetch():
        chrom = bam.get_reference_name(read.reference_id)
        if chrom in snp_positions:
            match_found = False
            match_all_positions = True
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                ref_pos_1_based = ref_pos + 1
                if ref_pos_1_based in snp_positions[chrom]:
                    match_found = True
                    alleles = snp_positions[chrom][ref_pos_1_based]
                    read_base = read.query_sequence[query_pos]
                    if filter_allele_mismatch and read_base not in alleles:
                        match_all_positions = False
                        break 
            if not match_all_positions: 
                print ('something went wrong')
            if match_found and match_all_positions:
                out_bam.write(read)
    
    bam.close()
    out_bam.close()

def main():
    parser = argparse.ArgumentParser(description="Filter reads that overlap with SNP positions")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("snp_file", help="SNP positions file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("--filter_allele_mismatch", action='store_true', help="Filter out reads that do not match the expected allele states")
    args = parser.parse_args()

    snp_positions = load_snp_positions(args.snp_file)
    filter_reads(args.bam_file, snp_positions, args.output_bam, args.filter_allele_mismatch)

if __name__ == "__main__":
    main()
