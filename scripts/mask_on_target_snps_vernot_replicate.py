import sys
import pysam

def get_snp_positions(snp_panel_file):
    snp_positions = {}
    with open(snp_panel_file) as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            chrom, pos1, *_ = line.strip().split('\t')
            chrom = f'chr{chrom}'  # Add 'chr' prefix to match BAM file naming convention
            pos1 = int(pos1)
            if chrom not in snp_positions:
                snp_positions[chrom] = set()
            snp_positions[chrom].add(pos1)
    return snp_positions

def mask_snps_in_bam(input_bam, snp_positions, output_bam):
    with pysam.AlignmentFile(input_bam, 'rb') as bam_in, pysam.AlignmentFile(output_bam, 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            seq = list(read.query_sequence)
            chrom = bam_in.get_reference_name(read.reference_id)
            if chrom in snp_positions:
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    ref_pos_1_based = ref_pos + 1
                    if ref_pos_1_based in snp_positions[chrom]:
                        if query_pos < len(seq):
                            seq[query_pos] = 'N'
                        else:
                            print(f"Index out of range: {query_pos} for read {read.query_name} with length {len(seq)}")
            read.query_sequence = ''.join(seq)
            bam_out.write(read)

if __name__ == "__main__":
    input_bam = sys.argv[1]
    snp_panel_file = sys.argv[2]
    output_bam = sys.argv[3]

    snp_positions = get_snp_positions(snp_panel_file)
    mask_snps_in_bam(input_bam, snp_positions, output_bam)