import sys
import pysam

def get_snp_positions(snp_panel_file):
    snp_positions = set()
    with open(snp_panel_file) as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            chrom, pos0, pos1 = line.strip().split('\t')[:3]
            chrom = f'chr{chrom}'  # Add 'chr' prefix to match BAM file naming convention
            snp_positions.add((chrom, int(pos0)))
    return snp_positions

def mask_snps_in_bam(input_bam, snp_positions, output_bam):
    with pysam.AlignmentFile(input_bam, 'rb') as bam_in, pysam.AlignmentFile(output_bam, 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            seq = list(read.query_sequence)
            chrom = bam_in.get_reference_name(read.reference_id)
            for pos in range(read.reference_start, read.reference_end):
                seq_pos = pos - read.reference_start
                if (chrom, pos) in snp_positions:
                    if seq_pos < len(seq):
                        seq[seq_pos] = 'N'
                    else:
                        print(f"Index out of range: {seq_pos} for read {read.query_name} with length {len(seq)}")
            read.query_sequence = ''.join(seq)
            bam_out.write(read)

if __name__ == "__main__":
    input_bam = sys.argv[1]
    snp_panel_file = sys.argv[2]
    output_bam = sys.argv[3]

    snp_positions = get_snp_positions(snp_panel_file)
    mask_snps_in_bam(input_bam, snp_positions, output_bam)
