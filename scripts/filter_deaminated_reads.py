import argparse
import pysam
from pyfaidx import Fasta
import csv
from scipy.stats import binom

def process_read(read, reference):
    # Extract the sequence of the read and its reference positions
    query_sequence = read.query_sequence
    positions = read.get_reference_positions(full_length=True)
    
    is_deaminated_read = False
    c_to_t_5 = 0
    total_reads_5 = 0
    c_to_t_3 = 0
    total_reads_3 = 0
    
    def is_cpg(positions, index, reference):
        """ Check if the position index is a CpG site in the reference """
        if index < len(positions) - 1 and positions[index] is not None and positions[index + 1] is not None:
            next_base = reference[read.reference_name][positions[index + 1]].seq
            return next_base == 'G'
        return False

    # Check the first base (5' end)
    if len(positions) > 0 and positions[0] is not None and not is_cpg(positions, 0, reference):
        ref_base_5 = reference[read.reference_name][positions[0]].seq
        read_base_5 = query_sequence[0]
        if ref_base_5 == 'C':
            total_reads_5 = 1
            if read_base_5 == 'T':
                c_to_t_5 = 1
                is_deaminated_read = True
        if (read_base_5 == 'A' and ref_base_5 == 'G'):  # already compare C to T one line before
            is_deaminated_read = True

    # Check the last base (3' end), it should be G to A
    if len(positions) > 0 and positions[-1] is not None and not is_cpg(positions, len(positions) - 1, reference):
        ref_base_3 = reference[read.reference_name][positions[-1]].seq
        read_base_3 = query_sequence[-1]
        if ref_base_3 == 'C':
            total_reads_3 = 1
            if read_base_3 == 'T':
                c_to_t_3 = 1
                is_deaminated_read = True
        if (read_base_3 == 'A' and ref_base_3 == 'G'):  # already compare C to T one line before 
            is_deaminated_read = True
    
    if is_deaminated_read: 
        return is_deaminated_read, c_to_t_5, total_reads_5, c_to_t_3, total_reads_3 
        
    # Check the remaining positions for deamination
    check_positions = positions[1:3] + positions[-3:-1] if len(positions) > 2 else []
    check_sequence = query_sequence[1:3] + query_sequence[-3:-1] if len(query_sequence) > 2 else []
    
    for i, (pos, base) in enumerate(zip(check_positions, check_sequence)):
        if pos is None or is_cpg(check_positions, i, reference):
            continue  # Skip positions that do not align to the reference or are CpG sites
        ref_base = reference[read.reference_name][pos].seq
        if (base == 'A' and ref_base == 'G') or (base == 'T' and ref_base == 'C'):
            is_deaminated_read = True
            break 

    return is_deaminated_read, c_to_t_5, total_reads_5, c_to_t_3, total_reads_3

def filter_reads(input_bam, output_deaminated, output_non_deaminated, reference_path, report_file, use_mapdamage):
    # Open the reference genome
    reference = Fasta(reference_path)
    
    c_to_t_5_total = 0
    total_reads_5_total = 0
    c_to_t_3_total = 0
    total_reads_3_total = 0
    
    # Open the input BAM file and create output BAM files
    with pysam.AlignmentFile(input_bam, "rb") as bamfile, \
         pysam.AlignmentFile(output_deaminated, "wb", template=bamfile) as deaminated_bam, \
         pysam.AlignmentFile(output_non_deaminated, "wb", template=bamfile) as non_deaminated_bam:
        
        # Iterate over each read in the BAM file. We only focus on forward reads here 
        for read in bamfile:
            if read.is_unmapped:
                continue  # Skip unmapped reads

            is_deaminated_read, c_to_t_5, total_reads_5, c_to_t_3, total_reads_3 = process_read(read, reference)
            c_to_t_5_total += c_to_t_5
            total_reads_5_total += total_reads_5
            c_to_t_3_total += c_to_t_3
            total_reads_3_total += total_reads_3

            if is_deaminated_read:
                deaminated_bam.write(read)
            else:
                non_deaminated_bam.write(read)
    
    if use_mapdamage:
        sample = report_file.split('_ct_report.csv')[0].split('temp/')[1] 
        mismatches = {
            "5p": {"C": 0, "C>T": 0},
            "3p": {"C": 0, "C>T": 0},
        }
        mapdamage_file = f"4_mapdamage_results/{sample}/misincorporation.txt"
        with open(mapdamage_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if int(row["Pos"]) == 1:
                    counts = mismatches[row["End"]]
                    for key in counts:
                        counts[key] += int(row[key])
        percent_ct_5 = mismatches["5p"]["C>T"] / mismatches["5p"]["C"] * 100 if mismatches["5p"]["C"] > 0 else 0 
        percent_ct_3 = mismatches["3p"]["C>T"] / mismatches["3p"]["C"] * 100 if mismatches["3p"]["C"] > 0 else 0 
        ci_5_lower = binom.ppf(0.025, mismatches["5p"]["C"], percent_ct_5 / 100) / mismatches["5p"]["C"] * 100 if mismatches["5p"]["C"] > 0 else 0
        ci_5_upper = binom.ppf(0.975, mismatches["5p"]["C"], percent_ct_5 / 100) / mismatches["5p"]["C"] * 100 if mismatches["5p"]["C"] > 0 else 0
        ci_3_lower = binom.ppf(0.025, mismatches["3p"]["C"], percent_ct_3 / 100) / mismatches["3p"]["C"] * 100 if mismatches["3p"]["C"] > 0 else 0
        ci_3_upper = binom.ppf(0.975, mismatches["3p"]["C"], percent_ct_3 / 100) / mismatches["3p"]["C"] * 100 if mismatches["3p"]["C"] > 0 else 0
        
    else: 
        percent_ct_5 = (c_to_t_5_total / total_reads_5_total) * 100 if total_reads_5_total > 0 else 0
        percent_ct_3 = (c_to_t_3_total / total_reads_3_total) * 100 if total_reads_3_total > 0 else 0
        ci_5_lower = binom.ppf(0.025, total_reads_5_total, percent_ct_5 / 100) / total_reads_5_total * 100 if total_reads_5_total > 0 else 0
        ci_5_upper = binom.ppf(0.975, total_reads_5_total, percent_ct_5 / 100) / total_reads_5_total * 100 if total_reads_5_total > 0 else 0
        ci_3_lower = binom.ppf(0.025, total_reads_3_total, percent_ct_3 / 100) / total_reads_3_total * 100 if total_reads_3_total > 0 else 0
        ci_3_upper = binom.ppf(0.975, total_reads_3_total, percent_ct_3 / 100) / total_reads_3_total * 100 if total_reads_3_total > 0 else 0

    # Format percentages with confidence intervals
    percent_ct_5_str = f"{percent_ct_5:.1f} ({ci_5_lower:.1f}-{ci_5_upper:.1f})"
    percent_ct_3_str = f"{percent_ct_3:.1f} ({ci_3_lower:.1f}-{ci_3_upper:.1f})"
    
    with open(report_file, 'w', newline='') as csvfile:
        fieldnames = [
            "5' C-to-T substitution frequency [%] (95% CI)",
            "3' C-to-T substitution frequency [%] (95% CI)"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({
            "5' C-to-T substitution frequency [%] (95% CI)": percent_ct_5_str,
            "3' C-to-T substitution frequency [%] (95% CI)": percent_ct_3_str
        })

def main():
    parser = argparse.ArgumentParser(description="Filter BAM reads based on deamination patterns.")
    parser.add_argument("input_bam", help="Path to the input BAM file")
    parser.add_argument("output_deaminated", help="Path to the output BAM file for deaminated reads")
    parser.add_argument("output_non_deaminated", help="Path to the output BAM file for non-deaminated reads")
    parser.add_argument("reference_path", help="Path to the reference genome FASTA file")
    parser.add_argument("report_file", help="Path to the output report CSV file")
    parser.add_argument("--use_mapdamage", help="Enable or disable mapDamage-like behavior", type=bool, default=False)

    args = parser.parse_args()

    filter_reads(args.input_bam, args.output_deaminated, args.output_non_deaminated, args.reference_path, args.report_file, args.use_mapdamage)

if __name__ == "__main__":
    main()