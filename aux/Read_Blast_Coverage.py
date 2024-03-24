from collections import defaultdict

def calculate_coverage(blast_file):
    target_ranges = defaultdict(list)

    with open(blast_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines

            fields = line.strip().split('\t')
            sstart, send, slen = map(int, fields[15:18])

            target_id = fields[1]  # Assuming target ID is in the second column

            target_ranges[target_id].append((sstart, send, slen))

    total_coverage_percentage = 0
    total_targets = 0

    # Calculate coverage percentages
    for target_id, ranges in target_ranges.items():
        full_range = combine_ranges(ranges)
        total_coverage = sum(end - start + 1 for start, end in full_range)
        total_length = max(end for _, end in full_range)
        coverage_percentage = (total_coverage / total_length) * 100

        print(f'Target: {target_id}, Coverage: {coverage_percentage:.2f}%')

        total_coverage_percentage += coverage_percentage
        total_targets += 1

    if total_targets > 0:
        average_coverage_percentage = total_coverage_percentage / total_targets
        print(f'Average Coverage: {average_coverage_percentage:.2f}%')
    else:
        print('No targets found.')

def combine_ranges(ranges):
    combined_range = []
    ranges.sort()  # Sort ranges by start position

    for start, end, _ in ranges:
        if not combined_range or start > combined_range[-1][1]:
            combined_range.append((start, end))
        else:
            combined_range[-1] = (combined_range[-1][0], max(end, combined_range[-1][1]))

    return combined_range



if __name__ == "__main__":
    blast_output_file = '../Genome_Processing/Mycoplasma/FrameRate/Mycoplasma_FrameRate_coding_DIAMOND_80_99'
    calculate_coverage(blast_output_file)

    # blast_output_file = '../Genome_Processing/Mycoplasma/FragGeneScan/Myco_ART_errFree_paired_DIAMOND_80_99'
    # calculate_coverage(blast_output_file)