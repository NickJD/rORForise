def visualize_reads_and_cds(cds_range, read_ranges_with_cds, line_length=100, gene_info="CDS Gene 686 - 1828"):
    max_position = max([cds_range[1]] + [end for _, _, _, end in read_ranges_with_cds])
    scale_factor = max_position / line_length

    # Function to scale positions
    def scale(pos):
        scaled_pos = int(pos / scale_factor)
        return min(scaled_pos, line_length - 1)

    # Initialize the visual representation of the CDS
    visual_cds = [' '] * line_length
    cds_start, cds_end = scale(cds_range[0]), scale(cds_range[1])
    for i in range(cds_start, cds_end + 1):
        visual_cds[i] = '-'
    visual_cds[cds_start] = '>'  # Mark CDS start
    visual_cds[cds_end] = '*'  # Mark CDS end
    visual_cds[cds_end] = '*'  # Add '*' at the end of the read CDS

    # Place gene info above the CDS line
    gene_info_visual = [' '] * line_length
    gene_info_start = max((cds_start + cds_end) // 2 - len(gene_info) // 2, 0)
    gene_info_end = min(gene_info_start + len(gene_info), line_length)
    gene_info_visual[gene_info_start:gene_info_end] = list(gene_info)[:line_length - gene_info_start]

    lines = [gene_info_visual, visual_cds]  # Start with the gene info line, followed by the CDS line
    read_info = {}
    read_id = 1

    # Define the function to check for read placement without overlap
    def can_place_read(existing_visual, s_start, s_end):
        return all(c == ' ' for c in existing_visual[s_start:s_end + 1])

    # Process each read
    for read_start, cds_start, cds_end, read_end in read_ranges_with_cds:
        s_read_start, s_cds_start, s_cds_end, s_read_end = map(scale, [read_start, cds_start, cds_end, read_end])

        new_visual = [' '] * line_length
        mid_point_read = (s_read_start + s_read_end) // 2
        for i in range(s_read_start, s_read_end + 1):
            new_visual[i] = '~' if i != mid_point_read else str(read_id)
        for i in range(max(s_read_start, s_cds_start), min(s_read_end, s_cds_end) + 1):
            new_visual[i] = '-' if i != mid_point_read else str(read_id)

        # Add '*' at the end of the read CDS
        new_visual[s_cds_end] = '*'
        new_visual[s_cds_start] = '>'

        read_info[
            read_id] = f"Read {read_id}: (Start: {read_start}, CDS Start: {cds_start}, CDS End: {cds_end}, End: {read_end})"
        read_id += 1

        # Attempt to place the read on an existing line
        placed = False
        for line in lines[2:]:  # Skip the first two lines which are for gene info and CDS
            if can_place_read(line, s_read_start, s_read_end):
                for i in range(s_read_start, s_read_end + 1):
                    line[i] = new_visual[i]
                placed = True
                break

        if not placed:
            lines.append(new_visual)

    # Convert line lists to strings for final output
    final_lines = [''.join(line) for line in lines]

    # Append read information
    final_lines.append("Read Information:")
    for id, info in sorted(read_info.items()):
        final_lines.append(info)

    return '\n'.join(final_lines)


# Example usage
cds_range = (686, 1828)
read_ranges_with_cds = [(650, 686, 798, 798), (700, 725, 845, 845), (850, 850, 1095, 1095), (1700, 1702, 1828, 1845)]

visualization = visualize_reads_and_cds(cds_range, read_ranges_with_cds, gene_info="Gene: 686 - 1828")
print(visualization)
