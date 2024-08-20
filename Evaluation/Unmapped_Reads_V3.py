import gzip
from collections import defaultdict

# Step 1: Parse BLAST output file to identify predictions with hits
predictions_with_hit = set()
with open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_DP_80', 'r') as blast_file:
    for line in blast_file:
        fields = line.strip().split('\t')
        predictions_with_hit.add(fields[0].replace('@', ''))  # Assuming the first column contains prediction IDs

# Step 2: Compare predictions with hits to all predictions to identify those without hits
predictions_all = defaultdict(str)
with open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined.faa', 'r') as predictions_file:
    sequence = ''
    for line in predictions_file:
        if line.startswith('>'):
            if sequence:  # If there's a sequence accumulated, add it to the dictionary
                predictions_all[prediction_id] = sequence
                sequence = ''  # Reset sequence
            prediction_id = line.strip().replace('>', '').replace('@', '')
        else:
            sequence += line.strip()  # Accumulate the sequence

    # Add the last sequence since there's no new header after it
    if sequence:
        predictions_all[prediction_id] = sequence

# Step 3: Extract IDs of predictions without hits
predictions_without_hit = {prediction_id: sequence for prediction_id, sequence in predictions_all.items() if
                           prediction_id not in predictions_with_hit}
count_unmapped = 0
# Step 4: Search for IDs in Eggnog mapper output file and extract COG category for unmapped reads
cog_categories_unmapped = defaultdict(int)
with gzip.open(
        '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations.gz',
        'rt') as eggnog_file:
    for line in eggnog_file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        prediction_id = fields[0][1:]  # Remove '@' from the beginning
        if prediction_id in predictions_without_hit:
            cogs = fields[6]  # Assuming COG category is in the 7th column
            for cog in cogs:
                if cog != 'Z' and cog != '-':
                    cog_categories_unmapped[cog] += 1

                count_unmapped += 1

# Step 5: Calculate the proportion of each COG category for unmapped reads
total_unmapped_predictions = sum(cog_categories_unmapped.values())
print("Number of Unmapped Predictions: " + str(total_unmapped_predictions))
unmapped_proportions = {cog: count / total_unmapped_predictions * 100 for cog, count in cog_categories_unmapped.items()}

# Step 6: Report the proportion of each COG category for unmapped reads
print("Proportion of COG categories for unmapped reads:")
for cog, proportion in unmapped_proportions.items():
    print(f"COG category {cog}: {proportion:.2f}%")

# Step 7: Collect COG assignments for mapped reads
cog_categories_mapped = defaultdict(int)
with gzip.open(
        '../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations.gz',
        'rt') as eggnog_file:
    for line in eggnog_file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        prediction_id = fields[0][1:]  # Remove '@' from the beginning
        if prediction_id in predictions_with_hit:
            cogs = fields[6]  # Assuming COG category is in the 7th column
            for cog in cogs:
                if cog != 'Z' and cog != '-':
                    cog_categories_mapped[cog] += 1

# Step 8: Calculate the proportion of each COG category for mapped reads
total_mapped_predictions = sum(cog_categories_mapped.values())
mapped_proportions = {cog: count / total_mapped_predictions * 100 for cog, count in cog_categories_mapped.items()}

# Report the proportion of mapped and unmapped reads with a COG assignment
total_mapped_reads = total_mapped_predictions + total_unmapped_predictions
proportion_mapped = total_mapped_predictions / len(predictions_with_hit) * 100
proportion_unmapped = total_unmapped_predictions / len(predictions_without_hit)  * 100

print(f"Proportion of mapped reads with a COG assignment: {proportion_mapped:.2f}%")
print(f"Proportion of unmapped reads with a COG assignment: {proportion_unmapped:.2f}%")
