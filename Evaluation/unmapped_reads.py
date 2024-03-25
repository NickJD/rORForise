import gzip
from collections import defaultdict

# Step 1: Parse BLAST output file to identify predictions with hits
predictions_with_hit = set()
with open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_DP_80', 'r') as blast_file:
    for line in blast_file:
        fields = line.strip().split('\t')
        predictions_with_hit.add(fields[0].replace('@',''))  # Assuming the first column contains prediction IDs

# Step 2: Compare predictions with hits to all predictions to identify those without hits
predictions_all = set()
with open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined.faa', 'r') as predictions_file:
    for line in predictions_file:
        if line.startswith('>'):
            prediction_id = line.strip().replace('>','').replace('@','')
            predictions_all.add(prediction_id)

# Step 3: Extract IDs of predictions without hits
predictions_without_hit = predictions_all - predictions_with_hit

# Step 4: Search for IDs in Eggnog mapper output file and extract COG category
cog_categories = defaultdict(int)
with gzip.open('../Genome_Processing/Staphylococcus_aureus_502A/FragGeneScan/FragGeneScan_ART_errFree_Combined_Emapper.emapper.annotations.gz', 'rt') as eggnog_file:
    for line in eggnog_file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        prediction_id = fields[0][1:]  # Remove '@' from the beginning
        if prediction_id in predictions_without_hit:
            cogs = fields[6]  # Assuming COG category is in the 7th column
            for cog in cogs:
                cog_categories[cog] += 1

# Step 5: Calculate the proportion of each COG category
total_unmapped_predictions = len(predictions_without_hit)
print("Number of Unmapped Predictions: " + str(total_unmapped_predictions))
proportions = {cog: count / total_unmapped_predictions * 100 for cog, count in cog_categories.items()}

# Step 6: Report the proportion of each COG category
for cog, proportion in proportions.items():
    print(f"COG category {cog}: {proportion:.2f}%")













