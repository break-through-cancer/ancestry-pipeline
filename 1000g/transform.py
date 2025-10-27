import pandas as pd

# Load your 1000 Genomes panel file
panel_file = "1000g/integrated_call_samples_v3.20130502.ALL.panel"
df = pd.read_csv(panel_file, sep='\t')

# Select only the columns needed for RFMix: sample and population
sample_map = df[['sample', 'pop']]

# Save as tab-delimited file
sample_map_file = "rfmix_sample_map.txt"
sample_map.to_csv(sample_map_file, sep='\t', index=False, header=False)

print(f"Sample map saved to {sample_map_file}")
