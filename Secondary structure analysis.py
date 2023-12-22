import random
import subprocess
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Get current working directory
os.getcwd()

# Insert path to bpRNA package
bpRNA_path = '/Path/to/bpRNA.pl'

def generate_sequences(template, count=1000):
    return [''.join(random.choice(['A', 'T', 'C', 'G']) if nucleotide == 'N' else nucleotide 
                    for nucleotide in template) for _ in range(count)]

def process_sequences(sequences):
    
    directory = os.getcwd()
    for i, sequence in enumerate(sequences):
        fa_file = f"{i+1}.fa"
        with open(fa_file, 'w') as file:
            file.write(f">{i+1}\n{sequence}\n")

        subprocess.run(f'RNAfold -T 22 -P DNA --noPS {fa_file} > {fa_file}.dbn.full', shell=True)
        os.remove(fa_file)

    subprocess.run(f'for i in *.dbn.full; do '
                   f'base=$(basename $i .dbn.full); '
                   f'head -n 3 $i | sed "s/ (.*//g" > ${{base}}.dbn; done', shell=True)

    subprocess.run(f'for i in *.dbn; do '
                   f'perl {bpRNA_path} $i; done', shell=True)

    subprocess.run(f'for i in *.st; do '
                   f'head -n 6 $i | tail -n 3 > ${{i%.st}}.input.st; done', shell=True)

    nucleotide_list, dot_bracket_list, bpRNA_list = [], [], []
    input_files = [file for file in os.listdir(directory) if file.endswith(".input.st")]

    for file in input_files:
        with open(os.path.join(directory, file), 'r') as f:
            lines = f.readlines()
            nucleotide_list.append(lines[0].strip().replace('U', 'T'))
            dot_bracket_list.append(lines[1].strip())
            bpRNA_list.append(lines[2].strip())

    with open(os.path.join(directory, 'nucleotide_sequences.txt'), 'w') as f:
        f.write('\n'.join(nucleotide_list))
    with open(os.path.join(directory, 'dot_bracket.txt'), 'w') as f:
        f.write('\n'.join(dot_bracket_list))
    with open(os.path.join(directory, 'bpRNA_structures.txt'), 'w') as f:
        f.write('\n'.join(bpRNA_list))

    return nucleotide_list, dot_bracket_list, bpRNA_list

def calculate_avg_structure_frequencies(bpRNA_list):
    # Initialize a dictionary to hold the total counts for each structure type
    structure_counts = {'E': 0, 'S': 0, 'H': 0, 'I': 0, 'B': 0, 'X': 0, 'M': 0, 'K': 0}

    # Sum the counts for each structure type across all sequences
    for seq in bpRNA_list:
        for char in structure_counts:
            structure_counts[char] += seq.count(char)

    # Calculate the average frequency for each structure type
    total_length = sum(len(seq) for seq in bpRNA_list)
    avg_frequencies = {char: count / total_length for char, count in structure_counts.items()}
    
    # Return the average frequencies
    return avg_frequencies

def calculate_avg_contiguous_lengths(bpRNA_list):
    # Initialize a dictionary to hold the total lengths for each structure type
    structure_lengths = {'E': [], 'S': [], 'H': [], 'I': [], 'B': [], 'X': [], 'M': [], 'K': []}

    # Iterate over each sequence in the bpRNA list
    for seq in bpRNA_list:
        # Find all contiguous stretches of each structure type
        for structure_type in structure_lengths.keys():
            matches = re.finditer(f"{structure_type}+", seq)
            for match in matches:
                structure_lengths[structure_type].append(len(match.group(0)))
    
    # Calculate the average length for each structure type
    avg_lengths = {s: (sum(lengths) / len(lengths) if lengths else 0) for s, lengths in structure_lengths.items()}
    
    return avg_lengths

# Export normalized counts and average contiguous lengths
def export_structural_data(bpRNA_list, filename):
    normalized_counts = calculate_avg_structure_frequencies(bpRNA_list)
    contiguous_lengths = calculate_avg_contiguous_lengths(bpRNA_list)
    data = pd.DataFrame({'Type': normalized_counts.keys(),
                         'Normalized_Count': normalized_counts.values(),
                         'Avg_Contiguous_Length': contiguous_lengths.values()})
    data.to_csv(filename, index=False)

# Export heatmap by position data
def export_heatmap_by_position_data(bpRNA_list, filename):
    max_length = max(len(seq) for seq in bpRNA_list)
    heatmap_data = pd.DataFrame(0, index=np.arange(max_length), columns=['E', 'S', 'H', 'I', 'B', 'X', 'M', 'K'])
    for seq in bpRNA_list:
        for position, element in enumerate(seq):
            if element in heatmap_data.columns:
                heatmap_data.loc[position, element] += 1
    heatmap_data.to_csv(filename)
    
    
def shannon_diversity_index(counts):
    """
    Calculate the Shannon Diversity Index for an array of counts.
    
    Parameters:
    counts (np.array): Array of counts for each category.
    
    Returns:
    float: Shannon Diversity Index.
    """
    proportions = counts / np.sum(counts)
    proportions = proportions[proportions > 0]  # Exclude zero counts
    return -np.sum(proportions * np.log(proportions))


def calculate_sdi_for_positions(bpRNA_list, structure_elements):
    max_length = max(len(seq) for seq in bpRNA_list)
    sdi_list = []
    for i in range(max_length):
        # Initialize counts for each structure element
        position_counts = {element: 0 for element in structure_elements}
        for seq in bpRNA_list:
            if i < len(seq):
                element = seq[i]
                if element in position_counts:
                    position_counts[element] += 1
                else:
                    # Handle unexpected elements by adding them to the dictionary
                    position_counts[element] = 1
        counts = np.array(list(position_counts.values()))
        sdi = shannon_diversity_index(counts)
        sdi_list.append(sdi)
    return sdi_list


# Templates
selex_template = "GATGTGAGTGTGTGACGAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCACAGAGAAGAAACAAGACC"
neomer_template = "TGTGTATAAGTCNNGAGGNNNGAATNNNAACGATCGGCGCCAACANNNCATTCNNNCAGANNTCTACTAGTCAC"

# Generate sequences
selex_sequences = generate_sequences(selex_template)
selex_results = process_sequences(selex_sequences)

subprocess.run('rm *fa*', shell=True)
subprocess.run('rm *st*', shell=True)
subprocess.run('rm *st', shell=True)

neomer_sequences = generate_sequences(neomer_template)
neomer_results = process_sequences(neomer_sequences)

subprocess.run('rm *fa*', shell=True)
subprocess.run('rm *st*', shell=True)
subprocess.run('rm *st', shell=True)

# Define the results
selex_bpRNA = selex_results[2]
neomer_bpRNA = neomer_results[2]

# Calculate average contiguous lengths
selex_lengths = calculate_avg_contiguous_lengths(selex_bpRNA)
neomer_lengths = calculate_avg_contiguous_lengths(neomer_bpRNA)

# Exporting data for SELEX and Neomer
export_structural_data(selex_bpRNA, 'selex_structural_data.csv')
export_structural_data(neomer_bpRNA, 'neomer_structural_data.csv')
export_heatmap_by_position_data(selex_bpRNA, 'selex_heatmap_by_position.csv')
export_heatmap_by_position_data(neomer_bpRNA, 'neomer_heatmap_by_position.csv')

print("Data export complete. Ready for R plotting.")

# bpRNA structure strings
structure_elements = 'BEIHMXS'

# Update the structure elements string if there are additional elements
selex_sdi = calculate_sdi_for_positions(selex_bpRNA, structure_elements)
neomer_sdi = calculate_sdi_for_positions(neomer_bpRNA, structure_elements)

# Find the length of the longest list
max_length = max(len(selex_sdi), len(neomer_sdi))

# Create a range for the positions based on the longest list
positions = range(1, max_length + 1)

# Pad the shorter list with NaN values to match the length of the longest list
selex_sdi += [np.nan] * (max_length - len(selex_sdi))
neomer_sdi += [np.nan] * (max_length - len(neomer_sdi))

# Now create the DataFrame
sdi_data = pd.DataFrame({
    'Position': positions,
    'SELEX_SDI': selex_sdi,
    'NEOMER_SDI': neomer_sdi
})

sdi_data.to_csv('./sdi_bpRNA_data.csv', index=False)

print("Shannon Diversity Index data exported for R plotting.")

# Calculate the average frequencies for each dataset
selex_frequencies = calculate_avg_structure_frequencies(selex_bpRNA)
neomer_frequencies = calculate_avg_structure_frequencies(neomer_bpRNA)

# Create dataframes from the frequency dictionaries
selex_df = pd.DataFrame(list(selex_frequencies.items()), columns=['Structure', 'Frequency'])
neomer_df = pd.DataFrame(list(neomer_frequencies.items()), columns=['Structure', 'Frequency'])

# Export the dataframes to CSV files
selex_df.to_csv('./selex_frequencies.csv', index=False)
neomer_df.to_csv('./neomer_frequencies.csv', index=False)

# Example usage with your data
selex_avg_lengths = calculate_avg_contiguous_lengths(selex_bpRNA)
neomer_avg_lengths = calculate_avg_contiguous_lengths(neomer_bpRNA)

# Convert to DataFrames
selex_lengths_df = pd.DataFrame(list(selex_avg_lengths.items()), columns=['Structure', 'AvgLength_Nucleotides'])
neomer_lengths_df = pd.DataFrame(list(neomer_avg_lengths.items()), columns=['Structure', 'AvgLength_Nucleotides'])

# Save to CSV
selex_lengths_df.to_csv('./selex_avg_contiguous_lengths_nucleotides.csv', index=False)
neomer_lengths_df.to_csv('./neomer_avg_contiguous_lengths_nucleotides.csv', index=False)