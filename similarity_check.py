import pandas as pd
import numpy as np
from Bio import pairwise2, SeqIO
from Bio.Align import substitution_matrices as matlist
import os
import warnings

warnings.filterwarnings("ignore")

matrix = matlist.load('PAM250')

def similarity_check_method_1(score_1, score_2, score_3):
    max1_3 = max(score_1, score_3)
    max_score = score_2/max1_3
    return max_score

def similarity_check_method_2(n, sequence1, sequence2):
    matches = sum(a==b for a,b in zip(sequence1, sequence2))
    return matches/n

def average_similarity(df1, df2):
    df1_length = len(df1)
    results_1, results_2 = [], []
    average_similarity_1, average_similarity_2 = 0, 0

    for i, sequence1 in df1.iterrows():
        sequence_ID = sequence1['ID']
        protein_sequence_1 = sequence1['Sequence']
        alignments_1 = pairwise2.align.globalds(protein_sequence_1, protein_sequence_1, matrix, -10, -0.5)
        score_1 = alignments_1[0].score

        sequence_similarity_1, sequence_similarity_2 = 0, 0

        for j, sequence2 in df2.iterrows():
                protein_sequence_2 = sequence2['Sequence']
                max_length_sequence = max(len(protein_sequence_1), len(protein_sequence_2))

                alignments_2 = pairwise2.align.globalds(protein_sequence_1, protein_sequence_2, matrix, -10, -0.5)
                alignments_3 = pairwise2.align.globalds(protein_sequence_2, protein_sequence_2, matrix, -10, -0.5)

                score_2 = alignments_2[0].score
                score_3 = alignments_3[0].score

                method1_score = similarity_check_method_1(score_1, score_2, score_3)
                method2_score = similarity_check_method_2(max_length_sequence, alignments_2[0][0], alignments_2[0][1])

                sequence_similarity_1 += method1_score
                results_1.append({'ID1': sequence1['ID'], 'ID2': sequence2['ID'], 'Similarity_1': method1_score})
                average_similarity_1 = sequence_similarity_1 / df1_length

                sequence_similarity_2 += method2_score
                results_2.append({'ID1': sequence1['ID'], 'ID2': sequence2['ID'], 'Similarity_2': method2_score})
                average_similarity_2 = sequence_similarity_2 / df1_length

                pd.DataFrame({
                    'ID': [sequence_ID],
                    'Average_Similarity_1': [average_similarity_1],
                    'Average_Similarity_2': [average_similarity_2]
                }).to_csv('average_similarity.csv', index=False, mode='a', header=False)    

    results_1_df = pd.DataFrame(results_1)
    results_2_df = pd.DataFrame(results_2)
    return sequence_ID, average_similarity_1, average_similarity_2, results_1_df, results_2_df    

# Load your data
data1_path = 'Swissprot_Train_Validation_dataset.csv'
data2_path = 'deeploc_data.fasta'

# Read the data from the CSV file
csv_data = pd.read_csv(data1_path)
# Select 'ACC' and 'Sequence' columns and rename 'ACC' to 'ID'
df1 = csv_data[['ACC', 'Sequence']].rename(columns={'ACC': 'ID'})

# Now, 'selected_columns' is a DataFrame with 'ID' and 'Sequence' columns.

# Read the data from the FASTA file
fasta_sequences = SeqIO.parse(data2_path, "fasta")

df2 = pd.DataFrame([(seq_record.id, str(seq_record.seq)) for seq_record in fasta_sequences], columns=['ID', 'Sequence'])

# Calculate average similarity and get all results dataframe 1 is compared with the dataframe 2 and dataframe 1s similarity check with the 2

sequence_ID, average_similarity_1, average_similarity_2, results_1, results_2  = average_similarity(df1, df2)

print(f'The overall average similarity is {average_similarity_2} and {average_similarity_2}')

# Save results to a CSV file
results_1.to_csv('similarity_results_1.csv', index=False)
results_2.to_csv('similarity_results_2.csv', index=False)