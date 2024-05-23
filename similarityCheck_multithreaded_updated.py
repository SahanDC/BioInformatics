import pandas as pd
from difflib import SequenceMatcher
from Bio import SeqIO
from multiprocessing import Pool
import os

from Bio import pairwise2
from Bio.Align import substitution_matrices as matlist

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

processed_ids = []
matrix = matlist.load('PAM250')

# Check if the CSV file already exists
if os.path.exists('average_similarity.csv'):
    processed_ids = pd.read_csv('average_similarity.csv')['ID'].tolist()
    sequence_ids = len(processed_ids)

def similarity_check_method_1(score_1, score_2, score_3):
    max1_3 = max(score_1, score_3)
    max_score = score_2/max1_3
    return max_score

def similarity_check_method_2(n, sequence1, sequence2):
    matches = sum(a==b for a,b in zip(sequence1, sequence2))
    return matches/n

def calculate_similarity(args):
    row1, df1 = args

    df_length = len(df1)
    sequence_similarity = 0
    seq_similarity_1 = 0
    seqw_similarity_2 = 0
    results_1 = []
    results_2 = []

    alignments_1 = pairwise2.align.globalds(row1['Sequence'], row1['Sequence'], matrix, -10, -0.5)
    score_1 = alignments_1[0].score

    for _, row2 in df1.iterrows():
        max_length_sequence = max(len(row1['Sequence']), len(row2['Sequence']))
        alignments_2 = pairwise2.align.globalds(row1['Sequence'], row2['Sequence'], matrix, -10, -0.5)
        alignments_3 = pairwise2.align.globalds(row2['Sequence'], row2['Sequence'], matrix, -10, -0.5)

        score_2 = alignments_2[0].score
        score_3 = alignments_3[0].score

        method1_score = similarity_check_method_1(score_1, score_2, score_3)
        method2_score = similarity_check_method_2(max_length_sequence, alignments_2[0][0], alignments_2[0][1])

        seq_similarity_1 += method1_score
        results_1.append({'ID1': row1['ID'], 'ID2': row2['ID'], 'Similarity': method1_score})
        average_similarity_1 = seq_similarity_1 / df_length

        seqw_similarity_2 += method2_score
        results_2.append({'ID1': row1['ID'], 'ID2': row2['ID'], 'Similarity': method2_score})
        average_similarity_2 = seqw_similarity_2 / df_length

    return row1["ID"],average_similarity_1, average_similarity_2, results_1, results_2
    
def calculate_average_similarity_and_return_all(df1, df2):
    total_similarity = 0
    all_results_1 = []
    all_results_2 = []
    with Pool() as pool:
        # for i, (id,average_similarity, results) in enumerate(pool.imap(calculate_similarity, [(row, df1) for _, row in df2.iterrows()])):
        args = [(row, df1) for _, row in df2.iterrows() if row['ID'] not in processed_ids or print(f'Sequence with ID {row["ID"]} has already been processed. Skipping...')]
        for i, (id, average_similarity_1,average_similarity_2, results_1, results_2) in enumerate(pool.imap(calculate_similarity, args)):
            print(f'Comparing sequence with ID {id} : Average similarity of sequence {sequence_ids+i+1} with all sequences in the second dataset: {average_similarity_1} and {average_similarity_2}')
            pd.DataFrame({'ID': [id], 'Average Similarity 1': [average_similarity_1], 'Average Similarity 2': [average_similarity_2]}).to_csv('average_similarity.csv', mode='a', header=False, index=False)
            total_similarity += average_similarity
            all_results_1.extend(results_1)
            all_results_2.extend(results_2)
    results_df_1 = pd.DataFrame(all_results_1)
    results_df_2 = pd.DataFrame(all_results_2)
    return total_similarity / len(df2), results_df_1, results_df_2

if __name__ == '__main__':
    # Load your data
    data1_path = 'Swissprot_Train_Validation_dataset.csv'
    data2_path = 'deeploc_data.fasta'

    # Read the data from the CSV file
    csv_data = pd.read_csv(data1_path)
    # Select 'ACC' and 'Sequence' columns and rename 'ACC' to 'ID'
    df1 = csv_data[['ACC', 'Sequence']].rename(columns={'ACC': 'ID'})

    # Read the data from the FASTA file
    fasta_sequences = SeqIO.parse(data2_path, "fasta")
    df2 = pd.DataFrame([(seq_record.id, str(seq_record.seq)) for seq_record in fasta_sequences], columns=['ID', 'Sequence'])

    # Calculate average similarity and get all results dataframe 1 is compared with the dataframe 2 and dataframe 1s similarity check with the 2
    average_similarity, similarity_results_1, similarity_results_2 = calculate_average_similarity_and_return_all(df1, df2)

    print(f'The overall average similarity is {average_similarity}')

    # Save results to a CSV file
    similarity_results_1.to_csv('similarity_results_!.csv', index=False)
    similarity_results_2.to_csv('similarity_results_2.csv', index=False)
