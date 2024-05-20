import pandas as pd
from difflib import SequenceMatcher
from Bio import SeqIO
from multiprocessing import Pool

def calculate_similarity(args):
    row1, df1 = args
    sequence_similarity = 0
    results = []
    for _, row2 in df1.iterrows():
        matcher = SequenceMatcher(None, row1['Sequence'], row2['Sequence'])
        similarity = matcher.ratio()
        sequence_similarity += similarity
        results.append({'ID1': row1['ID'], 'ID2': row2['ID'], 'Similarity': similarity})
    average_similarity = sequence_similarity / len(df1)
    return row1["ID"],average_similarity, results

def calculate_average_similarity_and_return_all(df1, df2):
    total_similarity = 0
    all_results = []
    with Pool() as pool:
        for i, (id,average_similarity, results) in enumerate(pool.imap(calculate_similarity, [(row, df1) for _, row in df2.iterrows()])):
            print(f'Comparing sequence with ID {id} : Average similarity of sequence {i+1} with all sequences in the second dataset: {average_similarity}')
            total_similarity += average_similarity
            all_results.extend(results)
    results_df = pd.DataFrame(all_results)
    return total_similarity / len(df2), results_df

if __name__ == '__main__':
    # Load your data
    data1_path = './data_files/TRAINING VALIDATION DEEPLOCK 2.0/Swissprot_Train_Validation_dataset.csv'
    data2_path = './data_files\deelloc1.0\deeploc_data.fasta'

    # Read the data from the CSV file
    csv_data = pd.read_csv(data1_path)
    # Select 'ACC' and 'Sequence' columns and rename 'ACC' to 'ID'
    df1 = csv_data[['ACC', 'Sequence']].rename(columns={'ACC': 'ID'})

    # Read the data from the FASTA file
    fasta_sequences = SeqIO.parse(data2_path, "fasta")
    df2 = pd.DataFrame([(seq_record.id, str(seq_record.seq)) for seq_record in fasta_sequences], columns=['ID', 'Sequence'])

    # Calculate average similarity and get all results dataframe 1 is compared with the dataframe 2 and dataframe 1s similarity check with the 2
    average_similarity, similarity_results = calculate_average_similarity_and_return_all(df1, df2)

    print(f'The overall average similarity is {average_similarity}')

    # Save results to a CSV file
    similarity_results.to_csv('similarity_results.csv', index=False)
