import pandas as pd
import numpy as np
import os
import re
import shutil
from datetime import date

NUM_DETECTORS = 4  # MS2, N, ORF1ab, S
CT_CUTOFF = 30  # Applied to N, ORF1ab only
TRANSFER_uL = 50

# Grab Ct values from PCR CSVs and put them in a dataframe
def collect_cts(pcr_results):
    pcr_dict = {}
    for filename in pcr_results:
        header_row = find_header(filename)
        batch = re.findall(r'\D(\d{5})\D', ' '+ filename +' ')[0]
        with open(filename, errors='ignore') as f:
            df = pd.read_csv(f, skiprows=header_row)
            df = clean_pcr_data(df)
            pcr_dict[batch] = df
        os.remove(filename)
    pcr_dict = dict(sorted(pcr_dict.items()))
    return pcr_dict

# Find number of rows to Ct value data
def find_header(filename):
    with open(filename, errors='ignore') as f:
        for row_num, line in enumerate(f):
            if line.startswith('Well'):
                return row_num

def clean_pcr_data(df):
    df = df.iloc[:, :5]  # First 5 columns contain necessary Ct data

    # Standardize column names and cell values between ABI and QuantStudio5
    df = df.rename(columns={'Well Position': 'Well', 
                            'Target Name': 'Detector', 
                            'CT': 'Ct'})
    df = df.drop(columns=['Task'])
    df = df.replace(to_replace='UNKNOWN', value='Unknown')
    df = df.replace(to_replace='Undetermined', value=0)
    df['Ct'] = pd.to_numeric(df['Ct'])
    return df

# Change sample order from A1, A2, A3... to A1, B1, C1...
def order_by_column(pcr_dict):
    for df in pcr_dict.values():
        df['Num'] = df.Well.str[1:]
        df['Num'] = pd.to_numeric(df['Num'])
        df['Alpha'] = df.Well.str[0]
        df.sort_values(by=['Num', 'Alpha'], inplace=True)
        df.drop(columns=['Num', 'Alpha'], inplace=True)
    return pcr_dict

# Format table as one row per sample and add columns for Ct values for each gene
def pivot_by_detector(pcr_dict):
    for batch, df in pcr_dict.items():
        df['Position'] = [i for i in range(1, df.shape[0]//NUM_DETECTORS+1) 
                          for _ in range(NUM_DETECTORS)]
        df = df.pivot(index=['Position', 'Well', 'Sample Name'], 
                      columns='Detector', values='Ct')
        df.reset_index(inplace=True)
        df = df.rename_axis(None, axis='columns')
        pcr_dict[batch] = df
    return pcr_dict

# Filter out samples with low viral detection
def filter_by_ct(pcr_dict):
    passed_dict = {}
    unsat_dict = {}
    for batch, df in pcr_dict.items():
        df["Pass"] = ((df['Sample Name'] != 'NC') 
                     & (df['Sample Name'] != 'PC')
                     & (df['N gene'] != 0)
                     & (df['N gene'] < CT_CUTOFF)
                     & (df['ORF1ab'] != 0)
                     & (df['ORF1ab'] < CT_CUTOFF))
        passed_dict[batch] = df[df['Pass'] == True]
        unsat_dict[batch] = df[df['Pass'] == False]
    return passed_dict, unsat_dict

# Read BatchSmash file into a dataframe and clean
def read_batch_smash(batch_smash_csv):
    batch_smash = pd.read_csv(batch_smash_csv)
    os.remove(batch_smash_csv)
    batch_smash = batch_smash.dropna(how='all')
    int_columns = ['Destination Location', 'Sample Location', 'Sample Name', 'Batch']
    batch_smash[int_columns] = batch_smash[int_columns].astype(int)
    batch_smash = batch_smash.reset_index(drop=True)
    return batch_smash

# Create a new BatchSmash dataframe based on the old BatchSmash
def create_batch_smash(old_batch_smash, passed_dict):
    num_to_transfer = old_batch_smash['Destination Location'].iloc[-1] - 1
    if num_to_transfer < 95:  # Don't transfer already complete plate
        new_batch_smash = old_batch_smash.tail(num_to_transfer)
    first_covseq = int(new_batch_smash['COVSEQ'].iloc[-1][-4:])
                                                         # last 4 digits is id
    # Append the new samples' data to the BatchSmash
    for batch, df in passed_dict.items():
        smash_df = pd.DataFrame(columns=['COVSEQ', 'Destination Location',
                                         'Sample Location', 'Sample Name', 'Batch'])
        smash_df['Sample Location'] = df['Position']
        smash_df['Sample Name'] = df['Sample Name']
        smash_df['Batch'] = int(batch)
        smash_df = smash_df.reset_index(drop=True)
        start_destination = new_batch_smash['Destination Location'].iloc[-1]
        smash_df['Destination Location'] = ((smash_df.index - 1 + start_destination) % 94) + 2
        new_batch_smash = pd.concat([new_batch_smash, smash_df])
    
    # Add COVSEQ ID column to the BatchSmash
    num_covseqs = len(new_batch_smash.index) // 94 + 1
    covseq_ids = ['COVSEQ_' + str(id) 
                  for id in range(first_covseq, first_covseq + num_covseqs)]
    covseq_repeat = np.repeat(covseq_ids, 94).tolist()
    covseq_column = covseq_repeat[:len(new_batch_smash.index)]
    new_batch_smash['COVSEQ'] = covseq_column

    return new_batch_smash, covseq_ids

def write_batch_smash(batch_smash):
    shutil.rmtree('Output', ignore_errors=True)
    os.mkdir('Output')
    new_batch_smash.to_csv(f'Output/BatchSmash {date.today()}.csv', index=False)

# Use new the BatchSmash file to create worklisting files used by the Tecan
def create_worklist(batch_smash, covseq_ids):
    for covseq_id in covseq_ids:
        worklist = batch_smash.loc[batch_smash['COVSEQ'] == covseq_id].copy()
        worklist.reset_index(inplace=True, drop=True)
        
        if len(worklist.index) == 94:
            # Change worklist dataframe from BatchSmash format to Tecan format
            batches = list(worklist['Batch'].unique())
            for batch in batches:
                source = f'Source{str(batches.index(batch) + 1)}'
                worklist.loc[worklist['Batch'] == batch, 'Sample Labware Label'] = source
            worklist = worklist.drop(columns=['COVSEQ', 'Batch', 'Sample Name'])
            worklist['Destination Labware Label'] = 'Destination'
            worklist['Volume'] = TRANSFER_uL

            # Rearrange columns to match Tecan format
            worklist = worklist[['Sample Labware Label',
                                 'Sample Location',
                                 'Destination Labware Label',
                                 'Destination Location',
                                 'Volume']]
            batches = [str(batch) for batch in batches]
            worklist.to_csv(f'Output/{covseq_id}_{"_".join(batches)}.csv', index=False)


def write_unsat(unsat_dict):
    for batch, df in unsat_dict.items():
        df = df.copy()
        df['Batch'] = int(batch)
        unsat_dict[batch] = df
    unsats = pd.concat([df for df in unsat_dict.values()])
    unsats.drop(df[(df['Sample Name'] == 'NC') | (df['Sample Name'] == 'PC')].index,
          inplace=True)
    unsats.to_csv(f'Output/UNSAT {date.today()}.csv', index=False)

if __name__ == '__main__':
    # Get samples that meet Ct cutoff criteria
    pcr_results = [f for f in os.listdir() if f.startswith('nCoV')]
    pcr_dict = collect_cts(pcr_results)
    pcr_dict = order_by_column(pcr_dict)
    pcr_dict = pivot_by_detector(pcr_dict)
    passed_dict, unsat_dict = filter_by_ct(pcr_dict)

    # Create BatchSmash file
    old_batch_smash = [f for f in os.listdir() if 'BatchSmash' in f][0]
    old_batch_smash = read_batch_smash(old_batch_smash)
    new_batch_smash, covseq_ids = create_batch_smash(old_batch_smash, passed_dict)
    write_batch_smash(new_batch_smash)

    # Create Worklisting files
    create_worklist(new_batch_smash, covseq_ids)

    # Write samples that do not meet sequencing criteria to UNSAT file
    write_unsat(unsat_dict)

