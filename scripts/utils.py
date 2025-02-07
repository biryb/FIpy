#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:57:44 2025

@author: imac9
"""

import numpy as np
import pandas as pd
from pyteomics import mzml
import time
import requests
import csv
from sklearn.cluster import DBSCAN


def infer_polarity(list_mzml_files):
    
    '''
    Find out if the data were acquired in negative or positive mode
    
    Parameters:
        list_mzml_files (list): List of .mzML file paths.

    '''
    mzml_obj =  mzml.read([x for x in list_mzml_files if 'mzML' in x] [0])
    spec= list(mzml_obj)
    if 'negative scan' in spec[0].keys():
        polarity = 'negative'
    elif 'negative scan' in spec[0].keys():
        polarity = 'positive'
    else:
        if 'NEG' in list_mzml_files[0].upper():
            polarity = 'negative'
        elif 'POS' in list_mzml_files[0].upper():
            polarity = 'positive'
        else:
            polarity='polarity was not inferred'
    return polarity

def find_apex_scan(mzml_obj):
    """
    Identify the scan with the highest sum of intensities in an mzML object.

    Parameters:
        mzml_obj (pyteomics.mzml.MzML): Parsed mzML object.

    Returns:
        int: Index of the scan with the highest total intensity.
    """
    max_intensity = 0
    apex_scan = -1

    for i, spectrum in enumerate(mzml_obj):
        if 'm/z array' in spectrum and 'intensity array' in spectrum:
            total_intensity = sum(spectrum['intensity array'])
            if total_intensity > max_intensity:
                max_intensity = total_intensity
                apex_scan = i
    
    return apex_scan

def read_mzml_data(mzml_file):
    """
    Reads m/z and intensity data from a given mzML file.
    """
    all_mzs, all_intensities = [], []
    with mzml.read(mzml_file) as reader:
        all_spectra = list(reader)
        for spectrum in all_spectra:
            all_mzs.extend(spectrum['m/z array'])
            all_intensities.extend(spectrum['intensity array'])
    return all_mzs, all_intensities

def cluster_mz_values(mz_values, intensities, tolerance):
    """
    Groups m/z values if they are within a specified tolerance, calculating the consensus m/z.
    """
    consensus_mzs = []
    current_cluster = []
    current_intensities = []
    dict_old2newmz = {}
    for mz, intensity in zip(mz_values, intensities):
        if not current_cluster:
            current_cluster.append(mz)
            current_intensities.append(intensity)
        else:
            if abs(mz - current_cluster[-1]) <= tolerance:
                current_cluster.append(mz)
                current_intensities.append(intensity)
            else:
                consensus_mzs.append(np.mean(current_cluster))
                dict_old2newmz[np.mean(current_cluster)] = current_cluster
                current_cluster = [mz]
                current_intensities = [intensity]

    # Add last cluster
    if current_cluster:
        consensus_mzs.append(np.average(current_cluster, weights=current_intensities))

    return consensus_mzs,dict_old2newmz

def find_consensus_spectra_within_file(list_mzml_files, tolerance=0.005):
    """
    Main function to combine m/z spectra from mzML files into consensus m/z values
    using simple merging based on a tolerance.
    """
    mz_dict = {}
    
    for mzml_file in list_mzml_files:
        if '.mzML' in mzml_file:
            print(f'Merging scans {mzml_file}')
            
            # Step 1: Read data from the mzML file
            all_mzs, all_intensities = read_mzml_data(mzml_file)

            # Step 2: Sort m/z values
            sorted_indices = np.argsort(all_mzs)
            sorted_mzs = np.array(all_mzs)[sorted_indices]
            sorted_intensities = np.array(all_intensities)[sorted_indices]
            # Step 3: Cluster m/z values based on tolerance
            consensus_mzs,dict_old2newmz = cluster_mz_values(sorted_mzs, sorted_intensities, tolerance)
            print(f'Went from {len(sorted_mzs)} ions to {len(consensus_mzs)} ions')

            # Step 4: Store the consensus m/z values
            mz_dict[mzml_file] = consensus_mzs
            
    return mz_dict,dict_old2newmz

def cluster_mz_values_across_files(mz_values, tolerance):
    """
    Groups m/z values from multiple files if they are within a specified tolerance,
    calculating the consensus m/z.
    
    :param mz_values: Sorted list of all m/z values from multiple files.
    :param tolerance: Maximum allowed difference for clustering.
    :return: List of final consensus m/z values.
    """
    consensus_mzs = []
    current_cluster = []

    for mz in mz_values:
        if not current_cluster:
            current_cluster.append(mz)
        else:
            if mz - current_cluster[-1] <= tolerance:
                current_cluster.append(mz)
            else:
                consensus_mzs.append(np.mean(current_cluster))
                current_cluster = [mz]

    # Add last cluster
    if current_cluster:
        consensus_mzs.append(np.mean(current_cluster))

    return consensus_mzs

def find_consensus_spectra_between_files(mz_dict, tolerance=0.005):
    """
    Main function to combine consensus m/z values across multiple files.
    
    :param mz_dict: Dictionary {file_name: [consensus_mzs]}.
    :param tolerance: Maximum allowed difference for clustering.
    :return: List of final consensus m/z values across files.
    """
    # Flatten all consensus m/z values from multiple files
    all_consensus_mzs = np.concatenate(list(mz_dict.values()))

    # Step 1: Sort m/z values
    sorted_mzs = np.sort(all_consensus_mzs)

    # Step 2: Cluster m/z values based on tolerance
    final_consensus_mzs = cluster_mz_values_across_files(sorted_mzs, tolerance)

    return final_consensus_mzs



def read_mzml_files(list_mzml_files):
    """
    Read multiple mzML files and extract apex surrounding scan data.

    Parameters:
        list_mzml_files (list): List of .mzML file paths.

    Returns:
        dict: A dictionary mapping filenames to extracted scan data.
    """
    dict_all_rawdata = {}
    for mzml_file in list_mzml_files:
        if 'mzML' in mzml_file:
            filename = mzml_file.split('/')[-1]
            dict_all_rawdata[filename] = extract_apex_surrounding_scans(mzml_file)
        
    return dict_all_rawdata


def merge_mz_int_file(dict_all_rawdata, tolerance):
    """
    Merge m/z values across multiple scans, grouping them within a tolerance.

    Parameters:
        dict_all_rawdata (dict): Input data in {'filename': [[mz1, int1], [mz2, int2], ...]} format.
        tolerance (float): Tolerance for grouping m/z values.

    Returns:
        dict: Dictionary mapping filenames to consensus m/z and summed intensities.
    """
    dict_mz2groupedmz = {}
    
    for filename, mz_int_list in dict_all_rawdata.items():
        print(f"Merging scans in file {filename}")

        # Flatten m/z and intensity arrays across all scans
        mzs = []
        intensities = []
        
        for mz_int_list_scan in mz_int_list:
            mzs.extend(mz_int_list_scan[0])
            intensities.extend(mz_int_list_scan[1])

        mzs = np.array(mzs)
        intensities = np.array(intensities)
        
        consensus_mzs = []
        consensus_ints = []
        
        for mz in mzs:
            mask = (mzs >= mz - tolerance) & (mzs <= mz + tolerance)
            grouped_mzs = mzs[mask]
            grouped_ints = intensities[mask]
            
            consensus_mzs.append(np.mean(grouped_mzs))
            consensus_ints.append(sum(grouped_ints))
        
        dict_mz2groupedmz[filename] = [consensus_mzs, consensus_ints]
    
    return dict_mz2groupedmz


def dict_to_df(dict_consensus_data):
    """
    Convert a dictionary of consensus m/z and intensity values into a DataFrame.

    Parameters:
        dict_consensus_data (dict): Dictionary mapping filenames to [m/z, intensities].

    Returns:
        pd.DataFrame: DataFrame with m/z values as indices and intensity values per file.
    """
    all_indices = set()
    for filename in dict_consensus_data:
        all_indices.update(dict_consensus_data[filename][0])
    
    # Create a mapping of all indices to their corresponding values
    values = {}
    for filename in dict_consensus_data:
        index_map = {index: None for index in all_indices}
        indices, vals = dict_consensus_data[filename]
        
        for idx, val in zip(indices, vals):
            index_map[idx] = val
        
        values[filename] = list(index_map.values())

    df = pd.DataFrame(values, index=sorted(all_indices))
    return df


def merge_files(df, tolerance):
    """
    Merge similar m/z values within a given tolerance and sum intensities.

    Parameters:
        df (pd.DataFrame): Input DataFrame with m/z values as indices.
        tolerance (float): Tolerance for merging m/z values.

    Returns:
        pd.DataFrame: Aggregated DataFrame with merged m/z values.
    """
    groups = (df.index.to_series().diff().abs() > tolerance).cumsum()
    dict_newind2mz = {}
    for value in groups.unique():
        rows = groups[groups==value]
        mz = rows.index
        new_mz = np.mean(mz)
        dict_newind2mz[value] = new_mz

    df_agg = df.groupby(groups).sum()
    df_agg.index = df_agg.index.map(dict_newind2mz)
    return df_agg


def filter_rare_ions(df, threshold=0.5):
    """
    Remove ions that appear in fewer samples than a given threshold.

    Parameters:
        df_merged (pd.DataFrame): Merged DataFrame with m/z values as indices.
        threshold (float): Minimum fraction of nonzero values required for retention.

    Returns:
        pd.DataFrame: Filtered DataFrame with rare ions removed.
    """
    df = df.replace(np.nan,0)
    df_filtered = df[df.apply(lambda row: (row <10).sum() / len(row) <= threshold, axis=1)]
    
    return df_filtered


def get_matching_mzs(list_mzs, tolerance, polarity):
    """
    Function to get m/z values from a publicly accessible Google Sheet based on the specified polarity.
    
    :param list_mzs: List of m/z values to check.
    :param polarity: 'positive' or 'negative' to select the appropriate m/z column.
    :param sheet_url: URL of the Google Sheets document (in CSV format).
    :return: List of matching m/z values from the Google Sheets.
    """
    # Load the Google Sheets data as a pandas DataFrame
    sheet_url = 'https://docs.google.com/spreadsheets/d/19FO85OiCjMch5Wy2OdVWB4hSUq15jaH2I1wTVdwN_So/export?format=csv'
    df = pd.read_csv(sheet_url)
    
    # Check for the correct polarity and select the appropriate column
    if polarity == 'negative':
        mz_column = 'mz_neg'  # Use 'mz_neg' column for negative polarity
    elif polarity == 'positive':
        mz_column = 'mz_pos'  # Use 'mz_pos' column for positive polarity
    else:
        raise ValueError("Polarity must be either 'positive' or 'negative'.")
    
    # Find matching m/z values in the list
    df_matched_data = pd.DataFrame()
    for mz in list_mzs:
        df_match = df[df[mz_column].between(mz-tolerance,mz+tolerance)]
        if df_match.shape[0]>0:
            df_match['mz_observed'] = mz
            df_matched_data = pd.concat([df_matched_data,df_match])

    df_matched_data['delta_mz'] = df_matched_data[mz_column]-df_matched_data.mz_observed
    average_error = abs(df_matched_data['delta_mz'].mean())
    return df_matched_data,average_error

def process_mzml_files_with_consensus_mzs(list_mzml_files, list_mzs, tolerance, apex_offset=8):
    """
    Process a list of mzML files to generate a gap-filled DataFrame of intensity values 
    for target m/z values across samples.
    
    For each file, the function:
      - Reads the file and finds the apex scan.
      - Processes a range of scans starting at the apex scan up to (apex + apex_offset).
      - For each spectrum in that range, it sums the intensities for m/z values that are 
        within ±tolerance of each target m/z from df_filtered.
      - Stores these summed intensities for each file.
      
    Parameters:
      list_mzml_files (list of str): Paths to mzML files.
      df_filtered (pd.DataFrame): DataFrame with index as target m/z values.
      tolerance (float): Allowed deviation when matching m/z values.
      apex_offset (int): Number of scans to process after the apex scan.
      
    Returns:
      pd.DataFrame: A gap-filled DataFrame with target m/z values as the index and 
                    filenames as columns, containing the summed intensities.
    """
    # Dictionary to store mapping for each file: {file_path: {target_mz: intensity, ...}, ...}
    dict_all_rawdata = {}
    
    # Process each file in the list
    for mzml_file in list_mzml_files:
        # Ensure the file is an mzML file based on its name
        if 'mzML' in mzml_file:
            print(f'Collecting data for {mzml_file}')
            # Extract filename from the full path
            filename = mzml_file.split('/')[-1]
            
            # Read the mzML file (first read to find the apex scan)
            mzml_obj = mzml.read(mzml_file)
            apex_scan = find_apex_scan(mzml_obj)
            
            # Re-read the mzML file so we can iterate over its spectra
            mzml_obj_2 = mzml.read(mzml_file)
            spectra = list(mzml_obj_2)
            
            # Define the scan range: from apex to (apex + apex_offset), ensuring we don't exceed the number of scans
            start = apex_scan-2
            end = min(len(spectra), apex_scan + apex_offset)
    
            # Dictionary for mapping target m/z to intensity for this file
            dict_mz2int = {}
            for i in range(start, end):
                spectrum = spectra[i]
                # Ensure both m/z and intensity arrays exist in the spectrum
                if 'm/z array' in spectrum and 'intensity array' in spectrum:
                    mzs = np.array(spectrum['m/z array'].tolist())
                    intensities = spectrum['intensity array'].tolist()
                    
                    for target_mz in list_mzs:
                        # Find indices where the m/z values are within the specified tolerance of the target
                        indices = np.where((mzs >= target_mz - tolerance) 
                        & (mzs <= target_mz + tolerance))[0].tolist()
                        if indices:
                            # Sum intensities for these indices
                            intensity_sum = sum([intensities[idx] for idx in indices])
                            dict_mz2int[target_mz] = intensity_sum
            
            # If any data was recorded, add it to the overall dictionary
            if dict_mz2int:
                dict_all_rawdata[mzml_file] = dict_mz2int
        else:
            print(f"{mzml_file} is not an mzML file")
    
    # If no files were processed, return an empty DataFrame
    if not dict_all_rawdata:
        print("No valid mzML data processed.")
        return pd.DataFrame()
    
    # Create a MultiIndex for the DataFrame (index = target m/z, columns = file paths)
    first_key = next(iter(dict_all_rawdata))
    gapfill_index = dict_all_rawdata[first_key].keys()
    df_gapfilled_data = pd.DataFrame(index=gapfill_index, columns=dict_all_rawdata.keys())

    # Convert the dict into a list of tuples (target_mz, file_path, intensity)
    data = [(target_mz, file_path, intensity)
            for file_path, mz2int in dict_all_rawdata.items()
            for target_mz, intensity in mz2int.items()]

    # Create a DataFrame directly from the list of tuples
    df_temp = pd.DataFrame(data, columns=['target_mz', 'file_path', 'intensity'])

    # Pivot to reshape the DataFrame and assign values to the final df_gapfilled_data
    df_gapfilled_data = df_temp.pivot(index='target_mz', columns='file_path', values='intensity')

 
    return df_gapfilled_data

def filter_high_variation_ions(df, threshold=0.20):
    # Find groups of replicates based on column names
    replicate_columns = {}
    
    for col in df.columns:
        sample_name = col.split('__')[0]  # Extract sample name (before the "__")
        if sample_name not in replicate_columns:
            replicate_columns[sample_name] = []
        replicate_columns[sample_name].append(col)
    
    # Store rows to drop
    df_rowcounter = pd.DataFrame(0, index=df.index,columns=['Counter'])
    # Process each sample's replicates
    for sample_name, columns in replicate_columns.items():
        if len(columns) > 1:  # Only process if there are multiple replicates
            for i, col_a in enumerate(columns):
                for col_b in columns[i+1:]:
                    # Compute mean safely, avoiding zero values
                    safe_mean = df[[col_a, col_b]].replace(0, np.nan).mean(axis=1)
                    
                    # Avoid division by zero issues
                    diff_percentage = (df[col_a] - df[col_b]).abs() / safe_mean
                    
                    # Replace NaNs with 0 (if all values were zero)
                    diff_percentage = diff_percentage.fillna(0)
                    diff_percentage = diff_percentage[diff_percentage<threshold]
                    df_rowcounter.loc[diff_percentage.index,'Counter'] +=1
                    
    n_samples = df.shape[0]
    df_rowcounter_sel = df_rowcounter[df_rowcounter.Counter>0.5*n_samples]
    # Drop identified rows
    df_filtered = df.drop(index=df_rowcounter_sel.index, errors='ignore')
        
    return df_filtered

def filter_mismatched_replicates(df):
    df = df.replace(np.nan,0)
    replicate_columns = {}
    
    for col in df.columns:
        sample_name = col.split('__')[0]  # Extract sample name (before the "__")
        if sample_name not in replicate_columns:
            replicate_columns[sample_name] = []
        replicate_columns[sample_name].append(col)
    inds_to_remove = []
    for i,row in df.iterrows():
        for samplegroup in replicate_columns:
            colnames = [x for x in replicate_columns[samplegroup]]
            if 0 in row[colnames].values:
                if not row[colnames].eq(0).all():
                    inds_to_remove.append(i)
    df = df.drop(inds_to_remove)
    return df
                

def load_annotation_file():
    url="https://docs.google.com/spreadsheets/d/15TVDmFFBHW4FWc8VfIlevmTx62t36c7BE7d01kIaVcQ/export?format=csv"
    df_annot = pd.read_csv(url)
    return df_annot

def annotate(df_data,df_annot,polarity,tolerance):
    df_annot['ion_ID'] = range(1,df_annot.shape[0]+1)
    masscol = 'monisotopic_molecular_weight'
    mass_proton = 1.007276466583
    if polarity=='negative':
        df_annot['mz'] = df_annot[masscol]-mass_proton
    elif polarity=='positive':
        df_annot['mz'] = df_annot[masscol]+mass_proton
    else:
        print('Polarity unknown, assuming negative')
        df_annot['mz'] = df_annot[masscol]-mass_proton
    df_data['annotation'] = ''
    for mz,row in df_data.iterrows():
        df_annot_mz = df_annot[df_annot['mz'].between(mz-tolerance,mz+tolerance)]
        if df_annot_mz.shape[0]>0:
            annots = '|'.join(df_annot_mz.name.astype(str))
            df_data.loc[mz,'annotation'] = annots
    return df_data
            