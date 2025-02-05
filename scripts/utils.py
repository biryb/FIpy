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


def infer_polarity(list_mzml_files):
    
    '''
    Find out if the data were acquired in negative or positive mode
    
    Parameters:
        list_mzml_files (list): List of .mzML file paths.

    '''
    mzml_obj =  mzml.read(list_mzml_files[0])
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


def extract_apex_surrounding_scans(mzml_file):
    """
    Extract m/z and intensity values for the three scans before and after the apex scan.

    Parameters:
        mzml_file (str): Path to the .mzML file.

    Returns:
        list: A list of [mzs, intensities] for the 3 scans before and after the apex scan.
    """
    if 'mzML' in mzml_file:
        print(f"Processing {mzml_file}")
        
        mzml_obj = mzml.read(mzml_file)
        apex_scan = find_apex_scan(mzml_obj)
        
        mzml_obj_2 = mzml.read(mzml_file)  # Re-opening to iterate again
        spectra = list(mzml_obj_2)
        
        start = apex_scan
        end = min(len(spectra), apex_scan + 6)
    
        data = []
        for i in range(start, end):
            spectrum = spectra[i]
            if 'm/z array' in spectrum and 'intensity array' in spectrum:
                mzs = spectrum['m/z array'].tolist()
                intensities = spectrum['intensity array'].tolist()
                data.append([mzs, intensities])
    else:
        print(f'{mzml_file} is not an .mzML file')
        data = []
    return data


def read_mzml(mzml_file):
    """
    Extract m/z and intensity values from an mzML file.

    Parameters:
        mzml_file (str): Path to the .mzML file.

    Returns:
        list: A list of [mzs, intensities] from all scans.
    """
    with mzml.read(mzml_file) as mzml_obj:
        spectra = list(mzml_obj)
        data = []
        
        for spectrum in spectra:
            if 'm/z array' in spectrum and 'intensity array' in spectrum:
                mzs = spectrum['m/z array'].tolist()
                intensities = spectrum['intensity array'].tolist()
                data.append([mzs, intensities])

    return data


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


def merge_mz_int_file(dict_all_rawdata, tolerance=0.001):
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


def merge_files(df, tolerance=0.001):
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


def filter_rare_ions(df_merged, threshold=0.5):
    """
    Remove ions that appear in fewer samples than a given threshold.

    Parameters:
        df_merged (pd.DataFrame): Merged DataFrame with m/z values as indices.
        threshold (float): Minimum fraction of nonzero values required for retention.

    Returns:
        pd.DataFrame: Filtered DataFrame with rare ions removed.
    """
    df_cleaned = df_merged[df_merged.apply(lambda row: (row == 0).sum() / len(row) <= threshold, axis=1)]
    
    return df_cleaned


def process_mzml_files_with_consensus_mzs(list_mzml_files, df_filtered, tolerance, apex_offset=6):
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
            start = apex_scan
            end = min(len(spectra), apex_scan + apex_offset)
    
            # Dictionary for mapping target m/z to intensity for this file
            dict_mz2int = {}
            for i in range(start, end):
                spectrum = spectra[i]
                # Ensure both m/z and intensity arrays exist in the spectrum
                if 'm/z array' in spectrum and 'intensity array' in spectrum:
                    mzs = np.array(spectrum['m/z array'].tolist())
                    intensities = spectrum['intensity array'].tolist()
                    
                    # Iterate through each target m/z in df_filtered.index
                    for target_mz in df_filtered.index:
                        # Find indices where the m/z values are within the specified tolerance of the target
                        indices = np.where((mzs >= target_mz - tolerance) & (mzs <= target_mz + tolerance))[0].tolist()
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
    
    print(f'Gapfilling')

    # Use the first file's dictionary keys (target m/z values) as the index for the final DataFrame
    first_key = next(iter(dict_all_rawdata))
    gapfill_index = dict_all_rawdata[first_key].keys()
    df_gapfilled_data = pd.DataFrame(index=gapfill_index, columns=dict_all_rawdata.keys())
    
    # Populate the DataFrame with the intensity data
    for file_path, mz2int in dict_all_rawdata.items():
        for target_mz, intensity in mz2int.items():
            df_gapfilled_data.loc[target_mz, file_path] = intensity
            
    return df_gapfilled_data

def load_annotation_data():
    # Open the CSV file as a text file from the package's 'data' directory.
    with pkg_resources.open_text('fipy.data', 'annotation.csv') as csv_file:
        # Read the CSV file into a Pandas DataFrame
        df = pd.read_csv(csv_file, encoding='latin1')
    return df


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
        df_annot_mz = df_annot[df_annot[masscol].between(mz-mztol,mz+mztol)]
        if df_annot_mz.shape[0]>0:
            annots = '|'.join(df_annot_mz.name)
            df_data.loc[mz,'annotation'] = annots
    return df_data
<<<<<<< HEAD
            
=======
>>>>>>> 6c8e9d29ca4d5dc57d006b94273c8818121db000
