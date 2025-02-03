#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:46:23 2025

@author: imac9
"""

from pyteomics import mzml
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import time
import pandas as pd
import random

#%%

dir_plots=r'/Volumes/bryback/_research_and_development/FIpy/plots'
#%%

def find_apex_scan(mzml_obj):
    """
    Find the scan with the highest sum of intensities in an mzML object.

    Parameters:
        mzml_obj (pyteomics.mzml.MzML): An mzML object parsed using Pyteomics.

    Returns:
        int: Scan index corresponding to the highest sum of intensities.
    """
    max_intensity = 0
    apex_scan = -1

    for i, spectrum in enumerate(mzml_obj):
        if 'm/z array' in spectrum and 'intensity array' in spectrum:
            total_intensity = sum(spectrum['intensity array'])
            if total_intensity > max_intensity:
                max_intensity = total_intensity
                apex_scan = i
    print(f'scan number {apex_scan} identified as apex')
    return apex_scan

def extract_apex_surrounding_scans(mzml_file):
    """
    Extract m/z and intensity values for 3 scans before and after the apex scan.

    Parameters:
        mzml_file (str): Path to the .mzML file.

    Returns:
        dict: A dictionary with the filename as the key and a list of 
              [mzs, intensities] for the 3 scans before and after the apex scan.
    """
    print(f'processing {mzml_file}')
    mzml_obj = mzml.read(mzml_file)
    apex_scan = find_apex_scan(mzml_obj)
    dict_rawdata = {}

    mzml_obj_2=mzml.read(mzml_file)
    spectra = list(mzml_obj_2)
    start = max(0, apex_scan - 3)
    end = min(len(spectra), apex_scan + 3)

    data = []
    for i in range(start, end):
        spectrum = spectra[i]
        if 'm/z array' in spectrum and 'intensity array' in spectrum:
            mzs = spectrum['m/z array'].tolist()
            intensities = spectrum['intensity array'].tolist()
            data.append([mzs, intensities])


    return data



def read_mzml(mzml_file):
    """
    Extract m/z and intensity values for 3 scans before and after the apex scan.

    Parameters:
        mzml_file (str): Path to the .mzML file.

    Returns:
        dict: A dictionary with the filename as the key and a list of 
              [mzs, intensities] for the 3 scans before and after the apex scan.
    """
    with mzml.read(mzml_file) as mzml_obj:
        dict_rawdata = {}

        with mzml.read(mzml_file) as mzml_obj_2:  # Re-open to allow re-iteration
            spectra = list(mzml_obj_2)

            data = []
            for spectrum in spectra:
                if 'm/z array' in spectrum and 'intensity array' in spectrum:
                    mzs = spectrum['m/z array'].tolist()
                    intensities = spectrum['intensity array'].tolist()
                    data.append([mzs, intensities])

    return data

def read_mzml_files(list_mzml_files):
    dict_all_rawdata = {}
    for mzml_file in list_mzml_files:
        filename = mzml_file.split('/')[-1]
        dict_all_rawdata[filename] = extract_apex_surrounding_scans(mzml_file)
    return dict_all_rawdata
#%%
dir_mzml = r'/Volumes/bryback/_research_and_development/FIpy/test_data/20230309135232_E00000014_NEG/mzML_files_centroided/'
os.chdir(dir_mzml)
list_mzml_files = os.listdir(dir_mzml)
dict_all_rawdata = read_mzml_files(list_mzml_files)333

#%%
dict_data = {}
random_files = random.sample(list_mzml_files,k=10)
random_files = [x for x in random_files if 'blank' not in x]
list_apex_scans = []
for file in random_files:
    with mzml.read(file) as mzml_obj:
        spectra = list(mzml_obj)
    with mzml.read(file) as mzml_obj:
        apex_scan = find_apex_scan(mzml_obj)
        list_apex_scans.append(apex_scan)
    with mzml.read(file) as mzml_obj:
        dict_data[file] = read_mzml(mzml_obj)


list_tics = []
for key in dict_data:
    dict_key = dict_data[key]
    i = 0
    apex_scan = list_apex_scans[i]
    start = apex_scan-3
    end = apex_scan+3
    for mz,ints in dict_key:
        total_int = sum(ints)
        if i == apex_scan:
            included_in_analysis = 'Yes;apex'
        elif start<=i<=end:
            included_in_analysis = 'Yes'
        else:
            included_in_analysis = 'No'
        list_tics.append([i,total_int,included_in_analysis])
        i+=1
        

df_tic = pd.DataFrame(list_tics,columns=['Scan number','TIC','datapoint_identity'])

sns.set(style="white", context="notebook", font_scale=1.5, rc={"axes.spines.top": True, "axes.spines.right": True})
fig, ax = plt.subplots(figsize=(10, 7), dpi=400)
fig = sns.scatterplot(x='Scan number',y='TIC',data=df_tic,edgecolor='k',hue='datapoint_identity', palette= {'Yes;apex':'red','Yes':'orange','No':'grey'},
                linewidth=1,s=75)
fig.legend(loc=1, fancybox=False, edgecolor='k', frameon=True,fontsize=15,title='Datapoint included')
plt.tight_layout()
plt.savefig("%s/tic_vs_scan.svg"%dir_plots, format="svg")


#%%

def merge_mz_int_file(data, tolerance=0.001):
    """
    Map m/z values across n vectors to consensus m/z within a tolerance.
    
    Parameters:
    - data (dict): Input dictionary in the format {'filename': [[mz1, int1], [mz2, int2], ...]}
    - tolerance (float): Tolerance for grouping m/z values.
    
    Returns:
    - dict: New dictionary with consensus m/z and summed intensities.
    """
    dict_mz2groupedmz = {}
    
    for filename,mz_int_list in dict_all_rawdata.items():
        print(f'Merging scans in file {filename}')
        # Create a long flat array of all m/zs and ints
        mzs = []
        ints = []
        for mz_int_list_scan in mz_int_list:
            mzs.extend(mz_int_list_scan[0])
            ints.extend(mz_int_list_scan[1])
        mzs = np.array(mzs)
        ints = np.array(ints)
        consensus_mzs = []
        consensus_ints = []
        for mz in mzs:
            grouped_mzs = mzs[(mzs >= (mz-tolerance)) & (mzs <= (mz+tolerance))]
            grouped_ints = ints[(mzs >= (mz-tolerance)) & (mzs <= (mz+tolerance))]
            consensus_mzs.append(np.mean(grouped_mzs))
            consensus_ints.append(sum(grouped_ints))
        dict_mz2groupedmz[filename] = [consensus_mzs,consensus_ints] 
    return dict_mz2groupedmz

start_time = time.time()
dict_consensus_data = merge_mz_int_file(dict_all_rawdata, tolerance=0.001)
elapsed_time = time.time() - start_time
print(f"The function took {elapsed_time:.4f} seconds to run.")

#%% Move from dicts to dataframes
# Step 1: Collect all the unique indices from all files

def dict_to_df(dict_consensus_data):
    all_indices = set()
    for filename in dict_consensus_data:
        all_indices.update(dict_consensus_data[filename][0])
    
    # Step 2: Create a dictionary to store the values for each file
    values = {}
    for filename in dict_consensus_data:
        # Create a dictionary where the keys are the full set of indices (from all files)
        # and the values are the corresponding values (fill missing indices with NaN)
        index_map = {index: None for index in all_indices}
        
        # Update the dictionary with actual values from the current file
        indices, vals = dict_consensus_data[filename]
        for idx, val in zip(indices, vals):
            index_map[idx] = val
        
        values[filename] = list(index_map.values())
    
    # Step 3: Create the DataFrame
    df = pd.DataFrame(values, index=sorted(all_indices))
    return df

df = dict_to_df(dict_consensus_data)

#%% Find consensus m/z's in the dataframe (inter-sample)
print('Merging files')
# Identify groups of rows where indices are within the tolerance
groups = (df.index.to_series().diff().abs() > tolerance).cumsum()
dict_newind2mz = {}
for value in groups.unique():
    rows = groups[groups==value]
    mz = rows.index
    new_mz = np.mean(mz)
    dict_newind2mz[value] = new_mz
    
# Aggregate by group, summing the values in each group
df_agg = df.groupby(groups).sum()
#%%

df_agg.index = df_agg.index.map(dict_newind2mz)
threshold = 0.5  # 50% zero values
df_cleaned = df_agg[df_agg.apply(lambda row: (row == 0).sum() / len(row) <= threshold, axis=1)]


#%%
df_cleaned['mz'] = df_cleaned.index.values
list_mz = [89.0244]
