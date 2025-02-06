#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:17:52 2025

@author: imac9
"""

import os
import time
import argparse  # Import argparse for argument parsing
from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import *
from . import utils

def main():
    parser = argparse.ArgumentParser(description='Process mzML files and perform analysis')
    parser.add_argument('dir_mzml', help='Directory containing mzML files')
    args = parser.parse_args()
    init_start_time = time.time()
    os.chdir(args.dir_mzml)
    list_mzml_files = os.listdir(args.dir_mzml)
    dir_analysis = args.dir_mzml
    polarity = infer_polarity(list_mzml_files)
    
    # Read in files and generate a dictionary with all the data
    print('Reading .mzML files')
    dict_all_rawdata = read_mzml_files(list_mzml_files)
    print('Merging .mzML files')
    start_time = time.time()
    dict_consensus_data = merge_mz_int_file(dict_all_rawdata, tolerance=0.005)
    elapsed_time = time.time() - start_time
    print(f"Merging m/z's within each file took {elapsed_time/60:.4f} minutes to run.")
    
    # Collect all the unique indices from all files
    df = dict_to_df(dict_consensus_data)
    
    # Find consensus m/z's in the dataframe (inter-sample)
    print('Merging files')
    start_time = time.time()
    
    # Identify groups of rows where indices are within the tolerance
    df_merged = merge_files(df)
    
    elapsed_time = time.time() - start_time
    print(f"Merging m/z's between files took {elapsed_time/60:.4f} minutes to run.")
    
    # Gap filling
    print('Gapfilling based on consensus m/z')
    start_time = time.time()
    #df_gapfilled = map_mz_to_consensus(list_mzml_files,df_merged,tolerance=0.005)
    df_gapfilled = process_mzml_files_with_consensus_mzs(list_mzml_files,df_merged,tolerance=0.005)
    elapsed_time = time.time() - start_time
    print(f"Gapfilling took {elapsed_time/60:.4f} minutes to run.")

    var_threshold = 0.2
    print(f'Filtering out ions that vary more than {var_threshold*100}% between replicates')
    df_filtered_ions = filter_high_variation_ions(df_gapfilled,threshold=var_threshold)
    print(f'Went from {df_gapfilled.shape[0]} to {df_filtered_ions.shape[0]} ions')
    prevalence_threshold = 0.5
    print(f'Filtering out ions present in less than {prevalence_threshold*100}% of the samples')
    df_filtered = filter_rare_ions(df_filtered_ions,threshold=prevalence_threshold)
    print(f'Went from {df_filtered_ions.shape[0]} to {df_filtered.shape[0]} ions')
    
    output_dir = os.path.join(dir_analysis, "fipy_output")
    # Create the folder if it doesn't exist
    try:  
        os.mkdir(output_dir)
        order_number = 1
    except OSError as error:  
        print("Folder exist, appending to existing folder")
        order_number = len(os.listdir(output_dir))
    print(f"Saving filtered raw data to {output_dir}")
    df_filtered.to_excel(os.path.join(output_dir, f"raw_output_{order_number}.xlsx"))
    df_annot = load_annotation_file()
    df_annotated_data = annotate(df_filtered,df_annot,polarity,tolerance=0.005)
    df_annotated_data.to_excel(os.path.join(output_dir, f"raw_annotated_output_{order_number}.xlsx"))
    
    print(f"Saving annotated data in: {output_dir}")        
    elapsed_time = time.time() - init_start_time
    print(f"Finished and took {elapsed_time/60:.4f} minutes to run.")
if __name__ == "__main__":
    main()
