#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:17:52 2025

@author: imac9
"""

import os
import time
import argparse  # Import argparse for argument parsing
from utils import *
from pathlib import Path
import pandas as pd

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process mzML files and perform analysis')
    parser.add_argument('dir_mzml', help='Directory containing mzML files')
    parser.add_argument('dir_analysis', help='Directory to save analysis results')

    # Parse the arguments
    args = parser.parse_args()

    # Now you can use args.dir_mzml and args.dir_analysis directly
    init_start_time = time.time()
    os.chdir(args.dir_mzml)
    list_mzml_files = os.listdir(args.dir_mzml)
    
    # Read in files and generate a dictionary with all the data
    print('Reading .mzML files')
    dict_all_rawdata = read_mzml_files(list_mzml_files)
    print('Merging .mzML files')
    start_time = time.time()
    dict_consensus_data = merge_mz_int_file(dict_all_rawdata, tolerance=0.001)
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
    print(f"Merging m/z's between files took {elapsed_time:.4f} seconds to run.")
    
    print('Filtering out ions present in less than 50% of the samples')
    df_filtered = filter_rare_ions(df_merged)
    
    output_dir = os.path.join(args.dir_analysis, "fipy_output")
    
    # Create the folder if it doesn't exist
    try:  
        os.mkdir(output_dir)
    except OSError as error:  
        print(error)

    print(f"Saving results in: {output_dir}")        
    df_filtered.to_excel(os.path.join(output_dir, "raw_output.xlsx"))
    elapsed_time = time.time() - init_start_time
    print(f"Finished and took {elapsed_time/60:.4f} minutes to run.")

if __name__ == "__main__":
    main()
