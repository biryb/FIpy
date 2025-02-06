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
import json

def main():
    tolerance = 0.005
    parser = argparse.ArgumentParser(description='Process mzML files and perform analysis')
    parser.add_argument('dir_mzml', help='Directory containing mzML files')
    args = parser.parse_args()
    init_start_time = time.time()
    os.chdir(args.dir_mzml)
    list_mzml_files = os.listdir(args.dir_mzml)
    dir_analysis = args.dir_mzml
    polarity = infer_polarity(list_mzml_files)
    print(f"Polarity is {polarity}")
    print('Merging m/z within each file')
    dict_all_rawdata,dict_old2newmz = find_consensus_spectra_within_file(list_mzml_files, tolerance=0.005)
    print('Merging m/z between files')
    list_mzs = find_consensus_spectra_between_files(dict_all_rawdata,tolerance=0.005)
    print(f'Found {len(list_mzs)} mzs')
    df_matched_data,average_error = get_matching_mzs(list_mzs,tolerance,polarity)
    print(df_matched_data)
    annot_tolerance = tolerance#max(average_error,0.001)
    print(f"Average deviation between expected and observed mz was {average_error}")
    # Gap filling
    print('Producing an intensity matrix consensus m/z')
    start_time = time.time()
    df_gapfilled = process_mzml_files_with_consensus_mzs(list_mzml_files,list_mzs,tolerance)
    elapsed_time = time.time() - start_time
    print(f"Producing intensity matrix took {elapsed_time/60:.4f} minutes to run.")
    var_threshold = 0.2
    print(f'Filtering out ions that vary more than {var_threshold*100}% between replicates')
    df_filtered_ions = filter_high_variation_ions(df_gapfilled,threshold=var_threshold)
    df_filtered_ions = filter_mismatched_replicates(df_filtered_ions)
    print(f'Went from {df_gapfilled.shape[0]} to {df_filtered_ions.shape[0]} ions')
    prevalence_threshold = 1
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
    with open(os.path.join(output_dir,"within_file_mz_mapping.json"), "w") as file:
        json.dump(dict_old2newmz, file, indent=4)  # `indent=4` makes it readable


    df_annot = load_annotation_file()
    df_annotated_data = annotate(df_filtered,df_annot,polarity,annot_tolerance)
    df_annotated_data.to_excel(os.path.join(output_dir, f"raw_annotated_output_{order_number}.xlsx"))
    
    print(f"Saving annotated data in: {output_dir}")        
    elapsed_time = time.time() - init_start_time
    print(f"Finished and took {elapsed_time/60:.4f} minutes to run.")

if __name__ == "__main__":
    main()
