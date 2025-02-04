# FIpy
Flow injection analysis processing workflow for untargeted metabolomics.

# What does FIpy do to my data?

1. It reads in the files and determines the highest intensity scans<br>
2. Within each file, it flattens the time dimension by merging m/z values from consecutive scans and summing up corresponding intensities<br>
3. It merges all the data into one large dataframe<br>
4. It merges m/z values in the dataframe based on similarity (rows with m/z 143.1221 and m/z 143.1220 would be combined)<br>
5. It filters out ions that were present in less than 50% of the files

# Data requirements
This script uses flow injection-MS1 data. It's been tested with 1 minute injections as in Fuhrer et al 2011 (https://pubs-acs-org.ezp-prod1.hul.harvard.edu/doi/10.1021/ac201267k), but could theoretically use any flow injection data. While developed for TOF's, there's no explicit parameters to exclude using it on Orbitrap data - however, it requires a high-resolution detector to produce meaningful outputs.
It uses the few scans that have the most ions in them as inferred from the TIC profile:

<img width="500" alt="image" src="https://github.com/user-attachments/assets/299fa61e-40c2-4a0c-a740-6efc9ac7e310" />

The input is a folder with .mzML files. Prior to running FIpy, convert your raw files to .mzML with the following settings

<img width="500" alt="image" src="https://github.com/user-attachments/assets/c90588ae-3b81-454f-8f84-51a3dd1add27" />

# Usage
Download FIpy script and run it in your IDE.

# Requirements
python 3.10.9<br>
pyteomics 4.6<br>
pandas 1.5.3<br>
numpy 1.23.5<br>
