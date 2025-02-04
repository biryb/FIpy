# FIpy
FIpy is a package for a processing and merging workflow for `mzML` files and can be used for untargeted mass spectrometry analysis of metabolites, lipids, and peptides.

# What does FIpy do to my data?

1. It reads in the files and determines the highest intensity scans<br>
2. Within each file, it flattens the time dimension by merging m/z values from consecutive scans and summing up corresponding intensities<br>
3. It merges all the data into one large dataframe<br>
4. It merges m/z values in the dataframe based on similarity (rows with m/z 143.1221 and m/z 143.1220 would be combined)<br>
5. It filters out ions that were present in less than 50% of the files
6. It writes the resulting dataframe into .xlsx file

# Data requirements
FIpy uses flow injection-MS1 data. It was developed for 1 minute injections as in Fuhrer et al 2011 (https://pubs-acs-org.ezp-prod1.hul.harvard.edu/doi/10.1021/ac201267k), but could theoretically use any flow injection data. While developed for TOF's, there's no explicit parameters to exclude using it on Orbitrap data - however, it requires a high-resolution detector to produce meaningful outputs.<br>

The input is a folder with .mzML files. Prior to running FIpy, convert your raw files to .mzML with the following settings

<img width="500" alt="image" src="https://github.com/user-attachments/assets/c90588ae-3b81-454f-8f84-51a3dd1add27" />

# Installation

You can install **FIpy** using `pip` either from GitHub or from a local directory.

## Install from GitHub

To install the package directly from GitHub, run the following command:
```bash
pip install git+https://github.com/biryb/FIpy.git
```

## Install from a Local Directory

If you'd like to install from your local repository, download the repository, navigate to the folder containing the `setup.py` file and run:
```bash
pip install .
```

# Usage

Command-Line Interface (CLI)<br>

The package provides command-line tools for processing your `mzML` files and analyzing them<br>

```bash
fipy <dir_mzml> <dir_analysis>
```

FIpy outputs an Excel file in dir_analysis

# Contact
Reach out to birgittaryback@gmail.com with questions and comments

