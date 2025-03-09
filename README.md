# DIpy
DIpy is a command line tool for a processing and merging workflow for `mzML` files and can be used for untargeted mass spectrometry analysis of metabolites, lipids, and peptides.

# What does DIpy do to my data?

1. It reads in the files and determines the highest intensity scans<br>
2. Within each file, it flattens the time dimension by merging m/z values from consecutive scans and summing up corresponding intensities<br>
3. Taking the {file : consensus m/z array} mapping, it looks for consensus m/z values between files <br>
4. It queries the .mzML files for the consensus m/z's and builds a pandas dataframe with the corresponding intensities<br>
4. The resulting consensus matrix is filtered three ways:
     - Ions that vary more than 20% between technical replicates in more than 50% of the samples
     - Ions that had 0 abundance in one technical replicate and above in the other
     - Ions that were present in less than 50% of the files
5. It annotates the ions based on accurate mass
6. It writes the resulting dataframe into .xlsx file

# Data requirements
DIpy uses flow injection-MS1 data. It was developed for TOF's, but there's no explicit parameters to exclude using it on Orbitrap data - however, it requires a high-resolution detector to produce meaningful outputs.<br>

If technical injection replicates are used (which is recommended), they should be named as uniquesampleidentifier__replicateID.mzML. DIpy uses this information to filter out ions that are very dissimilar between the technical replicates. It recognizes replicates from "__" specifically. <br>

The input is a folder with .mzML files. Prior to running DIpy, convert your raw files to .mzML with the following settings

<img width="500" alt="image" src="https://github.com/user-attachments/assets/c90588ae-3b81-454f-8f84-51a3dd1add27" />

# Installation

You can install **DIpy** using `pip` either from GitHub or from a local directory.

## Install from GitHub

To install the package directly from GitHub, run the following command:
```bash
pip install git+https://github.com/biryb/DIpy.git
```

## Install from a Local Directory

If you'd like to install from your local repository, download the repository, navigate to the folder containing the `setup.py` file and run:
```bash
pip install .
```

# Usage

It's recommended to use a virtual environment to run DIpy.<br>
If using Conda, open a terminal (Mac) or command line/powershell (Windows) and run the following commands:<br>

```bash
conda create --name DIpy python=3.10
conda activate DIpy
```

In Mac terminal, the prompt should appear as
```bash
(DIpy) id@macname DIpy %                                                                         
```

Then run
```bash
pip install git+https://github.com/biryb/DIpy.git
```
to install DIpy.<br>

That's it - you're ready to process samples. The basic usage is like this:

```bash
DIpy <dir_mzml>
```

DIpy creates a subfolder called "DIpy_output" in dir_mzml and saves results in that path as an Excel file.

# Contact
This is an actively managed project. Reach out to birgittaryback@gmail.com with questions and comments - especially if DIpy did not work on your system.
