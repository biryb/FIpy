# FIpy
Flow injection analysis processing workflow for untargeted metabolomics

# Data requirements
This script uses flow injection-MS data. It's been tested with 1 minute injections as in Fuhrer et al 2011 (), but could theoretically use any flow injection data. While developed for TOF's, there's no explicit parameters to exclude using it on Orbitrap data - however, it requires a high-resolution detector to produce meaningful outputs.
It uses the few scans that have the most ions in them as inferred from the TIC profile:

<img width="931" alt="image" src="https://github.com/user-attachments/assets/299fa61e-40c2-4a0c-a740-6efc9ac7e310" />


# Usage
Download FIpy script and run it in your IDE.

# Requirements
python 3.10.9<br>
pyteomics 4.6<br>
pandas 1.5.3<br>
numpy 1.23.5<br>
