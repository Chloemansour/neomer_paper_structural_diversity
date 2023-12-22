# A new method for the reproducible development of aptamers (Neomers) - Structural Diversity 
---
## Overview
This repository contains a python script and necessary data needed to generate the structural diversity data from the Nature Biotechnology submission titled, "A new method for the reproducible development of aptamers (Neomers)".

For further information on the Neomer method and result visualizations, please refer to our preprint: [here](https://biorxiv.org/cgi/content/short/2023.12.19.572437v1).

Copyright (c) 2023 NeoVentures Biotechnology Inc. and NeoVentures Biotechnology Europe SAS\
\
**Please note**: We do not have a license in this repository. The rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the software, data, and all other materials is fully restricted to developers from NeoVentures Biotechnology Inc. and NeoVentures Biotechnology Europe SAS. Permission for outside use can be granted upon request. 

## Setup
Code is written in python3. Please install the packages present in the requirements.txt file. You may use:
```
pip install -r requirements.txt
```
We use RNAfold and bpRNA programs to generate the the secondary predicted structures and structural annotations. Please make sure they are installed on your device prior to running the python script.

Please see the following instructions for download of:
1. [RNAfold package from ViennaRNA](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html#binary-packages)
2. [bpRNA](https://github.com/hendrixlab/bpRNA)

## Running the script
The python3 script *Secondary structure analysis.py* was used for the analysis.
\
On line 14, please change the following variable to the path to location of the bpRNA.pl file:
```
bpRNA_path = '/Path/to/bpRNA.pl'
```
Example outputs are provided in the Output data folder. 
To visualize the script outputs, please refer to the repository that contains the publication's visualization scripts: [here](https://github.com/Chloemansour/neomer_paper_visualizations).

## Contact
Please contact the corresponding author with any questions. 
