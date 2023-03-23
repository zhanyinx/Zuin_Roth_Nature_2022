[![Not Maintained](https://img.shields.io/badge/Maintenance%20Level-Not%20Maintained-yellow.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Code for [Zuin, Roth et al., Nature 2022](https://www.nature.com/articles/s41586-022-04570-y)

## Content

### HiC-Pro_parameters folder:

Within this folder, you find the HiC-Pro parameter file utilized for the analysis of the captureC data, as well as details regarding the target region in captureC.

### Nanopore folder

Within this folder, you find the Snakemake pipeline to analyse nanopore sequencing data. At the moment, the parameters are hard-coded at the beginning of the SnakeMake file. TODO: Introduce configuration file

### enhancer_insertion_mapping_individual_cell_lines

Within this folder, you find the R script to map Sanger sequencing data automatically (credits: Pia Mach)
