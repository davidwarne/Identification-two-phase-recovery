# R Code for Data Processing and Analysis for the Identification of Two-Phase Recovery Patterns using Reef Monitoring Data
This repository contains useful R functions and scripts to process coral cover data obtained from reef surveys to classify patterns in recovery trajectories following major disturbance events.

## Developer

David J. Warne (david.warne@qut.edu.au),
                School of Mathematical Sciences, 
                Faculty of Science, 
                Queensland Univeristy of Technology 
                
Google Scholar: (https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en)

Supported by [ARC Centre of Excellence for Mathematical and Statistical Frontiers](https://acems.org.au/home), and [QUT Centre for Data Science](https://research.qut.edu.au/qutcds/)

## Citation Information

This code is provided as supplementary information to the paper,

David J Warne, Kerryn A. Crossman, Wang Jin, Kerrie Mengersen, Kate Osborne, Matthew J. Simpson, Angus A. Thompson, Paul Wu, and Juan-C. Ortiz. Identification of two-phase recovery for interpretation of reef monitoring data. Journal of Applied Ecology, (https://doi.org/f10.1111/1365-2664.14039. 

## Data Location

The data used in this analysis can be obtained from,

Australian Institute of Marine Science (AIMS) [DOI:10.25845/ad6j-zm19](https://doi.org/10.25845/ad6j-zm19)

## Licensing
This source code is licensed under the GNU General Public License Version 3.
Copyright (C) 2021 David J. Warne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Contents

This folder contains a number of R functions and scripts that implement the analysis. 
```bash
The directory structure is as follows
|-- LTMPTools/ 
    |-- LTMPDataTools.R                                 Functions to process LTMP and MMP data
    |-- LTMPModellingTools.R                            Functions to fit models and analyse trajectories
|-- GBR_recovery_analysis/
    |-- setup_analysis_pipeline.R                       Script to define all parameters, directories and libraries used in analysis
    |-- script_main_init_transect_analysis_pipeline.R   Main script to run to repeat entire analysis
    |-- power_analysis_table.csv                        Supplementary Table showing power analysis results
    |-- data_processing/                                Support scripts for data processing
    |-- model_fitting/                                  Support scripts for model fitting
    |-- analysis/                                       Support scripts for analysis
    |-- plotting/                                       Scripts to draw all figures in the paper and supplementary material
```
## Usage

Follow these steps to run the analysis:

1. Load data from [here](https://doi.org/10.25845/ad6j-zm19)
2. Extract `data.zip` into folder Identification-two-phase-recovery/data/
3. Start RStudio
4. In RStudio set the working directory to Identification-two-phase-recovery/GBR_recovery_analysis
5. In the RStudio Console, run 
   `> source('script_main_init_transect_analysis_pipeline.R', echo=TRUE)` 
