# Source_Contributions_PM2.5_Mortality_Colombia
Analysis scripts supporting the study of PM$_{2.5}$ source contributions and pollution-related mortality in Colombia using WRF-Chem simulations and sectoral emission sensitivity experiments. This paper is submitted to Geohealth

This repository contains the scripts used to prepare emissions, process model outputs, and compute sectoral contributions to PM₂.₅ concentrations for the study “Source Contributions to PM₂.₅ Related Premature Mortality in Colombia.”

The workflow integrates global emission inventories, atmospheric chemistry simulations with WRF-Chem, and post-processing tools to quantify the contributions of major emission sectors to air pollution and related health impacts across Colombia.

## Emission Data Preparation
Scripts used to obtain and preprocess emission inventories:

#### EDGAR_download_submit.py
Automatically downloads global emission datasets from the EDGAR v8.1 inventory for multiple pollutants and sectors. 

#### EDGAR_processing_for_WRF_submit.py
Converts EDGAR emission files into WRF-Chem compatible NetCDF files, including unit conversions (e.g., tonnes → kg m⁻² s⁻¹) and grid formatting.

#### Subs_totals_submit.py
Groups EDGAR sectors into broader categories (e.g., agriculture, transport, industry) and computes sector-specific emissions by subtracting grouped sectors from the total inventory. 

#### fix_files_to_input_emiss_wrf_submit.py
Finalizes emission files for WRF-Chem by adjusting metadata, coordinate dimensions, and time variables.

## WRF-Chem Experiment Preparation

#### remove_fires_data_submit.py
Generates modified emission inputs by removing biomass-burning contributions, enabling sensitivity simulations without fire emissions. 

#### DL_emi_submit.inp
Configuration file used to prepare emissions for model runs.

#### namelist.input_submit
WRF-Chem simulation configuration.

### Model Output Processing
Scripts used to extract and process model outputs:

#### Extract_vars_submit.py
Extracts key variables such as PM₂.₅ and ozone from WRF-Chem output files and saves them as standalone datasets for analysis. 

#### Calculations_means_submit.ipynb
Computes spatial and temporal averages of pollutant concentrations.

#### Contribution_to_PM_results_submit.ipynb
Calculates sectoral contributions to PM₂.₅ concentrations from the sensitivity simulations.

#### Model_validation_MPI-C_submit.ipynb
Evaluates model performance against observations and reference datasets.

#### Test_substractions_submit.ipynb
Diagnostic notebook used to verify emission subtraction and sector attribution.
