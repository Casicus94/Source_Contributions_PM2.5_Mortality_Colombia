# Source_Contributions_PM2.5_Mortality_Colombia
Analysis scripts supporting the study of PM$_{2.5}$ source contributions and pollution-related mortality in Colombia using WRF-Chem simulations and sectoral emission sensitivity experiments. This paper has been accepted in GeoHealth

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


## Revision analysis (`revision_analysis/`)
Scripts added during the GeoHealth revision. Run from the project folder (paths are relative to it, or pass `--root`).

| Script | Purpose |
|---|---|
| `make_table1.py` | Table 1: mortality contributions by sector and department (`--method brendan` reproduces the submitted table, `--method corrected` the revised one). |
| `make_fig2.py`, `make_fig2d.py` | Figure 2: WRF-Chem PM₂.₅ evaluation against ACAG V6 and ground stations. |
| `make_figS1.py`, `make_figS2_S3.py` | Figures S1–S3: ACAG vs stations, and monthly evaluation. |
| `make_fig5_leading_sector.py` | Figure 5: leading sector of PM₂.₅-attributable mortality per department. |
| `Domains_figure_with_inset.py` | Figure 1: model domains with the Aburrá Valley inset. |
| `AWB_vs_FINN_check.py`, `AWB_FINN_overlap_and_masking.py`, `VIIRS_fire_climatology_check.py` | Checks on double counting between EDGAR agricultural waste burning and FINN fires, and on fire climatology (VIIRS). |
| `sensitivity_bias_R3M3.py`, `common.py` | Sensitivity of national sectoral mortality shares to model bias. |

`namelist.input_submit` was updated to the final configuration (`bio_emiss_opt = 3`).
