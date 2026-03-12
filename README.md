# Reef fish escape responses selectively match predator attack speeds

Code and data for the manuscript:

> Neven et al. *Reef fish escape responses selectively match predator attack speeds*. 

---

## Repository structure

```
├── data/
│   ├── filtered_observations.csv      		# Main analysis dataset (see below)
│   ├── all_observations.csv           		# All fish observations before filtering
│   ├── predator_attack_speeds.csv     		# Bar Jack (Caranx ruber) attack speeds
│   └── usable_videos.csv              		# Video metadata
│
├── notebooks/
│   ├── 00_metadata.ipynb              		# Summary counts reported in the Methods
│   ├── 01_video_processing.ipynb      		# Raw video data → filtered_observations.csv *
│   └── 09_figures.ipynb               		# Generate all paper figures
│
├── scripts/
│   ├── 02_random_variable.R           		# Test whether Dep_ID should be a random effect
│   ├── 03_univariable_analysis.R      		# Univariable GLMMs per predictor
│   ├── 04_multivariable_analysis.R    		# Multivariable GLMMs (alternative social context variables)
│   ├── 05_nonlinear_analysis_compare.R		# Compare linear / log / natural spline transformations
│   ├── 06_nonlinear_analysis_effect_sizes.R 	# Final model and effect sizes
│   ├── 07_test_species.R              		# Test species differences in escape response
│   └── 08_test_relationships_ns_lsp.R 		# Test two-way interactions (loom speed × context)
│
└── outputs/
    ├── S1_RandomEffect_Test.csv
    ├── S2_UnivariableNormalized.csv
    ├── S3_MultivariableNormalized_*.csv
    ├── S4_BestNonlinearForms_Context.csv
    ├── S5_Multivariable_WithNonLinear.csv
    ├── S6_Species_Comparison_Results.csv
    ├── S7_All_Interaction_LRT_Summary.csv
    ├── S8_CorrelationMatrix.csv
    └── figures/                        	# Figure 1–3 and S1–S4 (PNG, PDF, TIFF)
```

*`01_video_processing.ipynb` requires raw data not included in this repository (see note in the notebook).

---

## Reproducing the analysis

The dataset needed to reproduce the statistical analysis and figures is already included in `data/filtered_observations.csv`. Steps:

1. **Run R scripts 02–08** in order to reproduce all statistical results in `outputs/`.
2. **Run `notebooks/09_figures.ipynb`** to reproduce all figures in `outputs/figures/`.
3. **Run `notebooks/00_metadata.ipynb`** to reproduce summary counts reported in the Methods.

`notebooks/01_video_processing.ipynb` documents how `filtered_observations.csv` was derived from the raw stereo-video recordings. It cannot be run without the raw data, but is included for transparency.

---

## Analysis dataset

`data/filtered_observations.csv` contains one row per fish observation. It includes Brown Chromis (*Chromis multilineata*) and Bicolor Damselfish (*Stegastes partitus*) with a clear line of sight to the loom stimulus (viewing angle < 55°) and a measured body size < 500 mm. Response categories are first responders (FR) and non-responding fish that observed the stimulus in videos without other responding fish (NRL-NR).

---

## Software

**Python 3.13.6** (notebooks `00`, `01`, `09`):

```
numpy, pandas, opencv-python, scipy, matplotlib, seaborn, statsmodels
```

**R 4.5.1** (scripts `02`–`08`):

```
lme4, dplyr, readr, MuMIn, splines, tidyr, car, broom.mixed
```
