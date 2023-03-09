# Bayesian_ITC
Notebook for Bayesian analysis of isothermal titration calorimetry using a two-step binding model

## Manifest
The project includes several jupyter notebooks  for modeling of ITC data. A well-documented example notebook, `bayesian_itc_notebook.ipynb` as well as the module `itcfunctions.py` containing essential functions can be found in `the example_notebook/` folder. Additional models used in the work can be found in the model_notebooks folder. We have additionally included the twelve isotherms that were used in Estelle et al. in the isotherm_data, in both the un-processed form of .itc files that can be read by Origin or Nitpic, and csv files of integrated heats, injection volumes and initial concentrations.

## Requirements
- Python 3.4+
- Jupyter Notebook or Lab
- Numpy
- Scipy
- EMCEE
- Corner

## Usage

The example notebook contains everything needed to model calorimetry data two a two-step binding model, and can be run by running each notebook cell. Additional notebooks for other models can be found in the `model_notebooks\` folder.

### Setup options and parameters
The following variables are adjustable at the top of the notebook:

`conc_priors` - Boolean for whether concentrations are considered as model parameters. When set to False, model will run with fixed concentrations. When set to True, additional model parameters for each concentration of titrant and titrand will be included. \
\
`skipped_injections` - integer defining the number of injections ignored in modeling. Set to 1 by default to account for the first injection anomaly. \
\
`filename` - sets the name of the file that the run is saved to using EMCEE's built-in backend. \
\
`pt initial` and `lt initial` - Starting concentrations in the syringe and cell respectively. When modeling synthetic data, these are set up to be taken automatically from the `get_synthetic_itc()` function. When working with experimental data these need to be set manually. \
\
`n_walkers` - Number of individual MCMC chains. Should be at least 3x the number of parameters per EMCEE docs. \
\
`n_steps` - Number of steps taken by each walker. \
\
`pt_range` and `lt_range` - width of uniform prior for syringe and cell concentration respectively. The number listed represents a fraction of the stated value, e.g. 0.1 corresponds to a range of ±10%. 

### Synthetic Data
By default, the notebook is set up to run on synthetic data built in the `get_synthetic_itc()` function within `itcfunctions.py` Model parameters must be changed from within this function, which is called at the start of the notebook. 

### Experimental Data
For modeling experimental data, the line calling `get_synthetic_itc()` in the notebook should be commented out, and the line below calling `get_data()` should be uncommented. Experimental data should take the form of a two-column CSV file, with the integrated heat per injection on the left column, and the injection volume in the right. The units for each are microcalories and microliters respectively. The notebook will also need to be supplied with initial concentrations. For our published isotherms, these can be found in the 'integrated_heats_for_notebook' folder in the repository.
