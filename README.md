This repository contains the code to reproduce the model described in "A Hierarchical Model of Persistent and Transient Growth Variation Applied to Lake Superior Lake Trout" (Stebbins, Bence, Brenden,Hansen).

All R code for modeling is in the directory `R`. The scripts follow the workflow of the model:

* `evaluate_alternate_models.R` fits the growth model with every possible correlation structure of the von Bertalanffy parameters
* `best_fitting_model.R` runs the best fitting model and saves the parameter values.
* `persistent_transient_cases.R` simulates data with different levels of persistent and transient error and fits them to the best model
* `biphasic_sexspecific_models.R` fits a sex-specific version of the von Bertalanffy growth model and a biphasic growth model
* `simulation_code_for_HPCC.R` is the simulation code sent to the high power computing center (HPCC) for simulation testing
* `process_data.R` reads in `txt` files with original data and processes into .Rdata format
# TODO
* `backcalculation.R` documents the backcalculation methods

All TMB code (written in C++) is in the directory `TMB`.

All R code for processing the output of models is in the `analysis` directory.

Helper functions, including plotting scripts, are in the `helper` directory.

Output from models are in the `data-raw` directory.

If you are interested in the data files used in this model, please contact Michael J Hansen.

# TODO
1 - set up repo
2 - test full workflow (run models and fit)
3 - add plotting
4 - add analysis