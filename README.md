# metalearner_shap_surv
This repository contains key functions for paper titled "Evaluating meta-learners to analyze treatment heterogeneity in survival data: application to electronic health records of pediatric asthma care in COVID-19 pandemic".

## “simulated_data” folder contains exmple data for simulations:
### (1) balanced_setting1_traindata1.RData is training data for one simulation run under balanced observational study in scenario 1.
### (2) observational_10p1000n10000test_balanced_scenario1.RData is the test data for balanced observational study under scenario 1.
### (3) observational_10p1000n100test_kernelshap_balanced_scenario1.RData is a subsampled test data (to save computational load) for calculating SHAPley values on test data. 

## “simulation_functions” folder contains functions for simulations:
### (1) meta_learner_simulation_censor_survival_prob_kernelshap_ipcw.R: contains six meta-learner functions with SHAPley value computation.
### (2) weibull_true.R: contains benchmark parametric Weibull regressions models.
### (3) simulation_1run_observational_design_balanced.R: contains code examples for 1 simulation run.
