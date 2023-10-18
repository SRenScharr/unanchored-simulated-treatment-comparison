Simulation study code

The comprehensive simulation study of unanchored simulated treatment comparison appraoch developed in Ren et al. (2023) with 120 scenarios can be reproduced by firstly run Run.R to get the estimated treatment effects and then run Results.R to evaluate the estimated treatment effects.

1. All_functions: the functions needed in estimation
2. Model_input_SX_BB.csv: model input for the simulations and "X" indicates simulation scenarios
3. Output_functions: functions used to calculate bias, model SE, empirical SE, coverage
4. Results.R: R code to evaluate model results using bisa, model SE, empirical SE, coverage
5. Run.R: run the simulations with results saved in Output folder (mean in first column and SD in second column)
6. True correlation.R: calculate the Pearson correlation of binary varaibles simulated using the NORTA algorithm
7. True_trt_cal.R: calculate the true treament effect numerically, saved in Results folder


Example analysis code

Example analysis.R provides an example analysis code of implimenting the unanchored simulated treatment comparison approach using a simulated dataset. 
- The simulated dataset takes one Monte-Carlo repitition of scenario 1 as an example. 
- The data used in the analysis are saved in A_AgDSummary.csv and B_IPD.csv.