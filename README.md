# code_TRAC
This repository includes multiple R files that contain the functions for implementing all the point and variance estimation formulas presented in the paper [Tian Z, Chinchilli VM, Zhou S. TRAC: Robust Assessment of Inter-Rater Agreement and Reliability. Under Review.] For questions or comments about the code, please email Zizhong Tian at <zqt5121@psu.edu>.

functions_2raters.R = Functions for implementing the point estimation and interval estimation (SA, LRT, GOF methods) of TRAC measures and functions for the estimation of kappa, Phi, S, and AC1 methods used in the simulation experiments.

functions_Rraters.R = Functions for implementing the point and variance estimation of overall TRAC in the R-rater case and functions for calculating the true estimands of Conger's kappa and Gwet's AC1 when R=4.

examples_2raters.R = Source "functions_2raters.R" and run illustrative examples about the estimation of TRAC and other popular IRA methods in the 2x2 case.

examples_Rraters.R = Source "functions_Rraters.R" and illustrate the estimations of overall TRAC, Conger's kappa, and Gwet's AC1 in an R=4 example.
