# TVP IV-SVAR
## Replication files for "Nonparametric Time Varying IV-SVARs: Estimation and Inference" (Braun, Kapetanios, Marcellino 2025)
 
The code has been tested using Matlab 2022a and runtimes were recorded based on a standard Windows 10 laptop with an 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz. 
Please report any issues to robin.a.braun"at"frb.gov

The main empirical results (~8 minutes runtime) are replicated in "Empirical_Application.m". On a high-level, it executes the following code:
1) Loading in data for the oil-market model (dataset_oilmarket.xlsx) and industry-level industrial production indices (dataset_IP.csv) 
2) Run a pseudo-out of sample (conditional) forecast exercise as suggested in section (2.5), calling the function "choose_H_CF.m"
3) Estimate the baseline results calling the function "tvp_svar_iv_kernel.m"
4) Estimate a series of models discussed as the robustness exercises: 
-  	internal instrument VAR (calling "tvp_svar_internal_kernel.m")
-   two different bandwidths (calling "tvp_svar_iv_kernel.m")
5) Run industry-level models augmenting the baselinemodel with one industry at a time and create all the relevant Figures 

The code provided in "Minimal_working_example.m" is a compact version of the code and only replicates the results of the main model (Figure 7). The estimated runtime is 10 seconds.

The simulations (Figure 2 - 5) discussed in the main text are run in the file "Monte_carlo_simulations.m" (120 minutes runtime).

The supplementary Monte Carlo Results (Appendix D) are obtained by running the file "Appendix_Monte_Carlo_simulation.m" (30 minutes runtime).

The simple comparison between the kernel method, the path estimator and a Bayesian estimator (Appendix F) is obtained running "Appendix_Comparison_Estimators.m" (90 seconds runtime).



