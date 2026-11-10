##########################################################################
# "Identifying News Shocks from Forecasts"
#	Jonathan Adams & Philip Barrett, JAE 2026
# 	Replication manual. 
##########################################################################


##### Overview ####
This file outlines how to reproduce the main figure & tables in the above-mentioned paper.


#### Code ####
The main codes are listed below in the order in which they should be run:
1. robust_iterate.R		Creates estimates for the baseline and most robustness checks.
2. monetary_validation.R	Computes the monetary policy shock validation exercise.
3. var_decomp_compare.R		Compares the variance decomposition across specifications.
4. validate_via_sim_2025.R	Computes the variance decomposition for simulated data
In all files, you will need to set the output folder name, which currently defaults to "final_replication", as well as the working directory for the code.


#### Data ####
The main input data file are 
	data/final_replication_data.csv. 	The main empirical data file
	data/monetary_compare			A folder containing various monetary shocks
	data/simulated				A folder containing simulated model data
These are all read up by the code automatically and should not need changing.


#### Outputs ####
Outputs are typically saved in figures/[output_folder_name]/[specification_name], where the specification names are those looped over in robust_iterate.R.


##### Mapping of figures/tables to code ####

## Figures
Figure	Output file					Created in		Section	Line
2	irfs_surp_all_pc.pdf				robust_iterate.R	5.1	 720
3	mult_verify_epsilon_wide.pdf			robust_iterate.R	5.3	 851
4	monetary_compare_var_cum_1_lags_wide.pdf	monetary_validation.R	6.3	 355
5	irfs_coef_pol_v_scale_wide.pdf			robust_iterate.R	5.1b	 760
6	irfs_coef_dem_v_scale_wide			robust_iterate.R	5.1b	 744
	irfs_coef_sup_v_scale_wide			robust_iterate.R	5.1b	 752
7	consol_cf_compare_wide.pdf			robust_iterate.R	5.4	 963
8	mon_cf_compare_wide.pdf				robust_iterate.R	5.4	 950
10	var_decomp_compare_wide.pdf			var_decomp_compare.R	2	 101
11	aic.pdf						robust_iterate.R	3.1	 297
12	mult_verify_all_epsilon.pdf			robust_iterate.R	5.3	 870
13	monetary_compare_var_init_1_lags_wide.pdf	monetary_validation.R	6.2	 313
	monetary_compare_var_cum_4_lags_wide.pdf	monetary_validation.R	6.3	 357
14	data_cleaned_wide.pdf				robust_iterate.R	5.6	1052
15	shocks.pdf					robust_iterate.R	5.7	1071
16	irf_epsilon_decomp_pc50.pdf			robust_iterate.R	5.2	 821
17	cf_consol_y_wide.pdf				robust_iterate.R	5.4	 902
18	passive_cf_compare_wide.pdf			robust_iterate.R	5.4	 976
19	var_decomp_NK3_monetary_news_6shocks		validate_via_sim_2025.R 5.4	 429

## Tables
Table	Output file					Created in		Section	Line
3	var_decomp_pc50.tex				robust_iterate		6.1	1163
4	cf_irf_norm.tex					robust_iterate		5.4	1004
5	lit_multipliers_wide.tex			robust_iterate		3.8	579
