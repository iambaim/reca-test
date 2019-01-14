int predict(Input_predict *i_inPredict, Input_totcatch *i_inCatch, Input_cell *i_inCell, Data_COST *i_D_COST);
int run_predict(Data_age *i_D_age, Age_struct *i_age, Data_lin *i_D_lga, LW_struct *i_length, 
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC,
		Data_lin *i_D_wgl, LW_struct *i_weight, Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, 
		Data_lin *i_D_hsz, LW_struct *i_hsz, 
		Data_totcatch *i_D_totcatch, TC_struct *i_totcatch, TC_struct *i_totcatch_disc,
		Input_predict *i_inPredict);
int find_catch_at_age(Data_age *i_D_age, Age_struct *i_age, 
		      Data_lin *i_D_lga, LW_struct *i_length, Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		      Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC,
		      Data_lin *i_D_wgl, LW_struct *i_weight, Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
		      Data_lin *i_D_hsz, LW_struct *i_hsz,
		      Data_totcatch *i_D_totcatch,
		      int *i_nMC, int iter, int i_coastal_cod,
		      TC_struct *i_totcatch, double *i_mean_l, double *i_mean_w);
int calculate_means_age(double *mu_age,int c,Age_struct *i_age,Data_age *i_D_age,
			Data_totcatch *i_D_totcatch);
int calculate_means_lga(double *mu_lga,int c,LW_struct *i_length,Data_lin *i_D_lga,
			Data_totcatch *i_D_totcatch);
int calculate_means_wgl(double *mu_wgl,int c,LW_struct *i_weight,Data_lin *i_D_wgl,
			Data_totcatch *i_D_totcatch);
int calculate_means_hsz(double *mu_hsz,int c,LW_struct *i_hsz,Data_lin *i_D_hsz,
			Data_totcatch *i_D_totcatch);

