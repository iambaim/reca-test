int MCMC_model1_init(Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age,
		     Data_lin *i_D_lga, LW_struct *i_length, 
		     Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		     Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		      Data_COST *i_D_COST,int i_use_debug);
int MCMC_model1(Input_common *i_inCommon,Data_orig *i_D_orig, 
		Data_age *i_D_age, Age_struct *i_age, Age_struct *i_age_mean,
		Data_lin *i_D_lga, LW_struct *i_length, LW_struct *i_length_mean,
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC, LW_struct *i_length_CC_mean,
		Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		double *i_mod1_mean_inv_lik,double *i_age_mean_inv_lik,double *i_lga_mean_inv_lik,
		Data_COST *i_D_COST);
int MCMC_model1_it(int start_h,int i_force_acc,int i_len_only,
		   int i_it,int *o_acc_h,int i_it_tot,
		   Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age, 
		   Data_lin *i_D_lga, LW_struct *i_length, 
		   Data_lin *i_D_lga_CC, LW_struct *i_length_CC, 
		   Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		   Data_COST *i_D_COST, int i_save_sim, int i_write_alpha);
int MCMC_it_age(int start_h,int i_force_acc,int i_len_only,
		int i_it,int i_it_tot,int *o_acc_h,
		Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age,
		Data_lin *i_D_lga, LW_struct *i_length, 
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		Data_COST *i_D_COST,
		int save_sim, int i_nHaul,int i_write_alpha);
int MCMC_it_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h,
		int i_it,int i_nHaul,int i_age_error,int i_coastal_cod);
int MCMC_it_lga_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, Data_lin *i_D_lga, Data_g_a *i_D_g_a,
		     LW_struct *i_length, int i_start_h, int i_it, int i_nHaul, int i_age_error);
int MCMC_it_g_a(Data_orig *i_D_orig,Data_COST *i_D_COST,Data_age *i_D_age,Data_lin *i_D_lga,
		Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h,int i_nHaul,int i_it);
int MCMC_model2(Input_common *i_inCommon,
		Data_lin *i_D_wgl, LW_struct *i_weight, LW_struct *i_weight_mean,
		Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, LW_struct *i_weight_CC_mean,
		int i_coastal_cod,
		double *i_wgl_mean_inv_lik);
int MCMC_model2_it(int start_h, int it_tot, 		
		   Data_lin *i_D_wgl, LW_struct *i_weight, 
		   Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, 
		   int i_coastal_cod);
int MCMC_it_wgl(Data_lin *i_D_wgl,LW_struct *i_weight,int start_h);
int MCMC_it_lga_fixed(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int start_h,int i_it);
