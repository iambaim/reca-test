int initialize_model1(int i_seed, int *i_num_par,
		      Data_age *i_D_age, Age_struct **o_age, Age_struct **o_age_mean,
		      Data_lin *i_D_lga, LW_struct **o_length, LW_struct **o_length_mean,
		      Data_lin *i_D_lga_CC, LW_struct **o_length_CC, LW_struct **o_length_CC_mean,
		      Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, 
		      int i_age_errors, double *i_A2A, double i_delta_age,
		      int i_coastal_cod,int i_inc_hsz,int i_sim_ar,
		      double *i_pri_age_eff_mean, double *i_pri_age_eff_prec,
		      double *i_pri_age_prec_par, double *i_pri_age_ar,
		      double *i_pri_lga_eff_mean, double *i_pri_lga_eff_prec,
		      double *i_pri_lga_prec_par, double *i_pri_lga_ar,
		      int i_lga_fixed_model, double *i_lga_fixed_int, 
		      double *i_lga_fixed_slp, double *i_lga_fixed_tau,
		      double *i_lga_fixed_g_a_c,double *i_lga_fixed_g_a_theta, double *i_lga_fixed_g_a_gamma,
		      double *i_lga_fixed_int_CC, double *i_lga_fixed_slp_CC, double *i_lga_fixed_tau_CC,
		      double *i_lga_fixed_g_a_c_CC,double *i_lga_fixed_g_a_theta_CC, double *i_lga_fixed_g_a_gamma_CC,
		      int i_lga_cens_model, double *i_lga_cens_par);
int re_initialize_model1(Data_age *i_D_age, Age_struct *i_age, Age_struct *i_age_mean, 
			 Data_lin *i_D_lga, LW_struct *i_length, LW_struct *i_length_mean,
			 Data_lin *i_D_lga_CC, LW_struct *i_length_CC, LW_struct *i_length_CC_mean,
			 Data_g_a *i_D_g_a, int i_coastal_cod, int i_age_errors, double *i_A2A, int i_inc_hsz);
int initialize_model2(int i_seed, int *i_num_par,
		      Data_lin *i_D_wgl, LW_struct **o_weight, LW_struct **o_weight_mean,
		      Data_lin *i_D_wgl_CC, LW_struct **o_weight_CC, LW_struct **o_weight_CC_mean,
		      int i_sim_ar,
		      double *i_pri_wgl_eff_mean,double *i_pri_wgl_eff_prec,
		      double *i_pri_wgl_prec_par,double *i_pri_wgl_ar,
		      int i_wgl_fixed_model, double *i_wgl_fixed_int, 
		      double *i_wgl_fixed_slp, double *i_wgl_fixed_tau, 
		      double *i_wgl_fixed_int_CC, double *i_wgl_fixed_slp_CC, 
		      double *i_wgl_fixed_tau_CC, int i_coastal_cod);
int re_initialize_model2(Data_lin *i_D_wgl, LW_struct *i_weight, LW_struct *i_weight_mean, 
			 Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, LW_struct *i_weight_CC_mean, 
			 int i_coastal_cod);
int initialize_predict(Data_age *i_D_age, Age_struct **o_age, 
		       Data_lin *i_D_lga, LW_struct **o_length, Data_g_a *i_D_g_a,
		       Data_lin *i_D_lga_CC, LW_struct **o_length_CC, Data_g_a *i_D_g_a_CC,
		       Data_lin *i_D_wgl, LW_struct **o_weight, Data_lin *i_D_wgl_CC, LW_struct **o_weight_CC,
		       Data_lin *i_D_hsz, LW_struct **o_hsz,
		       int i_lga_cens_model, double *i_lga_cens_par, int i_coastal_cod);
int re_initialize_predict(Data_age *i_D_age, Age_struct *i_age, 
			  Data_lin *i_D_lga, LW_struct *i_length,
			  Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a, 
			  Data_lin *i_D_wgl, LW_struct *i_weight, 
			  Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod,
			  Data_lin *i_D_hsz, LW_struct *i_hsz);
int init_graph_model1(Data_age *i_D_age, Age_struct *i_age, Data_lin *i_D_lga, LW_struct *i_length,
		      Data_lin *i_D_lga_CC, LW_struct *i_length_CC, int i_coastal_cod, int i_constr);
int re_init_graph_model1(Data_age *i_D_age, Age_struct *i_age, LW_struct *i_length, 
			 LW_struct *i_length_CC, int i_coastal_cod);
int init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_coastal_cod,int i_inc_hsz,int i_sim_ar,
	     double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
	     double *i_pri_age_prec_par,double *i_pri_age_ar,
             Age_struct **o_age,Age_struct **o_age_mean);
int re_init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_inc_hsz,
		Age_struct **o_age,Age_struct **o_age_mean);
int alloc_age(Data_age *i_D_age,Age_struct **o_age);
int re_alloc_age(Data_age *i_D_age,Age_struct **o_age);
int init_lin(Data_lin *i_D_lin,int i_sim_ar,double *i_pri_lin_eff_mean,double *i_pri_lin_eff_prec,
	     double *i_pri_lin_prec_par,double *i_pri_lin_ar,int i_fixed_model,
	     LW_struct **o_lin,LW_struct **o_lin_mean);
int re_init_lin(Data_lin *i_D_lin,LW_struct **o_lin,LW_struct **o_lin_mean);
int alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin);
int re_alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin);
int init_lin_par(LW_struct *i_lin,Data_lin *i_D_lin);
