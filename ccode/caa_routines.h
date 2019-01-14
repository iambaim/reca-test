int makedata_age1(int i_nHaul,int i_nAges,int *i_a_vec,
		  Input_cov *i_cov, int *i_num_adj_area,int *i_adj_area,
		  Data_age **o_D_age,int i_fit);
int re_makedata_age1(int *i_num_adj_area,int *i_adj_area,Data_age **o_D_age,int i_fit);
int makedata_age2(Data_orig *i_D_orig,Data_age *i_D_age,int class_error);
int re_makedata_age2(Data_age *i_D_age);
int makedata_g_a(int i_age_ncat, int i_g_a_ncat, int i_nSeason, int *i_a_vec, int i_coastal_cod,
		 int i_g_a_model, double *i_par_init, int i_sample_c, int i_sample_theta,int i_sample_gamma,
		 Data_g_a **o_D_g_a);
int re_makedata_g_a(Data_g_a **o_g_a);
int makedata_lin1(int i_nHaul,double *i_haulweight,
		  Input_cov *i_int_cov,Input_cov *i_slp_cov,
		  int *i_num_adj_area,int *i_adj_area,
		  Data_lin **o_D_lin,int i_fit,int i_inc_hsz);
int re_makedata_lin1(int *i_num_adj_area,int *i_adj_area,Data_lin **o_D_lin,int i_fit);
int makedata_lga_suff(Data_lin *i_D_lga,Data_orig *i_D_orig,Data_g_a *i_D_g_a);
int re_makedata_lga_suff(Data_lin *i_D_lga);
int makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC,Data_orig *i_D_orig,
			 Data_g_a *i_D_g_a,Data_g_a *i_D_g_a_CC,int i_class_error);
int re_makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC);
int makedata_wgl_suff(Data_lin *i_D_lin,Data_orig *i_D_orig);
int re_makedata_wgl_suff(Data_lin *i_D_lin);
int makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC,Data_orig *i_D_orig);
int re_makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC);
int makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
                         Data_lin *i_D_lga,Data_l **o_D_l);
int re_makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
			    Data_lin *i_D_lga,Data_l **o_D_l);
int makedata_mcmc(Input_predict *i_inPredict, int i_nAges, int i_nlint, TC_struct **o_totcatch);
int re_makedata_mcmc(TC_struct **o_totcatch);
int makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,Data_lin *i_D_hsz,
                      Input_totcatch *i_inCatch,int i_inc_hsz,Data_totcatch **o_D_totcatch);
int re_makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,Data_lin *i_D_hsz,
			 Input_totcatch *i_inCatch,int i_inc_hsz,Data_totcatch **o_D_totcatch);
int compd(double *i_x,double *i_y);
int find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node);
int re_find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node);
int alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par);
int re_alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par);
double calc_eff(Data_cov *i_xcov,double **i_eff,int i_h);
double calc_eff_no_haul(Data_cov *i_xcov,double **i_eff,int i_h);
double calc_eff_suff_age(Data_cov *i_xcov,Eff_str *i_par,int i_h,int i_a);
double ldnorm(double x,double mu,double sigma,double logsigma);
int write_samples_model1(Data_age *i_D_age, Age_struct *i_age,
			 Data_lin *i_D_lga, LW_struct *i_length, Data_g_a *i_D_g_a,
			 Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a_CC, 
			 Data_CC *i_D_CC, int i_coastal_cod, 
			 int i_print_boat, int i_print_format, Data_COST *i_D_COST);
int write_samples_model2(Data_lin *i_D_wgl, LW_struct *i_weight,
			 Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
			 int i_coastal_cod, int i_print_boat, int i_print_format);
int write_it(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_print_boat);
int read_it(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_read_boat);
int write_it_ascii(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_print_boat);
int write_samples_haulsize(FILE *fp, Data_glm *i_glm, Eff_str *i_par);
int read_it_g_a(int i_it,Data_g_a *i_D_g_a);
int write_it_totcatch(FILE *fp,int i_ncat,int i_nlen,
		      TC_struct *i_totcatch,double *i_mean_l,double *i_mean_w);
int update_average_age(int i_n,Age_struct *i_age,Data_age *i_D_age,
                       Age_struct *x_age_mean);
int update_average_g_a(int i_n,Data_g_a *i_D_g_a);
int update_average_lin(int i_n,LW_struct *i_lin,Data_lin *i_D_lin,
                       LW_struct *x_lin_mean);
int update_mean(double *i_mean,double i_x,int i_n);
double scale_proposal(double x, double f, double *la);
void my_genmul(long n,double *p,long ncat,long *ix);
