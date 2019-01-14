int init_evaluate(int i_nHaul,int i_ncat);
int re_init_evaluate(int i_ncat);
int calc_lik_age(Age_struct *i_age,Data_age *i_D_age, double *o_loglik);
int Bayes_CV_model1(int i_it,Age_struct *i_age,Data_age *i_D_age,
		    LW_struct *i_lga,Data_lin *i_D_lga,
		    double *o_mean_inv_lik_mod1);
int Bayes_CV_age(int i_it,Age_struct *i_age,Data_age *i_D_age,
                 double *o_mean_inv_loglik_age);
double calc_lik_age_h(int i_h,Age_struct *i_age,Data_age *i_D_age);
int calc_lik_lin(LW_struct *i_lin,Data_lin *i_D_lin,double *o_loglik);
int Bayes_CV_lin(int i_it,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *o_mean_inv_loglik_lin);
double calc_lik_lin_h(int i_h,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *w_res);
int calc_resid_lga(Data_orig *i_D_orig,
                   Age_struct *i_age, Data_age *i_D_age,LW_struct *i_lin,Data_lin *i_D_lin,
		   Data_g_a *i_D_g_a,double *o_resid);
int calc_resid_wgl(double *i_totlength, double *i_totweight,int *i_nFishBoat,
                   LW_struct *i_weight,Data_lin *i_D_wgl,
		   double *o_resid);

