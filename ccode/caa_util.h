double lm_fit(int i_N,double *i_x,double *i_y,
	      double *o_beta0,double *o_beta1,double *o_ssq);
double lm_fit_suff(int i_ncat,double *i_g_a,double *i_sum_by_cat,double *i_sqsum_by_cat,double **i_suff,
		   double *o_beta0,double *o_beta1,double *o_ssq);
double dt(double x,double mu,double sigma,double df);
double dnorm(double x,double mu,double sigma);
double pnorm0(double x);
double pnorm(double x,double mu,double sigma);
