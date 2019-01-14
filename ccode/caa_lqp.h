int lqp_sample(int i_m,
	       double (*i_log_f)(double,double *),
	       double (*i_log_f1)(double,double *),
	       double (*i_log_f2)(double,double *),
               double *i_par_old,double i_x_old,double *i_par_new,
               double *w_knots,double **w_coef,double *w_prob,
               double *o_x_new,double *o_lacc);

