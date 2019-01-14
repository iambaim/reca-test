#include <R.h>
#include <Rdefines.h>

int makedata_COST(SEXP i_COSTList, Data_orig **o_D_orig, Data_COST **o_D_COST);
int re_makedata_COST(Data_orig **i_D_orig, Data_COST **i_D_COST);
int init_cens_COST(Data_orig *i_D_orig, Data_COST *i_D_COST);
int write_input_model1_COST(Data_orig *i_D_orig, Data_COST *i_D_COST);
int write_it_COST(FILE *fp,Data_COST *i_D_COST);
int read_it_COST(int i_it, Data_COST *i_D_COST);
int makedata_COST_predict(int i_nHaul,Data_COST **o_D_COST);
int sample_lambda_init_COST(Data_age *i_D_age,Data_COST *i_D_COST);
int sample_lambda_prior_COST(Data_COST *i_D_COST);
int resample_data_COST(Data_orig *i_D_orig, Data_COST *i_D_COST);
int make_suff_lga_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, Data_lin *i_D_lga,
		       Data_g_a *i_D_g_a, int i_start_h, int i_use_sim_age);
int suff_g_a_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, LW_struct *i_length, Data_age *i_D_age, 
		  Data_lin *i_D_lga, Data_g_a *i_D_g_a, int i_start_h, int i_nHaul, double **o_suff);
