int make_cell_constr_age(Data_age *i_D_age, Input_cov *i_cov);
int re_make_cell_constr(Data_glm *i_glm);
int make_cell_constr_lin(Data_lin *i_D_lin, Input_cov *i_int_cov, Input_cov *i_slp_cov);
int makedata_cell_dist_age(int i_num_cell_o, int i_num_cell_u,
			   int i_age_int_nC, double *i_age_int_E, double *i_age_int_C,
			   Age_struct *i_age, Data_age *i_D_age);
int makedata_cell_dist_lin(int *i_num_cell_o, int *i_num_cell_u,
			   int i_lin_int_nC, double *i_lin_int_E, double *i_lin_int_C,
			   int i_lin_slp_nC, double *i_lin_slp_E, double *i_lin_slp_C,
			   LW_struct *i_lin, Data_lin *i_D_lin);
int makedata_cell_dist_hsz(int i_num_cell_o, int i_num_cell_u,
			   int i_hsz_int_nC, double *i_hsz_int_E, double *i_hsz_int_C,
			   LW_struct *i_hsz, Data_lin *i_D_hsz);
int re_makedata_cell_dist_age(Age_struct *i_age,Data_age *i_D_age);
int re_makedata_cell_dist_lin(LW_struct *i_lin,Data_lin *i_D_lin);
int re_makedata_cell_dist_hsz(LW_struct *i_hsz,Data_lin *i_D_hsz);
int simulate_cell_effects(Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl,
			  LW_struct *i_hsz,Data_lin *i_D_hsz);
int simulate_cell_effects_CC(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_hsz,Data_lin *i_D_hsz,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC);
