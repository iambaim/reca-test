int write_input_model1(FILE *fp, Data_orig *i_D_orig, Input_common *i_inCommon, Input_age *i_inAge,
		       Input_lga *i_inLga, Input_prior *i_inPrior, Data_CC *i_D_CC);
int write_input_model2(FILE *fp, Data_orig *i_D_orig, Input_common *i_inCommon, Input_wgl *i_inWgl,
		       Input_prior *i_inPrior);
int write_input_predict(FILE *fp, Input_predict *i_inPredict, Input_age *i_inAge, Input_lga *i_inLga, 
			Input_wgl *i_inWgl, Input_totcatch *i_inCatch, Input_cell *i_inCell);
int write_inData(FILE *fp, Data_orig *i_D_orig, int i_inc_hsz);
int write_inCommon(FILE *fp, Input_common *i_inCommon);
int write_inAge(FILE *fp, Input_age *i_inAge, int i_nHaul, int i_pred);
int write_inLga(FILE *fp, Input_lga *i_inLga, int i_age_ncat, int i_nHaul, int i_pred, int i_it);
int write_inPrior(FILE *fp, Input_prior *i_inPrior, Input_age *i_inAge, Input_lga *i_inLga, 
		  Input_wgl *i_inWgl, int i_write_age, int i_write_lga, int i_write_wgl);
int write_inCC(FILE *fp, Data_CC *i_D_CC, int i_ncat);
int write_inWgl(FILE *fp, Input_wgl *i_inWgl, int i_nHaul, int i_pred, int i_it);
int write_inPredict(FILE *fp, Input_predict *i_inPredict, int i_nCell);
int write_inCell(FILE *fp, Input_cell *i_inCell, int i_cc);
int write_inCatch(FILE *fp, Input_totcatch *i_inCatch, int *i_age_n_cov, int *i_lga_n_cov, int *i_wgl_n_cov);
int write_input_cov(FILE *fp, Input_cov *i_inCov, int i_nAges, int i_nHaul, int i_pred);

int read_inData(FILE *fp, Data_orig **o_D_orig);
int re_read_inData(Data_orig **o_D_orig);
int read_inCommon(FILE *fp, Input_common **o_inCommon);
int re_read_inCommon(Input_common **o_inCommon);
