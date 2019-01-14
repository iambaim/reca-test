int alloc_objects_age_lga(Input_common **o_inCommon, Data_orig **o_D_orig, Input_age **o_inAge,
			  Input_lga **o_inLga, Data_CC **o_D_CC, Input_wgl **o_inHsz, Input_prior **o_inPrior);
int re_alloc_objects_age_lga(Input_common **o_inCommon,Data_orig **o_D_orig, Input_age **o_inAge,
			     Input_lga **o_inLga, Data_CC **o_D_CC, Input_wgl **o_inHsz, Input_prior **o_inPrior);
int alloc_objects_wgl(Input_common **o_inCommon, Data_orig **o_D_orig, Input_wgl **o_inWgl, 
		      Data_CC **o_D_CC, Input_prior **o_inPrior);
int re_alloc_objects_wgl(Input_common **o_inCommon,Data_orig **o_D_orig, 
			 Input_wgl **o_inWgl, Data_CC **o_D_CC, Input_prior **o_inPrior);
int readdata_common_ascii(char *i_filename, Input_common *i_inCommon, Data_orig *i_D_orig, Data_CC *i_D_CC);
int readdata_common(char *i_filename, Input_common *i_inCommon, Data_orig *i_D_orig, Data_CC *i_D_CC);
int readdata_orig(char *i_filename, Data_orig *i_D_orig);
int readdata_age_lga(char *i_filename, Data_orig *i_D_orig, Input_age *i_inAge, Input_lga *i_inLga);
int readdata_lga(char *i_filename, Data_orig *i_D_orig, Input_lga *i_inLga);
int readdata_prior(char *i_filename, Input_prior *i_inPrior);
int readdata_wgl(char *i_filename, Data_orig *i_D_orig, Input_wgl *i_inWgl);
int readdata_hsz(char *i_filename, Data_orig *i_D_orig, Input_wgl *i_inHsz);
int readdata_predict_ascii(char *i_filename, Input_predict **o_inPredict);
int readdata_predict(char *i_filename, Input_predict **o_inPredict);
int re_readdata_predict(Input_predict **o_inPredict);
int readdata_catch(char *i_filename, Input_totcatch **o_inCatch);
int re_readdata_catch(Input_totcatch **o_inCatch);
int readdata_cell(char *i_filename, Input_cell **o_inCell);
int re_readdata_cell(Input_cell **o_inCell);
int add_object_info_age_lga(Input_common *i_inCommon, Data_orig *i_D_orig, Input_age *i_inAge,
			    Input_lga *i_inLga, Input_wgl *i_inHsz);
int add_object_info_wgl(Input_common *i_inCommon, Data_orig *i_D_orig, Input_wgl *i_inWgl);
