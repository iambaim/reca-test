/*!
  \file caa_input.c
  \brief Containing routines for reading and writing different input data

  Only to be used for testing.
  \author Hanne Rognebakke
*/
#include "caa.h"
#include "caa_input.h"
#include "caa_read_write.h"

    

/*!
  \brief Writes input for model1 to file caa_input_model1.txt in the current directory. 
  \author Hanne Rognebakke
*/
int write_input_model1(FILE *fp, Data_orig *i_D_orig, Input_common *i_inCommon, Input_age *i_inAge,
		       Input_lga *i_inLga, Input_prior *i_inPrior, Data_CC *i_D_CC)
{
  int      err,it;

  fprintf(stderr,"write_input_model1: write_inData\n");
  err = write_inData(fp,i_D_orig,i_inCommon->inc_hsz);
  if(err)
    {
      write_warning("write_input_model1:Error calling write_inData\n");
      return(err);
    }

  fprintf(stderr,"write_input_model1: write_inCommon\n");
  err = write_inCommon(fp,i_inCommon);
  if(err)
    {
      write_warning("write_input_model1:Error calling write_inCommon\n");
      return(err);
    }

  fprintf(stderr,"write_input_model1: write_inAge\n");
  err = write_inAge(fp,i_inAge,i_D_orig->nHaul,0);
  if(err)
    {
      write_warning("write_input_model1:Error calling write_inAge\n");
      return(err);
    }

  fprintf(stderr,"write_input_model1: write_inLga\n");
  it = i_inCommon->burn_in + i_inCommon->num_it_inner*i_inCommon->num_it_outer;
  err = write_inLga(fp,i_inLga,i_inAge->nAges,i_D_orig->nHaul,0,it);
  if(err)
    {
      write_warning("write_input_model1:Error calling write_inLga\n");
      return(err);
    }

  fprintf(stderr,"write_input_model1: write_inPrior\n");
  err = write_inPrior(fp,i_inPrior,i_inAge,i_inLga,NULL,1,1,0);
  if(err)
    {
      write_warning("write_input_model1:Error calling write_inPrior\n");
      return(err);
    }

  if(i_D_orig->coastal_cod)
    {
      fprintf(stderr,"write_input_model1: write_inCC\n");
      err = write_inCC(fp,i_D_CC,(int)i_inAge->nAges/2);
      if(err)
	{
	  write_warning("write_input_model1:Error calling write_inCC\n");
	  return(err);
	}
    }
  


  return(0);
}		/* end of write_input_model1 */
    

/*!
  \brief Writes input for model2 to file caa_input_model2.txt in the current directory. 
  \author Hanne Rognebakke
*/
int write_input_model2(FILE *fp, Data_orig *i_D_orig, Input_common *i_inCommon, Input_wgl *i_inWgl,
		       Input_prior *i_inPrior)
{
  int      err,it;

  /* Common parameters */
  fprintf(fp,"mcmc_par=%d %d %d\n",i_inCommon->burn_in,i_inCommon->num_it_inner,
	  i_inCommon->num_it_outer);
  fprintf(fp,"seed=%d\n",i_inCommon->seed);
  fprintf(fp,"num_par=%d %d\n",i_inCommon->num_par2[0],i_inCommon->num_par2[1]);

  it = i_inCommon->burn_in + i_inCommon->num_it_inner*i_inCommon->num_it_outer;
  err = write_inWgl(fp,i_inWgl,i_D_orig->nHaul,0,it);
  if(err)
    {
      write_warning("write_input_model2:Error calling write_inWgl\n");
      return(err);
    }
  err = write_inPrior(fp,i_inPrior,NULL,NULL,i_inWgl,0,0,1);
  if(err)
    {
      write_warning("write_input_model2:Error calling write_inPrior\n");
      return(err);
    }

  return(0);
}		/* end of write_input_model2 */

    
/*!
  \brief Writes input for predict to file caa_input_predict.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik, Hanne Rognebakke
*/
int write_input_predict(FILE *fp, Input_predict *i_inPredict, Input_age *i_inAge, Input_lga *i_inLga, 
			Input_wgl *i_inWgl, Input_totcatch *i_inCatch, Input_cell *i_inCell)
{
  int  err=0;

  err = write_inPredict(fp,i_inPredict,i_inCatch->nCell);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inPredict\n");
      return(err);
    }
  err = write_inAge(fp,i_inAge,i_inPredict->nHaul,1);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inAge\n");
      return(err);
    }
  err = write_inLga(fp,i_inLga,i_inAge->nAges,i_inPredict->nHaul,1,0);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inLga\n");
      return(err);
    }
  err = write_inWgl(fp,i_inWgl,i_inPredict->nHaul,1,0);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inWgl\n");
      return(err);
    }
  err = write_inCell(fp,i_inCell,i_inPredict->coastal_cod);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inCell\n");
      return(err);
    }
  fprintf(stderr,"write_input_predict: write_inCatch - NB NB CHANGE!!\n");
  //  err = write_inCatch(fp,i_inCatch,i_inAge->cov->n_cov,i_inLga->int_cov->n_cov,i_inWgl->int_cov->n_cov);
  if(err)
    {
      write_warning("write_input_predict:Error calling write_inCatch\n");
      return(err);
    }


  return(err);
}               /* end of write_input_predict */



/*!
  \brief Write data in struct Data_orig to file
  \author Hanne Rognebakke
*/
int write_inData(FILE *fp, Data_orig *i_D_orig, int i_inc_hsz)
{
  int    i,h,f,cum_fish,ind_f,n;
  char   buffer[MAX_STR];

  fprintf(fp,"nFish %d\n",i_D_orig->nFish);
  fprintf(fp,"nBoat %d\n",i_D_orig->nBoat);
  fprintf(fp,"nHaul %d\n",i_D_orig->nHaul);
  fprintf(fp,"coastal_cod %d\n",i_D_orig->coastal_cod);
  if(i_inc_hsz)
    {
      fprintf(fp,"haul boat nFishBoat start_noAge num_noAge haulweight\n");
      for(i=0;i<i_D_orig->nHaul;i++)
	{
	  fprintf(fp,"%d %d %d %d %d %f\n",i,i_D_orig->boat[i],
		  i_D_orig->nFishBoat[i],i_D_orig->start_noAge[i],
		  i_D_orig->num_noAge[i],i_D_orig->haulweight[i]);
	}
    }
  else
    {
      fprintf(fp,"haul boat nFishBoat start_noAge num_noAge\n");
      for(i=0;i<i_D_orig->nHaul;i++)
	{
	  fprintf(fp,"%d %d %d %d %d\n",i,i_D_orig->boat[i],
		  i_D_orig->nFishBoat[i],i_D_orig->start_noAge[i],
		  i_D_orig->num_noAge[i]);
	}
    }
  for(i=0;i<i_D_orig->nHaul;i++)
    {
      if(i_D_orig->start_noAge[i]+i_D_orig->num_noAge[i]>i_D_orig->nFish)
	{
	  sprintf(buffer,"write_inData:start_noAge[%d](%d)+num_noAge[%d](%d)>nFish(%d) !\n",
		  i,i_D_orig->start_noAge[i],i,i_D_orig->num_noAge[i],i_D_orig->nFish);
	  write_warning(buffer);
	  return(1);
	}
    }
  int nReplength=0;
  int nMissAge=0;
  cum_fish = 0;
  n = 0;
  fprintf(fp,"coastal cod %d\n",i_D_orig->coastal_cod);
  fprintf(fp,"fish haul totage season totlength totweight replength");
  if(i_D_orig->coastal_cod)
    fprintf(fp," tottype\n");
  else
    fprintf(fp,"\n");
  
  for(h=0;h<i_D_orig->nHaul;h++)
    {
      /* print aged fish */
      for(f=0;f<(i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h]);f++)
	{
	  ind_f = cum_fish+f;
	  if(ind_f != n)
	    {
	      sprintf(buffer,"write_inData:Something is wrong: h=%d:ind_f=%d,nFish=%d\n",
		      h,ind_f,n);
	      write_warning(buffer);
	      return(1);
	    }
	  nReplength+= i_D_orig->replength[ind_f];
	  fprintf(fp,"%d %d %d %d %14.12f %d",ind_f,h,
		  i_D_orig->totage[ind_f],i_D_orig->season[ind_f],i_D_orig->totlength[ind_f],
		  i_D_orig->replength[ind_f]);
	  if(i_D_orig->coastal_cod)
	    fprintf(fp," %d\n",i_D_orig->tottype[ind_f]);
	  else
	    fprintf(fp,"\n");
	  n++;		  
	  if(i_D_orig->totage[ind_f]<0)
	    {
	      sprintf(buffer,"write_inData:Something is wrong: h=%d:totage[%d]=%d,num_noAge=%d,start_noAge=%d, age should be observed!\n",
		      h,ind_f,i_D_orig->totage[ind_f],i_D_orig->num_noAge[h],i_D_orig->start_noAge[h]);
	      write_warning(buffer);
	      return(1);
	    }
	}
      /* print fish with missing ages */
      for(f=0;f<i_D_orig->num_noAge[h];f++)
	{
	  ind_f = i_D_orig->start_noAge[h]+f;
	  if(ind_f != n)
	    {
	      sprintf(buffer,"write_inData:Something is wrong: h=%d:ind_f=%d,nFish=%d\n",
		      h,ind_f,n);
	      write_warning(buffer);
	      return(1);
	    }
	  fprintf(fp,"%d %d %d %d %14.12f %d",ind_f,h,
		  i_D_orig->totage[ind_f],i_D_orig->season[ind_f],i_D_orig->totlength[ind_f],
		  i_D_orig->replength[ind_f]);
	  nReplength+= i_D_orig->replength[ind_f];
	  nMissAge+= i_D_orig->replength[ind_f];
	  if(i_D_orig->coastal_cod)
	    fprintf(fp," %d\n",i_D_orig->tottype[ind_f]);
	  else
	    fprintf(fp,"\n");
	  n++;
	      if(i_D_orig->totage[ind_f]>-1)
		{
		  sprintf(buffer,"write_inData:Something is wrong: h=%d:totage[%d]=%d,num_noAge=%d, age should be missing!\n",
			  h,ind_f,i_D_orig->totage[ind_f],i_D_orig->num_noAge[h]);
		  write_warning(buffer);
		  return(1);
		}
	      if(i_D_orig->totlength[ind_f]<0)
		{
		  sprintf(buffer,"write_inData:Something is wrong: h=%d:totlength[%d]=%f, length should be observed!\n",
			  h,ind_f,i_D_orig->totlength[ind_f]);
		  write_warning(buffer);
		  return(1);
		}
	}
      cum_fish += i_D_orig->nFishBoat[h];
    } // end for(h=0;h<i_D_orig->nHaul;h++)
  if(n != i_D_orig->nFish)
    {
      sprintf(buffer,"write_inData:Something is wrong: totlength nFish=%d, input nFish=%d\n",n,i_D_orig->nFish);
      write_warning(buffer);
      return(1);
    }
  fprintf(fp,"nReplength %d\n",nReplength);
  fprintf(fp,"nMissAge %d\n",nMissAge);
  
  fprintf(fp,"n_int_len %d\n",i_D_orig->n_int_len);
  fprintf(fp,"int_len_lim int_len\n");
  for(i=0;i<i_D_orig->n_int_len;i++)
    fprintf(fp,"%f %f\n",i_D_orig->int_len_lim[i],i_D_orig->int_len[i]);

  return(0);
}               /* end of write_inData */


/*!
  \brief Write common parameters to file
  \author Hanne Rognebakke
*/
int write_inCommon(FILE *fp, Input_common *i_inCommon)
{
  int i;

  fprintf(fp,"seed %d\n",i_inCommon->seed);
  fprintf(fp,"burn_in %d\n",i_inCommon->burn_in);
  fprintf(fp,"num_it_inner %d\n",i_inCommon->num_it_inner);
  fprintf(fp,"num_it_outer %d\n",i_inCommon->num_it_outer);
  fprintf(fp,"num_par1");
  for(i=0;i<5;i++)
    fprintf(fp," %d",i_inCommon->num_par1[i]);
  fprintf(fp,"\n");
  fprintf(fp,"num_par2");
  for(i=0;i<2;i++)
    fprintf(fp," %d",i_inCommon->num_par2[i]);
  fprintf(fp,"\n");
  fprintf(fp,"constr %d\n",i_inCommon->constr);
  fprintf(fp,"sim_ar %d\n",i_inCommon->sim_ar);
  fprintf(fp,"use_debug %d\n",i_inCommon->use_debug);
  fprintf(fp,"filename_mcmc1 %s\n",i_inCommon->filename_mcmc1);
  fprintf(fp,"filename_mcmc2 %s\n",i_inCommon->filename_mcmc2);
  fprintf(fp,"filename_hsz_mcmc2 %s\n",i_inCommon->filename_hsz_mcmc2);
  fprintf(fp,"filename_hsz_it %s\n",i_inCommon->filename_hsz_it);
  fprintf(fp,"filename_hsz_hauleff %s\n",i_inCommon->filename_hsz_hauleff);
  fprintf(fp,"print_boat %d\n",i_inCommon->print_boat);
  fprintf(fp,"inc_hsz %d\n",i_inCommon->inc_hsz);

  return(0);
}               /* end of write_common_par */


/*!
  \brief Write age parameters to file
  \author Hanne Rognebakke
*/
int write_inAge(FILE *fp, Input_age *i_inAge, int i_nHaul, int i_pred)
{
  int i,err=0;
  
  /* Write specific age model parameters */
  fprintf(fp,"\n");
  fprintf(fp,"nAges %d\n",i_inAge->nAges);
  fprintf(fp,"avec");
  for(i=0;i<i_inAge->nAges;i++)
    fprintf(fp," %d",i_inAge->a_vec[i]);
  fprintf(fp,"\n");
  fprintf(fp,"errors %d\n",i_inAge->errors);
  if(i_inAge->errors)
    {
      fprintf(fp,"A2A");
      {
	for(i=0;i<(i_inAge->nAges*i_inAge->nAges);i++)
	  fprintf(fp," %f",i_inAge->A2A[i]);
      }
      fprintf(fp,"\n");
    }
  fprintf(fp,"delta_age %f\n",i_inAge->delta_age);
  /* Write covariate information */
  err = write_input_cov(fp, i_inAge->cov, i_inAge->nAges, i_nHaul, i_pred);
  if(err)
    {
      write_warning("write_inAge:Error calling write_input_cov\n");
      return(err);
    }

  return(err);
}               /* end of write_inAge */


/*!
  \brief Write lga parameters to file
  \author Hanne Rognebakke
*/
int write_inLga(FILE *fp, Input_lga *i_inLga, int i_age_ncat, int i_nHaul, int i_pred, int i_it)
{
  int    i;
  int    err=0;

  /* Write specific lga model parameters */
  fprintf(fp,"ga_model %d\n",i_inLga->g_a_model);
  fprintf(fp,"g_a_ncat %d\n",i_inLga->g_a_ncat);
  fprintf(fp,"g_a_nSeason %d\n",i_inLga->g_a_nSeason);
  if(i_inLga->g_a_model==1)
    {
      fprintf(fp,"g_a_par_init");
      for(i=0;i<3;i++)
	fprintf(fp,"%f",i_inLga->g_a_par_init[i]);
      fprintf(fp,"\n");
    }
  fprintf(fp,"fixed_model %d\n",i_inLga->fixed_model);
  if(i_inLga->fixed_model)
    {
      fprintf(fp,"fixed_int ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inLga->fixed_int[i]);
      fprintf(fp,"\nfixed_slp ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inLga->fixed_slp[i]);
      fprintf(fp,"\nfixed_tau ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inLga->fixed_tau[i]);
      fprintf(fp,"\n");
      if(i_inLga->g_a_model>0)
	{
	  fprintf(fp,"fixed_g_a_c ");
	  for(i=0;i<i_it;i++)
	    fprintf(fp,"%f ",i_inLga->fixed_g_a_c[i]);
	  fprintf(fp,"\nfixed_g_a_theta ");
	  for(i=0;i<i_it;i++)
	    fprintf(fp,"%f ",i_inLga->fixed_g_a_theta[i]);
	  fprintf(fp,"fixed_g_a_gamma ");
	  for(i=0;i<i_it;i++)
	    fprintf(fp,"%f ",i_inLga->fixed_g_a_gamma[i]);
	  fprintf(fp,"\n");
	}
    }
  fprintf(fp,"cens_model %d\n",i_inLga->cens_model);
  if(i_inLga->cens_model)
    {
      fprintf(fp,"cens_par");
      for(i=0;i<4;i++)
	fprintf(fp," %f",i_inLga->cens_par[i]);
      fprintf(fp,"\n");
    }

  /* Write covariate information */
  fprintf(fp,"covatiates intercept:\n");
  err = write_input_cov(fp, i_inLga->int_cov, 1, i_nHaul, i_pred);
  if(err)
    {
      printWarning("write_inLga:Error calling write_input_cov\n");
      return(err);
    }
  fprintf(fp,"covatiates slope:\n");
  err = write_input_cov(fp, i_inLga->slp_cov, 1, i_nHaul, i_pred);
  if(err)
    {
      printWarning("write_inLga:Error calling write_input_cov\n");
      return(err);
    }

  return(0);
}               /* end of write_inLga */


/*!
  \brief Write prior parameters to file
  \author Hanne Rognebakke
*/
int write_inPrior(FILE *fp, Input_prior *i_inPrior, Input_age *i_inAge, Input_lga *i_inLga, 
		  Input_wgl *i_inWgl, int i_write_age, int i_write_lga, int i_write_wgl)
{
  int i,j,k,ind,n,nfix,nfix_int,nfix_slp;

  if(i_write_age)
    {
      /* Write prior parameters for age model */
      fprintf(fp,"age_eff_mean");
      nfix = 0;
      ind = 0;
      for(i=0;i<i_inAge->cov->n_cov;i++)
	{
	  if(i_inAge->cov->random[i]==0)
	    {
	      nfix+=1;
	      for(j=0;j<i_inAge->cov->n_lev[i];j++)
		for(k=0;k<i_inAge->nAges;k++)
		  {
		    fprintf(fp," %f",i_inPrior->age_eff_mean[ind]);
		    ind++;
		  }
	    }
	}
      fprintf(fp,"\nNumber of fixed effects:%d\n",nfix);
      
      fprintf(fp,"age_eff_prec");
      for(i=0;i<nfix;i++)
	fprintf(fp," %f",i_inPrior->age_eff_prec[i]);
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inAge->cov->n_cov;i++)
	if(i_inAge->cov->random[i]==1)
	  n += 2;
      fprintf(fp,"age_prec_par");
      for(i=0;i<n;i++){
	fprintf(fp," %f",i_inPrior->age_prec_par[i]);
      }
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inAge->cov->n_cov;i++)
	if(i_inAge->cov->spatial[i]==1)
	  n += 2;
      fprintf(fp,"age_ar ");
      for(i=0;i<n;i++){
	fprintf(fp,"%f ",i_inPrior->age_ar[i]);
      }
      fprintf(fp,"\n");
    }

  if(i_write_lga)
    {
      /* Write prior parameters for lga model */
      n = 0;
      nfix_int = 0;
      for(i=0;i<i_inLga->int_cov->n_cov;i++){
	if(i_inLga->int_cov->random[i]==0){
	  nfix_int+=1;
	  n += i_inLga->int_cov->n_lev[i];
	}
      }
      nfix_slp = 0;
      for(i=0;i<i_inLga->slp_cov->n_cov;i++){
	if(i_inLga->slp_cov->random[i]==0){
	  nfix_slp+=1;       
	  n += i_inLga->slp_cov->n_lev[i];
	}
      }

      fprintf(fp,"lga_eff_mean");
      for(i=0;i<n;i++){
	fprintf(fp," %f",i_inPrior->lga_eff_mean[i]);
      }
      fprintf(fp,"\n");
      
      fprintf(fp,"lga_eff_prec");
      for(i=0;i<(nfix_int+nfix_slp);i++){
	fprintf(fp," %f",i_inPrior->lga_eff_prec[i]);
      }
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inLga->int_cov->n_cov;i++)
	if(i_inLga->int_cov->random[i]==1)
	  n += 2;
      for(i=0;i<i_inLga->slp_cov->n_cov;i++)
	if(i_inLga->slp_cov->random[i]==1)
	  n += 2;
      n += 2; // parameters for tau_obs
      fprintf(fp,"lga_prec_par");
      for(i=0;i<n;i++){
	fprintf(fp," %f",i_inPrior->lga_prec_par[i]);
      }
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inLga->int_cov->n_cov;i++)
	if(i_inLga->int_cov->spatial[i]>0)
	  n += 2;
      for(i=0;i<i_inLga->slp_cov->n_cov;i++)
	if(i_inLga->slp_cov->spatial[i]>0)
	  n += 2;
      fprintf(fp,"lga_ar ");
      for(i=0;i<n;i++)
	fprintf(fp,"%f ",i_inPrior->lga_ar[i]);
      fprintf(fp,"\n");
    }

   if(i_write_wgl)
    {
      nfix_int = 0;
      for(i=0;i<i_inWgl->int_cov->n_cov;i++)
	if(i_inWgl->int_cov->random[i]==0)
	  nfix_int+=1;
      nfix_slp = 0;
      for(i=0;i<i_inWgl->slp_cov->n_cov;i++)
	if(i_inWgl->slp_cov->random[i]==0)
	  nfix_slp+=1;       

      n = 0;
      for(i=0;i<i_inWgl->int_cov->n_cov;i++)
	n += i_inWgl->int_cov->n_lev[i]*nfix_int;
      for(i=0;i<i_inWgl->slp_cov->n_cov;i++)
	n += i_inWgl->slp_cov->n_lev[i]*nfix_slp;
      fprintf(fp,"wgl_eff_mean");
      for(i=0;i<n;i++)
	fprintf(fp," %f",i_inPrior->wgl_eff_mean[i]);
      fprintf(fp,"\n");
      
      fprintf(fp,"wgl_eff_prec");
      for(i=0;i<(nfix_int+nfix_slp);i++)
	fprintf(fp," %f",i_inPrior->wgl_eff_prec[i]);
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inWgl->int_cov->n_cov;i++)
	if(i_inWgl->int_cov->random[i]==1)
	  n += 2;
      for(i=0;i<i_inWgl->slp_cov->n_cov;i++)
	if(i_inWgl->slp_cov->random[i]==1)
	  n += 2;
      fprintf(fp,"wgl_prec_par");
      for(i=0;i<n;i++)
	fprintf(fp," %f",i_inPrior->wgl_prec_par[i]);
      fprintf(fp,"\n");
      
      n = 0;
      for(i=0;i<i_inWgl->int_cov->n_cov;i++)
	if(i_inWgl->int_cov->spatial[i]>0)
	  n += 2;
      for(i=0;i<i_inWgl->slp_cov->n_cov;i++)
	if(i_inWgl->slp_cov->spatial[i]>0)
	  n += 2;
      fprintf(fp,"wgl_ar ");
      for(i=0;i<n;i++)
	fprintf(fp,"%f ",i_inPrior->wgl_ar[i]);
      fprintf(fp,"\n");
    }

  return(0);
}               /* end of write_inPrior */


/*!
  \brief 
  \author Hanne Rognebakke
*/
int write_inCC(FILE *fp, Data_CC *i_D_CC, int i_ncat)
{
  int i;

  fprintf(fp,"class_error %d\n",i_D_CC->class_error);
  fprintf(fp,"ptype1_CC1");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype1_CC1[i]);
  fprintf(fp,"\nptype5_CC1");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype5_CC1[i]);
  fprintf(fp,"\nptype2_CC2");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype2_CC2[i]);
  fprintf(fp,"\nptype4_CC2");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype4_CC2[i]);
  fprintf(fp,"\nptype2_S4");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype2_S4[i]);
  fprintf(fp,"\nptype4_S4");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype4_S4[i]);
  fprintf(fp,"\nptype1_S5");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype1_S5[i]);
  fprintf(fp,"\nptype5_S5");
  for(i=0;i<i_ncat;i++)
    fprintf(fp," %f",i_D_CC->ptype5_S5[i]);
  fprintf(fp,"\n");

  return(0);
}               /* end of write_inCC */



/*!
  \brief Write wgl parameters to file
  \author Hanne Rognebakke
*/
int write_inWgl(FILE *fp, Input_wgl *i_inWgl, int i_nHaul, int i_pred, int i_it)
{
  int    i,err=0;

  fprintf(fp,"fixed_model %d\n",i_inWgl->fixed_model);
  if(i_inWgl->fixed_model)
    {
      fprintf(fp,"fixed_int ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inWgl->fixed_int[i]);
      fprintf(fp,"\nfixed_slp ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inWgl->fixed_slp[i]);
      fprintf(fp,"\nfixed_tau ");
      for(i=0;i<i_it;i++)
	fprintf(fp,"%f ",i_inWgl->fixed_tau[i]);
      fprintf(fp,"\n");
    }

  /* Write covariate information */
  fprintf(fp,"covatiates intercept:\n");
  err = write_input_cov(fp, i_inWgl->int_cov, 1, i_nHaul, i_pred);
  if(err)
    {
      printWarning("write_inWgl:Error calling write_input_cov\n");
      return(err);
    }
  fprintf(fp,"covatiates slope:\n");
  err = write_input_cov(fp, i_inWgl->slp_cov, 1, i_nHaul, i_pred);
  if(err)
    {
      printWarning("write_inWgl:Error calling write_input_cov\n");
      return(err);
    }
  return(0);
}               /* end of write_inWgl */


/*!
  \brief Write predict parameters to file
  \author Hanne Rognebakke
*/
int write_inPredict(FILE *fp, Input_predict *i_inPredict, int i_nCell)
{
  int i,err=0;

  fprintf(fp,"nMCMC %d\n",i_inPredict->nMCMC);
  fprintf(fp,"burnin %d\n",i_inPredict->burnin);
  fprintf(fp,"nHaul %d\n",i_inPredict->nHaul);
  fprintf(fp,"num_par1");
  for(i=0;i<5;i++)
    fprintf(fp," %d",i_inPredict->num_par1[i]);
  fprintf(fp,"\nnum_par2");
  for(i=0;i<2;i++)
    fprintf(fp," %d",i_inPredict->num_par2[i]);
  fprintf(fp,"\nn_MC");
  for(i=0;i<i_nCell;i++)
    fprintf(fp," %d",i_inPredict->n_MC[i]);
  fprintf(fp,"\nN_l_int %d\n",i_inPredict->N_l_int);
  fprintf(fp,"l_int");
  for(i=0;i<i_inPredict->N_l_int;i++)
    fprintf(fp," %f",i_inPredict->l_int[i]);
  fprintf(fp,"\n");
  fprintf(fp,"inc_hsz %d\n",i_inPredict->inc_hsz);
  fprintf(fp,"coastal_cod %d\n",i_inPredict->coastal_cod);
  fprintf(fp,"filename_mcmc1 %s\n",i_inPredict->filename_mcmc1);
  fprintf(fp,"filename_mcmc2 %s\n",i_inPredict->filename_mcmc2);
  fprintf(fp,"filename_hsz_mcmc2 %s\n",i_inPredict->filename_hsz_mcmc2);
  fprintf(fp,"filename_predict %s\n",i_inPredict->filename_predict);

  return(err);
}               /* end of write_inPredict */


/*!
  \brief Write predict cell parameters to file
  \author Hanne Rognebakke
*/
int write_inCell(FILE *fp, Input_cell *i_inCell, int i_cc)
{
  int    i,n;

  fprintf(fp,"num_cell_o");
  n = 5;
  if(i_cc)
    n = 9;
  for(i=0;i<n;i++)
    fprintf(fp," %d",i_inCell->num_cell_o[i]);
  fprintf(fp,"\nnum_cell_u");
  for(i=0;i<n;i++)
    fprintf(fp," %d",i_inCell->num_cell_u[i]);
  fprintf(fp,"\nage_int_nC %d\n",i_inCell->age_int_nC);

  return(0);
}               /* end of write_inCell */


/*!
  \brief Write total catch parameters to file
  \author Hanne Rognebakke
*/
int write_inCatch(FILE *fp, Input_totcatch *i_inCatch, 
		  int *i_age_n_cov, int *i_lga_n_cov, int *i_wgl_n_cov)
{
  int    i,c;

  fprintf(fp,"nCell %d\n",i_inCatch->nCell);
  fprintf(fp,"nFactors %d\n",i_inCatch->nFactors);

  fprintf(fp,"tot_fac_age_int");
  for(i=0;i<i_age_n_cov[0];i++)
    fprintf(fp," %d",i_inCatch->fac_age_int[i]);
  fprintf(fp,"\n");

  fprintf(fp,"tot_fac_lga_int ");
  for(i=0;i<i_lga_n_cov[0];i++)
    fprintf(fp," %d",i_inCatch->fac_lga_int[i]);
  fprintf(fp,"\n");
  fprintf(fp,"tot_fac_lga_slp ");
  for(i=0;i<i_lga_n_cov[1];i++)
    fprintf(fp," %d",i_inCatch->fac_lga_slp[i]);
  fprintf(fp,"\n");

  fprintf(fp,"tot_fac_wgl_int ");
  for(i=0;i<i_wgl_n_cov[0];i++)
    fprintf(fp," %d",i_inCatch->fac_wgl_int[i]);
  fprintf(fp,"\n");
  fprintf(fp,"tot_fac_wgl_slp ");
  for(i=0;i<i_wgl_n_cov[1];i++)
    fprintf(fp," %d",i_inCatch->fac_wgl_slp[i]);
  fprintf(fp,"\n");

  fprintf(fp,"cell ");
  for(i=0;i<i_inCatch->nFactors;i++)
    fprintf(fp,"factor%d ",i);
  fprintf(fp,"catch\n");
  for(c=0;c<i_inCatch->nCell;c++)
    {
      fprintf(fp,"%d ",c);
      for(i=0;i<i_inCatch->nFactors;i++)
	{
	  fprintf(fp,"%d ",i_inCatch->factors[i][c]);
	}
      fprintf(fp,"%lf\n",i_inCatch->catch[c]);
    }

  return(0);
}               /* end of write_inCatch */

/*!
  \brief Write parameters in struct Input_Cov to file

  \author Hanne Rognebakke
*/
int write_input_cov(FILE *fp, Input_cov *i_inCov, int i_nAges, int i_nHaul, int i_pred)
{
  int i,j;

  fprintf(fp,"n_cov");
  fprintf(fp," %d",i_inCov->n_cov);

  fprintf(fp,"\nn_lev");
  for(i=0;i<i_inCov->n_cov;i++)
    fprintf(fp," %d",i_inCov->n_lev[i]);
  fprintf(fp,"\nrandom");
  for(i=0;i<i_inCov->n_cov;i++)
    fprintf(fp," %d",i_inCov->random[i]);
  fprintf(fp,"\nspatial");
  for(i=0;i<i_inCov->n_cov;i++)
    fprintf(fp," %d",i_inCov->spatial[i]);
  fprintf(fp,"\ncontinuous");
  for(i=0;i<i_inCov->n_cov;i++)
    fprintf(fp," %d",i_inCov->continuous[i]);
  fprintf(fp,"\ninteraction");
  for(i=0;i<i_inCov->n_cov;i++)
    fprintf(fp," %d",i_inCov->interaction[i]);
  fprintf(fp,"\n");
  if(i_pred==0)
    {
      fprintf(fp,"c_cov");
      for(i=0;i<i_inCov->n_cov_i;i++)
	{
	  fprintf(fp,"cov %d: ",i);
	  for(j=0;j<i_nHaul;j++)
	    {
	      fprintf(fp," %d",i_inCov->c_cov_i[i][j]);
	    }
	  fprintf(fp,"\n");
	}
      for(i=0;i<i_inCov->n_cov_d;i++)
	{
	  fprintf(fp,"cov %d: ",i);
	  for(j=0;j<i_nHaul;j++)
	    {
	      fprintf(fp," %f",i_inCov->c_cov_d[i][j]);
	    }
	  fprintf(fp,"\n");
	}
    }
  if(i_pred==0)
    fprintf(fp,"write_input_cov: Write spatial information\n");
  if(i_pred==0)
    fprintf(fp,"write_input_cov: Write cell/interaction information\n");

  return(0);
}               /* end of write_input_cov */

