/*!
  \file caa_evaluate.c
  \brief Routines for model evaluations
  \author Geir Storvik

*/
#include "numrec.h"
#include "caa.h"
#include "caa_read_write.h"
#include "caa_evaluate.h"
#include "caa_routines.h"
#include "caa_init.h"
#include "caa_util.h"
#include "caa_sample_g_a.h"

static int                 s_N2=10;       
static double             *s_eps;
static double             *s_d_rep;
static double             *s_g;
static double             *s_entr;


/*!
  \author Geir Storvik
  \brief Initialization for evaluation routines

  At present only used for opening file for debugging.
*/
int init_evaluate(int i_nHaul,int i_ncat)
{
  int    b,i,j,k1,k2,n;
  double d_i,d_N;
  static int (*cmp)();

  #ifdef DEBUG_EVALUATE_FILE
  s_unit = fopen("caa_evaluate.txt","w");
  #endif

  #ifdef DEBUG_IND_KS
  s_ind_ks = fopen("caa_ind_ks.txt","w");
  #endif


  s_eps = CALLOC(i_nHaul*i_ncat+1,double);

  cmp = compd;
  n = i_nHaul*i_ncat;
  k1 = (int) ((double) n * 0.1);
  k2 = n-k1+1;
  s_d_rep = CALLOC(s_N2,double);
  s_g = CALLOC(s_N2,double);

  for(i=0;i<s_N2;i++)
    {
      for(j=0;j<n;j++)
	s_eps[j] = gennor(G_ZERO,G_ONE);
      qsort(s_eps,n,sizeof(double),cmp);
      s_d_rep[i] = fabs(fabs(s_eps[k2]) - fabs(s_eps[k1]));
    }
  qsort(s_d_rep,s_N2,sizeof(double),cmp);
  d_i = G_ONE;
  d_N = (double) s_N2;
  b = -1;
  while(s_d_rep[b+1]<(d_i*s_d_rep[s_N2-1]/d_N))
    b++;
  s_g[0] = b;
  for(i=0;i<s_N2;i++)
    {
      d_i = (double) (i+1);
      while(s_d_rep[b+1]<(d_i*s_d_rep[s_N2-1]/d_N)&&b<s_N2)
	b++;
      s_g[i] = b;
    }

  s_entr = CALLOC(i_nHaul,double);

  return(0);
}		/* end of init_evaluate */



/*!
  \author Geir Storvik
  \brief Re-initialize thing initialized in initialize_evaluate
*/
int re_init_evaluate(int i_ncat)
{
  #ifdef DEBUG_EVALUATE_FILE
  fclose(s_unit);
  #endif

  #ifdef DEBUG_IND_KS
  fclose(s_ind_ks);
  #endif

  FREE(s_eps);
  FREE(s_g);
  FREE(s_d_rep);
  FREE(s_entr);

   return(0);
}		/* end of re_init_evaluate */






/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates log-likelihood for age-model
*/
int calc_lik_age(Age_struct *i_age,Data_age *i_D_age,double *o_loglik)
{
  int  h;
  double sum;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  sum = G_ZERO;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    sum += calc_lik_age_h(h,i_age,i_D_age);
  
  (*o_loglik) = sum;

  return(0);
}		/* end of calc_lik_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the age and
     lga models

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_model1(int i_it,Age_struct *i_age,Data_age *i_D_age,
		    LW_struct *i_lga,Data_lin *i_D_lga,
		    double *o_mean_inv_lik_mod1)
{
  int  h,err=0;
  double loglik_age,loglik_lga,lik;
  double *res;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  res = CALLOC(i_D_lga->glm->nxcov,double);     // Free ok

  for(h=0;h<i_D_age->glm->nHaul;h++){
    loglik_age = calc_lik_age_h(h,i_age,i_D_age);
    loglik_lga = calc_lik_lin_h(h,i_lga,i_D_lga,res);
    lik = exp(loglik_age+loglik_lga);
    err = update_mean(&o_mean_inv_lik_mod1[h],G_ONE/lik,i_it);
    if(err){
      write_warning("Bayes_CV_model1:Error calling update_mean\n");
      return(err);
    }
    #ifdef DEBUG_EVALUATE_FILE
    fprintf(s_unit,"%d %d %lf\n",i_it,h,lik);
    #endif
  }
  
  FREE(res);

  return(0);
}		/* end of Bayes_CV_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the age model

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_age(int i_it,Age_struct *i_age,Data_age *i_D_age,
                 double *o_mean_inv_lik_age)
{
  int  err=0,h;
  double lik;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  for(h=0;h<i_D_age->glm->nHaul;h++){
    lik = exp(calc_lik_age_h(h,i_age,i_D_age));
    err = update_mean(&o_mean_inv_lik_age[h],G_ONE/lik,i_it);
    if(err){
      write_warning("Bayes_CV_age:Error calling update_mean\n");
      return(err);
    }
    #ifdef DEBUG_EVALUATE_FILE
    fprintf(s_unit,"%d %d %lf\n",i_it,h,lik);
    #endif
  }
  

  return(0);
}		/* end of Bayes_CV_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates the log-likelihood for the age model in a specific haul
*/
double calc_lik_age_h(int i_h,Age_struct *i_age,Data_age *i_D_age)
{
  int  a;
  double term1,term2,loglik, n_h;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  n_h = G_ZERO;
  term1 = G_ZERO;
  term2 = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      n_h += (double) i_D_age->Ages_fix[i_h][a];
      term1 += (double) i_D_age->Ages_fix[i_h][a] *i_age->alpha[i_h][a];
      term2 += exp(i_age->alpha[i_h][a]);
    }
  term2 = n_h * log(term2);
  loglik = term1-term2;
  
  return(loglik);
}		/* end of calc_lik_age_h */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates log-likelihood for the lga or wgl model
*/
int calc_lik_lin(LW_struct *i_lin,Data_lin *i_D_lin,double *o_loglik)
{
  int  h;
  double loglik, *res;

  res = CALLOC(i_D_lin->glm->nxcov,double);     // Free ok

  loglik = G_ZERO;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    loglik += calc_lik_lin_h(h,i_lin,i_D_lin,res);
  
  (*o_loglik) = loglik;

  FREE(res);

  return(0);
}		/* end of calc_lik_lin */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the lga or wgl model

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_lin(int i_it,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *o_mean_inv_lik_lin)
{
  int  err=0,h;
  double lik, *res;

  res = CALLOC(i_D_lin->glm->nxcov,double);     // Free ok

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_lin er en struct som har data estimater paa forskjellige 
     parametere */

  for(h=0;h<i_D_lin->glm->nHaul;h++){
    lik = exp(calc_lik_lin_h(h,i_lin,i_D_lin,res));
    err = update_mean(&o_mean_inv_lik_lin[h],G_ONE/lik,i_it);
    if(err){
      write_warning("Bayes_CV_lin:Error calling update_mean\n");
      return(err);
    }
  }
  

  FREE(res);

  return(0);
}		/* end of Bayes_CV_lin */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates the log-likelihood for the lga or wgl model in a specific haul
*/
double calc_lik_lin_h(int i_h,LW_struct *i_lin,Data_lin *i_D_lin,double *w_res)
{
  int  i,j;
  double term12,sum, N,l_h,loglik;
  Data_cov *xcov;

  /* LW_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_lin er en struct som har data estimater paa forskjellige 
     parametere */
  /* Se write_it_lin for forklaring på i_lin
     Se fish_sim.h for forklaring på Data_lin */
  N = (double) (i_D_lin->glm->suff[i_h][0][0]);

  /* Number of fish in haul times Number of nHaul */
  /* tau_cell,tau_area,tau_haul,tau_fish */
  term12 = -G_HALF * N*(log(G_TWO*G_PI)-log(i_lin->par->tau_obs));
  sum = G_ZERO;
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      xcov = i_D_lin->glm->xcov[i];
      w_res[i] = calc_eff(xcov,i_lin->par->eff[0][i],i_h) - 
	i_D_lin->glm->beta_hat[i_h][0][i];
    }

  l_h = i_D_lin->glm->ssq[i_h];
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    for(j=0;j<i_D_lin->glm->nxcov;j++)
    {
      l_h += i_D_lin->glm->suff[i_h][i][j] * w_res[i] * w_res[j];
    }
  sum += l_h;
  
  loglik = term12-G_HALF*sum*i_lin->par->tau_obs;

  return(loglik);
}		/* end of calc_lik_lin_h */


/*!
  \author Geir Storvik
  \brief Calculates residuals in lga model based on current simulated parameters
*/
int calc_resid_lga(Data_orig *i_D_orig,
                   Age_struct *i_age, Data_age *i_D_age,LW_struct *i_lin,Data_lin *i_D_lin,
		   Data_g_a *i_D_g_a,double *o_resid)
{
  int  a,f,h,i,ind;
  double *mu, *beta;
  Data_cov *xcov;

  beta = CALLOC(i_D_lin->glm->nxcov,double);
  mu = CALLOC(i_D_age->glm->ncat,double);

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      for(i=0;i<i_D_lin->glm->nxcov;i++)
	{
	  xcov = i_D_lin->glm->xcov[i];
	  beta[i] = calc_eff(xcov,i_lin->par->eff[0][i],h);
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]];
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  if(i_D_orig->totage[ind] > -1000 && i_D_orig->totlength[ind]> -1000.00)
	    o_resid[ind] = i_D_orig->totlength[ind]-mu[i_D_orig->totage[ind]-i_D_age->a_vec[0]];
	  else
	    o_resid[ind] = -99999.99;
	  ind++;
	}
    }
  #ifdef DEBUG_EVALUATE
  for(i=0;i<min(ind,100);i++)
    printf("i=%d,a=%d,res_lga=%lf\n",i,(int) i_D_orig->totage[i],o_resid[i]);
  #endif
  FREE(beta);
  FREE(mu);

  return(0);
}		/* end of calc_resid_lga */



/*!
  \author Geir Storvik
  \brief Calculates residuals in wgl model based on current simulated parameters
*/
int calc_resid_wgl(double *i_totlength, double *i_totweight,int *i_nFishBoat,
                   LW_struct *i_weight,Data_lin *i_D_wgl,
		   double *o_resid)
{
  int  f,h,i,ind;
  double mu, *beta;
  Data_cov *xcov;

  beta = CALLOC(i_D_wgl->glm->nxcov,double);

  ind = 0;
  for(h=0;h<i_D_wgl->glm->nHaul;h++)
    {
      for(i=0;i<i_D_wgl->glm->nxcov;i++)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  beta[i] = calc_eff(xcov,i_weight->par->eff[0][i],h);
	}
      for(f=0;f<i_nFishBoat[h];f++)
	{
	  if(i_totlength[ind]> -1000.00 && i_totweight[ind]> -1000.00)
	    {
	      mu = beta[0]+beta[1]*i_totlength[ind];
	      o_resid[ind] = i_totweight[ind]-mu;
	    }
	  else
	    o_resid[ind] = min(i_totlength[ind],i_totweight[ind]);
	  ind++;
	}
    }
  FREE(beta);

  return(0);
}		/* end of calc_resid_wgl */



