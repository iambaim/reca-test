#include "caa.h"
#include "caa_routines.h"






/*L:calc_summaries*
________________________________________________________________

		calc_summaries
________________________________________________________________

Name:		calc_summaries
Syntax:		
Description:    samples ages for hauls with missing ages.
                Assumes lengths are given by categories. For each
                length category posterior probabilites are calculated
                using prior and length, and ages are sampled from multinomial 
                distribution.
                Stored in ages per haul.
Side effects:   
Return value:   zero
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
int calc_summaries(int i_start_h,Age_struct *i_age,Data_age *i_D_age,Data_g_a *i_D_g_a,
		   LW_struct *i_length,Data_lin *i_D_lga,Data_l *i_D_l)
{
  int            a,h,i;
  double         mean,sum,mean_h,mean_l,mean_lga,N,N_h,mu[2];
  Data_cov      *xcov;
  double        *freq, *p;


  freq = CALLOC(i_D_age->glm->ncat,double);  // Free ok
  p= CALLOC(i_D_age->glm->ncat,double);      // Free ok

  N = 0;
  mean_l = 0.0;
  for(a=0;a<i_D_age->glm->ncat;a++)
    freq[a] = G_ZERO;
  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
    {
      N_h = i_D_lga->glm->suff[h][0][0];
      /* Calculate frequencies */
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
          mean = i_age->alpha[h][a];
          p[a] = exp(mean);
          sum += p[a];
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] /= sum;
	  freq[a] += p[a] * N_h;
	}
      /* Mean of length given age */
      for(i=0;i<i_D_lga->glm->nxcov;i++)
	{
	  xcov = i_D_lga->glm->xcov[i];
	  mu[i] = calc_eff(xcov,i_length->par->eff[0][i],h);
	}
      mean_h = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
          mean_lga = exp(mu[0]+mu[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]]+G_HALF/i_length->par->tau_obs);
          mean_h += p[a] * mean_lga;
	}
      mean_l += mean_h * N_h;
      N += N_h;
    }
  sum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    sum += freq[a];
  printf("Age frequencies\n");
  for(a=0;a<i_D_age->glm->ncat;a++)
    printf("Age[%d]=%8.2lf (%5.3lf)\n",a+2,freq[a],freq[a]/sum);

  mean_l /= N;

  printf("Mean_length = %lf\n",mean_l);

  FREE(freq);
  FREE(p);

  return(0);
}		/* end of calc_summaries */

