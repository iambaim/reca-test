#include "caa.h"



/*L:lm_fit*
________________________________________________________________

		lm_fit
________________________________________________________________

Name:		lm_fit
Syntax:		
Description:    Linear regression on y given x. Output is beta0,beta1,ssq
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double lm_fit(int i_N,double *i_x,double *i_y,
	      double *o_beta0,double *o_beta1,double *o_ssq)
{
  int    i;
  double n,n_1,res;
  double s_x,s_y,s_xx,s_yy,s_xy;

  s_x = s_y = s_xx = s_yy = s_xy = G_ZERO;

  for(i=0;i<i_N;i++)
    {
      n_1 = (double) i;
      n = (double) (i+1);
      s_x  = (s_x*n_1+i_x[i])/n;
      s_y  = (s_y*n_1+i_y[i])/n;
      s_xx = (s_xx*n_1+i_x[i]*i_x[i])/n;
      s_yy = (s_yy*n_1+i_y[i]*i_y[i])/n;
      s_xy = (s_xy*n_1+i_x[i]*i_y[i])/n;
    }
  if(fabs(s_xx-s_x*s_x)>0.0000000001)
      *o_beta1 = (s_xy-s_x*s_y)/(s_xx-s_x*s_x);
  else
      *o_beta1 = G_ZERO;
  *o_beta0 = s_y - *o_beta1*s_x;
  *o_ssq = G_ZERO;
  for(i=0;i<i_N;i++)
    {
      res = i_y[i] - *o_beta0 - *o_beta1*i_x[i];
      *o_ssq += res*res;
    }

  return(0);
}		/* end of lm_fit */



/*L:lm_fit_suff*
________________________________________________________________

		lm_fit_suff
________________________________________________________________

Name:		lm_fit_suff
Syntax:		
Description:    Linear regression on y given x. Use as input sufficient
                statistics. Output is beta0,beta1,ssq
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double lm_fit_suff(int i_ncat,double *i_g_a,double *i_sum_by_cat,double *i_sqsum_by_cat,double **i_suff,
		   double *o_beta0,double *o_beta1,double *o_ssq)
{
  int    a;
  double beta0,beta1,ssq;
  double s_x,s_y,s_xx,s_yy,s_xy;

  s_x = s_y = s_xx = s_yy = s_xy = G_ZERO;

  s_x = i_suff[0][1]/i_suff[0][0];
  s_xx = i_suff[1][1]/i_suff[0][0];
  for(a=0;a<i_ncat;a++)
    {
      s_y += i_sum_by_cat[a];
      s_xy += i_g_a[a]*i_sum_by_cat[a];
      s_yy += i_sqsum_by_cat[a];
    }
  s_y /= i_suff[0][0];
  s_xy /= i_suff[0][0];
  s_yy /= i_suff[0][0];

  if(fabs(s_xx-s_x*s_x)>0.0000000001)
      beta1 = (s_xy-s_x*s_y)/(s_xx-s_x*s_x);
  else
      beta1 = G_ZERO;
  beta0 = s_y - beta1*s_x;
  ssq = i_suff[0][0]*
        (s_yy+beta0*beta0+beta1*beta1*s_xx-
         G_TWO*beta0*s_y-G_TWO*beta1*s_xy+G_TWO*beta0*beta1*s_x);

  *o_beta0 = beta0;
  *o_beta1 = beta1;
  *o_ssq = ssq;

  return(0);
}		/* end of lm_fit_suff */



/*L:dnorm*
________________________________________________________________

		dnorm
________________________________________________________________

Name:		dnorm
Syntax:		
Description:    Normal density
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double dnorm(double x,double mu,double sigma)
{
  double res,d;

  res = x-mu;

  d = exp(-G_HALF*res*res/(sigma*sigma))/(sqrt(G_TWO*G_PI)*sigma);

  return(d);
}		/* end of dnorm */



/*L:pnorm*
________________________________________________________________

		pnorm
________________________________________________________________

Name:		pnorm
Syntax:		
Description:    Normal density
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double pnorm(double x,double mu,double sigma)
{
  double res,p;

  if(x>mu)
    {
      res = (x-mu)/(sqrt(G_TWO)*sigma);
      p = G_HALF+G_HALF*d_erff(res);
    }
  else
    {
      res = (mu-x)/(sqrt(G_TWO)*sigma);
      p = G_HALF-G_HALF*d_erff(res);
    }


  return(p);
}		/* end of pnorm */



/*L:pnorm0*
________________________________________________________________

		pnorm0
________________________________________________________________

Name:		pnorm0
Syntax:		
Description:    Cumulative distribution for standard Gaussian distribution
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double pnorm0(double x)
{
  double p;

  if(x>G_ZERO)
      p = G_HALF+G_HALF*d_erff(x/sqrt(G_TWO));
  else
      p = G_HALF-G_HALF*d_erff(-x/sqrt(G_TWO));


  return(p);
}		/* end of pnorm0 */



/*L:dt*
________________________________________________________________

		dt
________________________________________________________________

Name:		dt
Syntax:		
Description:    T density with df degrees of freedom
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double dt(double x,double mu,double sigma,double df)
{
  double res,d;

  res = (x-mu)/sigma;

  d = exp(-G_HALF*(df+G_ONE)*log(G_ONE+res*res/df))/sigma;

  return(d);
}		/* end of dt */


