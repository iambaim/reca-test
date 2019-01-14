#include "caa.h"
#include "caa_read_write.h"

#define  MAX(A,B)  ((A) > (B) ? (A) : (B))
#define  MIN(A,B)  ((A) < (B) ? (A) : (B))
static int lqp_optimize_f(double *x,double (*log_f)(double,double *),
                                    double (*log_f1)(double,double *),
			  double (*log_f2)(double,double *),double *i_par);
static int lqp_find_knots(int i_n,double i_x_opt,double (*log_f)(double,double *),
			  double (*log_f1)(double,double *),
			  double (*log_f2)(double,double *),double *i_par,
			  double *o_knots);
static int lqp_find_approx(int i_n,double *i_knots,double i_fopt,double *i_par,
			   double (*i_log_f)(double,double *),
			   double (*i_log_f1)(double,double *),
                           double **o_coef,double *o_prob,double *o_max_r);
static int lqp_sample_val(int i_m,double *i_knots,double i_fopt,double *i_par,
			  double (*i_log_f)(double,double *),
			  double **i_coef,double *i_prob,double i_max_r,
			  double *o_x,double *o_lacc);
static int lqp_lacc_val(int i_m,double *i_knots,double i_fopt,double *i_par,
			double (*i_log_f)(double,double *),
			double **i_coef,double *i_prob,double i_max_r,
			double i_x,double *o_lacc);
static double log_f_ratio(double i_x,double i_fopt,
                          double (*i_log_f)(double,double *),
                          double *i_par,double i_a,double i_b,double i_c);
static double log_f_ratio1(double x,double (*log_f1)(double,double *),
                           double *i_par,double i_a,double i_b);

/*C*

________________________________________________________________

        log_quadratic_proposal
        $Id: caa_lqp.c 1 2013-03-28 13:54:24Z hanne $
	Copyright 2003, 
________________________________________________________________
*/
int lqp_sample(int i_m,
	       double (*i_log_f)(double,double *),
	       double (*i_log_f1)(double,double *),
	       double (*i_log_f2)(double,double *),
               double *i_par_old,double i_x_old,double *i_par_new,
               double *w_knots,double **w_coef,double *w_prob,
               double *o_x_new,double *o_lacc)
{
  int  err;
  double mean,x;
  double max_r,lacc_new,lacc_old,fopt_new,fopt_old;

  mean = i_par_new[2];
  x = mean;
  err = lqp_optimize_f(&x,i_log_f,i_log_f1,i_log_f2,i_par_new);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_optimize\n");
      return(err);
    }
  err = lqp_find_knots(i_m,x,i_log_f,i_log_f1,i_log_f2,i_par_new,w_knots);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_find_knots\n");
      return(err);
    }
  /* Find spline approximation */      
  fopt_new = i_log_f(x,i_par_new);
  err = lqp_find_approx(i_m,w_knots,fopt_new,i_par_new,i_log_f,i_log_f1,
                        w_coef,w_prob,&max_r);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_find_approx\n");
      return(err);
    }
  /* Sampling */
  err = lqp_sample_val(i_m,w_knots,fopt_new,i_par_new,i_log_f,
		       w_coef,w_prob,max_r,&x,&lacc_new);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_sample_val\n");
      return(err);
    }
  *o_x_new = x;
  
  /* Finding reverse proposal distribution */
  mean = i_par_old[2];
  x = mean;
  err = lqp_optimize_f(&x,i_log_f,i_log_f1,i_log_f2,i_par_old);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_optimize\n");
      return(err);
    }
  err = lqp_find_knots(i_m,x,i_log_f,i_log_f1,i_log_f2,i_par_old,w_knots);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_find_knots\n");
      return(err);
    }
  /* Find spline approximation */      
  fopt_old = i_log_f(x,i_par_old);
  err = lqp_find_approx(i_m,w_knots,fopt_old,i_par_old,i_log_f,i_log_f1,
			w_coef,w_prob,&max_r);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_find_approx\n");
      return(err);
    }
  /* Finding lacc */
  err = lqp_lacc_val(i_m,w_knots,fopt_old,i_par_old,i_log_f,
		     w_coef,w_prob,max_r,i_x_old,&lacc_old);
  if(err)
    {
      write_warning("lqp_sample:Error calling lpq_lacc_val\n");
      return(err);
    }

  *o_lacc = lacc_new+fopt_new-lacc_old-fopt_old;

  return(0);
}



/*F:lqp_optimize_f*

________________________________________________________________

		lqp_optimize_f
________________________________________________________________

Name:		lqp_optimize_f
Syntax:		
Description: Optimize a function usin Newton's method
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_lqp.c 1 2013-03-28 13:54:24Z hanne $
________________________________________________________________
*/
static int lqp_optimize_f(double *x,double (*i_log_f)(double,double *),
			  double (*i_log_f1)(double,double *),
			  double (*i_log_f2)(double,double *),double *i_par)
{
  double  max_f,f,f1,f2,delta,x_new,r;

  max_f = i_log_f(*x,i_par);
  f1 = i_log_f1(*x,i_par);
  f2 = i_log_f2(*x,i_par);
  while(fabs(f1)>0.00000001)
    {
      delta = G_ONE;
      r = f1/f2;
      x_new = *x - delta*r;
      f = i_log_f(x_new,i_par);
      while(f<(max_f-0.0000001))
	{
          delta /= G_TWO;
	  x_new = *x - delta*r;
	  f = i_log_f(x_new,i_par);
	}
      *x = x_new;
      max_f = f;
      f1 = i_log_f1(*x,i_par);
      f2 = i_log_f2(*x,i_par);
    }
  
  return(0);
}		/* end of lqp_optimize_f */



/*F:lqp_find_knots*

________________________________________________________________

		lqp_find_knots
________________________________________________________________

Name:		lqp_find_knots
Syntax:		
Description: Find knots by first finding end points such that
             f(endpoint)=eps*f(optpoint)
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_lqp.c 1 2013-03-28 13:54:24Z hanne $
________________________________________________________________
*/
static int lqp_find_knots(int i_n,double i_x_opt,double (*i_log_f)(double,double *),
			  double (*i_log_f1)(double,double *),
			  double (*i_log_f2)(double,double *),double *i_par,
			  double *o_knots)
{
  int     i;
  double  eps,delta,f_opt,f_l,f_m,f_u,x_l,x_m=0,x_u;

  /* Finding lower end point */
  eps = log(0.0001);
  delta = 3.7/sqrt(-i_log_f2(i_x_opt,i_par));
  f_opt = i_log_f(i_x_opt,i_par);
  x_l = i_x_opt - delta;
  f_l = i_log_f(x_l,i_par);
  x_u = i_x_opt;
  while(f_l > (f_opt+eps))
    {
      x_u = x_l;
      delta *= G_TWO;
      x_l = i_x_opt - delta;
      f_l = i_log_f(x_l,i_par);
    }
  while(fabs(x_u-x_l)>0.0001)
    {
      x_m = G_HALF*(x_l+x_u);
      f_m = i_log_f(x_m,i_par);
      if(f_m > (f_opt+eps))
	x_u = x_m;
      else
        x_l = x_m;
    }
  o_knots[0] = x_m;

  /* Finding upper end point */
  delta = 3.7/sqrt(-i_log_f2(i_x_opt,i_par));
  x_l = i_x_opt;
  x_u = i_x_opt + delta;
  f_u = i_log_f(x_u,i_par);
  while(f_u > (f_opt+eps))
    {
      x_l = x_u;
      delta *= G_TWO;
      x_u = i_x_opt + delta;
      f_u = i_log_f(x_u,i_par);
    }
  while(fabs(x_u-x_l)>0.0001)
    {
      x_m = G_HALF*(x_l+x_u);
      f_m = i_log_f(x_m,i_par);
      if(f_m > (f_opt+eps))
	x_l = x_m;
      else
        x_u = x_m;
    }
  o_knots[i_n-1] = x_m;
  delta = (o_knots[i_n-1]-o_knots[0])/(double) (i_n-1);
  for(i=1;i<(i_n-1);i++)
    o_knots[i] = o_knots[i-1]+delta;

  return(0);
}		/* end of lqp_find_knots */



/*F:lqp_find_approx*

________________________________________________________________

		lqp_find_approx
________________________________________________________________

Name:		lqp_find_approx
Syntax:		
Description: For given knots, finds approximate function
    f_tilde(x) = \sum_k p_k s_k(x)I(x_{k} x<=x_{k+1})
    where s_k(x)=s_k^u(x)/\int s_k^u(x)
    while s_k^u(x)=p_{k,0}+p_{k,1}x
    and interpolates log f(x_k) and log f_{k+1}) 
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_lqp.c 1 2013-03-28 13:54:24Z hanne $
________________________________________________________________
*/
static int lqp_find_approx(int i_n,double *i_knots,double i_fopt,double *i_par,
			   double (*i_log_f)(double,double *),
			   double (*i_log_f1)(double,double *),
                           double **o_coef,double *o_prob,double *o_max_r)
{
  int     i;
  double  f_cur,f_prev,f1_cur,f1_prev,f_m;
  double  w1,w2,w3,k2_cur,k2_prev;
  double  a,a_u,a_l,b,c,max_r,x_l,x_m,x_h,f1;

  max_r = G_ZERO;
  /* First intervall from -infinity to knots[0] */
  f_cur = i_log_f(i_knots[0],i_par)-i_fopt;
  f1_cur = i_log_f1(i_knots[0],i_par);
  a = G_ZERO;
  b = f1_cur*0.9999; /* Makes f.tilde > f */
  c = f_cur - b*i_knots[0];
  o_prob[0] = exp(b*i_knots[0]+c)/b;
  o_coef[0][0] = a;
  o_coef[0][1] = b;
  o_coef[0][2] = c;
  k2_cur = i_knots[0]*i_knots[0];
  /* Intervals between knots[0] and knots[n-1] */
  for(i=1;i<i_n;i++)
    {
      f_prev = f_cur;
      f1_prev = f1_cur;
      f_cur = i_log_f(i_knots[i],i_par)-i_fopt;
      f1_cur = i_log_f1(i_knots[i],i_par);
      x_m = G_HALF*(i_knots[i-1]+i_knots[i]);
      f_m = i_log_f(x_m,i_par)-i_fopt;
      k2_prev = k2_cur;
      k2_cur = i_knots[i]*i_knots[i];
      /* Find parameters */
      a_u = (f_m-G_HALF*(f_prev+f_cur));
      a_l = x_m*x_m-G_HALF*(k2_prev+k2_cur);
      a = a_u/a_l;
      b = (f_cur-f_prev-a*(k2_cur-k2_prev))/
          (i_knots[i]-i_knots[i-1]);
      c = f_cur-(a*i_knots[i]+b)*i_knots[i];
      /* Find maximum between log f and log f.tilde */
      f1 = log_f_ratio1(i_knots[i-1],i_log_f1,i_par,a,b);
      if(f1<G_ZERO)
	{
          x_l = G_HALF*(i_knots[i-1]+i_knots[i]);
          x_h = i_knots[i];
          x_m = G_HALF*(x_l+x_h);
          while(fabs(x_h-x_l)>0.001)
	    {
              f1 = log_f_ratio1(x_m,i_log_f1,i_par,a,b);
              if(f1>G_ZERO)
		{
                  x_l = x_m;
		  x_m = G_HALF*(x_l+x_h);
		}
              else
		{
		  x_h = x_m;
		  x_m = G_HALF*(x_l+x_h);
		}
	    }
          max_r = MAX(max_r,log_f_ratio(x_m,i_fopt,i_log_f,i_par,a,b,c));
	}
      else
	{
          x_l = i_knots[i-1];
          x_h = G_HALF*(i_knots[i-1]+i_knots[i]);
          x_m = G_HALF*(x_l+x_h);
          while(fabs(x_h-x_l)>0.001)
	    {
              f1 = log_f_ratio1(x_m,i_log_f1,i_par,a,b);
              if(f1>G_ZERO)
		{
                  x_l = x_m;
		  x_m = G_HALF*(x_l+x_h);
		}
              else
		{
		  x_h = x_m;
		  x_m = G_HALF*(x_l+x_h);
		}
	    }
          max_r = MAX(max_r,log_f_ratio(x_m,i_fopt,i_log_f,i_par,a,b,c));
	}
      /* Calculate integral of s_k^u */ 
      w1 = G_HALF*sqrt(-G_PI/a)*exp(c-G_ONE_FOURTH*b*b/a);
      w2 = erfcc(sqrt(-a)*(i_knots[i-1]+b/(G_TWO*a)));
      w3 = erfcc(sqrt(-a)*(i_knots[i]+b/(G_TWO*a)));
      o_prob[i] = w1*(w2-w3);
      o_coef[i][0] = a;
      o_coef[i][1] = b;
      o_coef[i][2] = c;
   }
  /* Finally intervall from knots[n-1] to infinity */
  f_prev = f_cur;
  f1_prev = f1_cur;
  a = G_ZERO;
  b = f1_prev*0.9999;
  c = f_prev - b*i_knots[i_n-1];
  o_prob[i_n] = -exp(c)*exp(b*i_knots[i_n-1])/b;
  o_coef[i_n][0] = a;
  o_coef[i_n][1] = b;
  o_coef[i_n][2] = c;

  *o_max_r = max_r;
  return(0);
}		/* end of lqp_find_approx */



/*L:lqp_sample_val*
________________________________________________________________

		lqp_sample_val
________________________________________________________________

Name:		lqp_sample_val
Syntax:		
Description:
Side effects:   
Return value:   None
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int lqp_sample_val(int i_m,double *i_knots,double i_fopt,double *i_par,
			  double (*i_log_f)(double,double *),
			  double **i_coef,double *i_prob,double i_max_r,
			  double *o_x,double *o_lacc)
{
  int     i;
  double  u,p,a,b,c,beta,x,x_opt,max_r,e,r,f;
 
  /* Convert probabilities to cummulative ones */
  for(i=1;i<=i_m;i++)
    i_prob[i] = i_prob[i-1] + i_prob[i];
  /* Sample intervall */
  u = genunf(G_ZERO,i_prob[i_m]);
  i=0;
  while(u > i_prob[i])
    i++;
  /* Sample inside intervall */
  if(i==0)
    p = i_prob[i]/i_prob[i_m];
  else
      p = (i_prob[i]-i_prob[i-1])/i_prob[i_m];
  a = i_coef[i][0];
  b = i_coef[i][1];
  c = i_coef[i][2]-log(p);
  if(i==0)
    {
      /* Exact sampling */
      u = genunf(G_ZERO,G_ONE);
      x = i_knots[0]+log(u)/b;
    }
  else if(i<i_m)
    {
      /* Rejection sampling using exponential approximation */
      beta = a*(i_knots[i-1]+i_knots[i])+b;
      x_opt = G_HALF*(beta-b)/a;
      max_r = (a*x_opt+b-beta)*x_opt;
      e = exp(beta*(i_knots[i]-i_knots[i-1]))-G_ONE;
      u = genunf(G_ZERO,G_ONE);
      x = i_knots[i-1]+log(G_ONE+u*e)/beta;
      r = (a*x+b-beta)*x;
      while(genunf(G_ZERO,G_ONE)>exp(r-max_r))
	{
	  u = genunf(G_ZERO,G_ONE);
	  x = i_knots[i-1]+log(G_ONE+u*e)/beta;
	  r = (a*x+b-beta)*x;
	}
    }
  else
    {
      /* Exact sampling */
      u = genunf(G_ZERO,G_ONE);
      x = i_knots[i-1]+log(G_ONE+u)/b;
    }

  *o_x = x;
  c = i_coef[i][2];
  f = log_f_ratio(x,i_fopt,i_log_f,i_par,a,b,c);
  *o_lacc = f+log(i_prob[i_m]);
  return(0);
}		/* end of lqp_sample_val */



/*L:lqp_lacc_val*
________________________________________________________________

		lqp_lacc_val
________________________________________________________________

Name:		lqp_lacc_val
Syntax:		
Description:
Side effects:   
Return value:   None
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int lqp_lacc_val(int i_m,double *i_knots,double i_fopt,double *i_par,
			  double (*i_log_f)(double,double *),
			  double **i_coef,double *i_prob,double i_max_r,
			  double i_x,double *o_lacc)
{
  int     i,ind;
  double  a,b,c,f,sum_p;

  /* Convert to cummulative prob's */ 
  sum_p = i_prob[0];
  for(i=1;i<=i_m;i++)
    sum_p += i_prob[i];
  /* Find intervall for x */
  ind = 0;
  while(ind < i_m && i_x>i_knots[ind])
    ind++;
  a = i_coef[ind][0];
  b = i_coef[ind][1];
  c = i_coef[ind][2];
  f = log_f_ratio(i_x,i_fopt,i_log_f,i_par,a,b,c);
  *o_lacc = f+log(sum_p);
  return(0);
}		/* end of lqp_lacc_val */



/*L:log_f_ratio*
________________________________________________________________

		log_f_ratio
________________________________________________________________

Name:		log_f_ratio
Syntax:		
Description:
Side effects:   
Return value:   None
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double log_f_ratio(double i_x,double i_fopt,
                          double (*i_log_f)(double,double *),
                          double *i_par,double i_a,double i_b,double i_c)
{
  double logf;

  logf = i_log_f(i_x,i_par)-i_fopt-((i_a*i_x+i_b)*i_x+i_c);

  return(logf);
}		/* end of log_f_ratio */



/*L:log_f_ratio1*
________________________________________________________________

		log_f_ratio1
________________________________________________________________

Name:		log_f_ratio1
øSyntax:		
Description:
Side effects:   
Return value:   None
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double log_f_ratio1(double x,double (*i_log_f1)(double,double *),
                           double *i_par,double i_a,double i_b)
{
  double logf1;

  logf1 = i_log_f1(x,i_par)-G_TWO*i_a*x-i_b;

  return(logf1);
}		/* end of log_f_ratio1 */

