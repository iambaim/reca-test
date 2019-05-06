/* gdens.c
 * 
 * Copyright (C) 2001 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The author's contact information:
 *
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */


/*!
  \file gdens.c
  \brief Functions for constructing and working on log-linear and log-quadratic approximations to
  densities.
  
  \note This code could benefite from rewriting, such that each element is switched to a linear
   instead of quadratic it probability for the cell is low. The code is not that efficient nor
   robust as we would like. Linear extrapolation in the tails should be optional.
*/

#include <stddef.h>
#include <float.h>
#if !defined(__FreeBSD__)

#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "GMRFLib.h"
#include "GMRFLibP.h"

#define TWODIVSQRTPI    1.1283791670955125738961589031215452  /* 2/sqrt(pi) */

static const char RCSId[] = "$Id: gdens.c 1 2013-03-28 13:54:24Z hanne $";


/* 
   note: the representation is:  exp(c*x^2+b*x) 
 */

double erfi(double x)
{
    return erfi_(&x);
}
double lerf_diff(double x1, double x0)
{
    /* 
       compute log(erf(x1)-erf(x0)) safely, well, at least when x0 is high and x1-x0 is small, or
       -x0 is high and x1-x0 is small.

       we know that if x0 is high, then x1-x0 are close (hopefully), so if we're out of the limits
       for the 'erf' function, just before it breaks down numerically, we switch to an
       approximation, for high x0 and small (x1-x0).
       
    */
    double limit = 5.0, xmin, dx, value;
    
    if (x1 == x0)
	return -DBL_MAX;
    else
    {
	if (ABS(x1) < limit || ABS(x0) < limit)
	    return log((double)(erf(x1)-erf(x0)));
	else
	{
	    xmin  = MIN(x0, x1);
	    dx    = MAX(x0, x1) -xmin;
	    value = -SQR(xmin) + log((1.0-exp(-2.*dx*xmin))/(2.*xmin)) + log(TWODIVSQRTPI);
	}
    }
    return value;
}
double lerfi_diff(double x1, double x0)
{
    /* 
       same comment apply as for lerf_diff
     */
    double limit = 25.0, xmin, dx, value;

    if (x1 == x0)
    {
	return -DBL_MAX;
    }
    else
    {
	if (ABS(x1) < limit || ABS(x0) < limit)
	    return log(erfi(x1)-erfi(x0));
	else
	{
	    xmin  = MIN(x0, x1);
	    dx    = MAX(x0, x1) -xmin;
	    value = SQR(xmin) + log((exp(2.*dx*xmin)-1.0)/(2.*xmin)) + log(TWODIVSQRTPI);
	}
    }
    return value;
}
int GDensElmCompare(const void *e, const void *ee)
{
    GDensElm *elm, *eelm;

    elm  = (GDensElm *)e;
    eelm = (GDensElm *)ee;

    if (elm->lp > eelm->lp) return -1;
    if (elm->lp < eelm->lp) return  1;

    return 0;
}
int GDensFree(GDens *ptr)
{
    if (ptr)
    {
	FREE(ptr->elm);
	FREE(ptr);
    }
    return 0;
}
int GDens1Update(GDens *ptr)
{
    /* 
       update ptr
     */
    double lfmax, lpmax, lpsum, psum;
    int i, n;

    n = ptr->n;

    /* 
       this is just for scaling to prevent overflow/underflow
     */
    for(lfmax = ptr->elm[0].lfl, i=0;i<ptr->n;i++) 
    {
	lfmax = MAX(lfmax, ptr->elm[i].lfr);
	lfmax = MAX(lfmax, ptr->elm[i].lfl);
    }
    
    /* 
       compute various constants
     */
    for(i=0;i<n;i++) ptr->elm[i].b = (ptr->elm[i].lfr-ptr->elm[i].lfl)/ptr->elm[i].dx;

    for(i=0;i<n;i++)
    {
	if (ISZERO(ptr->elm[i].b))
	    ptr->elm[i].lp = ptr->elm[i].lfl -lfmax + log(ptr->elm[i].dx);
	else
	{
	    double arg = (exp(ptr->elm[i].dx*ptr->elm[i].b)-1.0)/ptr->elm[i].b;

	    if (ISZERO(arg) || arg < 0.)
		ptr->elm[i].lp = ptr->elm[i].lfl -lfmax + log(ptr->elm[i].dx);
	    else
		ptr->elm[i].lp = ptr->elm[i].lfl -lfmax + log(arg);
	}
    }

    for(lpmax = ptr->elm[0].lp, i=1;i<n;i++) lpmax = MAX(lpmax, ptr->elm[i].lp);
    for(i=0;i<n;i++) ptr->elm[i].lp -= lpmax;
    for(i=0, psum=0.0;i<n;i++) psum += exp(ptr->elm[i].lp);
    for(lpsum = log(psum), i=0;i<n;i++) ptr->elm[i].lp -= lpsum;
    
    /* 
       speedup the search for regions and sampling
     */
    qsort(ptr->elm, (size_t) ptr->n, sizeof(GDensElm), GDensElmCompare);

    return 0;
}
int GDens2Update(GDens *ptr)
{
    /* 
       update ptr
     */
    double lpmax, lpsum, psum;
    int i, n = ptr->n;
    
    /* 
       compute various constants
     */
    for(i=0;i<n;i++) 
    {
	double s,
	    x0 = ptr->elm[i].xl,  x1 = ptr->elm[i].xm,  x2 = ptr->elm[i].xr,
	    f0 = ptr->elm[i].lfl, f1 = ptr->elm[i].lfm, f2 = ptr->elm[i].lfr;

	ptr->elm[i].lfmax = f0;			  /* scale and shift */
	f1 -= f0; f2 -= f0; f0 = 0.0;
	x1 -= x0; x2 -= x0; x0 = 0.0;
	
	s = 1./((-SQR(x2) + x2*x1)*x1);
	ptr->elm[i].b = -s*(-f2*SQR(x1) +f1*SQR(x2));
	ptr->elm[i].c = s*(-f2*x1 +f1*x2);
    }
    
    /* 
       compute the integral over each cell, unnormalized
     */
    for(i=0;i<n;i++)
    {
	if (ISZERO(ptr->elm[i].c))
	{
	    if (ISZERO(ptr->elm[i].b))
		ptr->elm[i].lnormc = log(ptr->elm[i].xr-ptr->elm[i].xl);
	    else
	    {
		double tmp = (ptr->elm[i].xr-ptr->elm[i].xl)*ptr->elm[i].b;
		if (fabs(tmp) > FLT_EPSILON)
		{
		    /* 
		       all ok
		    */
		    ptr->elm[i].lnormc = log(fabs(exp((ptr->elm[i].xr-ptr->elm[i].xl)*ptr->elm[i].b) -1.0))
			- log(fabs(ptr->elm[i].b));
		}
		else
		{
		    /* 
		       this is the expansion for small 'tmp'
		    */
		    ptr->elm[i].lnormc = log(fabs(tmp)) - log(fabs(ptr->elm[i].b));
		}
	    }
	}
	else
	{
	    /* 
	       this can (and does) go wrong as the 'erf' and 'erfi' function may fail for high
	       values of the arguments. perhaps we need to intergrate directly
	       int(exp(x^2),x=low...high) ?
	     */
	    
	    double aa = 0.5*ptr->elm[i].b/sqrt(fabs(ptr->elm[i].c));
	    
	    if (ptr->elm[i].c <= 0.0)
		ptr->elm[i].lnormc = -0.25*SQR(ptr->elm[i].b)/ptr->elm[i].c +
		    log(1./(TWODIVSQRTPI*sqrt(-ptr->elm[i].c)))
		    + lerf_diff(sqrt(-ptr->elm[i].c)*(ptr->elm[i].xr-ptr->elm[i].xl) -aa, -aa);
	    else
		ptr->elm[i].lnormc = -0.25*SQR(ptr->elm[i].b)/ptr->elm[i].c +
		    log(1./(TWODIVSQRTPI*sqrt(ptr->elm[i].c)))
		    + lerfi_diff(sqrt(ptr->elm[i].c)*(ptr->elm[i].xr-ptr->elm[i].xl) +aa, aa);
	}
	if (ISINF(ptr->elm[i].lnormc))
	{
	    ptr->elm[i].lnormc = ptr->elm[i].b = ptr->elm[i].c = 0;
	    ptr->elm[i].lp     = ptr->elm[i].lfl = ptr->elm[i].lfm = ptr->elm[i].lfr = ptr->elm[i].lfmax = -FLT_MAX;
	}
	else
	{
	    ptr->elm[i].lp = ptr->elm[i].lnormc + ptr->elm[i].lfmax;
	}
    }

    for(lpmax = ptr->elm[0].lp, i=1;i<n;i++) lpmax = MAX(lpmax, ptr->elm[i].lp);
    for(i=0;i<n;i++) ptr->elm[i].lp -= lpmax;
    for(i=0, psum=0.0;i<n;i++) psum += exp(ptr->elm[i].lp);
    for(lpsum = log(psum), i=0;i<n;i++) ptr->elm[i].lp -= lpsum;
    
    /* 
       speedup the search for regions and sampling
     */
    qsort(ptr->elm, (size_t) ptr->n, sizeof(GDensElm), GDensElmCompare);

    if (0)
    {
	for(i=0;i<n;i++)
	{
	    printf("Gdens2 region %d\n", i);
	    printf("\txl %f xr %f\n", ptr->elm[i].xl, ptr->elm[i].xr);
	    printf("\tlnormc %f lp %f\n", ptr->elm[i].lnormc, ptr->elm[i].lp);
	}
    }

    return 0;
}
GDens *GDens1Init(double xmin, double xmax, int n, int n_min, GDensFunc_tp *lfunc, char *lfunc_arg)
{
    /* 
       arguments are:

       xmin : left  limit
       xmax : right limit
       n    : total number of interior points, which makes n+1 cells
       n_min: make first n_min interior points by equal division of (xmax-xmin), then
              n-n_min divisions by dividing the cell with the largest probability in half.
       lfunc: the function returning log(f(x)), unormalized

     */
    
    GDens *ptr;
    int i, nn;
    double xnew, f, df, dx;
    
    ptr       = Calloc(1, GDens);
    ptr->order= 1;
    ptr->xmin = MIN(xmin, xmax);
    ptr->xmax = MAX(xmin, xmax);
    ptr->n    = MAX(2, n_min);
    nn        = MAX(ptr->n, n);
    ptr->elm  = Calloc(nn+1, GDensElm);

    /* 
       make first equal spaced points
     */
    dx = (ptr->xmax -ptr->xmin)/((double)n_min);
    for(i=0;i<ptr->n;i++)
    {
	ptr->elm[i].x   = ptr->xmin + i*dx;
	ptr->elm[i].dx  = dx;
	ptr->elm[i].lfl = lfunc(ptr->elm[i].x, lfunc_arg);
    }
    for(i=0;i<ptr->n-1;i++) ptr->elm[i].lfr = ptr->elm[i+1].lfl;
    ptr->elm[ptr->n-1].lfr = lfunc(ptr->elm[ptr->n-1].x + ptr->elm[ptr->n-1].dx, lfunc_arg);

    GDens1Update(ptr);

    /* 
       decide now the rest of the points by dividing the cell with largest probability in half
     */
    for(i=n_min;i<n;i++)
    {
	/* 
	   split the first elm, and add one at the end
	 */
	if (ABS(ptr->elm[0].b) < FLT_EPSILON)
	    f = 0.5;
	else
	{
	    df = ptr->elm[0].b*ptr->elm[0].dx;
	    f  = log(0.5*(exp(df)+1.0))/df;
	}

	xnew                 = ptr->elm[0].x + ptr->elm[0].dx*f;
	ptr->elm[ptr->n].lfr = ptr->elm[0].lfr;
	ptr->elm[ptr->n].x   = xnew;
	ptr->elm[ptr->n].dx  = ptr->elm[0].dx*(1.-f);
	ptr->elm[0].lfr      = ptr->elm[ptr->n].lfl = lfunc(xnew, lfunc_arg);
	ptr->elm[0].dx      *= f;
	
	ptr->n++;
	GDens1Update(ptr);
    }

    return ptr;
}
GDens *GDens2Init(double xmin, double xmax, int n, int n_min, GDensFunc_tp *lfunc, char *lfunc_arg)
{
    /* 
       arguments are:

       xmin : left  limit
       xmax : right limit
       n    : total number of interior points, which makes n+1 cells
       n_min: make first n_min interior points by equal division of (xmax-xmin), then
              n-n_min divisions by dividing the cell with the largest probability in half.
       lfunc: the function returning log(f(x,lfunc_arg)), unormalized
     */
    
    GDens *ptr;
    int i, nn;
    double dx, df, f;
    
    ptr       = Calloc(1, GDens);
    ptr->order= 2;
    ptr->xmin = MIN(xmin, xmax);
    ptr->xmax = MAX(xmin, xmax);
    ptr->n    = MAX(1, n_min);
    nn        = MAX(ptr->n, n);
    ptr->elm  = Calloc(nn+1, GDensElm);

    /* 
       make first equal spaced points
     */
    dx = (ptr->xmax -ptr->xmin)/((double)n_min);
    for(i=0;i<ptr->n;i++)
    {
	ptr->elm[i].xl  = ptr->xmin + i*dx;
	ptr->elm[i].xr  = ptr->elm[i].xl + dx;
	ptr->elm[i].lfl = (i==0 ? lfunc(ptr->elm[i].xl, lfunc_arg) : ptr->elm[i-1].lfr);	
	ptr->elm[i].lfr = lfunc(ptr->elm[i].xr, lfunc_arg);

	/* 
	   now, try to decide a nice place to place the midpoint. currently, i try to split the
	   region into (approximately) equal integrals
	*/
	df = ptr->elm[i].lfr - ptr->elm[i].lfl;
	f  = (ISZERO(df) ? 0.5 : log(0.5*(exp(df)+1.0))/df);

	/* 
	   need also a more safe option than ISZERO
	 */
	if (ISZERO(f) || ISZERO(1.-f)) f = 0.5;
	if (f < FLT_EPSILON || 1.-f < FLT_EPSILON) f = 0.5;
	
	if (f == 0.5)
	{
	    /* 
	       this will emulate a log-linear spline
	    */
	    ptr->elm[i].xm = ptr->elm[i].xl + dx*f;
	    ptr->elm[i].lfm = 0.5*(ptr->elm[i].lfl + ptr->elm[i].lfr);
	}
	else
	{
	    ptr->elm[i].xm  = ptr->elm[i].xl + dx*f;
	    ptr->elm[i].lfm = lfunc(ptr->elm[i].xm, lfunc_arg);
	}
    }
    for(i=0;i<ptr->n-1;i++) ptr->elm[i].lfr = ptr->elm[i+1].lfl;
    ptr->elm[ptr->n-1].lfr = lfunc(ptr->elm[ptr->n-1].xr, lfunc_arg);

    GDens2Update(ptr);

    /* 
       decide now the rest of the points by dividing the cell with largest probability in half. try
       to split the region into (approximately) equal integrals
     */
    for(i=n_min;i<nn;i++)
    {
	/* 
	   split the first elm, and add one at the end
	 */
	ptr->elm[ptr->n].lfl = ptr->elm[0].lfm;
	ptr->elm[ptr->n].xl  = ptr->elm[0].xm;
	ptr->elm[ptr->n].lfr = ptr->elm[0].lfr;
	ptr->elm[ptr->n].xr  = ptr->elm[0].xr;
	df = ptr->elm[ptr->n].lfr - ptr->elm[ptr->n].lfl;
	f  = (ISZERO(df) ? 0.5 : log(0.5*(exp(df)+1.0))/df);
	if (ISZERO(f) || ISZERO(1.-f)) f = 0.5;
	if (f < FLT_EPSILON || 1.-f < FLT_EPSILON) f = 0.5;
	if (f == 0.5)
	{
	    /* 
	       this will emulate a log-linear spline
	    */
	    ptr->elm[ptr->n].xm  = ptr->elm[ptr->n].xl + f*(ptr->elm[ptr->n].xr -ptr->elm[ptr->n].xl);
	    ptr->elm[ptr->n].lfm = 0.5*(ptr->elm[ptr->n].lfl + ptr->elm[ptr->n].lfr);
	}
	else
	{
	    ptr->elm[ptr->n].xm  = ptr->elm[ptr->n].xl + f*(ptr->elm[ptr->n].xr -ptr->elm[ptr->n].xl);
	    ptr->elm[ptr->n].lfm = lfunc(ptr->elm[ptr->n].xm, lfunc_arg);
	}

	ptr->elm[0].lfr = ptr->elm[0].lfm;
	ptr->elm[0].xr  = ptr->elm[0].xm;
	df = ptr->elm[0].lfr - ptr->elm[0].lfl;
	f  = (ISZERO(df) ? 0.5 : log(0.5*(exp(df)+1.0))/df);
	if (ISZERO(f) || ISZERO(1.-f)) f = 0.5;
	if (f < FLT_EPSILON || 1.-f < FLT_EPSILON) f = 0.5;
	if (f == 0.5)
	{
	    /* 
	       this will emulate a log-linear spline
	    */
	    ptr->elm[0].xm  = ptr->elm[0].xl + f*(ptr->elm[0].xr -ptr->elm[0].xl);
	    ptr->elm[0].lfm = 0.5*(ptr->elm[0].lfl + ptr->elm[0].lfr);
	}
	else
	{
	    ptr->elm[0].xm  = ptr->elm[0].xl + f*(ptr->elm[0].xr -ptr->elm[0].xl);
	    ptr->elm[0].lfm = lfunc(ptr->elm[0].xm, lfunc_arg);
	}
	
	ptr->n++;
	GDens2Update(ptr);
    }

    return ptr;
}
double *cut_points(double fac, int n)
{
    /* 
       return the n+1 cutpoints between -fac...fac, that makes n regions of equal probability
    */

    static double fac_save = 0.0;
    static int    n_save   = 0;
    static double *cp      = NULL;

    double p_left, p, q, mean=0.0, sd=1.0, bound, x;
    int i, iwhich, status;

    if (fac == fac_save && n == n_save) return cp;

    fac_save = fac;
    n_save   = n;
    
    FREE(cp);
    cp = Calloc(n+1, double);

    iwhich = 1;
    cdfnor_(&iwhich, &p, &q, &fac, &mean, &sd, &status, &bound);
    p_left = q;

    cp[0] = -fac;
    cp[n] =  fac;
    for(i=1;i<n;i++)
    {
	iwhich = 2;
	p = p_left + i/((double)n)*(1.0 - 2.0*p_left);
	q = 1.0 -p;
	cdfnor_(&iwhich, &p, &q, &x, &mean, &sd, &status, &bound);
	cp[i] = x;
    }

    if (0)
    {
	for(i=0;i<n+1;i++) printf("cutpoint %d %f\n", i, cp[i]);
    }
    return cp;
}
GDens *GDens2InitNew(double fac, double mean, double stdev, int n, GDensFunc_tp *lfunc, char *lfunc_arg)
{
    /* 
       arguments are:

       limits are mean +- fac*stdev, and n-regions.

       lfunc: the function returning log(f(x,lfunc_arg)), unormalized
     */
    
    GDens *ptr;
    int i;
    double dx, df, f, *cpoints;
    
    fac = ABS(fac);

    ptr       = Calloc(1, GDens);
    ptr->xmin = mean -fac*stdev;
    ptr->xmax = mean +fac*stdev;
    ptr->n    = n;
    ptr->elm  = Calloc(n+1, GDensElm);


    cpoints = cut_points(fac, n);

    for(i=0;i<ptr->n;i++)
    {
	dx              = stdev*(cpoints[i+1]-cpoints[i]);
        ptr->elm[i].xl  = mean + cpoints[i]*stdev;
	ptr->elm[i].xr  = ptr->elm[i].xl + dx;
	ptr->elm[i].lfl = (i==0 ? lfunc(ptr->elm[i].xl, lfunc_arg) : ptr->elm[i-1].lfr);	
	ptr->elm[i].lfr = lfunc(ptr->elm[i].xr, lfunc_arg);

	/* 
	   now, try to decide a nice place to place the midpoint. currently, i try to split the
	   region into (approximately) equal integrals
	*/
	df = ptr->elm[i].lfr - ptr->elm[i].lfl;
	f  = (ISZERO(df) ? 0.5 : log(0.5*(exp(df)+1.0))/df);

	/* 
	   need also a more safe option than ISZERO
	 */
	if (ISZERO(f) || ISZERO(1.-f)) f = 0.5;
	if (f < FLT_EPSILON || 1.-f < FLT_EPSILON) f = 0.5;
	
	if (f == 0.5)
	{
	    /* 
	       this will emulate a log-linear spline
	    */
	    ptr->elm[i].xm = ptr->elm[i].xl + dx*f;
	    ptr->elm[i].lfm = 0.5*(ptr->elm[i].lfl + ptr->elm[i].lfr);
	}
	else
	{
	    ptr->elm[i].xm  = ptr->elm[i].xl + dx*f;
	    ptr->elm[i].lfm = lfunc(ptr->elm[i].xm, lfunc_arg);
	}
    }
    for(i=0;i<ptr->n-1;i++) ptr->elm[i].lfr = ptr->elm[i+1].lfl;
    ptr->elm[ptr->n-1].lfr = lfunc(ptr->elm[ptr->n-1].xr, lfunc_arg);

    GDens2Update(ptr);

    return ptr;
}
GDens *GDensInitNew(double fac, double mean, double stdev, int n, GDensFunc_tp *lfunc, char *lfunc_arg)
{
    return GDens2InitNew(fac, mean, stdev, n, lfunc, lfunc_arg);
}
GDens *GDensInit(double xmin, double xmax, int n, int n_min, int order, GDensFunc_tp *lfunc, char *lfunc_arg)
{
    if (!(order == 1 || order == 2))
    {
	GMRFLib_ERROR_NO_RETURN(GMRFLib_EPARAMETER);
	return NULL;
    }

    if (order == 1)
	return GDens1Init(xmin, xmax, n, n_min, lfunc, lfunc_arg);
    else
	return GDens2Init(xmin, xmax, n, n_min, lfunc, lfunc_arg);
}
int GDens1Eval(GDens *ptr, double *x, double *ldens, int flag)
{
    int i, r;
    double psum, p, xx;

    if (flag)					  /* sample, draw first segment then within the segment */
    {
	p = GMRFLib_uniform();
	r = ptr->n-1;
	for(i=0, psum=0.0;i<ptr->n-1;i++)          
	    if ((psum += exp(ptr->elm[i].lp)) > p)
	    {
		r = i; break;
	    }
	p = GMRFLib_uniform();
	if (ISZERO(ptr->elm[r].b))
	    *x = ptr->elm[r].x + p*ptr->elm[r].dx;
	else
	    *x = ptr->elm[r].x + log(1.0-p+p*exp(ptr->elm[r].b*ptr->elm[r].dx))/ptr->elm[r].b;
    }
    else					  /* evaluate the density, find the corresponding segment */
    {
	r = ptr->n-1;
	for(i=0;i<ptr->n-1;i++)          
	    if (*x >= ptr->elm[i].x && *x < ptr->elm[i].x+ptr->elm[i].dx) 
	    {
		r = i; break;
	    }
    }

    if (*x < ptr->xmin || *x > ptr->xmax)
    {
	fprintf(stdout, "\nGDens1: x=%f is not in [%f,%f], truncate!\n", *x, ptr->xmin, ptr->xmax);
	xx = MAX(ptr->xmin, MIN(ptr->xmax, *x));
    }
    else
	xx = *x;

    if(ISZERO(ptr->elm[r].b))
	*ldens = ptr->elm[r].lp - log(ptr->elm[r].dx);
    else
    {
	double arg = ptr->elm[r].b*exp(ptr->elm[r].b*(xx -ptr->elm[r].x))
		      /(exp(ptr->elm[r].b*ptr->elm[r].dx)-1.0);

	if (ISZERO(arg) || arg <= 0.0)
	    *ldens = ptr->elm[r].lp - log(ptr->elm[r].dx);
	else
	    *ldens = ptr->elm[r].lp + log(arg);
    }

    return 0;
}
int GDens2Eval(GDens *ptr, double *x, double *ldens, int flag)
{
    int i, r;
    double psum, p, xx, xlocalmax, lfmax, b, dx;

    if (flag)					  /* sample, draw first segment then within the segment */
    {
	p = GMRFLib_uniform();
	r = ptr->n-1;
	for(i=0, psum=0.0;i<ptr->n-1;i++)          
	    if ((psum += exp(ptr->elm[i].lp)) > p)
	    {
		r = i; break;
	    }
	/* 
	   here i use a rejection sampling algorithm with two proposals: a uniform and a linear in
	   log-scale. pick your choice ;-)
	 */
	if (1)
	{
	    /* 
	       uniform
	     */
	    dx = ptr->elm[r].xr-ptr->elm[r].xl;
	    if (ISZERO(ptr->elm[r].c))
		lfmax = MAX(ptr->elm[r].lfr, ptr->elm[r].lfl) -ptr->elm[r].lfmax;
	    else
	    {
		xlocalmax = -ptr->elm[r].b/(2.0*ptr->elm[r].c);
		xlocalmax = MAX(0.0, xlocalmax);
		xlocalmax = MIN(dx, xlocalmax);
		lfmax     = ptr->elm[r].c*SQR(xlocalmax)+ptr->elm[r].b*xlocalmax;
	    }
	    do
		xx = GMRFLib_uniform()*dx;
	    while(GMRFLib_uniform() > exp(ptr->elm[r].c*SQR(xx)+ptr->elm[r].b*xx-lfmax));
	    *x = xx + ptr->elm[r].xl;
	}
	else
	{
	    /* 
	       linear
	     */
	    dx = ptr->elm[r].xr -ptr->elm[r].xl;
	    b  = (ptr->elm[r].lfr -ptr->elm[r].lfl)/dx;

	    if (ISZERO(ptr->elm[r].c))
	    {
		p = GMRFLib_uniform();
		xx = log(1.0-p+p*exp(b*dx))/b;	  /* are on a straight line */
	    }
	    else
	    {
		xlocalmax = -(ptr->elm[r].b -b)/(2.0*ptr->elm[r].c);
		if (ptr->elm[r].c < 0.0)
		{
		    xlocalmax = MAX(0.0, xlocalmax);
		    xlocalmax = MIN(dx, xlocalmax);
		    lfmax     = ptr->elm[r].c*SQR(xlocalmax)+(ptr->elm[r].b -b)*xlocalmax;
		}
		else 
		{
		    double lf_left, lf_right;
		    lf_left = 0.0;
		    lf_right = ptr->elm[r].c*SQR(dx)+(ptr->elm[r].b -b)*dx;
		    lfmax = (lf_right > lf_left ? lf_right : lf_left);
		}

		do
		{
		    p = GMRFLib_uniform();
		    if (ISZERO(b))
			xx = dx*p;
		    else
			xx = log(1.0-p+p*exp(b*dx))/b;
		}
		while(GMRFLib_uniform() > exp(ptr->elm[r].c*SQR(xx)+(ptr->elm[r].b-b)*xx-lfmax));
	    }
	    *x = xx + ptr->elm[r].xl;
	}
    }
    else					  /* evaluate the density, find the corresponding segment */
    {
	r = ptr->n-1;
	for(i=0;i<ptr->n-1;i++)          
	    if (*x >= ptr->elm[i].xl && *x < ptr->elm[i].xr) 
	    {
		r = i; break;
	    }
    }

    if (*x < ptr->xmin || *x > ptr->xmax)
    {
	if (0) fprintf(stdout, "\nGDens2: x=%f is not in [%f,%f], truncate!\n", *x, ptr->xmin, ptr->xmax);
	xx = MAX(ptr->xmin, MIN(ptr->xmax, *x));  /* only happens in eval-mode */
    }
    else
    {
	xx = *x;
    }

    xx -= ptr->elm[r].xl;			  /* relative value within the cell */

    *ldens = ptr->elm[r].lp + ptr->elm[r].c*SQR(xx) + ptr->elm[r].b*xx -ptr->elm[r].lnormc;

    return 0;
}
int GDensEval(GDens *ptr, double *x, double *ldens, int flag)
{
    if (ptr->order == 1)
	return GDens1Eval(ptr, x, ldens, flag);
    else
	return GDens2Eval(ptr, x, ldens, flag);
}
int GDens1Sample(GDens *ptr, double *x, double *ldens)
{
    return GDens1Eval(ptr, x, ldens, 1);
}
int GDens2Sample(GDens *ptr, double *x, double *ldens)
{
    return GDens2Eval(ptr, x, ldens, 1);
}
int GDensSample(GDens *ptr,  double *x, double *ldens)
{
    if (ptr->order == 1)
	return GDens1Eval(ptr, x, ldens, 1);
    else
	return GDens2Eval(ptr, x, ldens, 1);
}
int GDens1LDens(GDens *ptr, double *x, double *ldens)
{
    return GDens1Eval(ptr, x, ldens, 0);
}
int GDens2LDens(GDens *ptr, double *x, double *ldens)
{
    return GDens2Eval(ptr, x, ldens, 0);
}
int GDensLDens(GDens *ptr,  double *x, double *ldens)
{
    if (ptr->order == 1)
	return GDens1Eval(ptr, x, ldens, 0);
    else
	return GDens2Eval(ptr, x, ldens, 0);
}
