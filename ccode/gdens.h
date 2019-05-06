/* gdens.h
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
 * RCSId: $Id: gdens.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file gdens.h
  \brief Typedefs and defines for \ref gdens.c
*/

#ifndef __GMRFLib_GDENS_H__
#define __GMRFLib_GDENS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#if !defined(__FreeBSD__)

#endif
#include <stdlib.h>

typedef struct
{
    double x, xl, xm, xr, dx, b, c, lp, lfl, lfr, lfm, lfmax, lnormc;
}
GDensElm;

typedef struct
{
    int n;
    int order;
    double xmin, xmax;
    GDensElm *elm;
}
GDens;

typedef double GDensFunc_tp(double, char *);

int cdfnor_(int *, double *, double *, double *, double *, double *, int *, double *);

double erf(double x);				  /* needed? */
double erfi(double x);
double lerf_diff(double x1, double x0);
double lerfi_diff(double x1, double x0);
int GDens1LDens(GDens *ptr, double *x, double *ldens);
int GDens1Sample(GDens *ptr, double *x, double *ldens);
int GDens1Eval(GDens *ptr, double *x, double *ldens, int flag);
GDens *GDens1Init(double xmin, double xmax, int n, int n_min,
		  GDensFunc_tp *lfunc, char *lfunc_arg);
int GDens1Update(GDens *ptr);



int GDens2LDens(GDens *ptr, double *x, double *ldens);
int GDens2Sample(GDens *ptr, double *x, double *ldens);
int GDens2Eval(GDens *ptr, double *x, double *ldens, int flag);
GDens *GDens2Init(double xmin, double xmax, int n, int n_min,
		  GDensFunc_tp *lfunc, char *lfunc_arg);
GDens *GDens2InitNew(double fac, double mean, double stdev, int n, GDensFunc_tp *lfunc, char *lfunc_arg);
int GDens2Update(GDens *ptr);

int GDensLDens(GDens *ptr, double *x, double *ldens);
int GDensSample(GDens *ptr, double *x, double *ldens);
int GDensEval(GDens *ptr, double *x, double *ldens, int flag);
GDens *GDensInit(double xmin, double xmax, int n, int n_min, int order,
		 GDensFunc_tp *lfunc, char *lfunc_arg);
GDens *GDensInitNew(double fac, double mean, double stdev, int n, GDensFunc_tp *lfunc, char *lfunc_arg);

int GDensFree(GDens *ptr);
int GDensElmCompare(const void *e, const void *ee);

__END_DECLS
#endif
