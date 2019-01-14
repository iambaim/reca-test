/* lapack-interfac.h
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
 * RCSId: $Id: lapack-interface.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file lapack-interface.h
  \brief Typedefs and defines for \ref lapack-interface.c
*/

#ifndef __GMRFLib_LAPACK_INTERFACE_H__
#define __GMRFLib_LAPACK_INTERFACE_H__

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
#include <malloc.h>
#endif
#include <stdlib.h>

/*!
  \brief Define BLAS_LEVEL2
*/
#define BLAS_LEVEL2 2
/*!
  \brief Define BLAS_LEVEL3
*/
#define BLAS_LEVEL3 3

int dpbtrf_(const char *, int *, int *, double *, int *, int *, int);
int dpbtf2_(const char *, int *, int *, double *, int *, int *, int);
int dpotrf_(const char *, int *, double *, int *, int *, int);
int dpotf2_(const char *, int *, double *, int *, int *, int);
int dtbsv_(const char *, const char *, const char *, int *, int *, double *, int *, double *, int *,
	   int, int, int);  
int dpotri_(const char *, int *, double *, int *, int *, int);
int dgemm_(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int
	   *, double *, double *, int *, int, int);  
int dgemv_(const char *, int *, int *, double *, double *, int *, double *, int *, double *, double
	   *, int *, int);  
double erfi_(double *);
int dpotrs_(const char *, int *, int *, double *, int *, double *, int *, int *, int);
int dchdc_(double *, int *, int *, double *, int *, int *, int *, double *);
int dtrmv_(const char *, const char *, const char *, int *, double *, int *, double *, int *, int, int, int);


int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps);
int GMRFLib_comp_posdef_inverse(double *matrix, int dim);
int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode);
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs);

__END_DECLS
#endif


