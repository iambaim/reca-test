/* lapack-interface.c
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
  \file lapack-interface.c
  \brief The interface towards the LAPACK routines written in fortran.
*/
static const char RCSId[] = "$Id: lapack-interface.c 1 2013-03-28 13:54:24Z hanne $";

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "GMRFLib.h"
#include "GMRFLibP.h"

int GMRFLib_comp_posdef_inverse(double *matrix, int dim)
{
    /* 
       overwrite a symmetric MATRIX with its inverse
    */
    int info=0, i, j;
    
    GMRFLib_ASSERT(GMRFLib_blas_level == BLAS_LEVEL2 || GMRFLib_blas_level == BLAS_LEVEL3, GMRFLib_EPARAMETER);
    switch(GMRFLib_blas_level)
    {
    case BLAS_LEVEL2:
	dpotf2_("L", &dim, matrix, &dim, &info, 1);
	break;
    case BLAS_LEVEL3:
	dpotrf_("L", &dim, matrix, &dim, &info, 1);
	break;
    }
    if (info) GMRFLib_ERROR(GMRFLib_ESINGMAT);

    dpotri_("L", &dim, matrix, &dim, &info, 1);
    if (info) GMRFLib_ERROR(GMRFLib_ESINGMAT);

    for(i=0;i<dim;i++)				  /* fill the U-part */
	for(j=i+1;j<dim;j++)
	    matrix[i+j*dim] = matrix[j+i*dim];

    return GMRFLib_SUCCESS;
}
int GMRFLib_comp_chol_semidef(double **chol, int **map, int *rank, double *matrix, int dim, double *logdet, double eps)
{
    /* 
       compute the ``cholesky factorisation'' for a positive semidefinite matrix. return malloc'ed
       factorization in chol, the malloc'ed mapping in map and the rank in *rank.

       if logdet, then compute the log determinant of the non-singular part

       eps is the smalles L[i,i] 
    */

    double *work, det;
    int job=1, info, i;
    
    GMRFLib_ASSERT(matrix, GMRFLib_EINVARG);
    GMRFLib_ASSERT(dim>0, GMRFLib_EINVARG);

    *chol = Calloc(SQR(dim), double); MEMCHK(*chol);
    *map  = Calloc(dim, int); MEMCHK(*map);
    work  = Calloc(dim, double); MEMCHK(work);

    memcpy(*chol, matrix, SQR(dim)*sizeof(double));

    dchdc_(*chol, &dim, &dim, work, *map, &job, &info, &eps);
    *rank = info;

    for(i=0;i<dim;i++) (*map)[i]--;		  /* convert to C index-ing */
    
    if (logdet)
    {
	for(det=0.0, i=0; i < *rank; i++) det += log((*chol)[i+i*dim]);
	*logdet = 2.0*det;
    }

    FREE(work);
    return GMRFLib_SUCCESS;
}
int GMRFLib_comp_chol_general(double **chol, double *matrix, int dim, double *logdet, int ecode)
{
    /* 
       return a malloc'ed cholesky factorisation of MATRIX in *chol and optional the
       log(determinant). if fail return `ecode'
    
    */
    int info=0, i, j;
    double *a = NULL, det;

    GMRFLib_ASSERT(chol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(dim >= 0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(GMRFLib_blas_level == BLAS_LEVEL2 || GMRFLib_blas_level == BLAS_LEVEL3, GMRFLib_EPARAMETER);
   
    if (!dim) 
    {
	*chol = NULL;
	return GMRFLib_SUCCESS;
    }
    
    a = Calloc(SQR(dim), double); MEMCHK(a);
    memcpy(a, matrix, SQR(dim)*sizeof(double));

    switch(GMRFLib_blas_level)
    {
    case BLAS_LEVEL2:
	dpotf2_("L", &dim, a, &dim, &info, 1);
	break;
    case BLAS_LEVEL3:
	dpotrf_("L", &dim, a, &dim, &info, 1);
	break;
    }
    
    if (info) GMRFLib_ERROR(ecode);

    if (logdet)
    {
	for(det=0.0, i=0;i<dim;i++) det += log(a[i+i*dim]);
	*logdet = 2.0*det;
    }

    for(i=0;i<dim;i++)				  /* set to zero the upper part */
	for(j=i+1;j<dim;j++)
	    a[i+j*dim] = 0.0;

    *chol = a;
    return GMRFLib_SUCCESS;
}
int GMRFLib_solveAxb_posdef(double *sol, double *chol, double *b, int dim, int nrhs)
{
    /* 
       solve Ax=b, where chol is lower Cholesky factor of A.
     */
    int info;
    
    GMRFLib_ASSERT(sol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(chol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(b, GMRFLib_EINVARG);
    GMRFLib_ASSERT(dim>=0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(nrhs>=0, GMRFLib_EINVARG);

    if (sol != b) memcpy(sol, b, dim*nrhs*sizeof(double));
    dpotrs_("L", &dim, &nrhs, chol, &dim, sol, &dim, &info, 1);

    if (info) GMRFLib_ERROR(GMRFLib_EPOSDEF);

    return GMRFLib_SUCCESS;
}
