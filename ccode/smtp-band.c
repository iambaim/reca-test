/* GMRFLib-smtp-band.c
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
  \file smtp-band.c
  \brief The implementation of the interface towards LAPACK's band solver.
*/

#include "GMRFLib.h"
#include "GMRFLibP.h"

static const char RCSId[] = "$Id: smtp-band.c 1 2013-03-28 13:54:24Z hanne $";

int GMRFLib_compute_reordering_BAND(int **remap, GMRFLib_graph_tp *graph)
{
    /* 
       compute the reordering from the graph using the routine in acm582.F
    */
    int i, j, lconnec, bandwidth, profile, error, space, ioptpro, worklen,
	*rstart, *connec, *degree, *work, simple;
    
    if (!graph || !graph->n) return GMRFLib_SUCCESS;
    
    /* 
       check if we have a simple solution --> no neigbours
    */
    for(i=0, simple=1; i<graph->n && simple; i++) simple = (graph->nnbs[i] > 0 ? 0 : 1);
    if (simple)
    {
	int *imap;
	imap = Calloc(graph->n, int); MEMCHK(imap);
	for(i=0;i<graph->n;i++) imap[i] = i;
	*remap = imap;
	return GMRFLib_SUCCESS;
    }

    /* 
       task I, reformat the graph to fit the fortran-routines
    */
    lconnec = 0;
    for(i=0;i<graph->n;i++) lconnec += graph->nnbs[i];
    if (lconnec == 0)
    {
	/* 
	   no connections in the graph, use the identity-map.
	*/
	*remap = Calloc(graph->n, int); MEMCHK(*remap);
	for(i=0;i<graph->n;i++) (*remap)[i] = i;
	return GMRFLib_SUCCESS;
    }
	
    connec = Calloc(lconnec, int); MEMCHK(connec);
    rstart = Calloc(graph->n, int);MEMCHK(rstart);
    degree = graph->nnbs;			  /* yes! */

    rstart[0] = 1;				  /* fortran indx'ing */
    for(i=1;i<graph->n;i++) rstart[i] = degree[i-1] + rstart[i-1];
    for(i=0;i<graph->n;i++)
	for(j=0;j<degree[i];j++)
	    connec[rstart[i]-1 + j] = graph->nbs[i][j]+1; /* fortran indx'ing */

    worklen = 6*graph->n+3;			  /* maximum over all graphs */
    work    = Calloc(worklen,int); MEMCHK(work);

    *remap = Calloc(graph->n, int); MEMCHK(*remap);
    for(i=0;i<graph->n;i++) (*remap)[i] = i+1; /* fortran indx'ing */
	
    ioptpro = 0;				  /* do bandwidth reduction */
    error   = space = 0;
    gpskca_(&graph->n, degree, rstart, connec, &ioptpro, &worklen,
	    *remap, work, &bandwidth, &profile, &error, &space);

    for(i=0;i<graph->n;i++) (*remap)[i]--;	  /* correct for fortran indx'ing */

    FREE(work);
    FREE(rstart);
    FREE(connec);

    if (error) GMRFLib_ERROR(GMRFLib_EREORDER);
    return GMRFLib_SUCCESS;
}
int GMRFLib_free_reordering_BAND(GMRFLib_sm_fact_tp *sm_fact)
{
    if (sm_fact) 
    {
	FREE(sm_fact->remap);
	sm_fact->bandwidth = 0;
    }
    return GMRFLib_SUCCESS;
}
int GMRFLib_build_sparse_matrix_BAND(double **bandmatrix, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, GMRFLib_graph_tp *graph, int *remap, int bandwidth)
{
#define BIDX(i,j) ((i)+(j)*nrow)		  /* band index'ing */

    /* 
       return a band-matrix BMATRIX in L-storage defining the precision matrix
    */

    int i, j, jj, ncol, nrow, node, nnode;

    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);
    GMRFLib_ASSERT(bandmatrix, GMRFLib_EINVARG);

    ncol        = graph->n;
    nrow        = bandwidth + 1;
    *bandmatrix = Calloc(ncol*nrow, double); MEMCHK(*bandmatrix);

    for(i=0;i<graph->n;i++)
    {
	node = remap[i];
	(*bandmatrix)[BIDX(0,node)] = (*Qfunc)(i, i, Qfunc_arg);

	for(j=0;j<graph->nnbs[i];j++)
	{
	    jj    = graph->nbs[i][j];
	    nnode = remap[jj];
	    if (nnode > node) (*bandmatrix)[BIDX(nnode-node, node)] = (*Qfunc)(i, jj, Qfunc_arg);
	}
    }
    return GMRFLib_SUCCESS;
#undef BIDX
}
int GMRFLib_factorise_sparse_matrix_BAND(double *band, GMRFLib_fact_info_tp *finfo, GMRFLib_graph_tp *graph, int bandwidth)
{
#define BIDX(i,j) ((i)+(j)*ldim)		  /* band index'ing */
    /* 
       compute the factorisation of 'band' and overwrite it with the Cholesky-factorisation
     */

    int error = 0, nband, ldim, i, j, k;

    GMRFLib_ASSERT((!(!band && (graph && graph->n > 0))), GMRFLib_EINVARG);
    if (!band && (!graph || (graph && graph->n == 0))) return GMRFLib_SUCCESS;
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    GMRFLib_ASSERT(GMRFLib_blas_level == BLAS_LEVEL2 || GMRFLib_blas_level == BLAS_LEVEL3, GMRFLib_EPARAMETER);
    
    nband = bandwidth;
    ldim  = bandwidth + 1;
    
    switch(GMRFLib_blas_level)
    {
    case BLAS_LEVEL2:
	dpbtf2_("L", &(graph->n), &nband, band, &ldim, &error, 1);
	break;
    case BLAS_LEVEL3:
	dpbtrf_("L", &(graph->n), &nband, band, &ldim, &error, 1);
	break;
    }
    if (error) GMRFLib_ERROR(GMRFLib_EPOSDEF);


    /* 
       provide some info about the factorization
    */
    for(i=0, k=0;i<graph->n;i++) k += graph->nnbs[i];

    finfo->n      = graph->n;			  /* size of Q */
    finfo->nnzero = k + graph->n;		  /* # non-zeros in Q */

    if (0)
    {
	/* 
	   this is wrong, as it will compute 0 and that is what matters
	*/
	for(i=0, k=0;i<graph->n;i++)
	    for(j=i+1;j < MIN(graph->n, i+nband+1); j++) k += !ISZERO(band[BIDX(j-i, i)]);
    }
    else
    {
	/* 
	   this is correct, as it will 'compute' 0's! 
	*/
	for(i=0, k=0;i<graph->n;i++) k += MIN(graph->n, i+nband+1) - (i+1);
    }
    
    finfo->nfillin = k - (finfo->nnzero - graph->n)/2; /* fillin in L */
    
    return GMRFLib_SUCCESS;
#undef BIDX
}
int GMRFLib_free_fact_sparse_matrix_BAND(double *bchol)
{
    FREE(bchol);
    return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth)
{
    /* 
       rhs in real world, bchol in mapped word

       solve L^Tx=rhs, rhs is overwritten by the solution
     */
    int nband, ldim, stride=1;

    GMRFLib_ASSERT(bchol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(rhs, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    
    nband = bandwidth;
    ldim  = nband+1;

    GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    dtbsv_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);

    GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
    return GMRFLib_SUCCESS;
}    
int GMRFLib_solve_llt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth)
{
    /* 
       rhs in real world, bchol in mapped word

       solve  Q x=rhs, where Q=L L^T
     */
    int nband, ldim, stride=1;

    GMRFLib_ASSERT(bchol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(rhs, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);

    nband = bandwidth;
    ldim  = nband+1;

    GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    dtbsv_("L", "N", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
    dtbsv_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
    GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

    return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth, 
			       int findx, int toindx, int remapped)
{
    /* 
       rhs in real world, bchol in mapped world.  solve L^Tx=b backward only from rhs[findx] up to
       rhs[toindx].  note that findx and toindx is in mapped world.  if remapped, do not
       remap/remap-back the rhs before solving.
       
     */
    int nband, ldim, stride=1, from, to;

    from  = findx+1;                              /* convert to fortran indxing */
    to    = toindx+1;                             /* convert to fortran indxing */
    nband = bandwidth;
    ldim  = nband+1;

    if (!remapped) GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    dtbsvspecial_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, &from, &to, 1, 1, 1);
    if (!remapped) GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

    return GMRFLib_SUCCESS;
}
int GMRFLib_comp_cond_meansd_BAND(double *cmean, double *csd, int indx, double *x, int remapped,
				  double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth) 
{
    /* 
       compute the conditonal mean and stdev for x[indx]|x[indx+1]...x[n-1] for the current value of
       x. if `remapped', then x is assumed to be remapped for possible (huge) speedup when used
       repeately.  note: indx is still the user world!!!

       example: approach 1,2,3 are equivalent (old style!)

       (*GMRFLib_uniform_init)(seed); set_stdgauss(x); memcpy(z, x, graph->n*sizeof(double));

       a1: gmrf_g_solve(x, bchol, graph);

       a2: for(i=graph->n-1;i>=0;i--){
               gmrf_g_compute_cmean_csd(&cmean, &csd, i, y, 0, bchol, graph);
               y[i] = cmean + csd*z[i];}

       a3: for(i=graph->n-1;i>=0;i--){
	       gmrf_g_compute_cmean_csd(&cmean, &csd, i, y, 1, bchol, graph);
	       y[graph->remap[i]] = cmean + csd*z[i];}
     */
    
    int nband, ldim, ii;

    GMRFLib_ASSERT(cmean, GMRFLib_EINVARG);
    GMRFLib_ASSERT(csd, GMRFLib_EINVARG);
    GMRFLib_ASSERT(x, GMRFLib_EINVARG);
    GMRFLib_ASSERT(bchol, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    GMRFLib_ASSERT(indx >= 0 && indx < graph->n, GMRFLib_EINVARG);
    
    nband = bandwidth;
    ldim  = nband+1;
    ii    = remap[indx] +1;			  /* fortran indxing */

    if (remapped)
	cmsd_(cmean, csd, &ii, &(graph->n), &nband, bchol, &ldim, x);
    else
    {
	GMRFLib_convert_to_mapped(x, NULL, graph, remap);
	cmsd_(cmean, csd, &ii, &(graph->n), &nband, bchol, &ldim, x);
	GMRFLib_convert_from_mapped(x, NULL, graph, remap);
    }
    return GMRFLib_SUCCESS;
}
int GMRFLib_log_determinant_BAND(double *logdet, double *bchol, GMRFLib_graph_tp *graph, int bandwidth)
{
    int ldim = bandwidth + 1, i;

    for(i=0, *logdet = 0.0; i<graph->n;i++) *logdet += log(bchol[i*ldim]);
    *logdet *= 2;
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_compute_Qinv_BAND(void)
{

    fprintf(stderr, "\n\n%s:%s:%1d: Qinv is not implemented for this solver. \n",
	    __FILE__, __GMRFLib_FuncName, __LINE__);
    fprintf(stderr, "%s:%s:%1d: Use  ``GMRFLib_smtp = GMRFLib_SMTP_TAUCS;'' instead \n",
	    __FILE__, __GMRFLib_FuncName, __LINE__);
    GMRFLib_ERROR(GMRFLib_ESMTP);
    return GMRFLib_ESMTP;
}

/* 
   from here on: undocumented features
*/
/* 
   from here is for internal use only. not documented
*/
int GMRFLib_bitmap_factorisation_BAND__internal(char *filename, double *band, GMRFLib_graph_tp *graph, int *remap, int bandwidth)
{
#define NBitsInByte 8
#define SETBIT(im,jm) { int ii; ii = (im)/NBitsInByte; \
                      GMRFLib_setbit(&bitmap[ ii+(jm)*m], NBitsInByte-1-((im)-ii*NBitsInByte)); }
#define BIDX(i,j) ((i)+(j)*ldim)		  /* band index'ing */
    
    int i, j, n = graph->n, m, nband, ldim;
    unsigned char *bitmap;
    FILE *fp;
    
    nband = bandwidth;
    ldim  = bandwidth + 1;

    m = n/NBitsInByte;
    if (m*NBitsInByte != n) m++;
    bitmap = Calloc(n*m, unsigned char); MEMCHK(bitmap);

    for(i=0;i<graph->n;i++)
	for(j=i;j<MIN(i+nband+1, graph->n);j++)
	    if (!ISZERO(band[BIDX(j-i, i)])) SETBIT(i, j);


    fp = fopen(filename, "w");
    if (fp)
    {
        fprintf(fp, "P4\n%1d %1d\n", n, n);
        for(i=0;i<n;i++) fwrite(&bitmap[i*m], m, 1, fp);
        fclose(fp);
    }
    else
        GMRFLib_ERROR(GMRFLib_EOPENFILE);

    return GMRFLib_SUCCESS;
#undef SETBIT
#undef NBitsInByte
#undef BIDX    
}
int GMRFLib_bitmap_factorisation_BAND(char *filename_body, double *band, GMRFLib_graph_tp *graph, int *remap, int bandwidth)
{
    /* 
       create a bitmap-file of the factorization
    */
    char filename[1024+1];

    sprintf(filename, "%s_L.pbm",(filename_body ? filename_body : "band_L"));
    GMRFLib_bitmap_factorisation_BAND__internal(filename, band, graph, remap, bandwidth);

    return GMRFLib_SUCCESS;
}
