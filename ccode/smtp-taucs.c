/* GMRFLib-smtp-taucs.c
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
  \file smtp-taucs.c
  \brief The interface towards the TAUCS-library
*/

#include <assert.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib.h"
#include "GMRFLibP.h"


static const char RCSId[] = "$Id: smtp-taucs.c 1 2013-03-28 13:54:24Z hanne $";


/* 
   how large should `nset' be before doing `memset' instead of a `for' loop.
*/
#define GMRFLib_NSET_LIMIT(nset, size, n)  MAX(10, (n)/10/(size))


/* 
   First some modified code from the TAUCS library. Checked to be ok version 2.0 and 2.2

   taucs_datatype is set to double in GMRFLib/taucs.h
*/

#if 0						  /* not in use */
taucs_ccs_matrix* my_taucs_dsupernodal_factor_to_ccs(void* vL)
{
    /* 
       this is to be called for a lower triangular double matrix only.

       it includes also zero terms as long as i>=j.
       
    */
    
    supernodal_factor_matrix* L = (supernodal_factor_matrix*) vL;
    taucs_ccs_matrix* C;
    int n,nnz;
    int i,j,ip,jp,sn,next;
    taucs_datatype v;
    int* len;

    n = L->n;

    len = (int*) malloc(n*sizeof(int));
    if (!len) return NULL;

    nnz = 0;

    for (sn=0; sn<L->n_sn; sn++) {
	for (jp=0; jp<(L->sn_size)[sn]; jp++)
	{
	    j      = (L->sn_struct)[sn][jp];
	    len[j] = 0;

	    for (ip=jp; ip<(L->sn_size)[sn]; ip++)
	    {
		i = (L->sn_struct)[sn][ ip ];
		if (i >= j)
		{
		    len[j] ++;
		    nnz ++;
		}
	    }
	    for (ip=(L->sn_size)[sn]; ip<(L->sn_up_size)[sn]; ip++)
	    {
		i = (L->sn_struct)[sn][ ip ];
		if (i >= j)
		{
		    len[j] ++;
		    nnz ++;
		}
	    }
	}
    }


    C = taucs_dccs_create(n,n,nnz);
    if (!C)
    {
	free(len);
	return NULL;
    }

#ifdef TAUCS_CORE_SINGLE
    C->flags = TAUCS_SINGLE;
#endif
#ifdef TAUCS_CORE_DOUBLE
    C->flags = TAUCS_DOUBLE;
#endif
#ifdef TAUCS_CORE_SCOMPLEX
    C->flags = TAUCS_SCOMPLEX;
#endif
#ifdef TAUCS_CORE_DCOMPLEX
    C->flags = TAUCS_DCOMPLEX;
#endif

    C->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;	  /* this was a bug in version 2.0 of taucs */

    (C->colptr)[0] = 0;
    for (j=1; j<=n; j++) (C->colptr)[j] = (C->colptr)[j-1] + len[j-1];

    free(len);

    for (sn=0; sn<L->n_sn; sn++)
    {
	for (jp=0; jp<(L->sn_size)[sn]; jp++)
	{
	    j    = (L->sn_struct)[sn][jp];
	    next = (C->colptr)[j];

	    for (ip=jp; ip<(L->sn_size)[sn]; ip++)
	    {
		i = (L->sn_struct)[sn][ ip ];
		v = (L->sn_blocks)[sn][ jp*(L->sn_blocks_ld)[sn] + ip ];

		if (i >= j)
		{
		    (C->rowind)[next]       = i;
		    (C->taucs_values)[next] = v;
		    next++;
		}
	    }
	    for (ip=(L->sn_size)[sn]; ip<(L->sn_up_size)[sn]; ip++)
	    {
		i = (L->sn_struct)[sn][ ip ];
		v = (L->up_blocks)[sn][ jp*(L->up_blocks_ld)[sn] + (ip-(L->sn_size)[sn]) ];

		if (i >= j)
		{
		    (C->rowind)[next]       = i;
		    (C->taucs_values)[next] = v;
		    next++;
		}
	    }
	}
    }
    return C;
}
#endif 

/* 
   copy a supernodal_factor_matrix
*/
supernodal_factor_matrix *my_taucs_supernodal_factor_matrix_copy(supernodal_factor_matrix *L)
{
#define CREATE(name,len,type) if (1)\
                              {\
                                  if (L->name)\
                                  {\
                                      LL->name = Calloc((len),type); \
                                      MEMCHK_RETVAL(LL->name,NULL); \
                                      memcpy(LL->name,L->name,(len)*sizeof(type));\
                                  }\
                                  else \
                                  {\
                                      LL->name = (type *)NULL;\
                                  }\
                              }
                             
    supernodal_factor_matrix *LL;
    int n, np, n_sn, i;
    
    if (!L) return NULL;
    LL        = Calloc(1, supernodal_factor_matrix);
    LL->flags = L->flags;
    LL->uplo  = L->uplo;
    LL->n     = L->n;
    LL->n_sn  = L->n_sn;
    
    n    = LL->n;
    np   = n+1;
    n_sn = LL->n_sn;
    
    CREATE(sn_size, np, int);
    CREATE(sn_up_size, np, int);
    CREATE(first_child, np, int);
    CREATE(next_child, np, int);
    CREATE(parent, np, int);

    CREATE(sn_blocks_ld, n_sn, int);
    CREATE(up_blocks_ld, n_sn, int);

    LL->sn_struct = Calloc(n, int*);  MEMCHK_RETVAL(LL->sn_struct, NULL);
    for(i=0;i<LL->n;i++)
    {
	CREATE(sn_struct[i], LL->sn_up_size[i], int);
    }
    
    LL->sn_blocks = Calloc(n_sn, double *); MEMCHK_RETVAL(LL->sn_blocks, NULL);
    for(i=0;i<LL->n_sn;i++)
    {
	CREATE(sn_blocks[i], SQR(LL->sn_size[i]), double);
    }

    LL->up_blocks = Calloc(n_sn, double *); MEMCHK_RETVAL(LL->up_blocks, NULL);
    for(i=0;i<LL->n_sn;i++)
    {
	CREATE(up_blocks[i], (LL->sn_up_size[i]- LL->sn_size[i])*(LL->sn_size)[i], double);
    }

    return LL;

#undef CREATE
}
/* 
   make a copy of a ccs-matrix
*/
taucs_ccs_matrix *my_taucs_dccs_copy(taucs_ccs_matrix *L, int flags)
{
    /* 
       copy a square matrix
    */
    
    taucs_ccs_matrix *LL;
    int n, nnz;
    
    GMRFLib_ASSERT_RETVAL(L->n == L->m, GMRFLib_EPARAMETER, NULL); /* must be square */

    if (!L) return NULL;

    n   = L->n;
    nnz = L->colptr[L->n];
    LL  = taucs_ccs_create(n, n, nnz, flags);
    
    memcpy(LL->colptr,   L->colptr,   (n+1)*sizeof(int));
    memcpy(LL->rowind,   L->rowind,   nnz  *sizeof(int));
    memcpy(LL->values.d, L->values.d, nnz  *sizeof(double));

    return LL;
}

int GMRFLib_compute_reordering_TAUCS_orig(int **remap, GMRFLib_graph_tp *graph)
{
    /* 
       this is the original version without treating global nodes spesifically.
    */
    
    int i, j, k, ic, ne, n, nnz, *perm = NULL, *iperm = NULL;
    taucs_ccs_matrix *Q;
    
    if (!graph || graph->n==0) return GMRFLib_SUCCESS;

    n = graph->n;
    for(i=0, nnz=n; i<n; i++) nnz += graph->nnbs[i];

    Q            = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
    Q->flags     = (TAUCS_PATTERN | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
    Q->colptr[0] = 0;
    
    for(i=0, ic=0; i<n; i++)
    {
	Q->rowind[ic++] = i;
	for(k=0, ne=1;k<graph->nnbs[i];k++)
	{
	    j = graph->nbs[i][k];
	    if (j > i) break;
	    Q->rowind[ic++] = j;
	    ne++;
	}
	Q->colptr[i+1] = Q->colptr[i]+ne;
    }
    
    if (1)
	taucs_ccs_order(Q, &perm, &iperm, "metis");	  /* use the metis library */
    else
    {
	FIXME("use identity reordering");
	taucs_ccs_order(Q, &perm, &iperm, "identity");	  /* use the identity mapping */
    }

    *remap = iperm;				  /* yes, this is correct */
    FREE(perm);
    taucs_ccs_free(Q);

    if (!*remap) GMRFLib_ERROR(GMRFLib_EREORDER);
    return GMRFLib_SUCCESS;
}
int GMRFLib_compute_reordering_TAUCS(int **remap, GMRFLib_graph_tp *graph)
{
    /* 
       new improved version which treats global nodes spesifically.
    */
    int i, j, k, ic, ne, n, ns, nnz, *perm = NULL, *iperm = NULL, limit, free_subgraph, *iperm_new = NULL, simple;
    char *fixed = NULL;
    taucs_ccs_matrix *Q;
    GMRFLib_graph_tp *subgraph;
    
    if (!graph || graph->n==0) return GMRFLib_SUCCESS;

    
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
       check if we have 'global' nodes
    */
    for(i=0, ne=0;i<graph->n;i++) ne = MAX(ne, graph->nnbs[i]);
    limit = MAX(50, graph->n/5.0);		  /* this is the limit for a 'global' node */
    if (ne > limit)
    {
	/* 
	   yes we have global nodes, make a new graph with these removed.
	*/
	fixed = Calloc(graph->n, char); MEMCHK(fixed);
	for(i=0;i<graph->n;i++) fixed[i] = (graph->nnbs[i] > limit);

	GMRFLib_compute_subgraph(&subgraph, graph, fixed);
	free_subgraph = 1;

	if (subgraph->n == 0)
	{
	    /* 
	       this is a weird event: abort treating the global nodes spesifically.
	    */
	    GMRFLib_free_graph(subgraph);
	    subgraph      = graph;
	    free_subgraph = 0;
	}
    }
    else
    {
	subgraph      = graph;
	free_subgraph = 0;
    }

    /* 
       continue with subgraph, which is the original graph minus global nodes
    */
    n = subgraph->n;
    for(i=0, nnz=n; i<n; i++) nnz += subgraph->nnbs[i];

    Q            = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
    Q->flags     = (TAUCS_PATTERN | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
    Q->colptr[0] = 0;
    
    for(i=0, ic=0; i<n; i++)
    {
	Q->rowind[ic++] = i;
	for(k=0, ne=1;k<subgraph->nnbs[i];k++)
	{
	    j = subgraph->nbs[i][k];
	    if (j > i) break;
	    Q->rowind[ic++] = j;
	    ne++;
	}
	Q->colptr[i+1] = Q->colptr[i]+ne;
    }
    
    if (1)
	taucs_ccs_order(Q, &perm, &iperm, "metis");	  /* use the metis library */
    else
    {
	FIXME("use identity reordering");
	taucs_ccs_order(Q, &perm, &iperm, "identity");	  /* use the identity mapping */
    }
    FREE(perm);
    taucs_ccs_free(Q);

    if (!free_subgraph)
    {
	/* 
	   no global nodes, then `iperm' is the reordering
	*/
	*remap = iperm;				  /* yes, this is correct */
    }
    else
    {
	/* 
	   global nodes, correct the reordering computed and add the global nodes at the end
	*/
	ns        = subgraph->n;
	n         = graph->n;
	iperm_new = Calloc(n, int); MEMCHK(iperm_new);

	for(i=0; i<ns; i++) iperm_new[i] = subgraph->mothergraph_idx[ iperm[i] ]; 
	for(i=0,j=ns; i<n; i++) if (fixed[i]) iperm_new[j++] = i;
	GMRFLib_ASSERT(j == n, GMRFLib_ESNH);	  /* just a check... */
	*remap = iperm_new;			  /* this is the reordering */

	FREE(iperm);				  /* this is no longer needed */
	GMRFLib_free_graph(subgraph);
    }

    if (!*remap) GMRFLib_ERROR(GMRFLib_EREORDER);
    return GMRFLib_SUCCESS;
}
int GMRFLib_free_reordering_TAUCS(GMRFLib_sm_fact_tp *sm_fact)
{
    if (sm_fact) FREE(sm_fact->remap);
    return GMRFLib_SUCCESS;
}
int GMRFLib_build_sparse_matrix_TAUCS(taucs_ccs_matrix **L, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg,
				      GMRFLib_graph_tp *graph, int *remap)
{
    int i, j, k, ic, ne, n, nnz, *perm = NULL, *iperm = NULL;
    taucs_ccs_matrix *Q;
    
    if (!graph || graph->n==0) 
    {
	*L = NULL;
	return GMRFLib_SUCCESS;
    }

    n = graph->n;
    
    for(i=0, nnz=n; i<n; i++) nnz += graph->nnbs[i];

    Q            = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE); GMRFLib_ASSERT(Q, GMRFLib_EMEMORY);
    Q->flags     = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
    Q->colptr[0] = 0;
    
    for(i=0, ic=0; i<n; i++)
    {
	Q->rowind[ic]   = i;
	Q->values.d[ic] = (*Qfunc)(i, i, Qfunc_arg);
	ic++;
	ne = 1;

	for(k=0;k<graph->nnbs[i];k++)
	{
	    j = graph->nbs[i][k];
	    if (j > i) break;

	    Q->rowind[ic]   = j;
	    Q->values.d[ic] = (*Qfunc)(i, j, Qfunc_arg);

	    ic++; ne++;
	}
	Q->colptr[i+1] = Q->colptr[i]+ne;
    }
    
    iperm = remap;			  /* yes, this is correct */
    perm  = Calloc(n, int); MEMCHK(perm);
    for(i=0;i<n;i++) perm[iperm[i]] = i;

    *L = taucs_ccs_permute_symmetrically(Q, perm, iperm); /* permute the matrix */

    taucs_ccs_free(Q);
    FREE(perm);

    return GMRFLib_SUCCESS;
}
int GMRFLib_factorise_sparse_matrix_TAUCS(taucs_ccs_matrix **L, supernodal_factor_matrix **symb_fact, GMRFLib_fact_info_tp *finfo)
{
    int flags, k, retval;
        
    if (!L) return GMRFLib_SUCCESS;

    /* 
       compute some info about the factorization
    */
    k = (*L)->colptr[(*L)->n] - (*L)->n;
    finfo->n      = (*L)->n;
    finfo->nnzero = 2*k + (*L)->n;

    flags = (*L)->flags;
    if (!*symb_fact) *symb_fact  = taucs_ccs_factor_llt_symbolic(*L);
    retval      = taucs_ccs_factor_llt_numeric(*L, *symb_fact); if (retval) GMRFLib_ERROR(GMRFLib_EPOSDEF);
    taucs_ccs_free(*L);
    *L          = taucs_supernodal_factor_to_ccs(*symb_fact);
    (*L)->flags = flags & ~TAUCS_SYMMETRIC;	  /* fixes a bug in ver 2.0 av TAUCS */
    taucs_supernodal_factor_free_numeric(*symb_fact); /* remove the numerics, preserve the symbolic */

    /* 
       some last info
    */
    k = (*L)->colptr[(*L)->n] - (*L)->n;
    finfo->nfillin = k - (finfo->nnzero - finfo->n)/2;

    return GMRFLib_SUCCESS;
}
int GMRFLib_factorise_sparse_matrix_TAUCS_OLD(taucs_ccs_matrix **L, GMRFLib_fact_info_tp *finfo)
{
    taucs_ccs_matrix *L_fact;
    int flags;
    int solver_option = 0;			  /* chose the solver here */
    int k;
    
    if (!L) return GMRFLib_SUCCESS;


    /* 
       compute some info about the factorization
    */
    k = (*L)->colptr[(*L)->n] - (*L)->n;
    finfo->n      = (*L)->n;
    finfo->nnzero = 2*k + (*L)->n;


    flags = (*L)->flags;

    switch(solver_option)
    {
    case 0:					  /* fastest */
	L_fact      = taucs_ccs_factor_llt_mf(*L); taucs_ccs_free(*L);
	if (!L_fact) GMRFLib_ERROR(GMRFLib_EPOSDEF);
	*L          = L_fact;
	L_fact      = taucs_supernodal_factor_to_ccs(*L); taucs_supernodal_factor_free(*L);
	*L          = L_fact;
	(*L)->flags = flags & ~TAUCS_SYMMETRIC;
	break;
	
    case 1:					  /* less fast, but less memory */
	L_fact      = taucs_ccs_factor_llt_ll(*L); taucs_ccs_free(*L);
	if (!L_fact) GMRFLib_ERROR(GMRFLib_EPOSDEF);
	*L          = L_fact;
	L_fact      = taucs_supernodal_factor_to_ccs(*L); taucs_supernodal_factor_free(*L);
	*L          = L_fact;
	(*L)->flags = flags & ~TAUCS_SYMMETRIC;
	break;

    case 2:					  /* slower */
	L_fact=taucs_ccs_factor_llt(*L, 0.0, 0);
	if (!L_fact) GMRFLib_ERROR(GMRFLib_EPOSDEF);
	taucs_ccs_free(*L);
	*L = L_fact;
	break;

    default:
	abort();
    }

    
    /* 
       some last info
    */
    k = (*L)->colptr[(*L)->n] - (*L)->n;
    finfo->nfillin = k - (finfo->nnzero - finfo->n)/2;

    return GMRFLib_SUCCESS;
}
int GMRFLib_free_fact_sparse_matrix_TAUCS(taucs_ccs_matrix *L, supernodal_factor_matrix *symb_fact)
{
    if (L) taucs_ccs_free(L);
    if (symb_fact) taucs_supernodal_factor_free(symb_fact);
    return GMRFLib_SUCCESS;
}
int GMRFLib_free_fact_sparse_matrix_TAUCS_OLD(taucs_ccs_matrix *L)
{
    if (L) taucs_ccs_free(L);
    return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap)
{
    double *b;
    
    GMRFLib_ASSERT(rhs, GMRFLib_EINVARG);
    GMRFLib_ASSERT(L, GMRFLib_EINVARG);

    GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    b = Calloc(graph->n, double); MEMCHK(b);
    memcpy(b, rhs, graph->n*sizeof(double));

    my_taucs_dccs_solve_lt(L, rhs, b);

    GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
    FREE(b);
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_solve_llt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap)
{
    double *b;
    
    GMRFLib_ASSERT(rhs, GMRFLib_EINVARG);
    GMRFLib_ASSERT(L, GMRFLib_EINVARG);

    GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    b = Calloc(graph->n, double); MEMCHK(b);
    memcpy(b, rhs, graph->n*sizeof(double));

    taucs_ccs_solve_llt(L, rhs, b);

    GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
    FREE(b);
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap, 
						 int findx, int toindx, int remapped)
{
    /* 
       rhs in real world, L in mapped world.  solve L^Tx=b backward only from rhs[findx] up to
       rhs[toindx].  note that findx and toindx is in mapped world.  if remapped, do not
       remap/remap-back the rhs before solving.
       
     */

    static double *b = NULL;
    static int len_b = 0;

    if (graph->n > len_b)
    {
	len_b = graph->n;
	b     = (double *)realloc(b, len_b*sizeof(double)); MEMCHK(b);
    }

    if (!remapped) GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
    memcpy(&b[toindx], &rhs[toindx], (graph->n-toindx)*sizeof(double));	/* this can be improved */

    my_taucs_dccs_solve_lt_special(L, rhs, b, findx, toindx); /* solve it */
    if (!remapped) GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

    return GMRFLib_SUCCESS;
}
int GMRFLib_comp_cond_meansd_TAUCS(double *cmean, double *csd, int indx, double *x, int remapped,
				   taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap)
{
    /* 
       compute the conditonal mean and stdev for x[indx]|x[indx+1]...x[n-1] for the current value of
       x. if `remapped', then x is assumed to be remapped for possible (huge) speedup when used
       repeately.  note: indx is still the user world!!!
     */
    int ii = remap[indx];

    if (remapped)
	my_taucs_cmsd(cmean, csd, ii, L, x);
    else
    {
	GMRFLib_convert_to_mapped(x, NULL, graph, remap);
	my_taucs_cmsd(cmean, csd, ii, L, x);
	GMRFLib_convert_from_mapped(x, NULL, graph, remap);
    }

    return GMRFLib_SUCCESS;
}
int GMRFLib_log_determinant_TAUCS(double *logdet, taucs_ccs_matrix *L)
{
    int i;

    *logdet = 0.0;
    for (i=0; i< L->n;i++) *logdet += log(L->values.d[(L->colptr)[i]]); 
    *logdet *= 2;

    return GMRFLib_SUCCESS;
}
int my_intcmp(const void *a, const void *b)
{
    int i = *(int *)a;
    int j = *(int *)b;
    if (i > j) return 1;
    if (i < j) return -1;
    return 0;
}
int GMRFLib_compute_Qinv_TAUCS(GMRFLib_problem_tp *problem, int storage)
{
    /* 
       strategy:

       try first to compute Qinv, if fail, then go back and add zero-terms to L until ok, and
       compute Qinv again.
    */

    int i, n;
    taucs_ccs_matrix *L  = NULL, *LL = NULL;	  /* to hold L matrices if new ones are built */
    map_ii **mis_elm = NULL;

    if (!problem) return GMRFLib_SUCCESS;
    n = problem->sub_sm_fact.L->n;

    /* 
       no-check, check-once or failsafe?
    */
    if ((storage & GMRFLib_QINV_NO_CHECK))
	GMRFLib_compute_Qinv_TAUCS_compute(problem, storage, NULL);
    else 
    {
	/* 
	   do some checking
	*/
	L = my_taucs_dccs_copy(problem->sub_sm_fact.L, TAUCS_DOUBLE|TAUCS_LOWER|TAUCS_TRIANGULAR);
	while ((mis_elm = GMRFLib_compute_Qinv_TAUCS_check(L)))
	{
	    LL = L; 
	    L  = GMRFLib_compute_Qinv_TAUCS_add_elements(LL, mis_elm);
	    taucs_ccs_free(LL);

	    for(i=0;i<n;i++) 
	    {
		map_ii_free(mis_elm[i]); FREE(mis_elm[i]);
	    }
	    FREE(mis_elm);

	    if ((storage & GMRFLib_QINV_CHECK_ONCE)) break; /* check once only */
	}
	GMRFLib_compute_Qinv_TAUCS_compute(problem, storage, L);
	taucs_ccs_free(L);
    }

    return GMRFLib_SUCCESS;
}
map_ii **GMRFLib_compute_Qinv_TAUCS_check(taucs_ccs_matrix *L)
{
    int i, j, k, jp, ii, kk, *nnbs, **nbs, *nnbsQ, nmissing=0, n, debug=0;
    map_ii **Qinv_L;				
    map_ii **missing;
    
    if (!L) return NULL;
    n = L->n;
    
    /* 
       construct a row-list of L_ij's including the diagonal
    */
    nnbs  = Calloc(n, int);   MEMCHK_RETVAL(nnbs,  NULL);
    nbs   = Calloc(n, int *); MEMCHK_RETVAL(nbs,   NULL);
    nnbsQ = Calloc(n, int);   MEMCHK_RETVAL(nnbsQ, NULL); /* number of elm in the Qinv_L[j] hash-table */
    
    for(i=0;i<n;i++)
	for(jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++) nnbs[ (L->rowind)[jp] ]++;

    for(i=0;i<n;i++)
    { 
	nbs[i]  = Calloc(nnbs[i], int); MEMCHK_RETVAL(nbs[i], NULL);
	nnbs[i] = 0;
    }
    
    for(j=0;j<n;j++)
    {						  
	for(jp = (L->colptr)[j]; jp < (L->colptr)[j+1]; jp++) /* including the diagonal */
	{
	    i                  = (L->rowind)[jp];
	    nbs[i][  nnbs[i] ] = j;		 
	    nnbs[i]++;
	    nnbsQ[MIN(i,j)]++;			  /* for the Qinv_L[] hash-table */
	}
    }
    for(i=0;i<n;i++) qsort(nbs[i], nnbs[i], sizeof(int), my_intcmp); /* is this needed ???? */
    
    Qinv_L = Calloc(n, map_ii *); MEMCHK_RETVAL(Qinv_L, NULL);
    for(i=0;i<n;i++) 
    {
	Qinv_L[i] = Calloc(1, map_ii); MEMCHK_RETVAL(Qinv_L[i], NULL);
	map_ii_init_hint(Qinv_L[i], nnbsQ[i]); 
	Qinv_L[i]->alwaysdefault = 0;		  /* return NULL if not there */
    }
    

    missing = Calloc(n, map_ii *); MEMCHK_RETVAL(missing, NULL);
    for(i=0;i<n;i++)
    {
	missing[i] = Calloc(1, map_ii); MEMCHK_RETVAL(missing[i], NULL);
	map_ii_init(missing[i]);
	missing[0]->alwaysdefault = 0;
    }

    if (1)
    {
	/* 
	   different versions
	*/
	if (1)
	{
	    /* 
	       fast version 1

	       not much to gain here.
	    */
	    char *Zj;
	    int *Zj_set, nset;
	    
	    Zj     = Calloc(n, char); MEMCHK_RETVAL(Zj, NULL);
	    Zj_set = Calloc(n, int);  MEMCHK_RETVAL(Zj_set, NULL);

	    for(j=n-1;j>=0;j--)
	    {
		nset = 0;
		for(k = -1; (k = map_ii_next(Qinv_L[j], k)) != -1;)
		{
		    kk             = Qinv_L[j]->contents[k].key;
		    Zj[kk]         = 1;
		    Zj_set[nset++] = kk;
		}

		for(ii = nnbs[j]-1; ii>=0; ii--)
		{
		    i = nbs[j][ii];
		    for(kk = L->colptr[i]+1; kk < L->colptr[i+1];kk++)
		    {
			k = L->rowind[kk];
			if (Zj[k] == 0 && L->values.d[kk] != 0.0)
			{
			    if (debug) printf("\tmissing %d %d\n", MIN(k, j), MAX(k, j));

			    map_ii_set(missing[MIN(k, j)], MAX(k, j), 1);
			    if (k < j) map_ii_set(Qinv_L[k], j, 1); /* also mark those who are missing */
			    Zj[k]          = 1;
			    Zj_set[nset++] = k;
			    nmissing++;
			}
		    }
		    Zj[i]          = 1;
		    Zj_set[nset++] = i;
		    map_ii_set(Qinv_L[i], j, 1); 
		}
		/*
		  if (j > 0) for(kk=0;kk<nset;kk++) Zj[Zj_set[kk]] = 0;
		*/
		
		if (j>0)				  /* not needed for j=0 */
		{
		    if (nset > GMRFLib_NSET_LIMIT(nset, sizeof(char), n))	  
			memset(Zj, 0, n*sizeof(char));
		    else
			for(kk=0;kk<nset;kk++) Zj[Zj_set[kk]] = 0; /* set those to zero */
		}
	    }
	    FREE(Zj); FREE(Zj_set);
	}
	else
	{
	    /* 
	       fast version 2

	       run almost as fast as the above, probably since the size of Zj is char and hence
	       small. but is to be preferred as it does not matter how large 'nset' is?

	    */
	    char *Zj;
	    
	    Zj = Calloc(n, char); MEMCHK_RETVAL(Zj, NULL);

	    for(j=n-1;j>=0;j--)
	    {
		for(k = -1; (k = map_ii_next(Qinv_L[j], k)) != -1;)
		    Zj[Qinv_L[j]->contents[k].key] = 1;

		for(ii = nnbs[j]-1; ii>=0; ii--)
		{
		    i = nbs[j][ii];
		    for(kk = L->colptr[i]+1; kk < L->colptr[i+1];kk++)
		    {
			k = L->rowind[kk];
			if (Zj[k] == 0 && L->values.d[kk] != 0.0)
			{
			    if (debug) printf("\tmissing %d %d\n", MIN(k, j), MAX(k, j));

			    map_ii_set(missing[MIN(k, j)], MAX(k, j), 1);
			    if (k < j) map_ii_set(Qinv_L[k], j, 1); /* also mark those who are missing */
			    Zj[k] = 1;
			    nmissing++;
			}
		    }
		    Zj[i] = 1;
		    map_ii_set(Qinv_L[i], j, 1); 
		}
		if (j>0) memset(Zj, 0, n*sizeof(char));
	    }
	    FREE(Zj);
	}
    }
    else
    {
	/* 
	   slow version, but don't delete!
	*/
	for(j=n-1;j>=0;j--)
	{
	    for(ii = nnbs[j]-1; ii>=0; ii--)
	    {
		i = nbs[j][ii];
		for(kk = L->colptr[i]+1; kk < L->colptr[i+1];kk++)
		{
		    k = L->rowind[kk];
		    if (!map_ii_ptr(Qinv_L[MIN(k, j)], MAX(k, j)))
		    {
			if (debug) printf("\tmissing %d %d\n", MIN(k, j), MAX(k, j));

			map_ii_set(missing[MIN(k, j)], MAX(k, j), 0);
			map_ii_set(Qinv_L[ MIN(k, j)], MAX(k, j), 0); /* also mark those who are missing */
			nmissing++;
		    }
		}
		map_ii_set(Qinv_L[MIN(i, j)], MAX(i, j), 0); 
	    }
	}
    }
    if (debug) P(nmissing);
    
    for(i=0;i<n;i++) 
    {
	FREE(nbs[i]); map_ii_free(Qinv_L[i]); FREE(Qinv_L[i]);
    }
    FREE(nbs); FREE(nnbs); FREE(Qinv_L);

    if (nmissing == 0)
    {
	/* 
	   free the missing-hash, and set it to NULL so that this function returns NULL
	*/
	for(i=0;i<n;i++) 
	{
	    map_ii_free(missing[i]); FREE(missing[i]);
	}
	FREE(missing);	
    }

    return missing;
}
taucs_ccs_matrix *GMRFLib_compute_Qinv_TAUCS_add_elements(taucs_ccs_matrix *L, map_ii **missing)
{
    int i, j, k, jp, jpp, n, nnz, nnz_new, nmissing, collen;
    taucs_ccs_matrix *LL;
    
    n   = L->n;
    /* 
       first count the number of missing terms
    */
    for(i=0, nmissing=0;i<n;i++)
	for(k = -1; (k = map_ii_next(missing[i], k)) != -1; )
	    nmissing++;

    nnz     = L->colptr[n];
    nnz_new = nnz + nmissing;

    LL = taucs_dccs_create(n, n, nnz_new);
    LL->flags = L->flags;			  /* maintain properties */

    LL->colptr[0] = 0;				  /* start at zero */
    for(j=0, jp=0;j<n;j++)
    {
	/* 
	   do each column j
	*/
	collen = L->colptr[j+1]-L->colptr[j];
	jpp    = L->colptr[j];
	
	memcpy(&(LL->rowind[jp]),   &(L->rowind[jpp]),   collen*sizeof(int));
	memcpy(&(LL->values.d[jp]), &(L->values.d[jpp]), collen*sizeof(double));

	jp += collen;
	
	/* 
	   now add new elms to this column
	*/
	for(k = -1; (k = map_ii_next(missing[j], k)) != -1; )
	{
	    LL->rowind[jp]   = missing[j]->contents[k].key;
	    LL->values.d[jp] = 0.0;
	    jp++;
	}
	LL->colptr[j+1] = jp;
    }
    
    if (0)
    {
        for(j=0;j<n;j++)
            for(jp=L->colptr[j];jp<L->colptr[j+1];jp++)
                printf("L[ %1d %1d ] = %.12f\n", L->rowind[jp], j, L->values.d[jp]);

        for(j=0;j<n;j++)
            for(jp=LL->colptr[j];jp<LL->colptr[j+1];jp++)
                printf("LL[ %1d %1d ] = %.12f\n", LL->rowind[jp], j, LL->values.d[jp]);
    }

    return LL;
}
int GMRFLib_compute_Qinv_TAUCS_compute(GMRFLib_problem_tp *problem, int storage, taucs_ccs_matrix *Lmatrix)
{
    /* 
       compute the elements in Qinv from the non-zero pattern of L (no checking). store them
       according to `storage': GMRFLib_QINV_ALL GMRFLib_QINV_NEIGB GMRFLib_QINV_DIAG
    */
    double *ptr, value, diag, *Zj;
    int i, j, k, jp, ii, kk, jj, iii, jjj, n,
	*nnbs=NULL, **nbs=NULL, *nnbsQ=NULL, *remove=NULL, nremove, *inv_remap=NULL,
	*Zj_set, nset;
    taucs_ccs_matrix *L;
    map_ii *mapping;
    map_id **Qinv_L, *q;

    L = (Lmatrix ? Lmatrix : problem->sub_sm_fact.L); /* chose matrix to use */
    n = L->n;
    
    /* 
       construct a row-list of L_ij's including the diagonal
    */
    nnbs  = Calloc(n, int);   MEMCHK(nnbs);
    nbs   = Calloc(n, int *); MEMCHK(nbs);
    nnbsQ = Calloc(n, int);   MEMCHK(nnbsQ);	  /* number of elm in the Qinv_L[j] hash-table */
    for(i=0;i<n;i++)
	for(jp = (L->colptr)[i]; jp < (L->colptr)[i+1]; jp++)
	    nnbs[ (L->rowind)[jp] ]++;

    for(i=0;i<n;i++)
    { 
	nbs[i]  = Calloc(nnbs[i], int); MEMCHK(nbs[i]);
	nnbs[i] = 0;
    }
    
    for(j=0;j<n;j++)
    {						  
	for(jp = (L->colptr)[j]; jp < (L->colptr)[j+1]; jp++) /* including the diagonal */
	{
	    i                  = (L->rowind)[jp];
	    nbs[i][  nnbs[i] ] = j;		 
	    nnbs[i]++;
	    nnbsQ[MIN(i,j)]++;			  /* for the Qinv_L[] hash-table */
	}
    }
    for(i=0;i<n;i++) qsort(nbs[i], nnbs[i], sizeof(int), my_intcmp);
    
    /* 
       setup the hash-table for storing Qinv_L
    */
    Qinv_L = Calloc(n, map_id *); MEMCHK(Qinv_L);
    for(i=0;i<n;i++) 
    {
	Qinv_L[i] = Calloc(1, map_id); MEMCHK(Qinv_L[i]);
	map_id_init_hint(Qinv_L[i], nnbsQ[i]); 
	Qinv_L[i]->alwaysdefault = 0;
    }

    
    Zj     = Calloc(n, double); MEMCHK(Zj);
    Zj_set = Calloc(n, int);    MEMCHK(Zj_set);

    for(j=n-1;j>=0;j--)
    {
	/* 
	   store those indices that are used and set only those to zero
	*/
	if (0) memset(Zj, 0, n*sizeof(double)); /* old and slow(er) solution */

	nset = 0;
	q    = Qinv_L[j];		  /* just to store the ptr */

	for(k = -1; (k = map_id_next(q, k)) != -1;)
	{
	    jj = q->contents[k].key;

	    Zj_set[nset++] = jj;
	    Zj[jj]         = q->contents[k].value;
	}
		
	for(ii = nnbs[j]-1; ii>=0; ii--)
	{
	    i     = nbs[j][ii];
	    diag  = L->values.d[L->colptr[i]];
	    value = (i == j ? 1./diag : 0.0);

	    for(kk = L->colptr[i]+1; kk < L->colptr[i+1];kk++) value -= L->values.d[kk] * Zj[L->rowind[kk]];

	    value         /= diag;
	    Zj[i]          = value;
	    Zj_set[nset++] = i;		

	    map_id_set(Qinv_L[i], j, value);
	}

	if (j>0)				  /* not needed for j=0 */
	{
	    if (nset > GMRFLib_NSET_LIMIT(nset, sizeof(double), n))	  
		memset(Zj, 0, n*sizeof(double));  /* faster if nset is large */
	    else
		for(kk=0;kk<nset;kk++) Zj[Zj_set[kk]] = 0.0; /* set those to zero */
	}
    }

    if (0)
    {
	/* 
	   keep this OLD version in the source
	*/
	for(j=n-1;j>=0;j--)
	{
	    for(ii = nnbs[j]-1; ii>=0; ii--)
	    {
		i     = nbs[j][ii];
		diag  = L->values.d[L->colptr[i]];
		value = (i == j ? 1./diag : 0.0);

		for(kk = L->colptr[i]+1; kk < L->colptr[i+1];kk++)
		{
		    k   = L->rowind[kk];
		    if ((ptr = map_id_ptr(Qinv_L[MIN(k, j)], MAX(k, j))))
			value -= L->values.d[kk] * *ptr;
		}

		value /= diag;
		map_id_set(Qinv_L[MIN(i, j)], MAX(i, j), value);
	    }
	}
    }

    /* 
       compute the mapping
    */
    inv_remap = Calloc(n, int); MEMCHK(inv_remap);
    for(k=0;k<n;k++) inv_remap[problem->sub_sm_fact.remap[k]] = k;


    /* 
       possible remove entries: options are GMRFLib_QINV_FULL GMRFLib_QINV_NEIGB GMRFLib_QINV_DIAG
    */
    if(storage & (GMRFLib_QINV_DIAG|GMRFLib_QINV_NEIGB))
    {
	remove = Calloc(n, int); MEMCHK(remove);
	for(i=0;i<n;i++)
	{
	    iii = inv_remap[i];
	    if (storage & GMRFLib_QINV_DIAG)
	    {
		for (k = -1, nremove=0; (k = map_id_next(Qinv_L[i], k)) != -1 ; )
		    if((j = Qinv_L[i]->contents[k].key) != i)
			remove[nremove++] = j;
	    }
	    else
	    {
		for (k = -1, nremove=0; (k = map_id_next(Qinv_L[i], k)) != -1; )
		{
		    j = Qinv_L[i]->contents[k].key;

		    if (j != i)
		    {
			jjj = inv_remap[j];
			if (!GMRFLib_is_neighb(iii, jjj, problem->sub_graph)) remove[nremove++] = j;
		    }
		}
	    }
	    for(k=0;k<nremove;k++) map_id_remove(Qinv_L[i], remove[k]);
	    map_id_adjustcapacity(Qinv_L[i]);
	}
    }
    
    /* 
       correct for constraints, if any. need `iremap' as the matrix terms, constr_m and qi_at_m, is
       in the sub_graph coordinates without reordering!

       not that this is correct for both hard and soft constraints, as the constr_m matrix contains
       the needed noise-term.
    */
    if (problem->sub_constr && problem->sub_constr->nc > 0)
    {
	for(i=0;i<n;i++)
	{
	    iii = inv_remap[i];
	    for(k=-1; (k = map_id_next(Qinv_L[i], k)) != -1 ; )
	    {
		j   = Qinv_L[i]->contents[k].key;
		jjj = inv_remap[j];
		
		map_id_get(Qinv_L[i], j, &value);

		for(kk=0;kk<problem->sub_constr->nc;kk++)
		    value -= problem->constr_m[iii+kk*n] * problem->qi_at_m[jjj+kk*n];

		map_id_set(Qinv_L[i], j, value);
	    }
	}
    }

    /* 
       done. store Qinv
    */
    problem->sub_inverse       = Calloc(1, GMRFLib_Qinv_tp); MEMCHK(problem->sub_inverse);
    problem->sub_inverse->Qinv = Qinv_L;
    
    /* 
       compute the mapping for lookup using GMRFLib_Qinv_get(). here, the user lookup using a global
       index, which is then transformed to the reordered sub_graph.
    */
    problem->sub_inverse->mapping = mapping = Calloc(1, map_ii); MEMCHK(mapping);
    map_ii_init_hint(mapping, n);
    mapping->alwaysdefault = 0;			  /* return NULL if not there */
    for(i=0;i<n;i++)
	map_ii_set(mapping, problem->sub_graph->mothergraph_idx[i], problem->sub_sm_fact.remap[i]);
    
    /* 
       cleanup
    */
    for(i=0;i<n;i++) FREE(nbs[i]); 
    FREE(nbs);FREE(nnbs);FREE(nnbsQ);FREE(inv_remap);FREE(remove);FREE(Zj);FREE(Zj_set);

    return GMRFLib_SUCCESS;
}
int my_taucs_dccs_solve_lt(void* vL, double* x, double* b)
{
    taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

    int i, j, jp;
    double  Aij, Aii;

    for (i=L->n-1; i>=0; i--)
    {
	for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++)
	{
	    j     = (L->rowind)[jp];
	    Aij   = (L->values.d)[jp] ;
	    b[i] -= x[j]*Aij;
	}

	jp   = (L->colptr)[i];
	j    = (L->rowind)[jp];
	Aii  = (L->values.d)[jp]; 
	x[i] = b[i] / Aii; 
    }

    return 0;
}
int my_taucs_dccs_solve_lt_special(void* vL, double* x, double* b, int from_idx, int to_idx)
{
    taucs_ccs_matrix* L = (taucs_ccs_matrix*) vL;

    int i, j, jp;
    double  Aij, Aii;

    for (i=from_idx; i>=to_idx; i--)
    {
	for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++)
	{
	    j     = (L->rowind)[jp];
	    Aij   = (L->values.d)[jp] ;
	    b[i] -= x[j]*Aij;
	}

	jp   = (L->colptr)[i];
	j    = (L->rowind)[jp];
	Aii  = (L->values.d)[jp]; 
	x[i] = b[i] / Aii; 
    }

    return 0;
}
int my_taucs_cmsd(double *cmean, double *csd, int idx, taucs_ccs_matrix *L, double *x)
{
    int j, jp;
    double  Aij, Aii, b;
    
    for (b=0.0, jp = (L->colptr)[idx]+1; jp < (L->colptr)[idx+1]; jp++)
    {
	j   = (L->rowind)[jp];
	Aij = (L->values.d)[jp] ;
	b  -= x[j]*Aij;
    }
	
    jp      = (L->colptr)[idx];
    Aii     = (L->values.d)[jp]; 
    *cmean  = b/Aii;
    *csd    = 1/Aii;

    return 0;
}
int my_taucs_check_flags(int flags)
{
#define CheckFLAGS(X) if (flags & X) printf(#X " is ON\n");if (!(flags & X)) printf(#X " is OFF\n")
    CheckFLAGS(TAUCS_INT);
    CheckFLAGS(TAUCS_DOUBLE);
    CheckFLAGS(TAUCS_SINGLE);
    CheckFLAGS(TAUCS_DCOMPLEX);
    CheckFLAGS(TAUCS_SCOMPLEX);
    CheckFLAGS(TAUCS_LOWER);
    CheckFLAGS(TAUCS_UPPER);
    CheckFLAGS(TAUCS_TRIANGULAR);
    CheckFLAGS(TAUCS_SYMMETRIC);
    CheckFLAGS(TAUCS_HERMITIAN);
    CheckFLAGS(TAUCS_PATTERN);
#undef CheckFLAGS
    return GMRFLib_SUCCESS;
}
int GMRFLib_bitmap_factorisation_TAUCS__internal(taucs_ccs_matrix *L, char *filename)
{
#define NBitsInByte 8
#define SETBIT(im,jm) { int ii; ii = (im)/NBitsInByte; \
                      GMRFLib_setbit(&bitmap[ ii+(jm)*m], NBitsInByte-1-((im)-ii*NBitsInByte)); }
    
    int i, j, jp, n = L->n, m;
    unsigned char *bitmap;
    FILE *fp;
    
    m = n/NBitsInByte;
    if (m*NBitsInByte != n) m++;
    bitmap = Calloc(n*m, unsigned char); MEMCHK(bitmap);

    for (i=0; i<n; i++)
    {
	SETBIT(i, i);
	for (jp = (L->colptr)[i]+1; jp < (L->colptr)[i+1]; jp++)
	{
	    j = (L->rowind)[jp];
	    SETBIT(i, j);
	}
    }
    
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
}
int GMRFLib_bitmap_factorisation_TAUCS(char *filename_body, taucs_ccs_matrix *L)
{
    /* 
       create a bitmap-file of the factorization
    */
    char filename[1024+1];

    sprintf(filename, "%s_L.pbm",(filename_body ? filename_body : "taucs_L"));
    GMRFLib_bitmap_factorisation_TAUCS__internal(L, filename);

    return GMRFLib_SUCCESS;
}

#undef GMRFLib_NSET_LIMIT
