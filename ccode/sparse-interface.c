/* GMRFLib-sparse-interface.c
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
  \file sparse-interface.c
  \brief Unified interface to the sparse-matrix libraries
*/

#include "GMRFLib.h"
#include "GMRFLibP.h"

static const char RCSId[] = "$Id: sparse-interface.c 1 2013-03-28 13:54:24Z hanne $";

#define SMTP_OK(s) (((s) == GMRFLib_SMTP_BAND) || ((s) == GMRFLib_SMTP_PROFILE) || ((s) == GMRFLib_SMTP_TAUCS))

/*!
  \brief Compute the reordering
*/
int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_compute_reordering_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph);
	break;
    }
    if (sm_fact->remap)				  /* need this still for the wa-routines. FIXME */
	sm_fact->bandwidth = GMRFLib_compute_bandwidth(graph, sm_fact->remap);
    else
	sm_fact->bandwidth = -1;
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Free the reordering
*/
int GMRFLib_free_reordering(GMRFLib_sm_fact_tp *sm_fact)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;

    if (sm_fact)
    {
        GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
        switch(GMRFLib_smtp)
        {
        case GMRFLib_SMTP_BAND:
            retval = GMRFLib_free_reordering_BAND(sm_fact);
            break;
        case GMRFLib_SMTP_PROFILE:
            retval = GMRFLib_free_reordering_PROFILE();
            break;
        case GMRFLib_SMTP_TAUCS:
            retval = GMRFLib_free_reordering_TAUCS(sm_fact);
            break;
        }
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*
  \brief Build a sparse matrix
*/
int GMRFLib_build_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact,
				GMRFLib_Qfunc_tp *Qfunc,
				char *Qfunc_arg,
				GMRFLib_graph_tp *graph)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_build_sparse_matrix_BAND(&(sm_fact->bchol), Qfunc, Qfunc_arg, graph, sm_fact->remap, sm_fact->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_build_sparse_matrix_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_build_sparse_matrix_TAUCS(&(sm_fact->L), Qfunc, Qfunc_arg, graph, sm_fact->remap);
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Factorise a sparse matrix
*/
int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{ 
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_factorise_sparse_matrix_BAND(sm_fact->bchol, &(sm_fact->finfo), graph, sm_fact->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_factorise_sparse_matrix_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_factorise_sparse_matrix_TAUCS(&(sm_fact->L), &(sm_fact->symb_fact), &(sm_fact->finfo));
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Free a factorisation of a sparse matrix
*/
int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    if (sm_fact) 
    {
	GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);

	switch(GMRFLib_smtp)
	{
	case GMRFLib_SMTP_BAND:
	    retval         = GMRFLib_free_fact_sparse_matrix_BAND(sm_fact->bchol);
	    sm_fact->bchol = NULL;
	    break;
	case GMRFLib_SMTP_PROFILE:
	    retval = GMRFLib_free_fact_sparse_matrix_PROFILE();
	    break;
	case GMRFLib_SMTP_TAUCS:
	    retval = GMRFLib_free_fact_sparse_matrix_TAUCS(sm_fact->L, sm_fact->symb_fact);
	    sm_fact->L         = NULL;
	    sm_fact->symb_fact = NULL;
	    break;
	default:
	    GMRFLib_ERROR(GMRFLib_ESNH);
	    break;
	}
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Solve \f$L^Tx=b\f$
*/
int GMRFLib_solve_lt_sparse_matrix(double *rhs, GMRFLib_sm_fact_tp *fact_tp,  GMRFLib_graph_tp *graph)
{
    /* 
       rhs in real world. solve L^Tx=rhs, rhs is overwritten by the solution
    */
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_solve_lt_sparse_matrix_BAND(rhs, fact_tp->bchol,  graph, fact_tp->remap, fact_tp->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_solve_lt_sparse_matrix_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_solve_lt_sparse_matrix_TAUCS(rhs, fact_tp->L, graph, fact_tp->remap);
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Solve \f$LL^Tx=b\f$  or \f$Qx=b\f$
*/
int GMRFLib_solve_llt_sparse_matrix(double *rhs, GMRFLib_sm_fact_tp *fact_tp,  GMRFLib_graph_tp *graph)
{
    /* 
       rhs in real world. solve  Q x=rhs, where Q=L L^T
    */
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;

    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_solve_llt_sparse_matrix_BAND(rhs, fact_tp->bchol,  graph, fact_tp->remap, fact_tp->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_solve_llt_sparse_matrix_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_solve_llt_sparse_matrix_TAUCS(rhs, fact_tp->L,  graph, fact_tp->remap);
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Solve \f$L^Tx=b\f$ for indices in an interval
*/
int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact,
			  GMRFLib_graph_tp *graph, 
			  int findx, int toindx, int remapped)
{
    /* 
       rhs in real world, bchol in mapped world. solve L^Tx=b backward only from rhs[findx] up to
       rhs[toindx]. note that findx and toindx is in mapped world. if remapped, do not
       remap/remap-back the rhs before solving.
       
       this routine is called to many times and the work is not that much, to justify
       GMRFLib_ENTER_ROUTINE;
     */
    int retval=GMRFLib_SUCCESS;

    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_solve_lt_sparse_matrix_special_BAND(rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_solve_lt_sparse_matrix_special_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_solve_lt_sparse_matrix_special_TAUCS(rhs, sm_fact->L, graph, sm_fact->remap, findx, toindx, remapped);
	break;
    }
    /* GMRFLib_LEAVE_ROUTINE; */
    return retval;
}
/*!
  \brief Compute the log determininant of \f$Q\f$
*/
int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;

    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_log_determinant_BAND(logdet, sm_fact->bchol, graph, sm_fact->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_log_determinant_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_log_determinant_TAUCS(logdet, sm_fact->L);
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
}
/*!
  \brief Compute conditional mean and standard deviation of \f$x[i]\f$ conditioned on {\f$x[j]\f$}
    for \f$j>i\f$
*/
int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped,
			     GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph) 
{
    int retval=GMRFLib_SUCCESS;
    /* 
       this routine is called to many times, and the work is minor, to justify
       GMRFLib_ENTER_ROUTINE;
    */

    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_comp_cond_meansd_BAND(cmean, csd, indx, x, remapped, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_comp_cond_meansd_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_comp_cond_meansd_TAUCS(cmean, csd, indx, x, remapped, sm_fact->L, graph, sm_fact->remap);
	break;
    }
    /* GMRFLib_LEAVE_ROUTINE; */
    return retval;
}
/*!
  \brief Produce a bitmap of the Cholesky triangle in the portable bitmap (pbm) format
*/
int GMRFLib_bitmap_factorisation(char *filename_body, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
    int retval=GMRFLib_SUCCESS;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_bitmap_factorisation_BAND(filename_body, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_bitmap_factorisation_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_bitmap_factorisation_TAUCS(filename_body, sm_fact->L);
	break;
    }
    return retval;
}
/*!
  \brief Wrapper for computing the (structural) inverse of \c Q.
*/
int GMRFLib_compute_Qinv(GMRFLib_problem_tp *problem, int storage)
{
    int retval=GMRFLib_SUCCESS;
    GMRFLib_ENTER_ROUTINE;
    
    GMRFLib_ASSERT(SMTP_OK(GMRFLib_smtp), GMRFLib_ESMTP);
    switch(GMRFLib_smtp)
    {
    case GMRFLib_SMTP_BAND:
	retval = GMRFLib_compute_Qinv_BAND();
	break;
    case GMRFLib_SMTP_PROFILE:
	retval = GMRFLib_compute_Qinv_PROFILE();
	break;
    case GMRFLib_SMTP_TAUCS:
	retval = GMRFLib_compute_Qinv_TAUCS(problem, storage);
	break;
    }
    GMRFLib_LEAVE_ROUTINE;
    return retval;
} 
