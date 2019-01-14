/* GMRFLib-sparse-interface.h
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
 * RCSId: $Id: sparse-interface.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file sparse-interface.h
  \brief Typedefs and defines for the interface to the sparse-matrix libraries.
*/

#ifndef __GMRFLib_SPARSE_INTERFACE_H__
#define __GMRFLib_SPARSE_INTERFACE_H__

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

/*!
  \brief Lapack's band-solver
*/
#define GMRFLib_SMTP_BAND     0 

/*!
  \brief An empty template
*/

#define GMRFLib_SMTP_PROFILE  1 
/*!
  \brief The TAUCS library
*/

#define GMRFLib_SMTP_TAUCS    2  


int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph);
int GMRFLib_free_reordering(GMRFLib_sm_fact_tp *sm_fact);
int GMRFLib_build_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, GMRFLib_graph_tp *graph);
int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph);
int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact);
int GMRFLib_solve_lt_sparse_matrix(double *rhs, GMRFLib_sm_fact_tp *fact_tp,  GMRFLib_graph_tp *graph);
int GMRFLib_solve_llt_sparse_matrix(double *rhs, GMRFLib_sm_fact_tp *fact_tp,  GMRFLib_graph_tp *graph);
int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx, int remapped);
int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph);
int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph);
int GMRFLib_compute_Qinv(GMRFLib_problem_tp *problem, int storage);


int GMRFLib_bitmap_factorisation(char *filename_body, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph);

__END_DECLS
#endif
