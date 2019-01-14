/* GMRFLib-smtp-taucs.h
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
 * RCSId: $Id: smtp-taucs.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file smtp-taucs.h
  \brief Typedefs and defines for the \ref smtp-taucs.c, which is the GMRFLib interface to the
  sparse matrix library TAUCS.
*/

#ifndef __GMRFLib_SMTP_TAUCS_H__
#define __GMRFLib_SMTP_TAUCS_H__

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

int GMRFLib_compute_reordering_TAUCS(int **remap, GMRFLib_graph_tp *graph);
int GMRFLib_free_reordering_TAUCS(GMRFLib_sm_fact_tp *sm_fact);
int GMRFLib_build_sparse_matrix_TAUCS(taucs_ccs_matrix **L, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, GMRFLib_graph_tp *graph, int *remap);
int GMRFLib_factorise_sparse_matrix_TAUCS_OLD(taucs_ccs_matrix **L, GMRFLib_fact_info_tp *finfo);
int GMRFLib_factorise_sparse_matrix_TAUCS(taucs_ccs_matrix **L, supernodal_factor_matrix **symb_fact, GMRFLib_fact_info_tp *finfo);
int GMRFLib_free_fact_sparse_matrix_TAUCS(taucs_ccs_matrix *L, supernodal_factor_matrix *symb_fact);
int GMRFLib_free_fact_sparse_matrix_TAUCS_OLD(taucs_ccs_matrix *L);
int GMRFLib_solve_lt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap);
int GMRFLib_solve_llt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap);
int GMRFLib_solve_lt_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap, int findx, int toindx, int remapped);
int GMRFLib_comp_cond_meansd_TAUCS(double *cmean, double *csd, int indx, double *x, int remapped, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap);
int GMRFLib_log_determinant_TAUCS(double *logdet, taucs_ccs_matrix *L);
int GMRFLib_compute_Qinv_validate_TAUCS(GMRFLib_problem_tp *problem, FILE *fp);
int GMRFLib_compute_Qinv_TAUCS(GMRFLib_problem_tp *problem, int storage);
int my_taucs_dccs_solve_lt(void* vL, double* x, double* b);
int my_taucs_check_flags(int flags);
int my_taucs_cmsd(double *cmean, double *csd, int idx, taucs_ccs_matrix *L, double *x);
int my_taucs_dccs_solve_lt_special(void* vL, double* x, double* b, int from_idx, int to_idx);
    
map_ii **GMRFLib_compute_Qinv_TAUCS_check(taucs_ccs_matrix *L);
taucs_ccs_matrix *GMRFLib_compute_Qinv_TAUCS_add_elements(taucs_ccs_matrix *L, map_ii **missing_elements);
int GMRFLib_compute_Qinv_TAUCS_compute(GMRFLib_problem_tp *problem, int storage, taucs_ccs_matrix *Lmatrix);

taucs_ccs_matrix *my_taucs_dccs_copy(taucs_ccs_matrix *L, int flags);
supernodal_factor_matrix *my_taucs_supernodal_factor_matrix_copy(supernodal_factor_matrix *L);

/* 
   internal functions here, not documented
*/
int GMRFLib_bitmap_factorisation_TAUCS__internal(taucs_ccs_matrix *L, char *filename);
int GMRFLib_bitmap_factorisation_TAUCS(char *filename_body, taucs_ccs_matrix *L);


__END_DECLS
#endif
