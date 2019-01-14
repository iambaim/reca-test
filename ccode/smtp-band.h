/* GMRFLib-smtp-band.h
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
 * RCSId: $Id: smtp-band.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file smtp-band.h
  \brief GMRFLib interface to the band-solver in the LAPACK library
*/

#ifndef __GMRFLib_SMTP_BAND_H__
#define __GMRFLib_SMTP_BAND_H__

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

void gpskca_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
void dtbsvspecial_(const char *, const char *, const char *, int *, int *, double *, int *, double *, int *, int *, int *, int, int, int);  
int cmsd_(double *, double *, int *, int *, int *, double *, int *, double *);

int GMRFLib_compute_reordering_BAND(int **remap, GMRFLib_graph_tp *graph);
int GMRFLib_free_reordering_BAND(GMRFLib_sm_fact_tp *sm_fact);
int GMRFLib_build_sparse_matrix_BAND(double **bandmatrix, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, GMRFLib_graph_tp *graph, int *remap, int bandwidth);
int GMRFLib_factorise_sparse_matrix_BAND(double *band,  GMRFLib_fact_info_tp *finfo, GMRFLib_graph_tp *graph, int bandwidth);
int GMRFLib_free_fact_sparse_matrix_BAND(double *bchol);
int GMRFLib_solve_lt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth);
int GMRFLib_solve_llt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth);
int GMRFLib_solve_lt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth, int findx, int toindx, int remapped);
int GMRFLib_log_determinant_BAND(double *logdet, double *bchol, GMRFLib_graph_tp *graph, int bandwidth);
int GMRFLib_comp_cond_meansd_BAND(double *cmean, double *csd, int indx, double *x, int remapped, double *bchol, GMRFLib_graph_tp *graph, int *remap, int bandwidth);
int GMRFLib_compute_Qinv_BAND(void);
int GMRFLib_bitmap_factorisation_BAND__iternal(char *filename, double *band, GMRFLib_graph_tp *graph, int *remap, int bandwidth);
int GMRFLib_bitmap_factorisation_BAND(char *filename_body, double *band, GMRFLib_graph_tp *graph, int *remap, int bandwidth);



__END_DECLS
#endif
