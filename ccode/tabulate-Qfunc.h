/* tabulate-Qfunc.h
 * 
 * Copyright (C) 2004 Havard Rue
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
 * RCSId: $Id: tabulate-Qfunc.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file tabulate-Qfunc.h
  \brief Typedefs for \ref tabulate-Qfunc.c
*/

#ifndef __GMRFLib_TABULATE_QFUNC_H__
#define __GMRFLib_TABULATE_QFUNC_H__


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

#include "hashP.h"


typedef struct
{
    int n;					  /* the size of the graph */
    map_id **values;				  /* hash-table for the values */
    double *scale;				  /* scaling */
}
GMRFLib_tabulate_Qfunc_arg_tp;

/*!
   \struct GMRFLib_tabulate_Qfunc_tp

   \brief The structure retured by \c GMRFLib_tabulate_Qfunc()
*/
typedef struct
{
    /*!
      \brief The Qfunction which returns the tabulated values
    */
    GMRFLib_Qfunc_tp *Qfunc;
    /*!
      \brief The arguments to GMRFLib_tabulate_Qfunc_tp::Qfunc
    */
    char *Qfunc_arg;		
}
GMRFLib_tabulate_Qfunc_tp;


double GMRFLib_tabulate_Qfunction(int node, int nnode, char *arg);
int GMRFLib_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp **tabulated_Qfunc, GMRFLib_graph_tp *graph,
			   GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, double *scale);
int GMRFLib_free_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp *tabulated_Qfunc);


__END_DECLS
#endif
