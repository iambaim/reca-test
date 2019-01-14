/* tabulate-Qfunc.c
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
 */

/*!
  \file tabulate-Qfunc.c
  \brief Tabulate a Qfunction. Useful for speeding up costly Qfunctions.

  These functions provide a tool for tabulate a Qfunction to gain speeup. On partiuclar example are
  those computed with \c GMRFLib_init_wa_problem(), which provides a convenient way to compute
  Qfunctions, but can sometimes be slow if the elements, perhaps apart from a common scale, does not
  change. For these cases then \c GMRFLib_tabulate_Qfunc provides a HUGE speedup.

  The functions are
  - \c GMRFLib_tabulate_Qfunc()  tabulate a Qfunction and return a \c GMRFLib_tabulate_Qfunc_tp object
  - \c GMRFLib_free_tabulate_Qfunc()  free's a \c GMRFLib_tabulate_Qfunc_tp object

  Example:
  \verbatim

  double scale;
  GMRFLib_wa_problem_tp *wa_problem = NULL;
  GMRFLib_tabulate_Qfunc_tp *tab = NULL;
  
  GMRFLib_init_wa_problem(&wa_problem, wagraph, waQfunc, NULL); // init the wa-problem, here scale=1
  GMRFLib_tabulate_Qfunc(&tab, wa_problem->graph, wa_problem->Qfunc, wa_problem->Qfunc_arg, &scale);

  scale = 2.0;					  
  for(i=0;i<wa_problem->graph->n;i++)
  {
      for(jj=0;jj<wa_problem->graph->nnbs[i];jj++)
      {
          j = wa_problem->graph->nbs[i][jj];
	  printf("Qfunc(%d,%d) = (wa) %f (tabulated) %f\n", i, j,
	       (*wa_problem->Qfunc)(i, j, wa_problem->Qfunc_arg), (*tab->Qfunc)(i, j, tab->Qfunc_arg)/scale);
       }
   }

  \endverbatim
*/

#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

static const char RCSId[] = "$Id: tabulate-Qfunc.c 1 2013-03-28 13:54:24Z hanne $";

#include "GMRFLib.h"
#include "GMRFLibP.h"

double GMRFLib_tabulate_Qfunction(int node, int nnode, char *arg)
{
    GMRFLib_tabulate_Qfunc_arg_tp *args;
    args = (GMRFLib_tabulate_Qfunc_arg_tp *) arg;
    return (*args->scale)*(*map_id_ptr(args->values[MIN(node, nnode)], MAX(node, nnode)));
}

/*!
  \brief Tabulate a Qfunction to gain speedup.

  \param[out] tabulate_Qfunc At output, \a tabulate_Qfunc is allocated as a pointer to a \c
  GMRFLib_tabulate_Qfunc_tp, holding a pointer to the Qfunction (GMRFLib_tabulate_Qfunc_tp::Qfunc)
  and its argument (GMRFLib_tabulate_Qfunc_tp::Qfunc_arg) which evaluate the tabulate Qfunction.

  \param[in] graph The graph

  \param[in] Qfunc  The Qfunction to be tabulate

  \param[in] Qfunc_arg The argument to the Qfunction to the tabulate

  \param[in] scale An optional argument which scales the Qfunction. If \c scale is not \c NULL, then
  GMRFLib_tabulate_Qfunc_tp::Qfunc returns \c Qfunc times *scale. If \c scale is \c NULL, then
  *scale set to 1.

*/
int GMRFLib_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp **tabulate_Qfunc, GMRFLib_graph_tp *graph,
			   GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg, double *scale)
{
    int i, j, k, count;
    static double one = 1.0;			  /* yes, this have to be static */
    GMRFLib_tabulate_Qfunc_arg_tp *arg;
    
    GMRFLib_ASSERT(tabulate_Qfunc, GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(graph, GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(Qfunc, GMRFLib_EPARAMETER);

    if (!scale) scale = &one;			  /* provide a default value */
    
    *tabulate_Qfunc              = Calloc(1, GMRFLib_tabulate_Qfunc_tp); MEMCHK(tabulate_Qfunc);
    (*tabulate_Qfunc)->Qfunc     = GMRFLib_tabulate_Qfunction; /* the Qfunction to use */
    arg                          = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp); MEMCHK(arg);
    (*tabulate_Qfunc)->Qfunc_arg = (char *)arg;
    
    arg->n      = graph->n;
    arg->scale  = scale;
    arg->values = Calloc(graph->n, map_id*); MEMCHK(arg->values);

    for(i=0;i<graph->n;i++)
    {
	arg->values[i] = Calloc(1, map_id); MEMCHK(arg->values[i]); /* allocate hash-table */

	for(j=0, count=1;j<graph->nnbs[i];j++)	  /* count the number of terms in the hash-table */
	    if (graph->nbs[i][j] > i) count++;

	map_id_init_hint(arg->values[i], count); /* init hash with count number of elms */

	map_id_set(arg->values[i], i, (*Qfunc)(i, i, Qfunc_arg)); /* diagonal */
	for(j=0, count=1;j<graph->nnbs[i];j++)
	{
	    k = graph->nbs[i][j];
	    if (k > i)				  /* store only those */
		map_id_set(arg->values[i], k, (*Qfunc)(i, k, Qfunc_arg));
	}
    }

    return GMRFLib_SUCCESS;
}
/*!
  \brief Free a tabulate_Qfunc-object

  \param[in] tabulate_Qfunc The object of type \c GMRFLib_tabulate_Qfunc_tp to be freed.
*/
int GMRFLib_free_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp *tabulate_Qfunc)
{
    int i;
    GMRFLib_tabulate_Qfunc_arg_tp *arg;

    if (!tabulate_Qfunc) return GMRFLib_SUCCESS;
    arg = (GMRFLib_tabulate_Qfunc_arg_tp *) tabulate_Qfunc->Qfunc_arg;
    for(i=0;i<arg->n;i++) 
    {
	map_id_free(arg->values[i]); FREE(arg->values[i]);
    }
    FREE(arg->values); FREE(arg); FREE(tabulate_Qfunc);

    return GMRFLib_SUCCESS;
}
