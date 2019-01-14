/* rw.c
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
  \file rw.c
  \brief Handy functions when using RW1, RW2, CRW1, CRW2 and approximate CRW2 models.
*/


static const char RCSId[] = "$Id: rw.c 1 2013-03-28 13:54:24Z hanne $";

#include "GMRFLib.h"
#include "GMRFLibP.h"


/*!
  \brief This function returns element Q_ij, with i=node and j=nnode, of the precision matrix for
  the RW defined in \c rwdef.

  This funcion allow also arguments ij that are not neighbours and can therefore also be used for
  pruning (using \ref GMRFLib_prune_graph()).

  \param[in] node   First node
  \param[in] nnode  Second node
  \param[in] rwdef  The definition of the RW1 or RW2

  \sa \ref GMRFLib_rwdef_tp, \ref GMRFLib_make_rw_graph, \ref GMRFLib_prune_graph
*/
double GMRFLib_rw(int node, int nnode, GMRFLib_rwdef_tp *rwdef)
{
    int imax, imin, idiff, edge;
    double prec;

    prec = (rwdef->prec ? *rwdef->prec : 1.0);
    
    GMRFLib_ASSERT(rwdef && rwdef->order >= 0 && rwdef->order <= 2, GMRFLib_EPARAMETER);

    if (rwdef->cyclic)
    {
	/* 
	   cyclic, simpler..
	*/

	idiff = MIN( ABS(node-nnode),  ABS(MIN(node, nnode)-MAX(node, nnode) + rwdef->n));

	switch(rwdef->order)
	{
	case 1:
	    switch(idiff)
	    {
	    case  0: return 2.0*prec;
	    case  1: return    -prec;
	    default: return 0.0;
	    }
	case 2:
	    switch(idiff)
	    {
	    case  0: return  6.0*prec;
	    case  1: return -4.0*prec;
	    case  2: return      prec;
	    default: return 0.0;
	    }
	}
    }
    else
    {
	imax  = MAX(node, nnode);
	imin  = MIN(node, nnode);
	idiff = imax - imin;
    
	if (idiff > rwdef->order) return 0.0;		  /* fast return */
    
	edge = rwdef->order;
    
	if (rwdef->order == 0 || (imax > edge && imax < rwdef->n -rwdef->order-1))
	{
	    /* 
	       internal node
	    */
	    switch(rwdef->order)
	    {
	    case 0: return prec*1.0; 
	    case 1: return prec*(idiff == 0 ? 2.0 : -1.0);
	    case 2: return prec*(idiff == 0 ? 6.0 : (idiff == 1 ? -4.0 : 1.0));
	    default: GMRFLib_ASSERT(0, GMRFLib_EINVARG);
	    }
	}
	else
	{
	    /* 
	       the edge.
	    */
	    if (imax > edge)
	    {
		/* 
		   map to the left egde
		*/
		imax = rwdef->n -1 -MIN(node, nnode);
		imin = rwdef->n -1 -MAX(node, nnode);
	    }
	    switch(rwdef->order)
	    {
	    case 1:
		return prec*(idiff == 1 ? -1.0 : (imax == 0 ? 1.0 : 2.0));
	    case 2:
		switch(idiff)
		{
		case 0:
		    switch(imax)
		    {
		    case 0: return prec*1.0;
		    case 1: return prec*5.0;
		    case 2: return prec*6.0;
		    default: GMRFLib_ASSERT(0, GMRFLib_EINVARG);
		    }
		case 1: return prec*(imin == 0 ? -2.0 : -4.0);
		case 2: return prec*1.0;
		default: GMRFLib_ASSERT(0, GMRFLib_EINVARG);
		}
	    default: GMRFLib_ASSERT(0, GMRFLib_EINVARG);
	    }
	}
    }
    GMRFLib_ASSERT(0, GMRFLib_EINVARG);

    return 0.0;
}
/*!
  \brief This function returns element Q_ij, with i=node and j=nnode, of the precision matrix for
  the CRW defined in \c rwdef.
 
  This funcion allow also arguments ij that are not neighbours and can therefore also be used for
  pruning (using \ref GMRFLib_prune_graph()).
 
  \param[in]  node  First node
  \param[in] nnode  Second node
  \param[in] crwdef The definition of the CRW1 or CRW2
 
  \sa \ref GMRFLib_crwdef_tp,  \ref GMRFLib_make_crw_graph,  \ref GMRFLib_prune_graph
*/
double GMRFLib_crw(int node, int nnode, GMRFLib_crwdef_tp *crwdef)
{
    /* 
       this is the continous version for order=1 or 2, which take into accont the positions
       consistently. the order=2 uses the augmentation with the velocity. be aware to set the
       pointers idelta, idelta2 and idelta3 to ZERO or compute them according to the description in
       the rw.h.

       TODO: make this also accept a cyclic argument. there should not be to much to change here. it
       is possible that it's easier if positions are present, then we can use the general expression
       only `defining' the work-entries correct.  in any case, then this function could be merged
       with GMRFLib_rw()? do we then need both?

       OOPS: if this is change to accept the cyclic case as well, then rememeber to change the
       make_graph routine as well!
    */
#define EVEN(num) (((num)/2)*2 == num)
#define TP_POS 0	
#define TP_VEL 1
#define SETUP_WORK_PTRS if (1){			\
	idelta  = &crwdef->work[2          ];	\
	idelta2 = &crwdef->work[2 +   (n+4)];	\
	idelta3 = &crwdef->work[2 + 2*(n+4)];	\
        isdelta = &crwdef->work[2 + 3*(n+4)];	\
        sidelta = &crwdef->work[2 + 4*(n+4)];}

    int i, use_pos, n, idiff, imin, imax, node_i=-1, node_tp=-1, nnode_i=-1, nnode_tp=-1, order;
    double prec;

    double *idelta = NULL, *idelta2 = NULL, *idelta3 = NULL, *isdelta = NULL, *sidelta = NULL;

    GMRFLib_ASSERT(crwdef && crwdef->order > 0 && crwdef->order <= 2, GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(crwdef->n > crwdef->order, GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(crwdef->layout == GMRFLib_CRW_LAYOUT_SIMPLE ||
		   crwdef->layout == GMRFLib_CRW_LAYOUT_PAIRS  ||
		   crwdef->layout == GMRFLib_CRW_LAYOUT_BLOCK, GMRFLib_EPARAMETER);

    prec    = (crwdef->prec ? *crwdef->prec : 1.0);
    n       = crwdef->n;
    use_pos = (crwdef->position ? 1 : 0);
    order   = crwdef->order;
    
    if (crwdef->position)			  
    {
	/* 
	   check `workspace', compute it if we don't have it.
	*/
	if (!crwdef->work)
	{
	    /* 
	       we do this simple; we compute the workspace for the union of all cases.
	    */
	    double *delta;

	    delta= Calloc(n-1, double); MEMCHK(delta);
	    for(i=0;i<n-1;i++) delta[i] = crwdef->position[i+1] - crwdef->position[i];
	     
	    crwdef->work = Calloc(5*(n+4), double); MEMCHK(crwdef->work);
	    SETUP_WORK_PTRS;

	    for(i=0;i<n-1;i++) 
	    {
		idelta[i]  = 1.0/delta[i];
		idelta2[i] = SQR(idelta[i]);
		idelta3[i] = idelta[i]*idelta2[i];
	    }
	    for(i=0;i<n-2;i++) 
	    {
		sidelta[i] = idelta[i]+idelta[i+1];
		isdelta[i] = 1.0/(delta[i]+delta[i+1]);
	    }
	    sidelta[n-2] = idelta[n-2];	  /* yes */
	    sidelta[-1]  = idelta[0];	  /* yes */

	    FREE(delta);
	}
	else
	    SETUP_WORK_PTRS;
    }

    if (order == 1)
    {
	/* 
	   first order model, this is the simple case with no agumentation 
	*/

	imin  = MIN(node, nnode);
	imax  = MAX(node, nnode);
    	idiff = imax - imin;
	if (idiff > order) return 0.0;		  /* fast return, nothing here */

	if (idiff == 0)
	{
	    if (use_pos)
	    {
		if (imin == 0)   return prec*idelta[0];
		if (imin == n-1) return prec*idelta[n-2];
		return prec*(idelta[imin-1] + idelta[imin]);
	    }
	    else
		return prec*((imin == 0 || imin == n-1) ? 1.0: 2.0);
	}
	else
	    return prec*(use_pos ? -idelta[imin] : -1.0);
    }

    /* 
       the second order model is a bit more involved

       first there is the switch beteen approximative/exact, then layout, the regular or irregular.
    */
    if (crwdef->layout == GMRFLib_CRW_LAYOUT_SIMPLE)
    {
	/* 
	   this is the approximative scheme with no augmentation.

	   do first the case with regular positions the irregular. the regular case is just a copy
	   from GMRFLib_rw()
	*/

	imax  = MAX(node, nnode);
	imin  = MIN(node, nnode);
	idiff = imax - imin;
	if (idiff > order) return 0.0;		  /* fast return, nothing here */
    
	if (!use_pos)
	{
	    /* 
	       regular
	    */
	    if ((imax > order && imax < n -order -1))
	    {
		/* 
		   internal node
		*/
		return prec*(idiff == 0 ? 6.0 : (idiff == 1 ? -4.0 : 1.0));
	    }
	    else
	    {
		if (imax > order)
		{
		    /* 
		       map to the left egde
		    */
		    imax = n -1 -MIN(node, nnode);
		    imin = n -1 -MAX(node, nnode);
		}
		switch(idiff)
		{
		case 0:
		    switch(imax)
		    {
		    case 0: return prec*1.0;
		    case 1: return prec*5.0;
		    case 2: return prec*6.0;
		    default: GMRFLib_ASSERT(0, GMRFLib_ESNH);
		    }
		case 1: return prec*(imin == 0 ? -2.0 : -4.0);
		case 2: return prec*1.0;
		default: GMRFLib_ASSERT(0, GMRFLib_ESNH);
		}
	    }
	}
	else
	{
	    /* 
	       irregular case
	    */
	    switch(idiff)
	    {
	    case 0:
		return  prec*2.0*(idelta2[imax-1]*isdelta[imax-2]
				  + idelta[imax-1]*idelta[imax]*sidelta[imax-1]
				  + idelta2[imax]*isdelta[imax]);
	    case 1:
		return -prec*2.0*idelta2[imax-1]*(idelta[imax-2] + idelta[imax]);
		
	    case 2:
		return  prec*2.0*idelta[imax-2]*idelta[imax-1]*isdelta[imax-2];

	    default:
		GMRFLib_ASSERT(0, GMRFLib_ESNH);
	    }
	}
    }
    else
    {
	/* 
	   this is the exact solution, which require augmentation
	*/
	
	switch(crwdef->layout)
	{
	    /* 
	       we need to make sure that node_i < nnode_i
	    */
	
	case GMRFLib_CRW_LAYOUT_PAIRS:
	    if (node < nnode)
	    {
		node_i   = node/2;
		nnode_i  = nnode/2;
		node_tp  = (EVEN(node)  ? TP_POS : TP_VEL);
		nnode_tp = (EVEN(nnode) ? TP_POS : TP_VEL);
	    }
	    else
	    {
		node_i   = nnode/2;
		nnode_i  = node/2;
		node_tp  = (EVEN(nnode) ? TP_POS : TP_VEL);
		nnode_tp = (EVEN(node)  ? TP_POS : TP_VEL);
	    }
	    break;
	case GMRFLib_CRW_LAYOUT_BLOCK:
	    if (node%n < nnode%n)			 
	    {
		node_i   = node  % n;
		nnode_i  = nnode % n;
		node_tp  = (node  < n ? TP_POS : TP_VEL);
		nnode_tp = (nnode < n ? TP_POS : TP_VEL);
	    }
	    else
	    {
		node_i   = nnode % n;
		nnode_i  = node  % n;
		node_tp  = (nnode < n ? TP_POS : TP_VEL);
		nnode_tp = (node  < n ? TP_POS : TP_VEL);
	    }
	    break;
	default:
	    GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	idiff = nnode_i - node_i;  assert(idiff >= 0);

	if(idiff > 1) return 0.0;			  /* nothing here */

	/* 
	   split in two cases to speedup
	*/
    
	if (use_pos)
	{
	    if (idiff == 0)
	    {
		if (node_tp == TP_POS)
		{
		    if (nnode_tp == TP_POS)
		    {					  /* TP_POS & TP_POS  */
			if (nnode_i == 0)   return prec*12.0*idelta3[0];
			if (nnode_i == n-1) return prec*12.0*idelta3[n-2];
			return prec*12.0*(idelta3[node_i] + idelta3[node_i-1]);
		    }
		    else
		    {					  /* TP_POS & TP_VEL */
			if (nnode_i == 0)   return prec*6.0*idelta2[0];
			if (nnode_i == n-1) return prec*(-6.0)*idelta2[n-2];
			return prec*6.0*(idelta2[node_i] -idelta2[node_i-1]);
		    }
		}
		else
		{
		    if (nnode_tp == TP_POS)
		    {					  /* TP_VEL & TP_POS  */
			if (nnode_i == 0)   return prec*6.0*idelta2[0];
			if (nnode_i == n-1) return prec*(-6.0*idelta2[n-2]);
			return prec*6.0*(idelta2[node_i] -idelta2[node_i-1]);
		    }
		    else
		    {					  /* TP_VEL & TP_VEL */
			if (nnode_i == 0)   return prec*4.0*idelta[0];
			if (nnode_i == n-1) return prec*4.0*idelta[n-2];
			return prec*4.0*(idelta[node_i] + idelta[node_i-1]);
		    }
		}
	    }
	    else
	    {
		if (node_tp == TP_POS)
		{
		    if (nnode_tp == TP_POS)
			return prec*(-12.0)*idelta3[node_i];
		    else
			return prec*6.0*idelta2[node_i];
		}
		else
		{					  
		    if (nnode_tp == TP_POS)
			return prec*(-6.0)*idelta2[node_i];
		    else
			return prec*2.0*idelta[node_i];
		}
	    }
	}
	else
	{
	    if (idiff == 0)
	    {
		if (node_tp == TP_POS)
		{
		    if (nnode_tp == TP_POS)
		    {					  /* TP_POS & TP_POS  */
			if (nnode_i == 0)   return prec*12.0;
			if (nnode_i == n-1) return prec*12.0;
			return prec*24.0;
		    }
		    else
		    {					  /* TP_POS & TP_VEL */
			if (nnode_i == 0)   return prec*6.0;
			if (nnode_i == n-1) return prec*(-6.0);
			return 0.0;
		    }
		}
		else
		{
		    if (nnode_tp == TP_POS)
		    {					  /* TP_VEL & TP_POS  */
			if (nnode_i == 0)   return prec*6.0;
			if (nnode_i == n-1) return prec*(-6.0);
			return 0.0;
		    }
		    else
		    {					  /* TP_VEL & TP_VEL */
			if (nnode_i == 0)   return prec*4.0;
			if (nnode_i == n-1) return prec*4.0;
			return prec*8.0;
		    }
		}
	    }
	    else
	    {
		if (node_tp == TP_POS)
		{
		    if (nnode_tp == TP_POS)
			return prec*(-12.0); /* POS & POS */
		    else
			return prec*6.0; /* POS & VEL */
		}
		else
		{					  
		    if (nnode_tp == TP_POS)
			return prec*(-6.0); /* VEL & POS */
		    else
			return prec*2.0; /* VEL & VEL */
		}
	    }
	}
	return 0.0;
    }

    GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);

#undef SETUP_WORK_PTRS
#undef EVEN
#undef TP_POS
#undef TP_VEL
}

/*
  Example for manual
 */

/*! \page ex_rw A worked out example smoothing a time-series data, using the routines in rw.c
  
Solve the same problem as in \ref ex_wa, now using the routines in rw.c

\par Program code:

\verbinclude example_doxygen_rw_1.txt

*/
