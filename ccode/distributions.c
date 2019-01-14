/* distributions.c
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
 *
 */

/*!
  \file distributions.c
  \brief Functions for statistical distributions used in GMRFLib
*/


#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include "GMRFLib.h"
#include "GMRFLibP.h"

static const char RCSId[] = "$Id: distributions.c 1 2013-03-28 13:54:24Z hanne $";

/* 
   functions in the dcdf-library that are used in this file
 */
void cdfchi_(int *which,double *p,double *q,double *x,double *df, int *status,double *bound);


/*!
  \brief The inverse cummulative distribution function for the Chi-square distribution

  \param[in] p The value such that \f$\mbox{Prob}(X \le x) = p\f$.
  \param[in] dof The degrees of freedom

  The returned value is \f$x\f$ such that  \f$\mbox{Prob}(X \le x) = p\f$.
*/
double GMRFLib_invChiSquare(double p, double dof)
{
    int which = 2,  status = 0;
    double q = 1.0 - p,  bound = 0.0, x;

    cdfchi_(&which, &p, &q, &x, &dof, &status, &bound);
    GMRFLib_ASSERT(status == 0, GMRFLib_EDCDFLIB);

    return x;
}
/*!
  \brief Return a sample from the standard normal distribution \f$N(0,1)\f$.
*/
double GMRFLib_stdnormal(void)
{
    return M_SQRT2*cos(2.0*M_PI*GMRFLib_uniform())*sqrt(-log(1.0-GMRFLib_uniform()));
}

    
    
    
