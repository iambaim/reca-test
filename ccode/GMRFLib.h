/* GMRFLib.h
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
 * RCSId: $Id: GMRFLib.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file GMRFLib.h
  \brief Include all include-files needed for using the GMRFLib in the correct order.
*/

#ifndef __GMRFLib_H__
#define __GMRFLib_H__

#define GMRFLib_VERSION_MAJOR    2
#define GMRFLib_VERSION_MINOR    2
#define GMRFLib_VERSION_REVISION 1
#define GMRFLib_VERSION          "2.2-1"

#include "taucs.h"
#include "compatibility.h"
#include "globals.h"
#include "error-handler.h"
#include "rw.h"
#include "graph.h"

#include "hash.h"
#include "tabulate-Qfunc.h"

#include "lapack-interface.h"
#include "problem-setup.h"
#include "optimize.h"
#include "gdens.h"
#include "hidden-approx.h"
#include "blockupdate.h"
#include "distributions.h"
#include "timer.h"
#include "wa.h"

#include "sparse-interface.h"
#include "smtp-band.h"
#include "smtp-profile.h"
#include "smtp-taucs.h"

#include "bitmap.h"			  /* needs both graph and problem and sparse */
#include "geo.h"

#endif
