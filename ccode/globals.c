/* globals.c
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
  \file globals.c
  \brief Set values of global variables.
*/


static const char RCSId[] = "$Id: globals.c 1 2013-03-28 13:54:24Z hanne $";

#define __GMRFLib_DONT_DEFINE_GLOBALS
#include "GMRFLib.h"
#include "GMRFLibP.h"
#undef __GMRFLib_DONT_DEFINE_GLOBALS

/*!
  \brief Define the sparse solver for the sparse-matrix computations.

  This variable defines the which method is used to factorise the sparse matrices. The current
  implementation includes
  - #GMRFLib_SMTP_BAND, using the band-matrix routines in \c LAPACK
  - #GMRFLib_SMTP_TAUCS, using the multifrontal supernodal factorisation in the \c TAUCS library.

  Default value is #GMRFLib_SMTP_TAUCS.\n\n
*/
int GMRFLib_smtp = GMRFLib_SMTP_TAUCS;


/*! 

  \brief Set the blas level in the Lapack routines using the band solver.

  When \c GMRFLib_smtp is set to #GMRFLib_SMTP_BAND, then there is an additional option of chosing
  the blas level for the routines used.  The two options are
  - level 2 (#BLAS_LEVEL2),  and
  - level 3 (#BLAS_LEVEL3)
  
  Use level 2 for smaller problems and level 3 for larger problems. Default is #BLAS_LEVEL3.

  \remark For larger problems you may want to use the TAUCS library instead.
 */
int GMRFLib_blas_level = BLAS_LEVEL3;

/*!
  \brief Use internal lookup-tables or not.

  GMRFLib use in the file wa.c internal lookup-tables to speed up the calculations (default).  This
  feature can be turned off by setting \c GMRFLib_use_wa_table_lookup() to \c #GMRFLib_FALSE.\n\n
*/
int GMRFLib_use_wa_table_lookup = GMRFLib_TRUE;

/*! 
  \brief Verify graph after reading.

  GMRFLib can automatically verify a graph after reading it from file. By default this feature is
  off, but can be turned on setting \c GMRFLib_verify_graph_read_from_disc() to \c #GMRFLib_TRUE \n\n
*/
int GMRFLib_verify_graph_read_from_disc = GMRFLib_FALSE;


/*!
  \brief Collect timing statistics

  GMRFLib collects by default various internal statistics of the CPU demanding functions.  This
  statistics can be extracted using the functions \c GMRFLib_timer_full_report() or
  \c GMRFLib_timer_report(). This feature can be turned off by setting 
  \c GMRFLib_collect_timer_statistics()
  to \c #GMRFLib_FALSE. \n\n
*/
int GMRFLib_collect_timer_statistics     = GMRFLib_TRUE;

/*!
  \brief The CPU timing function

  Writing a function that (correctly) returns the used CPU time from a fixed reference, is tricky
  and often OS dependent. You may override the GMRFLib's implementation by letting \c GMRFLib_get_ctime()
  point to your function. \n\n
*/
GMRFLib_get_ctime_tp  *GMRFLib_get_ctime = GMRFLib_ctime_default; 

/* 
   default uniform(0,1) generator and support utilities. these routines are defined in random.c
 */
GMRFLib_uniform_tp           ran;
GMRFLib_uniform_init_tp      ran_init;
GMRFLib_uniform_getstate_tp  ran_getstate;
GMRFLib_uniform_setstate_tp  ran_setstate;

GMRFLib_uniform_tp          *GMRFLib_uniform          = ran;
GMRFLib_uniform_init_tp     *GMRFLib_uniform_init     = ran_init;
GMRFLib_uniform_getstate_tp *GMRFLib_uniform_getstate = ran_getstate;
GMRFLib_uniform_setstate_tp *GMRFLib_uniform_setstate = ran_setstate;
