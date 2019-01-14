/* compatibility.h
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
 * RCSId: $Id: compatibility.h 1 2013-03-28 13:54:24Z hanne $
 *
 */


/*!
  \file compatibility.h
  \brief Include defines to ensure backwards-compatbility with version of GMRFLib less than 2.0.

  This file defines `old names' to the `new name' to ensure backwards-compatbility with version of
  GMRFLib less or equal to 2.0-alpha. The mapping is as follows:

  \verbatim
  #define GMRFLib_KEEP_bchol      GMRFLib_KEEP_chol
  #define GMRFLib_create_lattice  GMRFLib_make_lattice_graph
  #define GMRFLib_create_graph    GMRFLib_make_empty_graph
  #define GMRFLib_create_constr   GMRFLib_make_empty_constr

  #define GMRFLib_hidden_approx   GMRFLib_init_problem_hidden
  #define GMRFLib_blockupdate2    GMRFLib_blockupdate_hidden
  \endverbatim  

  Please use the new names.
*/

                                                                                                    
#ifndef __GMRFLib_COMPATBILITY_H__
#define __GMRFLib_COMPATBILITY_H__
 
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
 
/* 
   versions < 2.0
*/
#define GMRFLib_KEEP_bchol      GMRFLib_KEEP_chol
#define GMRFLib_create_lattice  GMRFLib_make_lattice_graph
#define GMRFLib_create_graph    GMRFLib_make_empty_graph
#define GMRFLib_create_constr   GMRFLib_make_empty_constr

/*
  version >= 2.0
*/
#define GMRFLib_hidden_approx GMRFLib_init_problem_hidden
#define GMRFLib_blockupdate2  GMRFLib_blockupdate_hidden

__END_DECLS
#endif
