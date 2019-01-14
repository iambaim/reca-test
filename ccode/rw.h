/* rw.h
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
 * RCSId: $Id: rw.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file rw.h
  \brief Typedefs used to define RW1, RW2, CRW1 and CRW2 models.
*/

#ifndef __GMRFLib_RW_H__
#define __GMRFLib_RW_H__

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
  \struct GMRFLib_rwdef_tp rw.h
  \brief The typedef used to define RW1 and RW2 models with regular locations.

  This defines the RW1 and RW2 models with regular locations, which is an argument to \ref
  GMRFLib_rw(), to compute an entry in the presicion matrix.  \ref GMRFLib_make_rw_graph() can be
  used to create the appropirate graph.

  \sa \ref GMRFLib_crwdef_tp, \ref GMRFLib_rw(),  \ref GMRFLib_make_rw_graph()
*/
typedef struct
{
    /*! 
      \brief The size of the graph or the number of locations
    */
    int n;					

    /*!
      \brief The order or the random walk. Must be 1 or 2.
    */
    int order;			

    /*!
      \brief A flag for using cyclic boundary conditions or not.

      Use cyclic boundary conditions if cyclic=\c TRUE, otherwise do not use cyclic boundary
      conditions.
    */
    int cyclic;	

    /*!
      \brief A possible pointer to precision.

      If prec!=\c NULL, then this pointer points to the precision of the RW, otherwise, the precision
      is assumed to be equal to 1
     */
    double *prec;			
}
GMRFLib_rwdef_tp;

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], ..., z[n-1].
*/
#define GMRFLib_CRW_LAYOUT_SIMPLE 0

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], z'[0], z[1], z'[1]....
*/
#define GMRFLib_CRW_LAYOUT_PAIRS  1

/*!
  \brief Define the layout for the CRW2 model: Layout z[0], ..., z[n-1], z'[0], ..., z'[n-1]
*/
#define GMRFLib_CRW_LAYOUT_BLOCK  2


/*!
  \struct GMRFLib_crwdef_tp rw.h

  \brief The typedef used to define (C)RW1 and CRW2 models with irregular locations.

  This defines the (C)RW1, CRW2 and approximate CRW2 models with irregular locations, which is an
  argument to to \ref GMRFLib_crw(), to compute an entry in the presicion matrix and \ref
  GMRFLib_make_crw_graph() which create the appropirate graph.

  \sa GMRFLib_rwdef_tp, GMRFLib_crw, GMRFLib_make_crw_graph
*/
typedef struct
{
    /*! 
      \brief The size of the graph or the number of locations
    */
    int n;					

    /*!
      \brief The order or the random walk. Must be 1 or 2.
    */
    int order;			

    /*!
      \brief A possible pointer to precision.

      If prec!=\c NULL, then this pointer points to the precision of the CRW, otherwise, the precision
      is assumed to be equal to 1
     */
    double *prec;			

    /*! 
       \brief An array with the positions for x[0]...x[n-1]

       \c position[i] is the position to x[i], for i=0...n-1.  <em> It is assumed, and not checked
       for, that the positions are increasing </em>, i.e. position[0] < position[1] < ... <
       position[n-1].  If position is \c NULL, then the positions are assumed to be 0, 1, ..., n-1.
    */
    double *position;
    
    /*!
      \brief Define the memory layout for the CRW2 model.
      
      This argument defines the memory layout for the CRW2 model. If it is equal to \c
      #GMRFLib_CRW_LAYOUT_SIMPLE, then the approximative model is used where no augmentation with
      velocities are needed. If it is equal to \c #GMRFLib_CRW_LAYOUT_BLOCK or \c
      #GMRFLib_CRW_LAYOUT_PAIRS, then the exact augmented solution is used.

      The exact solution of the CRW2 model augment the state-vector z with its derivatives z' (or
      velocities). Both vectors are stored in the general array x which is of length 2*n. The term
      ij in the presicion matrix, and also the graph itself, depends then on how z and z' are stored
      in x.

      If \c layout = #GMRFLib_CRW_LAYOUT_PAIRS then <em>x = (z[0], z'[0], z[1], z'[1], ...., z[n-1],
      z'[n-1])</em>.

      If \c layout = #GMRFLib_CRW_LAYOUT_BLOCK then <em>x = (z[0], z[1], ...., z[n-1], z'[0], z'[1],
      ..., z'[n-1])</em>.

      If \c layout = #GMRFLib_CRW_LAYOUT_SIMPLE then <em>x = (z[0], z[1], ...., z[n-1])</em>.

    */
    int layout;	


    /*!
    \brief Work array for the \c GMRFLib_crw() function.

    The array \c work is a working array for the \c GMRFLib_crw() function. It must be initialised
    to \c NULL, and then \c GMRFLib_rw() allocate the space needed and calculate the contents. This
    array can be free'ed using \c free().

    The workarray depends ONLY of the positions, and does NOT depend of the precision, order or the
    layout.

    */
    double *work;
}
GMRFLib_crwdef_tp;

double GMRFLib_rw( int node, int nnode, GMRFLib_rwdef_tp  *rwdef);
double GMRFLib_crw(int node, int nnode, GMRFLib_crwdef_tp *crwdef);

__END_DECLS
#endif
