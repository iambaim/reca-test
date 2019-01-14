/* GMRFLibP.h
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
 * RCSId: $Id: GMRFLibP.h 1 2013-03-28 13:54:24Z hanne $
 *
 */

/*!
  \file GMRFLibP.h
  \brief Internal include-file for the GMRFLib source.
*/
  
#ifndef __GMRFLibP_H__
#define __GMRFLibP_H__

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


#include <float.h>
#include <math.h>
#include <assert.h>
#include <stddef.h>
#if !defined(__FreeBSD__)
#include <malloc.h> 
#endif
#include <stdlib.h>
#include <stdio.h>

#define ABS(x)    ((x) >= 0.0 ? (x) : (-(x)))
#define Calloc(n,type) ((n) ? (type *)calloc((size_t)(n),sizeof(type)) : (type *)NULL)
#define FIXME(msg) printf("\n%s:%1d:%s: FIXME [%s]\n",  __FILE__, __LINE__, __GMRFLib_FuncName,(msg?msg:""))
#define FIXME1(msg) if (1) { static int first=1; if (first) { first=0; FIXME(msg); }}
#define FREE(ptr) {if(ptr)free((void *)(ptr));ptr=NULL;}
#define ISINF(x) ((x) > DBL_MAX ? 1:0)
#define ISZERO(x) ((x) == 0.0)
#define LEGAL(i,n) ((i) >= 0 && (i) < (n))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MEMCHK(ptr) GMRFLib_ASSERT(ptr,GMRFLib_EMEMORY)
#define MEMCHK_RETVAL(ptr,retval) GMRFLib_ASSERT_RETVAL(ptr,GMRFLib_EMEMORY,retval)
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MOD(i,n) (((i)+(n))%(n))
#define Malloc(n,type) ((n) ? (type *)malloc((size_t)(n)*sizeof(type)) : (type *)NULL)
#define P(x) printf("line[%1d] " #x "=[%.10f]\n",__LINE__,(double)x)
#define PP(msg,pt) printf("%d: %s ptr " #pt " = 0x%x\n",__LINE__,msg,pt)
#define SQR(x)    ((x)*(x))

/* from /usr/include/assert.h. use __GMRFLib_FuncName to define name of current function.

   Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__' which contains the
   name of the function currently being defined.  This is broken in G++ before version 2.6.  C9x has
   a similar variable called __func__, but prefer the GCC one since it demangles C++ function names.
*/

#ifndef __GNUC_PRERQ
#  if defined __GNUC__ && defined __GNUC_MINOR__
#   define __GNUC_PREREQ(maj, min) ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#  else
#   define __GNUC_PREREQ(maj, min) 0
#  endif
#endif

#if defined __GNUC__
# if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#    define  __GMRFLib_FuncName   __PRETTY_FUNCTION__
#  else
#    if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#      define __GMRFLib_FuncName  __func__
#    else
#      define __GMRFLib_FuncName  ((const char *) "(function-name unavailable)")
#    endif
#  endif
#else
#  if defined __func__
#    define __GMRFLib_FuncName __func__
#  else
#    define __GMRFLib_FuncName ((const char *) "(function-name unavailable)")
#  endif
#endif


/* 
   parts taken from /usr/include/tcl.h
 */
#ifdef __STRING
#  define __GMRFLib_STRINGIFY(x) __STRING(x)
#else
#  ifdef _MSC_VER
#    define __GMRFLib_STRINGIFY(x) #x
# else
#    ifdef RESOURCE_INCLUDED
#      define __GMRFLib_STRINGIFY(x) #x
#    else
#      ifdef __STDC__
#        define __GMRFLib_STRINGIFY(x) #x
#      else
#        define __GMRFLib_STRINGIFY(x) "x"
#      endif
#    endif
#  endif
#endif


/* 
   should not be needed, but the compiler complain...
 */
#ifndef M_E
# define M_E            2.7182818284590452354   /* e */
#endif
#ifndef M_LOG2E
# define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif
#ifndef M_LOG10E
# define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif
#ifndef M_LN2
# define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif
#ifndef M_LN10
# define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif
#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
# define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_PI_4
# define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif
#ifndef M_1_PI
# define M_1_PI         0.31830988618379067154  /* 1/pi */
#endif
#ifndef M_2_PI
# define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif
#ifndef M_2_SQRTPI
# define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif
#ifndef M_SQRT2
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif
#ifndef M_SQRT1_2
# define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

#ifdef isnan
#  define GMRFLib_ISNAN(a) isnan(a)
#else
#  if defined(NAN)
#    define GMRFLib_ISNAN(a) (((a)==NAN)?1:0)
#  else
#    define GMRFLib_ISNAN(a) 0
#  endif
#endif
#ifdef isinf
#  define GMRFLib_ISINF(a) isinf(a)
#else
#  if defined(INFINITY)
#    define GMRFLib_ISINF(a) (((a)==INFINITY)?1:0)
#  else
#    define GMRFLib_ISINF(a) 0
#  endif
#endif

__END_DECLS
#endif
