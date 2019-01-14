/* timer.c
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
  \file timer.c
  \brief Functions to measure (CPU) time and display statistics on the time spent in various
  routines and number times they are called.
 
  The normal usage is to use \c GMRFLib_ctime() to get the CPU time used since a fixed reference,
  and \c GMRFLib_timer_full_report() to display statistics of the most computational demanding
  routines in \c GMRFLib (if  \c GMRFLib_collect_timer_statistics = \c TRUE).
 
  \sa GMRFLib_collect_timer_statistics
*/
                                                                                                                  
                                                                                                                  
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
 
#include "GMRFLib.h"
#include "GMRFLibP.h"
#include "hashP.h"
 
 
static const char RCSId[] = "$Id: timer.c 1 2013-03-28 13:54:24Z hanne $";
 
static map_strvp GMRFLib_timer_hashtable;
static int G_first = 1;
static int G_debug = 0;
 
/* 
   choose default timer according to arch
*/
#if defined(__linux__) || defined(__sun__) || defined(__FreeBSD__)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
double GMRFLib_ctime_default(void)
{
  struct rusage a;
  getrusage(RUSAGE_SELF,&a);
  return (double)
    (double) ((a.ru_utime).tv_sec +(a.ru_stime).tv_sec ) +
    (double) ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec) * 1.0e-6;
}
#else
double GMRFLib_ctime_default(void)
{
    double taucs_ctime();
    return taucs_ctime();
}
#endif
                                                                                                                  
/*!
  \brief Return the CPU time used since a fixed reference.
                                                                                                                  
  This function use the function defined in the global variable
  \c GMRFLib_get_ctime(), and if this
  is \c NULL, then the function \c clock() is used.
                                                                                                                  
  \sa GMRFLib_get_ctime
*/
double GMRFLib_ctime(void)
{
    /*
       return CPU used since any fixed reference
    */
    if (GMRFLib_get_ctime)
        return (*GMRFLib_get_ctime)();
    else
        return (double) (clock()/CLOCKS_PER_SEC);
}
                                                                                                                  
/*
  here are functions for storing and reporting timing information
*/
int GMRFLib_timer_compare(const void *a,  const void *b)
{
    GMRFLib_timer_hashval_tp *aa, *bb;
                                                                                                                  
    aa = (GMRFLib_timer_hashval_tp *) (((map_strvp_element *)a)->value);
    bb = (GMRFLib_timer_hashval_tp *) (((map_strvp_element *)b)->value);
                                                                                                                  
                                                                                                                  
    /*
       sort by name if both are empty
    */
    if (!aa->ntimes && !bb->ntimes) return strcmp(aa->name, bb->name);
                                                                                                                  
    /*
       othewise, by average time spent
    */
    if (aa->ntimes && !bb->ntimes) return -1;
    if (bb->ntimes && !aa->ntimes) return 1;
    return (aa->ctime_acc/aa->ntimes >  bb->ctime_acc/bb->ntimes ? -1 : 1);
}
int GMRFLib_timer_enter(const char *name)
{
    GMRFLib_timer_hashval_tp *p;
    void **vpp;
                                                                                                                  
    if (G_first)
    {
        G_first = 0;
        map_strvp_init_hint(&GMRFLib_timer_hashtable, 30); /* about the number of elmements in the hash-table */
    }
                                                                                                                  
    if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable, (char *)name)))
    {
        p = (GMRFLib_timer_hashval_tp *) *vpp;
    }
    else
    {
        p       = Calloc(1, GMRFLib_timer_hashval_tp);
        p->name = strdup(name);
        map_strvp_set(&GMRFLib_timer_hashtable, (char *)name, (void *)p);
    }
    
    /* 
       if ctime_ref > 0.0, then this routine is already initialized. in this case, we keep the first.
    */
    if (p->ctime_ref <= 0.0) p->ctime_ref = GMRFLib_ctime();
    
    if (G_debug) printf("%s: name %s ctime_ref %f\n", __GMRFLib_FuncName, name, p->ctime_ref);
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_timer_leave(const char *name)
{
    GMRFLib_timer_hashval_tp *p;
    void **vpp;
    double used;
                                                                                                                  
    if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable, (char *)name)))
        p = (GMRFLib_timer_hashval_tp *) *vpp;
    else
        GMRFLib_ASSERT(vpp != NULL, GMRFLib_ESNH);

    if (p->ctime_ref < 0.0)
    {
	/* 
	   this is an ``illegal instruction''. _timer_leave is called without a corresponding call
	   to _timer_enter. this happens, with purpose, with some of the `_intern' routines.
	*/
	return GMRFLib_SUCCESS;
    }
    used           = GMRFLib_ctime() - p->ctime_ref;
    used           = MAX(0.0, used);              /* yes */
    p->ctime_acc  += used;
    p->ctime_acc2 += SQR(used);
    if (p->ntimes)
    {
        p->ctime_min = MIN(p->ctime_min, used);
        p->ctime_max = MAX(p->ctime_max, used);
    }
    else
        p->ctime_min = p->ctime_max = used;
                                                                                                                  
    p->ctime_ref = -1.0;			  /* flag it specially */
    p->ntimes++;
                                                                                                                  
    if (G_debug) printf("%s: name %s ctime_ref %f used %f\n", __GMRFLib_FuncName, name, p->ctime_ref, used);
    
    return GMRFLib_SUCCESS;
}
const char *GMRFLib_timer_strip_store(const char *name)
{
    /* 
       if name ends with `_intern', then remove it from name.

       this routine returns a ptr to a static storage containing the modified name, use 'strdup'
       alternatively.
    */
#define LEN_MAX 1024
    int idx, len_name, len_suffix;
    char *suffix = "_store";
    static char nm[LEN_MAX+1];
    
    if (!name) return name;

    len_name   = strlen(name); GMRFLib_ASSERT_RETVAL(len_name <= LEN_MAX, GMRFLib_ESNH, name);
    len_suffix = strlen(suffix);
    
    strcpy(nm, name);
    idx = len_name - len_suffix;

    if (idx > 0  && (strcmp(&nm[idx],  suffix) == 0)) /* yes, i want '> 0' */
    {
	nm[idx] = '\0';
	return (const char *)nm;
    }
    else
	return name;
#undef LEN_MAX    
}
int GMRFLib_timer_print_entry(FILE *ffp, GMRFLib_timer_hashval_tp *p)
{
    fprintf(ffp, "%-40s %10.6f %6d %10.6f %10.6f %10.6f %10.6f\n",
            p->name, (p->ntimes ? p->ctime_acc/p->ntimes : 0.0),
            (int)p->ntimes, p->ctime_acc,
            (p->ntimes ? sqrt(MAX(0.0, p->ctime_acc2/p->ntimes - SQR(p->ctime_acc/p->ntimes))) : 0.0),
            p->ctime_min, p->ctime_max);
                                                                                                                  
    return GMRFLib_SUCCESS;
}
/*!
  \brief Write statistics collected for the function named \c name to \c fp.
                                                                                                                  
  \param[in] fp Pointer to (an already open) file the report is written to.
  \param[in] name The name of the function for which statistics is to be reported. If \c name =
  \c NULL, then the statistics for all functions are displayed
                                                                                                                  
  \sa GMRFLib_timer_full_report
*/
int GMRFLib_timer_report(FILE *fp, const char *name)
{
                                                                                                                  
    GMRFLib_timer_hashval_tp *p;
    void **vpp;
    FILE *ffp;
    char *sep = "------------------------------------------------------------------------------------------------------";
                                                                                                                  
    if (G_first) return GMRFLib_SUCCESS;
    ffp = (fp ? fp : stdout);
                                                                                                                  
    fprintf(ffp, "\n\nGMRFLib report on time usage\n%-40s %10s %6s %10s %10s %10s %10s\n%s\n",
            "Routine", "Mean", "N", "Total", "Stdev", "Min", "Max", sep);
    if (name)
    {
        if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable, (char *)name)))
        {
            p = (GMRFLib_timer_hashval_tp *) *vpp;
            GMRFLib_timer_print_entry(ffp, p);
        }
    }
    else
    {
        map_strvp_element *all;
        long count, i;
                                                                                                                  
        map_strvp_getall(&GMRFLib_timer_hashtable, &all, &count);
                                                                                                                  
        qsort(all, count, sizeof(map_strvp_element), GMRFLib_timer_compare);
                                                                                                                  
        for(i=0;i<count;i++) GMRFLib_timer_print_entry(ffp, (GMRFLib_timer_hashval_tp *) (all[i].value));
        FREE(all);
    }
    fprintf(ffp, "%s\n", sep);
                                                                                                                  
    return GMRFLib_SUCCESS;
}
/*!
  \brief Write  statistics collected for all functions to \c fp.
                                                                                                                  
  This is the routine to use to print internal statistics of number of times each routine is called
  and its CPU usage.
                                                                                                                  
  \remark This function is simply a wrapper for \c GMRFLib_timer_report(fp,NULL).
                                                                                                                  
  \sa  GMRFLib_collect_timer_statistics
*/
int GMRFLib_timer_full_report(FILE *fp)
{
    if (G_first) return GMRFLib_SUCCESS;
    return GMRFLib_timer_report(fp, NULL);
}

