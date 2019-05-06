/*!
  \file random.c
  \brief RNG functions based on the MT19937 random generator. GMRFLib use this library by default.

  The routines in random.c are most easily available to the user as the following function pointers:

  - GMRFLib_uniform() used to return a Uniform(0,1) variable.
  - GMRFLib_uniform_init() to initialize the RNG
  - GMRFLib_uniform_getstate() which returns in a malloc'ed array the state
    of the RNG
  - GMRFLib_uniform_setstate() which set the state in the RNG.

  \sa GMRFLib_uniform, GMRFLib_uniform_init, GMRFLib_uniform_getstate, GMRFLib_uniform_setstate
  
  The details are as follows:
  
  A C-program for MT19937: Real number version([0,1)-interval) (1998/4/6) srandom() generates one
  pseudorandom real number (double) which is uniformly distributed on [0,1)-interval, for each
  call. srandom(seed) set initial values to the working area of 624 words. Before ran(),
  srandom(seed) must be called once. (seed is any 32-bit integer except for 0).  Integer generator
  is obtained by modifying two lines.  Coded by Takuji Nishimura, considering the suggestions by
  Topher Cooper and Marc Rieffel in July-Aug. 1997.

  Copyright for the MT19937 code [which is modified by hrue@math.ntnu.no]:
  \verbatim
  This library is free software; you can redistribute it and/or modify it under the terms of the GNU
  Library General Public License as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.  This library is distributed in the hope that
  it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public License for more details.
  You should have received a copy of the GNU Library General Public License along with this library;
  if not, write to the Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA * 02111-1307
  USA

  Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When you use this, send an email to:
  matumoto@math.keio.ac.jp with an appropriate reference to your work.

  REFERENCE M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally Equidistributed
  Uniform Pseudo-Random Number Generator", ACM Transactions on Modeling and Computer Simulation,
  Vol. 8, No. 1, January 1998, pp 3--30.
\endverbatim

  Example of usage:
\verbatim
  int main()
  { 
    int i, j, seed;
    char *state;
    
    printf("seed?\n");
    fscanf(stdin, "%d", &seed);
    printf("seed %d\n", seed);
    ran_init(seed);
    for (j=0; j<10; j++) 
        printf("%10.8f\n", ran());
    printf("\n");

    ran_init(seed);
    for (j=0; j<10; j++) 
        printf("%10.8f\n", ran());


    printf("get state\n");
    state = ran_getstate();
    for(i=0;i<2*N;i++) if (i==0 || i == 2*N-1) printf(" %.12f", ran());
    printf("\n");
    printf("set state\n");
    ran_setstate(state); free(state);
    for(i=0;i<2*N;i++) if (i==0 || i == 2*N-1) printf(" %.12f", ran());
    printf("\n");

    return 0;
  }
\endverbatim

*/

#include <string.h>
#include <stdio.h>
#include <time.h>
#if !defined(__FreeBSD__)

#endif
#include <stdlib.h>
#include <stddef.h>

#ifdef MY_MALLOC
#include "my_malloc.h"
#endif
#ifndef RLIB__
#include "random.h"
#endif

/*
  Period parameters
*/  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df			  /* constant vector a */
#define UPPER_MASK 0x80000000			  /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff			  /* least significant r bits */

/*
  Tempering parameters
*/   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned int mt[N];			  /* the array for the state vector  */
static int mti=N+1;				  /* mti==N+1 means mt[N] is not initialized */
static double ran_count = 0.0;			  /* count each call to ran() */

char *ran_getstate_expert(size_t *len)
{
    char *state;

    *len = sizeof(unsigned int)*N+sizeof(int);
    state = (char *)calloc((unsigned) *len, sizeof(char));
    memcpy(state, &mti, sizeof(int));
    memcpy(&state[sizeof(int)], mt, *len -sizeof(int));

    return state;
}
/*!
  \brief Return the state of the RNG in a malloced char array.
  \sa ran_setstate
*/
char *ran_getstate(void)
{
    size_t len;
    return ran_getstate_expert(&len);
}
/*!
  \brief Set the state of the RNG
  \sa ran_getstate
*/
int ran_setstate(char *state)
{
    size_t len = sizeof(unsigned int)*N+sizeof(int);
 
    memcpy(&mti, state, sizeof(int));
    memcpy(mt, &state[sizeof(int)], len-sizeof(int));

    return 0;
}
/*!
  \brief Initialize the RNG
*/
int ran_init(unsigned int seed)
{
    /*
      setting initial seeds to mt[N] using the generator Line 25 of Table 1 in [KNUTH 1981, The Art
      of Computer Programming Vol. 2 (2nd Ed.), pp102]
    */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
    return 0;
}
int ran_init_(unsigned int seed)
{
    return ran_init(seed);
}
/*!
  \brief Return a Uniform(0,1) random number
*/
double ran(void)
{
    unsigned int y;
    static unsigned int mag01[2]={0x0, MATRIX_A};

    ran_count++;				  /* global counter */

    if (mti >= N)
    {
	/* 
	   generate N words at one time
	*/
        int kk;

        if (mti == N+1)   /* if ran_init() has not been called, */
            ran_init(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++)
	{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++)
	{
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    /* 
       double: [0,1)-interva
       return y for integer generation
     */
    return ((double)y * 2.3283064365386963e-10 );
}
double ran_(void)
{
    return ran();
}
double ran_counter(void)
{
    return ran_count;
}
#ifdef TEST
int main()
{ 
    int i, j, seed;
    char *state;
    
    printf("seed?\n");
    fscanf(stdin, "%d", &seed);
    printf("seed %d\n", seed);
    ran_init(seed);
    for (j=0; j<10; j++) 
        printf("%10.8f\n", ran());
    printf("\n");

    ran_init(seed);
    for (j=0; j<10; j++) 
        printf("%10.8f\n", ran());


    /* 
       test 'state'
     */

    printf("get state\n");
    state = ran_getstate();
    for(i=0;i<2*N;i++) if (i==0 || i == 2*N-1) printf(" %.12f", ran());
    printf("\n");
    printf("set state\n");
    ran_setstate(state); free(state);
    for(i=0;i<2*N;i++) if (i==0 || i == 2*N-1) printf(" %.12f", ran());
    printf("\n");

    return 0;
}
#endif
