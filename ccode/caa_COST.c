/*!
  \file caa_COST.c
  \brief Routines for COST project
*/
#include <R.h>
#include <Rdefines.h>
#include "caa.h"
#include "caa_Rextensions.h"
#include "caa_read_write.h"
#include "caa_routines.h"

    

/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing 

  Makes a struct of type Data_orig (see caa.h for definition)

  Space allocated in this routine is reallocated in re_makedata_COST
*/
int makedata_COST(SEXP i_COSTList, Data_orig **o_D_orig, Data_COST **o_D_COST)
{
  Data_orig     *D_orig;
  Data_COST     *D_COST;
  Data_obs      *D_obs;
  Data_mland    *D_mland;
  int            i,f,h,n,s,t,N_sim,N_orig;
  int            l_int,n_trip,n_fish,N_int,nHaul,nSize;
  int            ind,ind_alk,ind_fish,ind_fish_l,ind_haul,ind_size,ind_orig,ind_t;
  long          *lengths,*lengths2;
  double         l,psum;
  double        *P_l,*p,*int_len,*cens_mu,*cens_tau,*cens_pri;
  SEXP           elmt = R_NilValue;

  FILE          *caa_debug;
  #ifdef DEBUG_COST
  caa_debug = fopen("caa_debug_COST.txt","w");
  #endif

  /* Allocating space for COST object */
  D_COST = CALLOC(1,Data_COST);

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_par")))
    D_COST->num_var = INTEGER_VALUE(elmt); 
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "resample.l")))
    D_COST->resample = INTEGER_VALUE(elmt); 


  /* Observer data */
  D_obs = CALLOC(1,Data_obs);

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "n_trip_obs")))
    D_obs->n_trip = INTEGER_VALUE(elmt); // number of trips with observer data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_trip_obs")))
    D_obs->num_trip = INTEGER_POINTER(AS_INTEGER(elmt)); // number of hauls pr trip 

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_haul_disc")))
    D_obs->num_haul_disc = INTEGER_POINTER(AS_INTEGER(elmt)); // number of length-measured discarded fish pr haul

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "season_obs")))
    D_obs->season = INTEGER_POINTER(AS_INTEGER(elmt)); // observed month

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "l_disc")))
    D_obs->l_disc = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length categories for discard samples

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "lfreq_disc")))
    D_obs->lfreq_disc = INTEGER_POINTER(AS_INTEGER(elmt)); // number at length for discards

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "haulsize_disc")))
    D_obs->haulsize_disc = NUMERIC_POINTER(AS_NUMERIC(elmt)); // number of discards in haul

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "sampsize_disc")))
    D_obs->sampsize_disc = NUMERIC_POINTER(AS_NUMERIC(elmt)); // number of discards sampled

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_alk_disc")))
    D_obs->num_alk = INTEGER_POINTER(AS_INTEGER(elmt)); // number of discard age-length data within trip

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_l_disc")))
    D_obs->alk_l = NUMERIC_POINTER(AS_NUMERIC(elmt)); // lengths for discard age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_a_disc")))
    D_obs->alk_a = INTEGER_POINTER(AS_INTEGER(elmt)); // ages for discard age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_lfreq_disc")))
    D_obs->alk_lfreq = INTEGER_POINTER(AS_INTEGER(elmt)); // numbers at length for discard age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_trip_land")))
    D_obs->num_trip_land = INTEGER_POINTER(AS_INTEGER(elmt)); // number of size classes pr trip with landings

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_size_land")))
    D_obs->num_size_land = INTEGER_POINTER(AS_INTEGER(elmt)); // number of measured landed fish pr size class

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "l_land")))
    D_obs->l_land = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length categories for landing samples

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "lfreq_land")))
    D_obs->lfreq_land = INTEGER_POINTER(AS_INTEGER(elmt)); // number at length for landings

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "totsize_land")))
    D_obs->totsize_land = NUMERIC_POINTER(AS_NUMERIC(elmt)); // total weight landed in size class

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "sampsize_land")))
    D_obs->sampsize_land = NUMERIC_POINTER(AS_NUMERIC(elmt)); // weight of landings sampled for lengths in size class

  /* Market landing data */
  D_mland = CALLOC(1,Data_mland);
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "n_trip_mland")))
    D_mland->n_trip = INTEGER_VALUE(elmt); // number of trips with market landing data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_trip_mland")))
    D_mland->num_trip = INTEGER_POINTER(AS_INTEGER(elmt)); // number of size classes pr trip with market landings

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "season_mland")))
    D_mland->season = INTEGER_POINTER(AS_INTEGER(elmt)); // observed month

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_alk_mland")))
    D_mland->num_alk = INTEGER_POINTER(AS_INTEGER(elmt)); // number of market landing age-length data within trip

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_l_mland")))
    D_mland->alk_l = NUMERIC_POINTER(AS_NUMERIC(elmt)); // lengths for market landing age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_a_mland")))
    D_mland->alk_a = INTEGER_POINTER(AS_INTEGER(elmt)); // ages for market landing age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "alk_lfreq_mland")))
    D_mland->alk_lfreq = INTEGER_POINTER(AS_INTEGER(elmt)); // numbers at length for market landing age-length data

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_size_mland")))
    D_mland->num_size = INTEGER_POINTER(AS_INTEGER(elmt)); // number of measured market landing fish pr size class

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "l_mland")))
    D_mland->l = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length categories for market landing samples

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "lfreq_mland")))
    D_mland->lfreq = INTEGER_POINTER(AS_INTEGER(elmt)); // number at length for market landings

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "totsize_mland")))
    D_mland->totsize = NUMERIC_POINTER(AS_NUMERIC(elmt)); // total weight for market landing in size class

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "sampsize_mland")))
    D_mland->sampsize = NUMERIC_POINTER(AS_NUMERIC(elmt)); // weight of market landings sampled for lengths in size class

  /* Allocating space for censoring parameters */
  D_COST->cens = CALLOC(1,cens_struct);
  D_COST->cens->ncat = D_obs->n_trip+D_mland->n_trip;
  D_COST->cens->r = CALLOC(D_COST->cens->ncat,double);
  D_COST->cens->mu = CALLOC(3,double);
  D_COST->cens->tau = CALLOC(3,double);
  cens_mu = NUMERIC_POINTER(getListElement(i_COSTList, "cens.mu"));
  cens_tau = NUMERIC_POINTER(getListElement(i_COSTList, "cens.tau"));
  cens_pri = NUMERIC_POINTER(getListElement(i_COSTList, "cens.pri"));
  for(i=0;i<3;i++)
    D_COST->cens->mu[i] = cens_mu[i];
  for(i=0;i<3;i++)
    D_COST->cens->tau[i] = cens_tau[i];
  D_COST->cens->mu_prior_mean = cens_pri[0];
  D_COST->cens->mu_prior_prec = cens_pri[1];
  D_COST->cens->a_prior = cens_pri[2];
  D_COST->cens->b_prior = cens_pri[3];


  /* Allocating space for 'original' parameters */

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "num_fish")))
    n_fish = INTEGER_VALUE(elmt);

  n_trip = D_obs->n_trip+D_mland->n_trip;
  D_orig = CALLOC(1,Data_orig);
  D_orig->nFishBoat = CALLOC(n_trip,int); // Free ok
  D_orig->totage = CALLOC(n_fish,int);  // Free ok 
  D_orig->totlength = CALLOC(n_fish,double); // Free ok
  D_orig->totweight = CALLOC(n_fish,double); // Free ok
  D_orig->replength = CALLOC(n_fish,int);  // Free ok
  D_orig->discard = CALLOC(n_fish,int);  // Free ok
  D_orig->landed = CALLOC(n_fish,int);  // Free ok
  D_orig->start_noAge = CALLOC(n_trip,int); // Free ok
  D_orig->start_Age = CALLOC(n_trip,int); // Free ok
  D_orig->num_noAge = CALLOC(n_trip,int);  // Free ok
  D_orig->haulweight = CALLOC(n_trip,double); // Free ok
  D_orig->season = CALLOC(n_trip,int);   // Free ok
  D_orig->n_discard = CALLOC(n_trip,int);   // Free ok
  D_orig->n_landed = CALLOC(n_trip,int);   // Free ok
  D_orig->coastal_cod = 0;
  D_orig->nHaul = n_trip;
  D_orig->nFish = n_fish;

  if(!Rf_isNull(elmt = getListElement(i_COSTList, "n_int_len")))
    D_orig->n_int_len = INTEGER_VALUE(elmt); // number of intervals for length
  N_int = D_orig->n_int_len;
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "int_len_lim")))
    D_orig->int_len_lim = NUMERIC_POINTER(AS_NUMERIC(elmt)); // lower limits of length-intervals
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "int_len_vec")))
    D_orig->int_len = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length value for intervals

  for(i=0;i<n_fish;i++)
    D_orig->totweight[i] = 0.0;

  lengths = CALLOC(N_int,long);      // Free ok
  P_l = CALLOC(N_int,double);      // Free ok
  p = CALLOC(2,double);      // Free ok
  lengths2 = CALLOC(2,long);      // Free ok


  printf("\nStart simulate total lengths for observer data\n");

  /* Simulate total lengths for observer data */
  ind_fish = 0;
  ind_fish_l = 0;
  ind_haul = 0;
  ind_size = 0;
  ind_alk = 0;
  ind_orig = 0;
  ind = 0;
  for(t=0;t<D_obs->n_trip;t++)
    {
       /* Discard data */
      D_orig->start_noAge[t] = ind_orig + D_obs->num_alk[t];
      D_orig->start_Age[t] = ind_orig;
      D_orig->num_noAge[t] = N_int;
      D_orig->nFishBoat[t] = D_obs->num_alk[t]+N_int;
      D_orig->season[t] = D_obs->season[t];
      D_orig->n_discard[t] = 0;
      D_orig->n_landed[t] = 0;
      N_sim = 0;
      N_orig = 0;
      ind_orig = D_orig->start_noAge[t];
      for(f=0;f<N_int;f++)
	{
	  D_orig->totage[ind_orig] = -99999;
	  D_orig->totlength[ind_orig] = D_orig->int_len[f];
	  D_orig->replength[ind_orig] = 0;
	  D_orig->discard[ind_orig] = 0;
	  D_orig->landed[ind_orig] = 0;
	  ind_orig++;
	}
      ind_orig = D_orig->start_noAge[t];
      for(h=0;h<D_obs->num_trip[t];h++)
	{
	  if(D_obs->num_haul_disc[ind_haul]>0)
	    {
	      nHaul = 0;
	      for(i=0;i<N_int;i++)
		P_l[i] = 0.0;
	      for(f=0;f<D_obs->num_haul_disc[ind_haul];f++)
		{
		  l = D_obs->l_disc[ind_fish];
		  l_int = 0;
		  while(l > D_orig->int_len_lim[l_int])
		    l_int++;
		  P_l[l_int] += D_obs->lfreq_disc[ind_fish];
		  D_orig->replength[ind_orig+l_int] += D_obs->lfreq_disc[ind_fish];
		  D_orig->discard[ind_orig+l_int] += D_obs->lfreq_disc[ind_fish];
		  nHaul += D_obs->lfreq_disc[ind_fish];
		  ind_fish++;
		}
	      // convert to probabilities
	      for(i=0;i<N_int;i++)
		P_l[i] /= nHaul;
	      // number of fish to be simulated
	      if(nHaul==0)
		n=0;
	      else
		n = (int) nHaul*(D_obs->haulsize_disc[ind_haul]/D_obs->sampsize_disc[ind_haul]-1);
	      N_sim += nHaul+n;
	      N_orig += nHaul;
	      my_genmul(n,P_l,N_int,lengths);
	      for(i=0;i<N_int;i++)
		{
		  D_orig->replength[ind_orig+i] += (int) lengths[i];
		  D_orig->discard[ind_orig+i] += (int) lengths[i];
		}
	    }
	  ind_haul++;
	}
      /* Landing data */
      for(s=0;s<D_obs->num_trip_land[t];s++)
	{
	  //	  if(D_obs->num_size_land[ind_size]==0)
	  nSize = 0;
	  for(i=0;i<N_int;i++)
	    P_l[i] = 0.0;
	  for(f=0;f<D_obs->num_size_land[ind_size];f++)
	    {
	      l = D_obs->l_land[ind_fish_l];
	      l_int = 0;
	      while(l > D_orig->int_len_lim[l_int])
		l_int++;
	      P_l[l_int] += D_obs->lfreq_land[ind_fish_l];
	      D_orig->replength[ind_orig+l_int] += D_obs->lfreq_land[ind_fish_l];
	      D_orig->landed[ind_orig+l_int] += D_obs->lfreq_land[ind_fish_l];
	      nSize += D_obs->lfreq_land[ind_fish_l];
	      ind_fish_l++;
	    }
	  // convert to probabilities
	  for(i=0;i<N_int;i++)
	    P_l[i] /= nSize;
	  // number of fish to be simulated
	  n = nSize*(D_obs->totsize_land[ind_size]/D_obs->sampsize_land[ind_size]-1);
	  N_sim += nSize+n;
	  N_orig += nSize;
          my_genmul(n,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      D_orig->replength[ind_orig+i] += (int) lengths[i];
	      D_orig->landed[ind_orig+i] += (int) lengths[i];
	    }
	  ind_size++; 
	}
      
      //scale to original size
      if(N_orig>0)
	{
	  psum = 0;
	  for(i=0;i<N_int;i++)
	    {
	      P_l[i] = D_orig->replength[ind_orig+i]*N_orig/N_sim;
	      psum += P_l[i];
	    }
	  for(i=0;i<N_int;i++)
	    P_l[i] /= psum;
	  my_genmul(N_orig,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      if(lengths[i]==0)
		{
		  D_orig->replength[ind_orig+i] = 0;
		  D_orig->discard[ind_orig+i] = 0;
		  D_orig->landed[ind_orig+i] = 0;
		}
	      else
		{
		  p[0] = (double) D_orig->discard[ind_orig+i]/(D_orig->discard[ind_orig+i]+D_orig->landed[ind_orig+i]);
		  p[1] = (double) D_orig->landed[ind_orig+i]/(D_orig->discard[ind_orig+i]+D_orig->landed[ind_orig+i]);
		  D_orig->replength[ind_orig+i] = (int) lengths[i];
		  my_genmul(D_orig->replength[ind_orig+i],p,2,lengths2);
		  D_orig->discard[ind_orig+i] = (int) lengths2[0];
		  D_orig->landed[ind_orig+i] = (int) lengths2[1];
		}
	    }
	}
      // put the age-length data into D_orig object
      n = 0;
      for(f=0;f<D_obs->num_alk[t];f++)
	{
	  D_orig->totage[ind] = D_obs->alk_a[ind_alk];
	  D_orig->totlength[ind] = D_obs->alk_l[ind_alk];
	  D_orig->replength[ind] = D_obs->alk_lfreq[ind_alk];
	  D_orig->discard[ind] = D_obs->alk_lfreq[ind_alk];
	  // remove length count for lengths with missing ages
	  l_int = 0;
	  while(D_obs->alk_l[ind_alk] > D_orig->int_len_lim[l_int])
	    l_int++;
	  D_orig->replength[ind_orig+l_int] -= D_obs->alk_lfreq[ind_alk];
	  D_orig->discard[ind_orig+l_int] -= D_obs->alk_lfreq[ind_alk];
	  if(D_orig->discard[ind_orig+l_int]<0) // remove age data from sample
	    {
	      n += 1;
	      D_orig->replength[ind] = D_obs->alk_lfreq[ind_alk]-1;
	      D_orig->discard[ind] = D_obs->alk_lfreq[ind_alk]-1;
	      D_orig->replength[ind_orig+l_int] = D_orig->landed[ind_orig+l_int];
	      D_orig->discard[ind_orig+l_int] = 0;
	    }
	  ind_alk++;
	  ind++;
	}
      printf("trip=%d,number of age-length data removed=%d\n",t,n);
      for(i=D_orig->start_noAge[t];i<D_orig->start_noAge[t]+D_orig->num_noAge[t];i++)
	{
	  D_orig->n_landed[t] += D_orig->landed[i];
	  D_orig->n_discard[t] += D_orig->discard[i];
	}
      ind += N_int;
      ind_orig += N_int;
    }

  #ifdef DEBUG_COST
  n=0;
  for(t=0;t<D_obs->n_trip;t++)
    {
      fprintf(caa_debug,"t=%d,nFishBoat=%d,start_noAge=%d,num_noAge=%d\n",
	      t,D_orig->nFishBoat[t],D_orig->start_noAge[t],D_orig->num_noAge[t]);
      n += D_orig->nFishBoat[t];
    }
  fprintf(caa_debug,"n=%d,totage[i],totlength[i],replength[i]:\n",n);
  n=0;
  for(i=0;i<n_fish;i++)
    {
      fprintf(caa_debug,"i=%d,%d,%f,%d\n",i,D_orig->totage[i],
	      exp(D_orig->totlength[i]),D_orig->replength[i]);
      n += D_orig->replength[i];
    }
  fprintf(caa_debug,"n=%d\n",n);
  #endif

  printf("\nStart simulate total lengths for market landing data\n");
  ind_fish = 0;
  ind_size = 0;
  ind_alk = 0;
  ind_t = D_obs->n_trip;
  for(t=0;t<D_mland->n_trip;t++)
    {
      D_orig->start_noAge[ind_t] = ind_orig + D_mland->num_alk[t];
      D_orig->start_Age[ind_t] = ind_orig;
      D_orig->num_noAge[ind_t] = N_int;
      D_orig->nFishBoat[ind_t] = D_mland->num_alk[t]+N_int;
      D_orig->season[ind_t] = D_mland->season[t];
      D_orig->n_discard[ind_t] = 0;
      D_orig->n_landed[ind_t] = 0;
      N_sim = 0;
      N_orig = 0;
      ind_orig = D_orig->start_noAge[ind_t];
      for(f=0;f<N_int;f++)
	{
	  D_orig->totage[ind_orig] = -99999;
	  D_orig->totlength[ind_orig] = D_orig->int_len[f];
	  D_orig->replength[ind_orig] = 0;
	  D_orig->discard[ind_orig] = 0;
	  D_orig->landed[ind_orig] = 0;
	  ind_orig++;
	}
      ind_orig = D_orig->start_noAge[ind_t];
      for(s=0;s<D_mland->num_trip[t];s++)
	{
	  nSize = 0;
	  for(i=0;i<N_int;i++)
	    P_l[i] = 0.0;
	  for(f=0;f<D_mland->num_size[ind_size];f++)
	    {
	      l = D_mland->l[ind_fish];
	      l_int = 0;
	      while(l > D_orig->int_len_lim[l_int])
	      	l_int++;
	      P_l[l_int] += D_mland->lfreq[ind_fish];
	      D_orig->replength[ind_orig+l_int] += D_mland->lfreq[ind_fish];
	      D_orig->landed[ind_orig+l_int] += D_mland->lfreq[ind_fish];
	      nSize += D_mland->lfreq[ind_fish];
	      ind_fish++;
	    }
	  // convert to probabilities
	  for(i=0;i<N_int;i++)
	    P_l[i] /= nSize;
	  // number of fish to be simulated
	  n = nSize*(D_mland->totsize[ind_size]/D_mland->sampsize[ind_size]-1);
	  N_sim += nSize+n;
	  N_orig += nSize;
          my_genmul(n,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      D_orig->replength[ind_orig+i] += (int) lengths[i];
	      D_orig->landed[ind_orig+i] += (int) lengths[i];
	    }
	  ind_size++; 
	}
      //scale to original size
      if(N_orig>0)
	{
	  psum = 0;
	  for(i=0;i<N_int;i++)
	    {
	      P_l[i] = D_orig->replength[ind_orig+i]*N_orig/N_sim;
	      psum += P_l[i];
	    }
	  for(i=0;i<N_int;i++)
	    P_l[i] /= psum;
	  my_genmul(N_orig,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      if(lengths[i]==0)
		{
		  D_orig->replength[ind_orig+i] = 0;
		  D_orig->landed[ind_orig+i] = 0;
		}
	      else
		{
		  D_orig->replength[ind_orig+i] = (int) lengths[i];
		  D_orig->landed[ind_orig+i] = (int) lengths[i];
		}
	    }
	}
      // put the age-length data into D_orig object
      n = 0;
      for(f=0;f<D_mland->num_alk[t];f++)
	{
	  D_orig->totage[ind] = D_mland->alk_a[ind_alk];
	  D_orig->totlength[ind] = D_mland->alk_l[ind_alk];
	  D_orig->replength[ind] = D_mland->alk_lfreq[ind_alk];
	  D_orig->landed[ind] = D_mland->alk_lfreq[ind_alk];
	  // remove length count for lengths with missing ages
	  l_int = 0;
	  while(D_mland->alk_l[ind_alk] > D_orig->int_len_lim[l_int])
	    l_int++;
	  D_orig->replength[ind_orig+l_int] -= D_mland->alk_lfreq[ind_alk];
	  D_orig->landed[ind_orig+l_int] -= D_mland->alk_lfreq[ind_alk];
	  if(D_orig->landed[ind_orig+l_int]<0)// remove age data from sample
	    {
	      n -= D_orig->landed[ind_orig+l_int];
	      D_orig->replength[ind] = D_mland->alk_lfreq[ind_alk]+D_orig->landed[ind_orig+l_int];
	      D_orig->landed[ind] = D_mland->alk_lfreq[ind_alk]+D_orig->landed[ind_orig+l_int];
	      D_orig->replength[ind_orig+l_int] = D_orig->discard[ind_orig+l_int];
	      D_orig->landed[ind_orig+l_int] = 0;
	    }
	  ind_alk++;
	  ind++;
	}
      printf("trip=%d,number of age-length data removed=%d\n",t,n);
      for(i=D_orig->start_noAge[t];i<D_orig->start_noAge[t]+D_orig->num_noAge[t];i++)
	{
	  D_orig->n_landed[ind_t] += D_orig->landed[i];
	}
      ind += N_int;
      ind_orig += N_int; 
      ind_t++;
    }
  printf("\n");

  printf("Allocating space and initalize simulated discards for market landing data\n");
  /* Allocating space and initalize simulated discards for market landing data */
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "n_int_len_disc")))
    N_int = INTEGER_VALUE(elmt); // number of intervals for length
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "int_len_vec_disc")))
    int_len = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length value for intervals
  if(!Rf_isNull(elmt = getListElement(i_COSTList, "int_len_lim_disc")))
    D_mland->int_len_lim = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length value for intervals
  n_fish = (N_int)*D_mland->n_trip;
  D_mland->N_int_disc = N_int;
  D_mland->l_disc = CALLOC(n_fish,double); //Free ok
  D_mland->lfreq_disc = CALLOC(n_fish,int); //Free ok
  ind = 0;
  for(t=0;t<D_mland->n_trip;t++)
    {
      for(f=0;f<N_int;f++)
	{
	  D_mland->l_disc[ind] = int_len[f];
	  D_mland->lfreq_disc[ind] = 0;
	  ind++;
	}
    }
  D_mland->lambda = CALLOC(D_mland->n_trip,double); //Free ok
  printf("\n");

  #ifdef DEBUG_COST
  fclose(caa_debug);
  #endif

  FREE(lengths);
  FREE(P_l);
  FREE(lengths2);
  FREE(p);

  D_COST->obs = D_obs;
  D_COST->mland = D_mland;

  *o_D_orig = D_orig;
  *o_D_COST = D_COST;

  return(0);
}		/* end of makedata_COST */

    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in makedata_COST
*/
int re_makedata_COST(Data_orig **i_D_orig, Data_COST **i_D_COST)
{
  Data_orig   *D_orig;
  Data_COST   *D_COST;

  D_orig = *i_D_orig;
  D_COST = *i_D_COST;

  FREE(D_COST->cens->r);
  FREE(D_COST->cens->mu);
  FREE(D_COST->cens->tau);
  FREE(D_COST->cens);

  FREE(D_COST->mland->l_disc);
  FREE(D_COST->mland->lfreq_disc);
  FREE(D_COST->mland->lambda);

  FREE(D_COST->obs);
  FREE(D_COST->mland);
  FREE(D_COST);

  FREE(D_orig->nFishBoat);
  FREE(D_orig->totage);
  FREE(D_orig->totlength);
  FREE(D_orig->totweight);
  FREE(D_orig->replength);
  FREE(D_orig->discard);
  FREE(D_orig->landed);
  FREE(D_orig->start_noAge);
  FREE(D_orig->start_Age);
  FREE(D_orig->num_noAge);
  FREE(D_orig->haulweight);
  FREE(D_orig->season);
  FREE(D_orig->n_discard);
  FREE(D_orig->n_landed);
  FREE(D_orig);
  
  return(0);
}		/* end of re_makedata_COST */


    

/*!
  \author Hanne Rognebakke
  \brief Initialize cens_struct in COST object
*/
int init_cens_COST(Data_orig *i_D_orig, Data_COST *i_D_COST)
{
  int i,f,h;
  int f_end;
  int sum;
  double quant,r,l;

  i_D_COST->cens->k = i_D_COST->cens->mu[0];
  i_D_COST->cens->m = i_D_COST->cens->mu[1];
  for(h=0;h<i_D_COST->cens->ncat;h++)
    {
      f = i_D_orig->start_noAge[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      sum = 0;
      if(i_D_orig->n_landed[h]>0)
	{
	  quant = 0.05*i_D_orig->n_landed[h];
	  l = 0;
	  while(sum<quant && f<f_end)
	    {
	      if(i_D_orig->landed[f]>0 && l<0.01)
		l = i_D_orig->totlength[f];
	      sum += i_D_orig->landed[f];
	      f++;
	    }
	  if(f == f_end)
	    {
	      r = i_D_COST->cens->mu[2];
	    }
	  else
	    r = i_D_orig->totlength[f-1];
	}
      else if(i_D_orig->n_discard[h]>0)
	{
	  quant = 0.95*i_D_orig->n_discard[h];
	  l = 0;
	  while(sum<quant && f<f_end)
	    {
	      if(i_D_orig->discard[f]>0)
		l = i_D_orig->totlength[f];
	      sum += i_D_orig->discard[f];
	      f++;
	    }
	  if(f == f_end)
	    if(sum>0)
	      r = l;
	    else
	      r = i_D_COST->cens->mu[2];
	  else
	    r = i_D_orig->totlength[f];
	}
      else
	r = i_D_COST->cens->mu[2];
      i_D_COST->cens->r[h] = r;
    }
  
  if(i_D_COST->cens->k<1.0)
    i_D_COST->cens->mu[0] = log(i_D_COST->cens->k)-log(1-i_D_COST->cens->k);
  else
    i_D_COST->cens->mu[0] = 100.0;
  if(0)
    {
      printf("mu[0]=%f,mu[1]=%f,mu[2]=%f\n",
	     i_D_COST->cens->mu[0],i_D_COST->cens->mu[1],i_D_COST->cens->mu[2]);
      printf("tau[0]=%f,tau[1]=%f,tau[2]=%f\n",
	     i_D_COST->cens->tau[0],i_D_COST->cens->tau[1],i_D_COST->cens->tau[2]);
      printf("mu_mean=%f,mu_prec=%f,a=%f,b=%f\n",
	     i_D_COST->cens->mu_prior_mean,i_D_COST->cens->mu_prior_prec,
	     i_D_COST->cens->a_prior,i_D_COST->cens->b_prior);
    }

  return(0);
}		/* end of init_cens_COST */


    

/*!
  \author Hanne Rognebakke
  \brief Writes COST data after using the routine makedata_COST

  Only to be used for testing.
*/
int write_input_model1_COST(Data_orig *i_D_orig, Data_COST *i_D_COST)
{
  SEXP      elmt = R_NilValue;
  int       a,h,i,nBoatsObs,nBoatsMl,nFishObs,nFishMl,n,nFish;
  FILE     *caa_input;

  caa_input = fopen("caa_input_model1_COST.txt","w");
  

  fprintf(caa_input,"Observer data\n");
  nBoatsObs = i_D_COST->obs->n_trip;
  fprintf(caa_input,"Number of trips=%d\n",nBoatsObs);
  nFishObs = 0;
  nFish = 0;
  for(h=0;h<nBoatsObs;h++)
    {
      fprintf(caa_input,"i=%d,nFishBoat=%d,start_Age=%d,start_noAge=%d,num_noAge=%d,season=%d,ndisc=%d,nland=%d\n",
	      h,i_D_orig->nFishBoat[h],i_D_orig->start_Age[h],i_D_orig->start_noAge[h],
	      i_D_orig->num_noAge[h],i_D_orig->season[h],
	      i_D_orig->n_discard[h],i_D_orig->n_landed[h]);
      nFishObs += i_D_orig->nFishBoat[h];
      nFish += i_D_orig->n_landed[h];
    }
  fprintf(caa_input,"n_landed_tot=%d\n",nFish);

  fprintf(caa_input,"n=%d,trip[i],totage[i],totlength[i],replength[i],discard[i],landed[i]:\n",
	  nFishObs);
  h = 0;
  n = i_D_orig->nFishBoat[0]-1;
  for(i=0;i<nFishObs;i++)
    {
      fprintf(caa_input,"%d,%d,%d,%f,%d,%d,%d\n",i,h,
	      i_D_orig->totage[i],i_D_orig->totlength[i],i_D_orig->replength[i],
	      i_D_orig->discard[i],i_D_orig->landed[i]);
      if(i==n)
	{
	  h++;
	  n += i_D_orig->nFishBoat[h];
	}
    }


  fprintf(caa_input,"Market landing data\n");
  nBoatsMl = i_D_COST->mland->n_trip;
  fprintf(caa_input,"Number of trips=%d\n",nBoatsMl);
  nFishMl = 0;
  nFish = 0;
  for(h=nBoatsObs;h<nBoatsObs+nBoatsMl;h++)
    {
      fprintf(caa_input,"i=%d,nFishBoat=%d,start_Age=%d,start_noAge=%d,num_noAge=%d,season=%d,ndisc=%d,nland=%d\n",
	      h,i_D_orig->nFishBoat[h],i_D_orig->start_Age[h],i_D_orig->start_noAge[h],
	      i_D_orig->num_noAge[h],i_D_orig->season[h],
	      i_D_orig->n_discard[h],i_D_orig->n_landed[h]);
      nFishMl += i_D_orig->nFishBoat[h];
      nFish += i_D_orig->n_landed[h];
    }
  fprintf(caa_input,"n_landed_tot=%d\n",nFish);
  
  fprintf(caa_input,"n=%d,trip[i],totage[i],totlength[i],replength[i],discard[i],landed[i]:\n",
	  nFishMl);
  h = nBoatsObs;
  for(i=nFishObs;i<nFishObs+nFishMl;i++)
    {
      fprintf(caa_input,"%d,%d,%d,%f,%d,%d,%d\n",i,h,i_D_orig->totage[i],
	      i_D_orig->totlength[i],i_D_orig->replength[i],
	      i_D_orig->discard[i],i_D_orig->landed[i]);
      if(i==n && i<nFishObs+nFishMl-1)
	{
	  h++;
	  n += i_D_orig->nFishBoat[h];
	}
    }

  fclose(caa_input);

  return(0);
}		/* end of write_input_model1_COST */



/*!
  \author Hanne Rognebakke
  \brief  Writes the MCMC samples for COST specific parameters to the mcmc_COST-vector

  All parameters are stored sequentially
*/
int write_it_COST(FILE *fp,Data_COST *i_D_COST)
{
  int  i,n,ind;
  
  ind = 0;
  printf("write_it_COST: write iteration to binary file!!!\n");

  n = 0;
  for(i=0;i<3;i++)
    {
      i_D_COST->mcmc[ind+n] = i_D_COST->cens->mu[i];
      i_D_COST->mcmc[ind+n+1] = i_D_COST->cens->tau[i];
      n += 2;
    }

  ind += n;

  i_D_COST->mcmc[ind] = i_D_COST->cens->k;
  i_D_COST->mcmc[ind+1] = i_D_COST->cens->m;

  ind += 2;
  for(i=0;i<i_D_COST->cens->ncat;i++)
    {
      i_D_COST->mcmc[ind+i] = i_D_COST->cens->r[i];
    }

  return(0);
}		/* end of write_it_COST */


/*!
  \author Hanne Rognebakke
  \brief  Picks out for iteration it from mcmc_COST-vector

  All parameters are stored sequentially
*/
int read_it_COST(int i_it, Data_COST *i_D_COST)
{
  int  i,n,ind;

  ind = i_it*i_D_COST->num_var;
  n = 0;
  for(i=0;i<3;i++)
    {
      i_D_COST->cens->mu[i] = i_D_COST->mcmc[ind+n];
      i_D_COST->cens->tau[i] = i_D_COST->mcmc[ind+n+1];
      n += 2;
    }

  ind += n;

  i_D_COST->cens->k = i_D_COST->mcmc[ind];
  i_D_COST->cens->m = i_D_COST->mcmc[ind+1];

  ind += 2;
  for(i=0;i<i_D_COST->cens->ncat;i++)
    {
      i_D_COST->cens->r[i] = i_D_COST->mcmc[ind+i];
    }

  return (0);
}		/* end of read_it_COST */

int makedata_COST_predict(int i_nHaul,Data_COST **o_D_COST)
{
  Data_COST     *D_COST;


  /* Allocating space for COST object */
  D_COST = CALLOC(1,Data_COST);

  D_COST->obs = NULL;
  D_COST->mland = NULL;  

  /* Allocating space for censoring parameters */
  D_COST->cens = CALLOC(1,cens_struct);
  D_COST->cens->ncat = i_nHaul;
  D_COST->cens->r = CALLOC(D_COST->cens->ncat,double);
  D_COST->cens->mu = CALLOC(3,double);
  D_COST->cens->tau = CALLOC(3,double);

  *o_D_COST = D_COST;

  return(0);
}		/* end of makedata_COST_predict */

int sample_lambda_init_COST(Data_age *i_D_age,Data_COST *i_D_COST)
{
  int a,h,T,N_land;
  double sumLambda,sumsqLambda,meanLambda,varLambda;

  T = i_D_COST->mland->n_trip;
  sumLambda = G_ZERO;
  for(h=0;h<T;h++)
    {
      N_land = 0;
      for(a=0;a<i_D_age->glm->ncat;a++)
	N_land += i_D_age->Ages_land[h][a];

      i_D_COST->mland->lambda[h] = 2*N_land;
      sumLambda += i_D_COST->mland->lambda[h];
    }
  meanLambda = sumLambda/T;
  sumsqLambda = G_ZERO;
  for(h=0;h<T;h++)
    sumsqLambda += (i_D_COST->mland->lambda[h]-meanLambda)*(i_D_COST->mland->lambda[h]-meanLambda);
  varLambda = sumsqLambda/(T-1);

  i_D_COST->mland->c = meanLambda*meanLambda/varLambda;
  i_D_COST->mland->d = meanLambda/varLambda;

  return(0);
}		/* end of sample_lambda_COST */

int sample_lambda_prior_COST(Data_COST *i_D_COST)
{
  int i,h,T,num;
  double sumLambda,sumLogLambda;
  double c_new,c_old,log_new,log_old,accProb,u;
  double e = 0.001;
  double f = 0.001;
  double fac = 0.01;

  sumLambda = G_ZERO;
  sumLogLambda = G_ZERO;
  for(h=0;h<i_D_COST->mland->n_trip;h++)
    {
      sumLambda += i_D_COST->mland->lambda[h];
      sumLogLambda += log(i_D_COST->mland->lambda[h]);
    }
  T = i_D_COST->mland->n_trip;

  num=1000;
  c_old = i_D_COST->mland->c;
  for(i=0;i<num;i++)
    {
      c_new = scale_proposal(c_old,fac,NULL);
      
      log_new = -T*log(exp(gammln(c_new))) + (c_new-1)*sumLogLambda
	         +log(exp(gammln(e+c_new*T))) - (e+c_new*T)*log(f+sumLambda);
      log_old = -T*log(exp(gammln(c_old))) + (c_old-1)*sumLogLambda
	         +log(exp(gammln(e+c_old*T))) - (e+c_old*T)*log(f+sumLambda);
      accProb = log_new - log_old;
      u = genunf(G_ZERO,G_ONE);
      if(accProb > -1.0e32 && accProb < 1.0e32 && log(u) < accProb)
	c_old = c_new;
    }
  i_D_COST->mland->c = c_old;

  i_D_COST->mland->d = gengam(f+sumLambda,e+i_D_COST->mland->c*T);

  return(0);
}		/* end of sample_lambda_prior_COST */

    

/*!
  \author Hanne Rognebakke
  \brief */
int resample_data_COST(Data_orig *i_D_orig, Data_COST *i_D_COST)
{
  int            i,f,h,n,s,t;
  int            l_int,N_int,nHaul,nSize,N_sim,N_orig;
  int            ind,ind_alk,ind_fish,ind_fish_l,ind_haul,ind_size,ind_orig,ind_t;
  long          *lengths,*lengths2;
  double         l,psum;
  double        *P_l,*p;

  N_int = i_D_orig->n_int_len;

  lengths = CALLOC(N_int,long);      // Free ok
  P_l = CALLOC(N_int,double);      // Free ok
  p = CALLOC(2,double);      // Free ok
  lengths2 = CALLOC(2,long);      // Free ok

  /* Simulate total lengths for observer data */
  ind_fish = 0;
  ind_fish_l = 0;
  ind_haul = 0;
  ind_size = 0;
  ind_alk = 0;
  ind_orig = 0;
  ind = 0;
  for(t=0;t<i_D_COST->obs->n_trip;t++)
    {
      /* Discard data */
      N_sim = 0;
      N_orig = 0;
      ind_orig = i_D_orig->start_noAge[t];
      i_D_orig->n_discard[t] = 0;
      i_D_orig->n_landed[t] = 0;
      for(f=0;f<N_int;f++)
	{
	  i_D_orig->replength[ind_orig] = 0;
	  i_D_orig->discard[ind_orig] = 0;
	  i_D_orig->landed[ind_orig] = 0;
	  ind_orig++;
	}
      ind_orig = i_D_orig->start_noAge[t];
      for(h=0;h<i_D_COST->obs->num_trip[t];h++)
	{
	  if(i_D_COST->obs->num_haul_disc[ind_haul]>0)
	    {
	      nHaul = 0;
	      for(i=0;i<N_int;i++)
		P_l[i] = 0.0;
	      for(f=0;f<i_D_COST->obs->num_haul_disc[ind_haul];f++)
		{
		  l = i_D_COST->obs->l_disc[ind_fish];
		  l_int = 0;
		  while(l > i_D_orig->int_len_lim[l_int])
		    l_int++;
		  P_l[l_int] += i_D_COST->obs->lfreq_disc[ind_fish];
		  i_D_orig->replength[ind_orig+l_int] += i_D_COST->obs->lfreq_disc[ind_fish];
		  i_D_orig->discard[ind_orig+l_int] += i_D_COST->obs->lfreq_disc[ind_fish];
		  nHaul += i_D_COST->obs->lfreq_disc[ind_fish];
		  ind_fish++;
		}
	      // convert to probabilities
	      for(i=0;i<N_int;i++)
		P_l[i] /= nHaul;
	      // number of fish to be simulated
	      if(nHaul==0)
		n=0;
	      else
		n = (int) nHaul*(i_D_COST->obs->haulsize_disc[ind_haul]/i_D_COST->obs->sampsize_disc[ind_haul]-1);
	      N_sim += nHaul+n;
	      N_orig += nHaul;
	      my_genmul(n,P_l,N_int,lengths);
	      for(i=0;i<N_int;i++)
		{
		  i_D_orig->replength[ind_orig+i] += (int) lengths[i];
		  i_D_orig->discard[ind_orig+i] += (int) lengths[i];
		}
	    }
	  ind_haul++;
	}
      /* Landing data */
      for(s=0;s<i_D_COST->obs->num_trip_land[t];s++)
	{
	  nSize = 0;
	  for(i=0;i<N_int;i++)
	    P_l[i] = 0.0;
	  for(f=0;f<i_D_COST->obs->num_size_land[ind_size];f++)
	    {
	      l = i_D_COST->obs->l_land[ind_fish_l];
	      l_int = 0;
	      while(l > i_D_orig->int_len_lim[l_int])
		l_int++;
	      P_l[l_int] += i_D_COST->obs->lfreq_land[ind_fish_l];
	      i_D_orig->replength[ind_orig+l_int] += i_D_COST->obs->lfreq_land[ind_fish_l];
	      i_D_orig->landed[ind_orig+l_int] += i_D_COST->obs->lfreq_land[ind_fish_l];
	      nSize += i_D_COST->obs->lfreq_land[ind_fish_l];
	      ind_fish_l++;
	    }
	  // convert to probabilities
	  for(i=0;i<N_int;i++)
	    P_l[i] /= nSize;
	  // number of fish to be simulated
	  n = nSize*(i_D_COST->obs->totsize_land[ind_size]/i_D_COST->obs->sampsize_land[ind_size]-1);
	  N_sim += nSize+n;
	  N_orig += nSize;
          my_genmul(n,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      i_D_orig->replength[ind_orig+i] += (int) lengths[i];
	      i_D_orig->landed[ind_orig+i] += (int) lengths[i];
	    }
	  ind_size++; 
	}
      //scale to original size
      if(N_orig>0)
	{
	  psum = 0;
	  for(i=0;i<N_int;i++)
	    {
	      P_l[i] = i_D_orig->replength[ind_orig+i]*N_orig/N_sim;
	      psum += P_l[i];
	    }
	  for(i=0;i<N_int;i++)
	    P_l[i] /= psum;
	  my_genmul(N_orig,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      if(lengths[i]==0)
		{
		  i_D_orig->replength[ind_orig+i] = 0;
		  i_D_orig->discard[ind_orig+i] = 0;
		  i_D_orig->landed[ind_orig+i] = 0;
		}
	      else
		{
		  p[0] = (double) i_D_orig->discard[ind_orig+i]/(i_D_orig->discard[ind_orig+i]+i_D_orig->landed[ind_orig+i]);
		  p[1] = (double) i_D_orig->landed[ind_orig+i]/(i_D_orig->discard[ind_orig+i]+i_D_orig->landed[ind_orig+i]);
		  i_D_orig->replength[ind_orig+i] = (int) lengths[i];
		  my_genmul(i_D_orig->replength[ind_orig+i],p,2,lengths2);
		  i_D_orig->discard[ind_orig+i] = (int) lengths2[0];
		  i_D_orig->landed[ind_orig+i] = (int) lengths2[1];
		}
	    }
	}
      // put the age-length data into D_orig object
      // remove length count for lengths with missing ages
      n = 0;
      for(f=0;f<i_D_COST->obs->num_alk[t];f++)
	{
	  i_D_orig->totage[ind] = i_D_COST->obs->alk_a[ind_alk];
	  i_D_orig->totlength[ind] = i_D_COST->obs->alk_l[ind_alk];
	  i_D_orig->replength[ind] = i_D_COST->obs->alk_lfreq[ind_alk];
	  i_D_orig->discard[ind] = i_D_COST->obs->alk_lfreq[ind_alk];
	  l_int = 0;
	  while(i_D_COST->obs->alk_l[ind_alk] > i_D_orig->int_len_lim[l_int])
	    l_int++;
	  i_D_orig->replength[ind_orig+l_int] -= i_D_COST->obs->alk_lfreq[ind_alk];
	  i_D_orig->discard[ind_orig+l_int] -= i_D_COST->obs->alk_lfreq[ind_alk];
	  if(i_D_orig->discard[ind_orig+l_int]<0)
	    {
	      n += 1;
	      i_D_orig->replength[ind] = i_D_COST->obs->alk_lfreq[ind_alk]-1;
	      i_D_orig->discard[ind] = i_D_COST->obs->alk_lfreq[ind_alk]-1;
	      i_D_orig->replength[ind_orig+l_int] = i_D_orig->landed[ind_orig+l_int];
	      i_D_orig->discard[ind_orig+l_int] = 0;
	    }
	  ind_alk++;
	  ind++;
	}
      for(i=i_D_orig->start_noAge[t];i<i_D_orig->start_noAge[t]+i_D_orig->num_noAge[t];i++)
	{
	  i_D_orig->n_landed[t] += i_D_orig->landed[i];
	  i_D_orig->n_discard[t] += i_D_orig->discard[i];
	}
      ind += N_int;
      ind_orig += N_int;
    }//end for(t=0;t<i_D_COST->obs->n_trip;t++)


  /* Simulate total lengths for market landing data */

  ind_fish = 0;
  ind_size = 0;
  ind_alk = 0;
  ind_t = i_D_COST->obs->n_trip;
  for(t=0;t<i_D_COST->mland->n_trip;t++)
    {
      i_D_orig->n_discard[ind_t] = 0;
      i_D_orig->n_landed[ind_t] = 0;
      ind_orig = i_D_orig->start_noAge[ind_t];
      N_sim = 0;
      N_orig = 0;
      for(f=0;f<N_int;f++)
	{
	  i_D_orig->replength[ind_orig] = 0;
	  i_D_orig->discard[ind_orig] = 0;
	  i_D_orig->landed[ind_orig] = 0;
	  ind_orig++;
	}
      ind_orig = i_D_orig->start_noAge[ind_t];
      for(s=0;s<i_D_COST->mland->num_trip[t];s++)
	{
	  nSize = 0;
	  for(i=0;i<N_int;i++)
	    P_l[i] = 0.0;
	  for(f=0;f<i_D_COST->mland->num_size[ind_size];f++)
	    {
	      l = i_D_COST->mland->l[ind_fish];
	      l_int = 0;
	      while(l > i_D_orig->int_len_lim[l_int])
	      	l_int++;
	      P_l[l_int] += i_D_COST->mland->lfreq[ind_fish];
	      i_D_orig->replength[ind_orig+l_int] += i_D_COST->mland->lfreq[ind_fish];
	      i_D_orig->landed[ind_orig+l_int] += i_D_COST->mland->lfreq[ind_fish];
	      nSize += i_D_COST->mland->lfreq[ind_fish];
	      ind_fish++;
	    }
	  // convert to probabilities
	  for(i=0;i<N_int;i++)
	    P_l[i] /= nSize;
	  // number of fish to be simulated
	  n = nSize*(i_D_COST->mland->totsize[ind_size]/i_D_COST->mland->sampsize[ind_size]-1);
	  N_sim += nSize+n;
	  N_orig += nSize;
          my_genmul(n,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      i_D_orig->replength[ind_orig+i] += (int) lengths[i];
	      i_D_orig->landed[ind_orig+i] += (int) lengths[i];
	    }
	  ind_size++; 
	}
      //scale to original size
      if(N_orig>0)
	{
	  psum = 0;
	  for(i=0;i<N_int;i++)
	    {
	      P_l[i] = i_D_orig->replength[ind_orig+i]*N_orig/N_sim;
	      psum += P_l[i];
	    }
	  for(i=0;i<N_int;i++)
	    P_l[i] /= psum;
	  my_genmul(N_orig,P_l,N_int,lengths);
	  for(i=0;i<N_int;i++)
	    {
	      if(lengths[i]==0)
		{
		  i_D_orig->replength[ind_orig+i] = 0;
		  i_D_orig->landed[ind_orig+i] = 0;
		}
	      else
		{
		  i_D_orig->replength[ind_orig+i] = (int) lengths[i];
		  i_D_orig->landed[ind_orig+i] = (int) lengths[i];
		}
	    }
	}
      // remove length count for lengths with missing ages
      n = 0;
      for(f=0;f<i_D_COST->mland->num_alk[t];f++)
	{
	  i_D_orig->totage[ind] = i_D_COST->mland->alk_a[ind_alk];
	  i_D_orig->totlength[ind] = i_D_COST->mland->alk_l[ind_alk];
	  i_D_orig->replength[ind] = i_D_COST->mland->alk_lfreq[ind_alk];
	  i_D_orig->landed[ind] = i_D_COST->mland->alk_lfreq[ind_alk];
	  l_int = 0;
	  while(i_D_COST->mland->alk_l[ind_alk] > i_D_orig->int_len_lim[l_int])
	    l_int++;
	  i_D_orig->replength[ind_orig+l_int] -= i_D_COST->mland->alk_lfreq[ind_alk];
	  i_D_orig->landed[ind_orig+l_int] -= i_D_COST->mland->alk_lfreq[ind_alk];
	  if(i_D_orig->landed[ind_orig+l_int]<0)// remove age data from sample
	    {
	      n -= i_D_orig->landed[ind_orig+l_int];
	      i_D_orig->replength[ind] = i_D_COST->mland->alk_lfreq[ind_alk]+i_D_orig->landed[ind_orig+l_int];
	      i_D_orig->landed[ind] = i_D_COST->mland->alk_lfreq[ind_alk]+i_D_orig->landed[ind_orig+l_int];
	      i_D_orig->replength[ind_orig+l_int] = i_D_orig->discard[ind_orig+l_int];
	      i_D_orig->landed[ind_orig+l_int] = 0;
	    }
	  ind_alk++;
	  ind++;
	}
      for(i=i_D_orig->start_noAge[t];i<i_D_orig->start_noAge[t]+i_D_orig->num_noAge[t];i++)
	{
	  i_D_orig->n_landed[ind_t] += i_D_orig->landed[i];
	}
      ind += N_int;
      ind_orig += N_int; 
      ind_t++;
    }//end for(t=0;t<i_D_COST->mland->n_trip;t++)

  FREE(lengths);
  FREE(P_l);
  FREE(lengths2);
  FREE(p);

  return(0);
}		/* end of makedata_COST */


/*!
  \author Geir Storvik/Hanne Rognebakke
  \brief Calculate sufficient statistics for length given age model

  The sufficient statistics are conditional on the \f$g(a)\f$ function.
  If i_use_sim_age=1, the data were ages are simulated is included, else
  if 0, only observed age data is used.

  Sufficient statistics are for each haul
  - the number of fish per age
  - the sum of $g(a)$
  - the square sum of \f$g(a)\f$
  - Least squares intercept, slope and sum of squares
*/
int make_suff_lga_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, Data_lin *i_D_lga,
		       Data_g_a *i_D_g_a, int i_start_h, int i_use_sim_age)
{
  int      a,h,sum_tot=0;
  double   g_a,maxAges,N_h,N_h_a,sum_g,sum_g2,sum_l,sum_l2,sum_gl,beta0,beta1,ssq;
  double  *Ages,*sum_by_cat,*sqsum_by_cat;

  Ages = CALLOC(i_D_g_a->ncat,double);
  sum_by_cat = CALLOC(i_D_g_a->ncat,double);
  sqsum_by_cat = CALLOC(i_D_g_a->ncat,double);


  for(h=i_start_h;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  Ages[a] = G_ZERO;
	  sum_by_cat[a] = G_ZERO; 
	  sqsum_by_cat[a] = G_ZERO;
	}
      
      if(h<i_D_COST->obs->n_trip) // if observer data 
	{
	  if(i_D_orig->start_noAge[h]>i_D_orig->start_Age[h]) // if observed ages in given haul
	    {
	      for(a=0;a<i_D_g_a->ncat;a++)
		{
		  Ages[a] = i_D_lga->Ages[h][a];
		  sum_tot += i_D_lga->Ages[h][a];
		  sum_by_cat[a] = i_D_lga->sum_by_cat[h][a];
		  sqsum_by_cat[a] = i_D_lga->sqsum_by_cat[h][a];
		}
	    }
	}
      else // market landing data
	{
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      if(i_D_lga->Ages_fix[h][a]>0) // if observed ages in given haul and age group
		{
		  Ages[a] = i_D_lga->Ages[h][a];
		  sum_tot += i_D_lga->Ages[h][a];
		  sum_by_cat[a] = i_D_lga->sum_by_cat[h][a];
		  sqsum_by_cat[a] = i_D_lga->sqsum_by_cat[h][a];
		}
	    }
	}
      
      sum_g = G_ZERO;
      sum_g2 = G_ZERO;
      sum_gl = G_ZERO;
      sum_l = G_ZERO;
      sum_l2 = G_ZERO;
      maxAges = G_ZERO;
      N_h = G_ZERO;
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  g_a = i_D_g_a->g_a[a];
          N_h_a = Ages[a];
          sum_g +=  N_h_a * g_a;
          maxAges = max(maxAges,Ages[a]);
          sum_g2 +=  N_h_a * g_a * g_a;
          sum_gl += g_a * sum_by_cat[a];
	  sum_l += sum_by_cat[a];
	  sum_l2 += sqsum_by_cat[a];
          N_h += N_h_a;
	}
      if(N_h < 0.0001)
	{
	  /* No data */
          beta0 = G_ZERO;
          beta1 = G_ZERO;
          ssq = G_ZERO;
	}
      else if(fabs(maxAges-N_h)<0.001)
	{ /* Only one age group sampled */
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g2-sum_g)<0.0001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g*sum_g-N_h*sum_g2)<0.0000000000001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else
	{
	  beta1 = (sum_l*sum_g-N_h*sum_gl)/(sum_g*sum_g-N_h*sum_g2);
	  beta0 = (sum_l-beta1*sum_g)/N_h;
	  ssq = sum_l2+beta0*beta0*N_h+beta1*beta1*sum_g2-
	    G_TWO*(beta0*sum_l+beta1*sum_gl-beta0*beta1*sum_g)+0.00000001;
	}
      if(ssq < G_ZERO || !(beta0 > -9999999.99 && beta0 < 999999999.99) ||
                         !(beta0 > -9999999.99 && beta0 < 999999999.99))
	{
	  printf("h=%d,beta0=%lf,beta1=%lf,ssq=%lf\n",h,beta0,beta1,ssq);
          for(a=0;a<i_D_g_a->ncat;a++)
	    printf("%f %f %f %f\n",Ages[a],g_a,sum_by_cat[a],sqsum_by_cat[a]);
	  write_warning("make_suff_lga:Something is wrong\n");
	  return(1);
	}
      i_D_lga->glm->beta_hat[h][0][1] = beta1;
      i_D_lga->glm->beta_hat[h][0][0] = beta0;
      i_D_lga->glm->ssq[h] = ssq;
      i_D_lga->glm->suff[h][0][0] = N_h;
      i_D_lga->glm->suff[h][0][1] = sum_g;
      i_D_lga->glm->suff[h][1][0] = sum_g;
      i_D_lga->glm->suff[h][1][1] = sum_g2;
      if(0)
	printf("h=%d, beta0=%12.10f, beta1=%12.10f, suff00=%10.8f, suff01=%10.8f, suff11=%10.8f\n",
	       h,beta0,beta1,N_h,sum_g,sum_g2);

      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = N_h * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff[h][0][1]*i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][2] = N_h * i_D_lga->haulweight[h]*i_D_lga->haulweight[h];
	  
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff[h][0][2];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff[h][1][2];
	  i_D_lga->glm->beta_hat[h][0][2] = G_ZERO;
	}
    }
  FREE(Ages);
  FREE(sum_by_cat);
  FREE(sqsum_by_cat);

  return(0);
}               /* end of make_suff_lga_COST */

int suff_g_a_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, LW_struct *i_length, Data_age *i_D_age, 
		  Data_lin *i_D_lga, Data_g_a *i_D_g_a, int i_start_h, int i_nHaul, double **o_suff)
{
  int         a,h,i;
  double      sum_l,sum_l2,N_h_a,weight;
  double     *beta;
  Data_glm   *glm;

  glm = i_D_lga->glm;
  beta = CALLOC(glm->nxcov,double); //Free ok


  if(glm->nxcov==2)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  o_suff[2][a] = G_ZERO;
	  for(h=i_start_h;h<i_nHaul;h++)
	    {
	      if(((h<i_D_COST->obs->n_trip) && (i_D_orig->start_noAge[h]>i_D_orig->start_Age[h])) || (i_D_lga->Ages_fix[h][a]>0)) 
		// if observer data and observed ages in given haul or if market landing data and observed ages in given haul and age group
		{
		  sum_l = i_D_lga->sum_by_cat[h][a];
		  sum_l2 = i_D_lga->sqsum_by_cat[h][a];
		  N_h_a = (double) i_D_lga->Ages[h][a];
		  for(i=0;i<glm->nxcov;i++)
		    beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
		  o_suff[0][a] += sum_l2 - 2*beta[0]*sum_l + N_h_a*beta[0]*beta[0];
		  o_suff[1][a] += -2*(sum_l*beta[1] - N_h_a*beta[0]*beta[1]);
		  o_suff[2][a] += beta[1]*beta[1]*N_h_a;
		}
	    }
	}
    }
  else if(glm->nxcov==3)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  o_suff[2][a] = G_ZERO;
	  for(h=i_start_h;h<i_nHaul;h++)
	    {
	      sum_l = i_D_lga->sum_by_cat_fix[h][a];
	      sum_l2 = i_D_lga->sqsum_by_cat_fix[h][a];
	      weight = i_D_lga->haulweight[h];
	      N_h_a = (double) i_D_lga->Ages_fix[h][a];
	      for(i=0;i<glm->nxcov;i++)
		beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	      o_suff[0][1] += sum_l2 - 2*beta[0]*sum_l + N_h_a*beta[0]*beta[0] 
		              - 2*beta[2]*weight*sum_l + 2*beta[0]*beta[2]*weight
                              + beta[2]*weight*beta[2]*weight;
	      o_suff[1][a] += -2*(sum_l*beta[1] - N_h_a*beta[0]*beta[1]
                              -N_h_a*beta[2]*weight*beta[1]);
	      o_suff[2][a] += beta[1]*beta[1]*N_h_a;
	    }
	}
    }

  FREE(beta);

  return(0);
}               /* end of suff_g_a */
