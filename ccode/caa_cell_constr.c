/*!
  \file caa_routines.c
  \brief Routines for simulating cell effects
*/
#include "caa.h"
#include "caa_read_write.h"
#include "caa_cell_constr.h"
#include "caa_chol.h"
#include "caa_utl.h"

/*!
  \brief Makes the cell structure for simulating cell effects

  Memory allocated in this routine is reallocated in ::re_make_cell_constr
  \author Geir Storvik
*/
int make_cell_constr_age(Data_age *i_D_age, Input_cov *i_cov)
{
  int         j,k,ind,ncell;
  double     *constr_cell, *Sigma_cell;
  Data_cov   *xcov;

    
  xcov = i_D_age->glm->xcov[0];
  #ifdef DEBUG_PROG
  fprintf(stderr,"age icell=%d,ncat=%d\n",xcov->icell,i_D_age->glm->ncat);
  #endif
  if(xcov->icell>=0)
    {
      ncell = xcov->n_fac[xcov->icell]*i_D_age->glm->ncat;
      xcov->Q_cell = CALLOC2_d(ncell,ncell);
      xcov->Qchol_cell = CALLOC2_d(ncell,ncell);
      xcov->cell_vec = CALLOC(ncell,double);
      Sigma_cell = i_cov->Sigma_cell;
      ind=0;
      for(j=0;j<ncell;j++)
	{
	  for(k=0;k<ncell;k++)
	    {
	      xcov->Q_cell[j][k] = Sigma_cell[ind];
	      ind++;
	    }
	}
      #ifdef DEBUG_CELL
      FILE *unit;
      unit = fopen("Q_cell_age.dat","w");
      for(j=0;j<ncell;j++)
	{
	  for(k=0;k<ncell;k++)
	    {
	      fprintf(unit,"%f ",xcov->Q_cell[j][k]);
	    }
	  fprintf(unit,"\n");
	}
      fclose(unit);
      #endif

      xcov->n_constr_cell = i_cov->nconstr_cell;
      constr_cell = i_cov->constr_cell;

      xcov->constr_cell = CALLOC2_d(xcov->n_constr_cell,ncell);
      ind = 0;
      for(j=0;j<xcov->n_constr_cell;j++)
	{
	  for(k=0;k<ncell;k++)
	    {
	      xcov->constr_cell[j][k] = constr_cell[ind];
	      ind++;
	    }
	}
      #ifdef DEBUG_CELL
      unit = fopen("constr_cell_age.dat","w");
      for(j=0;j<xcov->n_constr_cell;j++)
	{
	  for(k=0;k<ncell;k++)
	    fprintf(unit,"%lf ",xcov->constr_cell[j][k]);
	  fprintf(unit,"\n");
	}
      fclose(unit);
      #endif
    }
  return(0);
}               /* end of make_cell_constr_age */

/*!
  \brief Reallocate memory allocated in ::make_cell_constr_age
  \author Hanne Rognebakke
*/
int re_make_cell_constr(Data_glm *i_glm)
{
  int    i,ncell;
  Data_cov *xcov;

  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->icell>=0)
	{
	  ncell = xcov->n_fac[xcov->icell]*i_glm->ncat;
	  FREE2_d(xcov->Q_cell,ncell);
	  FREE2_d(xcov->Qchol_cell,ncell);
	  FREE(xcov->cell_vec);
	  FREE2_d(xcov->constr_cell,xcov->n_constr_cell);
	}
    }

  return(0);
}               /* end of re_make_cell_constr_age */

/*!
  \brief Makes the cell structure for simulating cell effects

  Memory allocated in this routine is reallocated in ::re_make_cell_constr
  \author Geir Storvik
*/
int make_cell_constr_lin(Data_lin *i_D_lin, Input_cov *i_int_cov, Input_cov *i_slp_cov)
{
  int         i,j,k,ind,ncell;
  double     *constr_cell;
  double     *Sigma_cell;
  Data_cov   *xcov;
  Input_cov  *input_cov;
  
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      if(i==0)
	input_cov = i_int_cov;
      else if(i==1)
	input_cov = i_slp_cov;
      else
	{
	  write_warning("make_cell_constr_lin: Only 2 nxcov possible\n");
	  return(1);
	}
      xcov = i_D_lin->glm->xcov[i];
      if(xcov->icell>=0)
	{
	  ncell = xcov->n_fac[xcov->icell];
          #ifdef DEBUG_PROG
	  fprintf(stderr,"lin icell=%d,ncell=%d,nconstr=%d\n",
		  xcov->icell,ncell,input_cov->nconstr_cell);
          #endif
	  xcov->Q_cell = CALLOC2_d(ncell,ncell);
	  xcov->Qchol_cell = CALLOC2_d(ncell,ncell);
	  xcov->cell_vec = CALLOC(ncell,double);
	  Sigma_cell = input_cov->Sigma_cell;
	  ind=0;
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  xcov->Q_cell[j][k] = Sigma_cell[ind];
		  ind++;
		}
	    }

          #ifdef DEBUG_CELL
	  FILE *unit;
	  unit = fopen("Q_cell_lin.dat","w");
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%f ",xcov->Q_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
          #endif
	  xcov->n_constr_cell = input_cov->nconstr_cell;
	  constr_cell = input_cov->constr_cell;
	  xcov->constr_cell = CALLOC2_d(xcov->n_constr_cell,ncell);
	  ind = 0;
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  xcov->constr_cell[j][k] = constr_cell[ind];
		  ind++;
		}
	    }
          #ifdef DEBUG_CELL
	  unit = fopen("constr_cell_lin.dat","w");
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%lf ",xcov->constr_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
          #endif
	}
    }
  return(0);
}

/*!
  \author Geir Storvik
  \brief  Put age cell distributions into structs 
  This routine is called from ::predict
*/  
int makedata_cell_dist_age(int i_num_cell_o, int i_num_cell_u,
			   int i_age_int_nC, double *i_age_int_E, double *i_age_int_C,
			   Age_struct *i_age, Data_age *i_D_age)
{
  int       j,k;
  double   *E=NULL,*C=NULL;
  Data_cov *xcov;

  i_age->par->cell = CALLOC(i_D_age->glm->nxcov,double *);

  xcov = i_D_age->glm->xcov[0];
  xcov->icell = i_D_age->glm->xcov[0]->icell;
  #ifdef DEBUG_PROG
  fprintf(stderr,"makedata_cell_dist_age:icell=%d\n",xcov->icell);
  #endif
  if(xcov->icell>=0)
    {
      xcov->cell_dist = CALLOC(1,Cell_dist);
      xcov->cell_dist->n_o = i_num_cell_o;
      xcov->cell_dist->n_u = i_num_cell_u;
      i_age->par->cell[0] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
      if(xcov->cell_dist->n_u>0)
	{
	  xcov->cell_dist->n_C = i_age_int_nC;
	  E = i_age_int_E;
	  C = i_age_int_C;
	  
	  xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	  for(j=0;j<xcov->cell_dist->n_u;j++)
	    for(k=0;k<xcov->cell_dist->n_o;k++)
	      {
		xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      }
	  
	  xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	  for(j=0;j<xcov->cell_dist->n_u;j++)
	    for(k=0;k<xcov->cell_dist->n_C;k++)
	      xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	}
    }

  return(0);
}

/*!
  \author Geir Storvik
  \brief  Put cell distributions into structs
*/  
int makedata_cell_dist_lin(int *i_num_cell_o, int *i_num_cell_u,
			   int i_lin_int_nC, double *i_lin_int_E, double *i_lin_int_C,
			   int i_lin_slp_nC, double *i_lin_slp_E, double *i_lin_slp_C,
			   LW_struct *i_lin, Data_lin *i_D_lin)
{
  int   i,j,k,ind;
  double *E=NULL,*C=NULL;
  Data_cov *xcov;

  i_lin->par->cell = CALLOC(i_D_lin->glm->nxcov,double *);

  ind = 0;
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      xcov = i_D_lin->glm->xcov[i];
      xcov->icell = i_D_lin->glm->xcov[i]->icell;
      #ifdef DEBUG_PROG
      fprintf(stderr,"makedata_cell_dist_lin:icell=%d\n",xcov->icell);
      #endif
      if(xcov->icell>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = i_num_cell_o[ind];
	  xcov->cell_dist->n_u = i_num_cell_u[ind];
	  i_lin->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = i_lin_int_nC;
		  E = i_lin_int_E;
		  C = i_lin_int_C;
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = i_lin_slp_nC;
		  E = i_lin_slp_E;
		  C = i_lin_slp_C;
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  {
		    //printf("j=%d,k=%d,ind=%d\n",j,k,j*xcov->cell_dist->n_o+k);
		    xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
		    //printf("E[%d]=%f\n",j*xcov->cell_dist->n_o+k,E[j*xcov->cell_dist->n_o+k]);
		  }
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }
  return(0);
}

/*!
  \author Hanne Rognebakke
  \brief  Put haulsize cell distributions into structs
*/  
int makedata_cell_dist_hsz(int i_num_cell_o, int i_num_cell_u,
			   int i_hsz_int_nC, double *i_hsz_int_E, double *i_hsz_int_C,
			   LW_struct *i_hsz, Data_lin *i_D_hsz)
{
  int       j,k;
  double   *E=NULL,*C=NULL;
  Data_cov *xcov;

  i_hsz->par->cell = CALLOC(i_D_hsz->glm->nxcov,double *);

  xcov = i_D_hsz->glm->xcov[0];
  xcov->icell = i_D_hsz->glm->xcov[0]->icell;
  fprintf(stderr,"makedata_cell_dist_hsz:icell=%d\n",xcov->icell);
  if(xcov->icell>=0)
    {
      xcov->cell_dist = CALLOC(1,Cell_dist);
      xcov->cell_dist->n_o = i_num_cell_o;
      xcov->cell_dist->n_u = i_num_cell_u;
      i_hsz->par->cell[0] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
      if(xcov->cell_dist->n_u>0)
	{
	  xcov->cell_dist->n_C = i_hsz_int_nC;
	  E = i_hsz_int_E;
	  C = i_hsz_int_C;
	  
	  xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	  for(j=0;j<xcov->cell_dist->n_u;j++)
	    for(k=0;k<xcov->cell_dist->n_o;k++)
	      {
		xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      }
	  
	  xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	  for(j=0;j<xcov->cell_dist->n_u;j++)
	    for(k=0;k<xcov->cell_dist->n_C;k++)
	      xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	}
    }

  return(0);
}

int re_makedata_cell_dist_age(Age_struct *i_age,Data_age *i_D_age)
{
  Data_cov *xcov;

  xcov = i_D_age->glm->xcov[0];
  if(xcov->icell>=0 && xcov->cell_dist->n_u>0)
    {
      FREE(i_age->par->cell[0]);
      FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
      FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
      FREE(xcov->cell_dist);
    }
  FREE(i_age->par->cell);

  return(0);
}

int re_makedata_cell_dist_lin(LW_struct *i_lin,Data_lin *i_D_lin)
{
  int    i;
  Data_cov *xcov;

  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      xcov = i_D_lin->glm->xcov[i];
      if(xcov->icell>=0 && xcov->cell_dist->n_u>0)
	{
	  FREE(i_lin->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_lin->par->cell);

  return(0);
}

int re_makedata_cell_dist_hsz(LW_struct *i_hsz,Data_lin *i_D_hsz)
{
  Data_cov *xcov;

  xcov = i_D_hsz->glm->xcov[0];
  if(xcov->icell>=0 && xcov->cell_dist->n_u>0)
    {
      FREE(i_hsz->par->cell[0]);
      FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
      FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
      FREE(xcov->cell_dist);
    }
  FREE(i_hsz->par->cell);

  return(0);
}


/*!
  \author Geir Storvik
  \brief  Simulate unobserved cell effects
  This routine is called from ::predict
*/  
int simulate_cell_effects(Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl,
			  LW_struct *i_hsz,Data_lin *i_D_hsz)
{
  
  int       a,i,j,k,n,ncat;
  double   *eps, *cell, sd;
  Data_cov *xcov;
  
  n = 1;
  if(i_D_age->glm->xcov[0]->icell>=0)
    n = MAX(n,i_D_age->glm->xcov[0]->cell_dist->n_C);
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_lga->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_wgl->glm->xcov[i]->cell_dist->n_C);
    }
  if(i_D_age->glm->xcov[0]->ihaulsize>0) // Haulsize included in model
    {
      if(i_D_hsz->glm->xcov[0]->icell>=0)
	n = MAX(n,i_D_hsz->glm->xcov[0]->cell_dist->n_C);
    }
  eps = CALLOC(n,double);
  
  ncat = i_D_age->glm->ncat;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_age->glm->xcov[i];
	  cell = i_age->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_age->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(a=0;a<ncat;a++)
	    {
	      for(k=0;k<n;k++)
		cell[a*n+k] = i_age->par->eff[a][i][xcov->icell][k];
	    }
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[ncat*n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[ncat*n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[ncat*n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  if(i_D_age->glm->xcov[0]->ihaulsize>0) // Haulsize included in model
    {
      i = 0;
      if(i_D_hsz->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_hsz->glm->xcov[i];
	  cell = i_hsz->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_hsz->par->tau[i][i_D_hsz->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_hsz->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga->glm->xcov[i];
	  cell = i_length->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  cell = i_weight->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  FREE(eps);

  return(0);
}

/*!
  \author Geir Storvik
  \brief  Simulate unobserved cell effects
*/  
int simulate_cell_effects_CC(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_hsz,Data_lin *i_D_hsz,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC)
{
  
  int       a,i,j,k,n,ncat;
  double   *eps, *cell, sd;
  Data_cov *xcov;
  
  n = 1;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_age->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_lga->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_wgl->glm->xcov[i]->cell_dist->n_C);
    }
  if(i_D_age->glm->xcov[0]->ihaulsize>0) // Haulsize included in model
    {
      if(i_D_hsz->glm->xcov[0]->icell>=0)
	n = MAX(n,i_D_hsz->glm->xcov[0]->cell_dist->n_C);
    }
  eps = CALLOC(n,double);
  
  /* age */
  ncat = i_D_age->glm->ncat;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_age->glm->xcov[i];
	  cell = i_age->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_age->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(a=0;a<ncat;a++)
	    {
	      for(k=0;k<n;k++)
		{
		  cell[a*n+k] = i_age->par->eff[a][i][xcov->icell][k];
		}
	    }
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[ncat*n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[ncat*n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[ncat*n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  /* hsz */
  if(i_D_age->glm->xcov[0]->ihaulsize>0) // Haulsize included in model
    {
      i = 0;
      if(i_D_hsz->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_hsz->glm->xcov[i];
	  cell = i_hsz->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_hsz->par->tau[i][i_D_hsz->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_hsz->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  /* lga */
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga->glm->xcov[i];
	  cell = i_length->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* lga - coastal cod */
  for(i=0;i<i_D_lga_CC->glm->nxcov;i++)
    {
      if(i_D_lga_CC->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga_CC->glm->xcov[i];
	  cell = i_length_CC->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length_CC->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length_CC->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* wgl */
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  cell = i_weight->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* wgl - coastal cod */
  for(i=0;i<i_D_wgl_CC->glm->nxcov;i++)
    {
      if(i_D_wgl_CC->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl_CC->glm->xcov[i];
	  cell = i_weight_CC->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight_CC->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight_CC->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  FREE(eps);

  return(0);
}
