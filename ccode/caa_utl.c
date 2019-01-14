#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "caa.h"
#include "caa_read_write.h"

int **CALLOC2_i(int i_nrow,int i_ncol)
{
  int i;
  int **m;	

  m = CALLOC(i_nrow,int *);
  for(i=0;i<i_nrow;i++)
    m[i] = CALLOC(i_ncol,int);

  return m;
}				/* end of CALLOC2_i */


void FREE2_i(int**i_matrix,int i_nrow)
{
  int  i;
  for(i=0;i<i_nrow;i++)
     FREE(i_matrix[i]);
  FREE(i_matrix);
  return;
}				/* end of FREE2_i */

double **CALLOC2_d(int i_nrow,int i_ncol)
{
  int i;
  double **m;		

  m = CALLOC(i_nrow,double *);
  for(i=0;i<i_nrow;i++)
    m[i] = CALLOC(i_ncol,double);

  return m;
}				/* end of CALLOC2_d */


void FREE2_d(double**i_matrix,int i_nrow)
{
  int  i;
  for(i=0;i<i_nrow;i++)
     FREE(i_matrix[i]);
  FREE(i_matrix);
  return;
}				/* end of FREE2_d */

