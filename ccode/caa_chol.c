/*!
  \file caa_chol.c
  \brief Contains various routines related to cholesky decomposition
  \author Geir Storvik

  Some of the routimes are small modifications of the Numrec library.
*/
#include "caa.h"
#include "caa_chol.h"

/*!
  \author Geir Storvik
  \brief Caclulate Cholesky decomposition of a positive-definite symmetric matrix

  On input, only the upper triangle of the matrix \f$a\f$a is needed. 
  The Cholesky factor \f$L\f$ is returned in the lower triangle of a. 

  This is a modification of the choldc routine in Numerical Recipes in C changed
  so that the matrix starts at index 0 and such that the diagonal of \f$L\f$ is stored
  in a and not as a separate vector.
*/
int choldc0(double **a, int n)
{
	void nrerror(char error_text[]);
	int i,j,k;
	double sum;

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					return(1);
				a[i][i]=sqrt(sum);
			} else a[j][i]=sum/a[i][i];
		}
	}
        return(0);
}

/*!
  \author Geir Storvik
  \brief Solves the set of \f$n\f$ linear equations \f$\bf A\cdot \bf x=\bf b\f$.

  Assumes \f$\bf A\f$ is positive definite and symmetric of dimension n. 
  Further assumes that
  the routine ::choldc0 is called first so that the Cholesky decomposition of \f$\bf A\f$
  is stored in the lower triangle of a. 

  The solution vector is returned in x.

  This is a modification of the cholsl routine in Numerical Recipes in C changed
  so that the matrix starts at index 0 and such that the diagonal of \f$\bf L\f$ 
  is stored in a and not as a separate vector.
*/  
void cholsl0(double **a, int n, double b[], double x[])
{
  double *w;

  w = CALLOC(n,double);
  cholll0(a,n,b,w);
  chollTl0(a,n,w,x);
  FREE(w);
}
/*!
  \author Geir Storvik
  \brief Solves the set of \f$n\f$ linear equations \f$\bf A\cdot \bf x=\bf b\f$.

  Assumes \f$\bf A\f$ is positive definite and symmetric of dimension n. 
  Further assumes that
  the routine ::choldc0 is called first so that the Cholesky decomposition of \f$\bf A\f$
  is stored in the lower triangle of a. 

  The solution vector is returned in x.

  This is a modification of the cholsl routine in Numerical Recipes in C changed
  so that the matrix starts at index 0 and such that the diagonal of \f$\bf L\f$ 
  is stored in a and not as a separate vector.
*/  
void cholsl0_old(double **a, int n, double b[], double x[])
{
	int i,k;
	double sum;

	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
		x[i]=sum/a[i][i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
		x[i]=sum/a[i][i];
	}
}


/*!
  \author Geir Storvik
  \brief Calculates \f$L^{-1}b\f$ where \f$L\f$ is the Cholesky of \f$A\f$

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The solution vector is returned in x.
*/
void cholll0(double **a, int n, double b[], double x[])
{
  int i,k;
  double sum;
  
  for(i=0;i<n;i++)
    {
      sum=b[i];
      for(k=0;k<i;k++)
	sum -= a[i][k]*x[k];
      x[i]=sum/a[i][i];
    }
}
/*!
  \author Geir Storvik
  \brief Calculates \f$L^{-1}b\f$ where \f$L\f$ is the Cholesky of \f$A\f$

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The solution vector is returned in x.
*/
void cholll0_old(double **a, int n, double b[], double x[])
{
	int i,k;
	double sum;

	for (i=n-1;i>=0;i--) {
		for (sum=b[i],k=i+1;k<n;k++) 
		  sum -= a[k][i]*x[k];
		x[i]=sum/a[i][i];
	}
}

/*!
  \author Geir Storvik
  \brief Inverts a symmetric positive matrix

  Also calculates the determinant of the matrix.

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The inverse matrix is returned in a_inv.
*/
int cholinv0(double **a,int n,double **a_inv,double *det)
{
 int     i;
 double *b;

 b = (double *) Mmatrix_1d(0,n-1,sizeof(double),1);  // Free ok
 *det = G_ZERO;
 for(i=0;i<n;i++)
   *det += log(a[i][i]);
 *det = exp(G_TWO*(*det));
 b[0] = G_ONE;
 cholsl0(a,n,b,a_inv[0]);
 for(i=1;i<n;i++)
   {
    b[i-1] = G_ZERO;
    b[i] = G_ONE;
    cholsl0(a,n,b,a_inv[i]);
   }
 b = (double *) Fmatrix_1d(&b[0]);
 return(0);
}

/*!
  \author Geir Storvik
  \brief Calculates \f$(L^T)^{-1}b\f$ where \f$L\f$ is the Cholesky of \f$A\f$

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that   the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The solution vector is returned in x.
*/
void chollTl0(double **a, int n, double b[], double x[])
{
  int i,k;
  double sum;
  
  for (i=n-1;i>=0;i--) 
    {
      sum = b[i];
      for(k=i+1;k<n;k++)
	sum -= a[k][i]*x[k];
      x[i]=sum/a[i][i];
    }
}

/*!
  \author Geir Storvik
  \brief Calculates \f$\bf b^T\bf A \bf b\f$ using Cholesky decomposition related to 
  the matrix in Q_cell.
*/
double symssq0(double **a, int n, double b[])
{
  int    i,j;
  double ssq=G_ZERO;

  for(i=0;i<n;i++)
    {
      ssq += b[i]*a[i][i]*b[i];
      for(j=0;j<i;j++)
        ssq += 2*b[i]*a[i][j]*b[j];
    }
  return(ssq);
}

/*!
  \author Geir Storvik
  \brief Calculates \f$\bf b^T\bf A \bf b\f$ using Cholesky decomposition

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The sum of squares is the return value of the function.
*/
double cholssq0(double **a, int n, double b[])
{
  int    i,j;
  double ssq=G_ZERO,lx;

  for(i=0;i<n;i++)
    {
      lx = G_ZERO;
      for(j=i;j<n;j++)
        lx += b[j]*a[j][i];
      ssq += lx*lx;
    }
  return(ssq);
} 

/*!
  \author Geir Storvik
  \brief Calculates \f$\bf b^T\bf A^{-1} \bf b\f$ using Cholesky decomposition

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the Cholesky 
  decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The sum of squares is the return value of the function.
*/
double cholssqinv0(double **a, int n, double b[])
{
  int    i;
  double ssq=G_ZERO, *w;

  w = CALLOC(n,double);
  cholll0(a,n,b,w);
  for(i=0;i<n;i++)
    ssq += w[i]*w[i];
  FREE(w);
  return(ssq);
} 

/*!
  \author Geir Storvik
  \brief Calculates the log-determinant of a positive definite matrix

  Assumes \f$\bf A=\bf L\bf L^T\f$ is positive definite and symmetric. 
  Further assumes that the routine ::choldc0 is called first so that the 
  Cholesky decomposition of \f$\bf A\f$ is stored in the lower triangle of a. 

  The log-determinant is the return value of the function.
*/
double chollogdet0(double **a, int n)
{
  int    i;
  double ldet=G_ZERO;

  for(i=0;i<n;i++)
      ldet += G_TWO*log(a[i][i]);
  return(ldet);
} 

