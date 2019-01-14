/*!
  \file caa_chol.h
  \brief Declaring external functions defined in caa_chol.c
  \author Geir Storvik
*/
int cholinv0(double **a,int n,double **a_inv,double *det);
int choldc0(double **a, int n);
void cholsl0(double **a, int n, double b[], double x[]);
void cholsl0_old(double **a, int n, double b[], double x[]);
void cholll0(double **a, int n, double b[], double x[]);
void cholll0_old(double **a, int n, double b[], double x[]);
void chollTl0(double **a, int n, double b[], double x[]);
double symssq0(double **a, int n, double b[]);
double cholssq0(double **a, int n, double b[]);
double cholssqinv0(double **a, int n, double b[]);
double chollogdet0(double **a, int n);

