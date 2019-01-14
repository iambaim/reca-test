float *vector();
float **matrix();
float **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
float **submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();
void gauleg(double, double, double *, double *, int);
void gauher(double x[], double w[], int n);
void choldc(double **a, int n, double p[]);
int powell(double p[], double **xi, int n, double ftol, int *iter, 
            double *fret,double (*func)(double []));
void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk);
double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	     double *xmin);
void bessik(double x, double xnu, double *ri, double *rk, double *rip,
	    double *rkp);
double gammln(double xx);
double d_erff(double x);
double erfcc(double x);
double rtflsp(double (*func)(double), double x1, double x2, double xacc);
void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	   double *prob,long *ind);
double probks(double alam);
void sort(unsigned long n, double arr[]);


