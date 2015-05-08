#include "mex.h"


void  subspace_iter_sv(double*, double*, double*,
                       int*, int*, int*, int*,
                       double*,
                       int*,
                       double*, double*,
                       int*, int*, int*, int* );


int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */
int IONE=1;

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  double *a, *u, *s, *res_sv, *z;
  int    *it_sv;
  double percentage, eps;
  int    p, m, n, l, i, j, version, maxit, lwork, it, n_sv, ierr, lm;

  if(nrhs!=6) 
    mexErrMsgIdAndTxt( "MATLAB:subspace_iter_sv:invalidNumInputs",
                       "Seven inputs required.");
  if(nlhs!=2) 
    mexErrMsgIdAndTxt( "MATLAB:subspace_iter_sv:invalidNumOutputs",
                       "Two outputs required.");

  /* get a pointer to the matrix */
  a = mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  /*  get the scalar input l */
  l = mxGetScalar(prhs[1]);
  lm = l*m;
  
  /*  get the scalar input p */
  p = mxGetScalar(prhs[2]);

  /*  get the scalar input m */
  percentage = mxGetScalar(prhs[3]);

  /*  get the scalar input eps */
  eps = mxGetScalar(prhs[4]);

  /*  get the scalar input maxit */
  maxit = mxGetScalar(prhs[5]);
  
  u      = (double*)malloc(l*m*sizeof(double));
  s      = (double*)malloc(l*sizeof(double));
  res_sv = (double*)malloc(l*sizeof(double));
  it_sv  = (int*)malloc(l*sizeof(int));
  ISEED[0]=0; ISEED[1]=0; ISEED[2]=0; ISEED[3]=1; 
  dlarnv_(&IONE, ISEED, &lm, u);

  subspace_iter_sv_(a, u, s, &m, &n, &l, &p,
                &percentage, &maxit, &eps, 
                res_sv, it_sv, &it, &n_sv, &ierr );

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n_sv, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);

  for(j=0; j<n_sv; j++)
    for(i=0; i<m; i++)
      z[j*m+i] = u[j*m+i];


  /*  set the output pointer to the output matrix */
  plhs[1] = mxCreateDoubleMatrix( (mwSize)n_sv, (mwSize)1, mxREAL);
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[1]);

  for(j=0; j<n_sv; j++)
      z[j] = s[j];

  free(u);
  free(s);
  free(res_sv);
  free(it_sv);

}
