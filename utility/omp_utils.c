#include "omp_utils.h"
#include <ctype.h>
#include <math.h>


const char FULL_GAMMA_STR[] = "full";
const char SPARSE_GAMMA_STR[] = "sparse";

/* quicksort, public-domain C implementation by Darel Rex Finley. */
/* modification: sorts the array data[] as well, according to the values in the array vals[] */

#define  MAX_LEVELS  300

void quicksort(mwIndex vals[], double data[], mwIndex n) {
  
  long piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
  double datapiv;
  
  beg[0]=0;
  end[0]=n;
  
  while (i>=0) {
    
    L=beg[i]; 
    R=end[i]-1;
    
    if (L<R) {
      
      piv=vals[L];
      datapiv=data[L];
      
      while (L<R) 
      {
        while (vals[R]>=piv && L<R) 
          R--;
        if (L<R) {
          vals[L]=vals[R];
          data[L++]=data[R];
        }
        
        while (vals[L]<=piv && L<R) 
          L++;
        if (L<R) { 
          vals[R]=vals[L];
          data[R--]=data[L];
        }
      }
      
      vals[L]=piv;
      data[L]=datapiv;
      
      beg[i+1]=L+1;
      end[i+1]=end[i];
      end[i++]=L;
      
      if (end[i]-beg[i] > end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
      }
    }
    else {
      i--;
    }
  }
}


/* find maximum of absolute values */

mwIndex maxabs(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double absval, maxval = SQR(*c);   /* use square which is quicker than absolute value */

  for (k=1; k<m; ++k) {
    absval = SQR(c[k]);
    if (absval > maxval) {
      maxval = absval;
      maxid = k;
    }
  }
  return maxid;
}


/* compute y := alpha*x + y */

void vec_sum(double alpha, double x[], double y[], mwSize n)
{
  mwIndex i;

  for (i=0; i<n; ++i) {
    y[i] += alpha*x[i];
  }
}


/* compute y := alpha*A*x */

void mat_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m)
{
  mwIndex i, j, i_n;
  double *Ax;

  Ax = mxCalloc(n,sizeof(double));

  for (i=0; i<m; ++i) {
    i_n = i*n;
    for (j=0; j<n; ++j) {
      Ax[j] += A[i_n+j] * x[i];
    }
  }

  for (j=0; j<n; ++j) {
    y[j] = alpha*Ax[j];
  }

  mxFree(Ax);
}

/* find maximum of vector */

mwIndex maxpos(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double val, maxval = *c;

  for (k=1; k<m; ++k) {
    val = c[k];
    if (val > maxval) {
      maxval = val;
      maxid = k;
    }
  }
  return maxid;
}


/* solve L*x = b */

void backsubst_L(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= L[j*n+i]*x[j];
    }
    x[i] = rhs/L[i*n+i];
  }
}


/* solve L'*x = b */

void backsubst_Lt(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= L[(i-1)*n+j]*x[j];
    }
    x[i-1] = rhs/L[(i-1)*n+i-1];
  }
}


/* solve U*x = b */

void backsubst_U(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= U[j*n+i-1]*x[j];
    }
    x[i-1] = rhs/U[(i-1)*n+i-1];
  }
}


/* solve U'*x = b */

void backsubst_Ut(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= U[i*n+j]*x[j];
    }
    x[i] = rhs/U[i*n+i];
  }
}


/* back substitution solver */

void backsubst(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  if (tolower(ul) == 'u') {
    backsubst_U(A, b, x, n, k);
  }
  else if (tolower(ul) == 'l') {
    backsubst_L(A, b, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be ''U'' or ''L''");
  }
}


/* solve equation set using cholesky decomposition */

void cholsolve(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  double *tmp;

  tmp = mxMalloc(k*sizeof(double));

  if (tolower(ul) == 'l') {
    backsubst_L(A, b, tmp, n, k);
    backsubst_Lt(A, tmp, x, n, k);
  }
  else if (tolower(ul) == 'u') {
    backsubst_Ut(A, b, tmp, n, k);
    backsubst_U(A, tmp, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be either ''U'' or ''L''");
  }

  mxFree(tmp);
}


/* perform a permutation assignment y := x(ind(1:k)) */

void vec_assign(double y[], double x[], mwIndex ind[], mwSize k)
{
  mwIndex i;

  for (i=0; i<k; ++i)
    y[i] = x[ind[i]];
}



