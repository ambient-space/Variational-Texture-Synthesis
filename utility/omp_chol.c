/**************************************************************************
 *
 * File name: ompmex.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 18.8.2009
 *
 *************************************************************************/

#include "omp_utils.h"
#include <math.h>
#include <string.h>


/* Input Arguments */


#define IN_DtX        prhs[0]
#define IN_G          prhs[1]
#define IN_T          prhs[2]

/* Output Arguments */

#define	GAMMA_OUT     plhs[0]


/***************************************************************************************/




mxArray * omp_chol_core(double DtX[], double G[], size_t m, size_t L,int T)
{
    mxArray *Gamma;
    size_t i, j, signum, pos, *ind, gamma_count, k,  *gammaIr, *gammaJc;
    int *selected_atoms;
    double *alpha, *r, *Lchol, *c, *Gsub, *Dsub, sum, *gammaPr, *tempvec1, *tempvec2;
            
    Gamma = mxCreateSparse(m, L, L*T, mxREAL);
    gammaPr = mxGetPr(Gamma);
    gammaIr = mxGetIr(Gamma);
    gammaJc = mxGetJc(Gamma);
    gamma_count = 0;
    gammaJc[0] = 0;
    
    /*** helper arrays ***/
    
    alpha = (double*)malloc(m*sizeof(double));        /* contains D'*residual */
    ind = (size_t*)malloc(T*sizeof(size_t));        /* indices of selected atoms */
    selected_atoms = (int*)malloc(m*sizeof(int));     /* binary array with 1's for selected atoms */
    c = (double*)calloc(T,sizeof(double));            /* orthogonal projection result */
    
    /* Cholesky decomposition of D_I'*D_I */
    Lchol = (double*)calloc(T*T,sizeof(double));
    
    /* temporary vectors for various computations */
    tempvec1 = (double*)calloc(m,sizeof(double));
    tempvec2 = (double*)malloc(m*sizeof(double));
    
    /* matrix containing G(:,ind) - the columns of G corresponding to the selected atoms, in order of selection */
    Gsub = (double*)malloc(m*T*sizeof(double));
    
    
    /**********************   perform omp for each signal   **********************/
    
    for (signum=0; signum<L; ++signum) {
        
        /* initialize residual norm and deltaprev for error-omp */
        
        
        if (T>0) {
            
            /* initialize alpha := DtX */
            
            memcpy(alpha, DtX + m*signum, m*sizeof(double));
            
            /* mark all atoms as unselected */
            
            for (i=0; i<m; ++i) {
                selected_atoms[i] = 0;
            }
            
        }
        
        /* main loop */
        
        i=0;
        while (i<T) {
            
            
            /* index of next atom */
            
            pos = maxabs(alpha, m);
            
            /* stop criterion: selected same atom twice, or inner product too small */
            
            if (selected_atoms[pos] || alpha[pos]*alpha[pos]<1e-14) {
                
                break;
            }
            
            
            /* mark selected atom */
            ind[i] = pos;
            selected_atoms[pos] = 1;
            
            /* append column to Gsub or Dsub */
            
            
            memcpy(Gsub+i*m, G+pos*m, m*sizeof(double));
            
            /*** Cholesky update ***/
            
            if (i==0) {
                *Lchol = 1;
            }
            else {
                
                /* incremental Cholesky decomposition: compute next row of Lchol */
                
                vec_assign(tempvec1, Gsub+i*m, ind, i);          /* extract tempvec1 := Gsub(ind,i) */
                
                backsubst('L', Lchol, tempvec1, tempvec2, T, i);   /* compute tempvec2 = Lchol \ tempvec1 */
                
                for (j=0; j<i; ++j) {                              /* write tempvec2 to end of Lchol */
                    Lchol[j*T+i] = tempvec2[j];
                }
                
                /* compute Lchol(i,i) */
                sum = 0;
                
                for (j=0; j<i; ++j) {         /* compute sum of squares of last row without Lchol(i,i) */
                    sum += SQR(Lchol[j*T+i]);
                }
                if ( (1-sum) <= 1e-14 ) {     /* Lchol(i,i) is zero => selected atoms are dependent */
                    break;
                }
                
                Lchol[i*T+i] = sqrt(1-sum);
            }
            
            
            
            i++;
            /* perform orthogonal projection and compute sparse coefficients */
            
            vec_assign(tempvec1, DtX + m*signum, ind, i);   /* extract tempvec1 = DtX(ind) */
            
            cholsolve('L', Lchol, tempvec1, c, T, i); /* solve LL'c = tempvec1 for c */
            
            /* update alpha = D'*residual */
            mat_vec(1, Gsub, c, tempvec1, m, i); /* compute tempvec1 := Gsub*c */
            
            memcpy(alpha, DtX + m*signum, m*sizeof(double));    /* set alpha = D'*x */
            
            vec_sum(-1, tempvec1, alpha, m);  /* compute alpha := alpha - tempvec1 */
            
            
        }
        
        /*** generate output vector gamma ***/
        
        /* sort the coefs by index before writing them to gamma */
        quicksort(ind,c,i);

        /* append coefs to gamma and update the indices */
        for (j=0; j<i; ++j) {
            gammaPr[gamma_count] = c[j];
            gammaIr[gamma_count] = ind[j];
            gamma_count++;
        }
        gammaJc[signum+1] = gammaJc[signum] + i;
    }
    
    /* free memory */
    
    free(Gsub);
    free(tempvec2);
    free(tempvec1);
    free(Lchol);
    free(c);
    free(selected_atoms);
    free(ind);
    free(alpha);
    
    return Gamma;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

{
    double *DtX, *G, *tmp_gamma;
    int gmode, profile, T;
    mwSize m, n, L;   /* D is n x m , X is n x L, DtX is m x L */
    
    
    /* get parameters */
    
    DtX = G = 0;
    
    DtX = mxGetPr(IN_DtX);
    
    G = mxGetPr(IN_G);
    
    T = (int)(mxGetScalar(IN_T)+1e-2);

    m = mxGetM(IN_DtX);
    L = mxGetN(IN_DtX);
    
    if (mxGetN(IN_G)!=mxGetM(IN_G)) {
        mexErrMsgTxt("G must be a square matrix.");
    }
    if (mxGetN(IN_G) != m) {
        mexErrMsgTxt("DtX and G have incompatible sizes.");
    }
    
    GAMMA_OUT = omp_chol_core(DtX,G,m, L, T);
    
    return;
}

