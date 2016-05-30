/* Ryan Webster, 2015 */


#include "mex.h"
#include <string.h>


/* Input Arguments */

#define	X_IN	 prhs[0]
#define sample_ids_IN  prhs[1]
#define blocksize_IN   prhs[2]
#define im_sz_IN   prhs[3]

/* Output Arguments */

#define	im_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[])
        
{
    double *sample_ids, *im, *X, *weights_pr;
    mxArray *weights;
    mwSize num_samples, sz[3];
    mwIndex i,j, k,c, l, im_offset, X_offset, blocksize,num_pels,data_sz;
    
    blocksize = mxGetPr(blocksize_IN)[0];
    
    for (k=0;k<3;k++){
        sz[k] = (unsigned int) mxGetPr(im_sz_IN)[k];
    }

    num_samples = (mxGetDimensions(sample_ids_IN))[0];
    
    im_OUT = mxCreateNumericArray(3, sz, mxDOUBLE_CLASS, mxREAL);
    weights = mxCreateNumericArray(3,sz, mxDOUBLE_CLASS, mxREAL);
    
    sample_ids = mxGetPr(sample_ids_IN);
    im = mxGetPr(im_OUT);
    weights_pr = (double *)mxGetPr(weights);
    X = mxGetPr(X_IN);
    
    num_pels = blocksize*blocksize*sz[2];
    data_sz = blocksize*sizeof(double);

    /* iterate over samples, numel(X)>numel(x), so we respect its linear memory*/
    for (j = 0; j < num_samples; j++) {
        /* copy blocks from each channel*/
        for (c = 0; c < sz[2];c++){
            /* copy each col seperately*/
            for (k = 0; k < blocksize; k++) {
                im_offset = sample_ids[j] - 1 + k*sz[0] + c*sz[0]*sz[1]; //-1 for MATLAB, <3
                X_offset = num_pels*j + blocksize*k + blocksize*blocksize*c;
                for (l=0; l<blocksize; l++) {
                    im[im_offset + l] = im[im_offset+l] + X[X_offset+l];
                    weights_pr[im_offset + l] = weights_pr[im_offset+l] + 1;
                }
            } 
        }
    }

    /* avg. in image domain */
    for (c = 0; c < sz[2];c++){
        for (j = 0; j < sz[1]; j++) {
            im_offset = j*sz[0] + c*sz[0]*sz[1];
            for (i= 0; i < sz[0]; i++) {
                if (weights_pr[im_offset + i] > 0) {
                    im[im_offset + i] = im[im_offset + i]/weights_pr[im_offset + i];
                }
            }
        }
    }
    return;
}

