/*Ryan Webster, 2015
 *
 *
 */
#include "mex.h"
#include <string.h>


/* Input Arguments */

#define	im_IN	 prhs[0]
#define sample_ids_IN  prhs[1]
#define blocksize_IN   prhs[2]
#define	im_sz_IN	 prhs[3]

/* Output Arguments */

#define	X_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[], 
		             int nrhs, const mxArray*prhs[])
     
{ 
    double *sample_ids,*im, *X,*sz;
    mwSize num_samples;
    mwIndex i,j, k, l, ind, blocksize,num_pels,data_sz;
    
    blocksize = mxGetPr(blocksize_IN)[0];
    sz = mxGetPr(im_sz_IN);
    num_samples = (mxGetDimensions(sample_ids_IN))[0];

    X_OUT = mxCreateDoubleMatrix(blocksize*blocksize*sz[2],num_samples, mxREAL);

    
    sample_ids = mxGetPr(sample_ids_IN);
    im = mxGetPr(im_IN);
    X = mxGetPr(X_OUT);

    
    num_pels = blocksize*blocksize*sz[2];
    data_sz = blocksize*sizeof(double);
    
    
    /* iterate over all blocks */
    for (j = 0; j < num_samples; j++) {
        /* copy single block */
        for (i = 0; i < sz[2]; i++){
            for (k = 0; k < blocksize; k++) {
                ind = sample_ids[j] - 1 + i*sz[0]*sz[1] + k*sz[0];
                memcpy(X + num_pels*j + blocksize*blocksize*i + blocksize*k , im + ind, data_sz);
            }
        }
    }

    return;
}

