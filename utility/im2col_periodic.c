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
    size_t num_samples, im_sz[3];
    size_t blocksize,num_pels;
    
    size_t im_id,im_c_id, X_id, X_k_id, X_c_id, X_j_id, k, c, j, i, j_k,i_k, j_per, i_per;
    size_t num_channels;
    size_t num_block_pels = blocksize*blocksize*num_channels;
    
    blocksize = mxGetPr(blocksize_IN)[0];
    sz = mxGetPr(im_sz_IN);
    num_samples = (mxGetDimensions(sample_ids_IN))[0];
    
    im_sz[0] = (size_t) sz[0];
    im_sz[1] = (size_t) sz[1];
    im_sz[2] = (size_t) sz[2];
    num_channels = (size_t)im_sz[2];
    X_OUT = mxCreateDoubleMatrix(blocksize*blocksize*im_sz[2],num_samples, mxREAL);

    sample_ids = mxGetPr(sample_ids_IN);
    im = mxGetPr(im_IN);
    X = mxGetPr(X_OUT);

    for ( k = 0; k < num_samples; k++) {
    	j_k = ((size_t)(sample_ids[k]-1))/im_sz[0];
    	i_k = ((size_t)(sample_ids[k]-1)) % im_sz[0];

    	
    	X_k_id = blocksize*blocksize*num_channels*k;
         //for each channel
        for ( c = 0; c < num_channels;c++){
			//copy block into image, column wise 
			X_c_id = blocksize*blocksize*c;
			im_c_id = im_sz[1]*im_sz[0]*c;
            for ( j = 0; j < blocksize; j++) {
            	j_per = (j_k + j) % im_sz[1];
            	X_j_id = j*blocksize;
            	for ( i = 0; i < blocksize; i++) {
            		i_per = (i_k + i) % im_sz[0];
            		X_id = X_k_id + X_c_id + X_j_id + i;
		            im_id = im_c_id + j_per*im_sz[0] + i_per;
		            X[X_id] = im[im_id];
                }
            }
        }
    }
    
    return;
}

