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
    double *sample_ids, *im, *X, *weights,*sz;
    
    
    size_t im_id,im_c_id, X_id, X_k_id, X_c_id,
            X_j_id, k, c, i, j, j_k,i_k, j_per, i_per,blocksize,num_channels,num_block_pels;
    
    mxArray *weights_mat;
    size_t num_samples;
    mwSize im_sz[3];
    
    blocksize = mxGetPr(blocksize_IN)[0];
    sz = mxGetPr(im_sz_IN);
    im_sz[0] = (unsigned int) sz[0];
    im_sz[1] = (unsigned int) sz[1];
    im_sz[2] = (unsigned int) sz[2];

    num_samples = (mxGetDimensions(sample_ids_IN))[0];
    
    im_OUT = mxCreateNumericArray(3, im_sz, mxDOUBLE_CLASS, mxREAL);
    weights_mat = mxCreateNumericArray(3,im_sz, mxDOUBLE_CLASS, mxREAL);
    
    sample_ids = mxGetPr(sample_ids_IN);
    im = mxGetPr(im_OUT);
    weights = mxGetPr(weights_mat);
    X = mxGetPr(X_IN);
    
	num_channels = im_sz[2];
    num_block_pels = blocksize*blocksize*num_channels;

    
     //iterate over all block indices 
    for ( k = 0; k < num_samples; k++) {
    	j_k = ((size_t)sample_ids[k]-1)/im_sz[0];
    	i_k = ((size_t)sample_ids[k]-1) % im_sz[0];
        
//         printf("(size_t)(sample_ids[k]-1) = %u \n", (size_t)(sample_ids[k]-1));
//     	printf("k,i_k,j_k = %u, %u, %u \n", k,i_k,j_k);
//     	printf("im_sz = %u, %u, %u \n", im_sz[0],im_sz[1],im_sz[2]);
    	X_k_id = num_block_pels*k;
         //for each rgb channel (1/num_chanels of the column vector) 
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
//                     printf("X_id, im_id = %u, %u \n", X_id,im_id);
		            im[im_id] = im[im_id] + X[X_id];
		            weights[im_id] = weights[im_id] + 1;
                }
            }
        }
    }
    
//     printf("averageing\n");
    /* avg. in image domain */
    for (c = 0; c < num_channels; c++){
		for (j = 0; j < im_sz[1]; j++) {
			im_id = j*im_sz[0] + c*im_sz[1]*im_sz[0];
		
			  for (i = 0; i < im_sz[0]; i++) {
				    if (weights[im_id + i] > 0) {
				        im[im_id + i] = im[im_id + i]/weights[im_id + i];
				    }
			  }
		}
    }
    return;
}

