#include <cstring>
#include <cassert>
#include "mex.h"
#include "class_handle.hpp"
#include "activeSetTruncTV1D_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    ActiveSetTruncTV1D *tv;

    int nThreads, length;

    // New
    if (!strcmp("new", cmd)) 
    {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");

        nThreads  = *(mxGetPr(prhs[1]));
        length    = *(mxGetPr(prhs[2]));

        tv = new ActiveSetTruncTV1D(nThreads, length);
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<ActiveSetTruncTV1D>(tv);
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
	mexErrMsgTxt("Second input should be a class instance handle.");

   
    // Delete
    if (!strcmp("delete", cmd)) 
    {
        // Destroy the C++ object
        destroyObject<ActiveSetTruncTV1D>(prhs[1]);
        
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    ActiveSetTruncTV1D *TV_instance= convertMat2Ptr<ActiveSetTruncTV1D>(prhs[1]);
   
    // Matlab Command XXX: [y s A_n iter] = optimize( handle, u, W, epsilon, A, max(A(:)));
    // Matlab Command XXX: [y s A_n iter] = truncatedTVInterface('opt', handle, u, W, epsilon, A, max(A(:)));
    // Call the various class methods
    // optimize 
    // Matlab gives transpose
    if (!strcmp("opt", cmd)) 
    {
        // Check parameters
        if (nrhs < 2)
            mexErrMsgTxt("Optimize: Unexpected arguments.");

        int length   = mxGetM(prhs[2]);
        int nThreads = mxGetN(prhs[2]);

        assert(length   == mxGetM(prhs[3]) + 1);
        assert(nThreads == mxGetN(prhs[3]) );

        double *u  = mxGetPr(prhs[2]);
        double *W  = mxGetPr(prhs[3]);
        double eps = mxGetScalar(prhs[4]);
        int    *A  = (int *)mxGetData(prhs[5]);
        int    np  = mxGetScalar(prhs[6]);

        //plhs[0] = mxCreateDoubleMatrix(length, nThreads, mxREAL);
        //plhs[1] = mxCreateDoubleMatrix(length, nThreads, mxREAL);

        //double *y = mxGetPr(plhs[0]);
        //double *s = mxGetPr(plhs[1]);
        double *y = new double(nThreads*length);
        double *s = new double(nThreads*length);
        
        tv = new ActiveSetTruncTV1D(nThreads, length);

        //TV_instance->initialize(y, s, W, u, A, np, eps);
        tv->initialize(y, s, W, u, A, np, eps);

        //int iter = TV_instance->optimize();
        int iter = tv->optimize();

        for(int i=0; i < nThreads; ++i)
        {
            for(int j=0; j < length; ++j)
            {
                mexPrintf("%f\t", y[i*length + j]);
            }
            mexPrintf("\n");
        }
        mexPrintf("\n");

        for(int i=0; i < nThreads; ++i)
        {
            for(int j=0; j < length; ++j)
            {
                mexPrintf("%f\t", s[i*length + j]);
            }
            mexPrintf("\n");
        }
        mexPrintf("\n");


        mexPrintf("nThreads = %d; length = %d\n", nThreads, length);

        delete tv;
        //mexPrintf("Number of Iterations  = %d \n", iter);

        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
