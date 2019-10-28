#ifndef _TRUNC_DUAL_FRAME_PROJECTOR_H_
#define _TRUNC_DUAL_FRAME_PROJECTOR_H_

#include <cstddef>
#include "truncFrameProjector.h"

class TruncDualFrameProjector: public TruncFrameProjector
{
    public:
        TruncDualFrameProjector(const size_t nThreads, const size_t w, const size_t h);
        ~TruncDualFrameProjector();
        
        void parallelProject(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u, double epsilon); 
        void parallelProjectS(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u, double epsilon); 
    
    private:
        double **t;
};

#endif
