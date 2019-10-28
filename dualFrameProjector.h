#ifndef _DUAL_FRAME_PROJECTOR_H_
#define _DUAL_FRAME_PROJECTOR_H_

#include <cstddef>
#include "frameProjector.h"

class DualFrameProjector: public FrameProjector
{
    public:
        DualFrameProjector(const size_t nThreads, const size_t w, const size_t h);
        ~DualFrameProjector();
        
        void parallelProject(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u); 
        void parallelProjectS(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u); 
    
    private:
        double **t;
};

#endif
