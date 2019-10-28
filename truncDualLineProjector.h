#ifndef _TRUNC_DUAL_LINE_PROJECTOR_H_
#define _TRUNC_DUAL_LINE_PROJECTOR_H_

#include <cstddef>
#include "truncLineProjector.h"

class TruncDualLineProjector: public TruncLineProjector
{
    public:
        TruncDualLineProjector(size_t nThreads, const size_t *n);
        ~TruncDualLineProjector();
        
        void parallelProject(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u, double epsilon); 
    
    private:
        double **t;
};

#endif
