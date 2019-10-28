#ifndef _DUAL_LINE_PROJECTOR_ALT_H_
#define _DUAL_LINE_PROJECTOR_ALT_H_

#include <cstddef>
#include "lineProjector.h"

class DualLineProjectorAlt: public LineProjector
{
    public:
        DualLineProjectorAlt(size_t nThreads, const size_t *n);
        ~DualLineProjectorAlt();
        
        void parallelProject(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u); 
	void parallelReflect(double **y, double **s, const double *const*tp, const double *const*u, const double *const*W);
    
    private:
        double **t;
};

#endif
