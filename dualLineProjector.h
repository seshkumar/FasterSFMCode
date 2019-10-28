#ifndef _DUAL_LINE_PROJECTOR_H_
#define _DUAL_LINE_PROJECTOR_H_

#include <cstddef>
#include "lineProjector.h"

class DualLineProjector: public LineProjector
{
    public:
        DualLineProjector(size_t nThreads, const size_t *n);
        ~DualLineProjector();
        
        void parallelProject(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u); 
        void parallelProjectS(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u); 
	void parallelReflect(double **y, double **s, const double *const*tp, const double *const*u, const double *const*W);
    
    private:
        double **t;
};

#endif
