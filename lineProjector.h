#ifndef _LINE_PROJECTOR_H_
#define _LINE_PROJECTOR_H_

#include <cstddef>
#include "activeSetChain_base.h"

class LineProjector
{
    public:
        LineProjector(const size_t nThreads, const size_t *n);
        ~LineProjector();
        void parallelProject(double **y, double **s, const double *const*W, const double *const*u);
        size_t project(double *y, double *s, const double *W, const double *u, size_t thread);
        size_t *iter;

    protected:
        const size_t nThreads;
        size_t *n; 


        ActiveSetChain_base **lineTV;
};


#endif
