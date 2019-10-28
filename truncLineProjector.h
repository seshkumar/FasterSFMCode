#ifndef _TRUNC_LINE_PROJECTOR_H_
#define _TRUNC_LINE_PROJECTOR_H_

#include <cstddef>
#include "activeSetTruncChain.h"

class TruncLineProjector
{
    public:
        TruncLineProjector(const size_t nThreads, const size_t *n);
        ~TruncLineProjector();
        void parallelProject(double **y, double **s, const double *const*W, const double *const*u, double epsilon);
        size_t project(double *y, double *s, const double *W, const double *u, double epsilon, size_t thread);
        size_t *iter;

    protected:
        const size_t nThreads;
        size_t *n; 


        ActiveSetTruncChain **lineTV;
};


#endif
