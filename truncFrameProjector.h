#ifndef _TRUNC_FRAME_PROJECTOR_H_
#define _TRUNC_FRAME_PROJECTOR_H_

#include <cstddef>
#include "activeSetTrunc2D.h"

class TruncFrameProjector
{
    public:
        TruncFrameProjector(const size_t nThreads, const size_t w, const size_t h);
        ~TruncFrameProjector();
        void parallelProject(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, double epsilon);
        size_t project(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, double epsilon, size_t thread);
        size_t *iter;

    protected:
        const size_t nThreads;
        const size_t w;
        const size_t h;

        ActiveSetTrunc2D **frameTV;
};


#endif
