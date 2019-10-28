#ifndef _FRAME_PROJECTOR_H_
#define _FRAME_PROJECTOR_H_

#include <cstddef>
#include "activeSet2D_base.h"

class FrameProjector
{
    public:
        FrameProjector(const size_t nThreads, const size_t w, const size_t h);
        ~FrameProjector();
        void parallelProject(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u);
        size_t project(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, size_t thread);
        size_t *iter;

    protected:
        const size_t nThreads;
        const size_t w;
        const size_t h;

        ActiveSet2D_base **frameTV;
};


#endif
