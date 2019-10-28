#include <iostream>
#include "truncFrameProjector.h"

TruncFrameProjector::TruncFrameProjector(size_t nThreads, const size_t w, const size_t h)
              :nThreads(nThreads), w(w), h(h)
{
    iter   = new size_t[nThreads];

    frameTV = new ActiveSetTrunc2D*[nThreads];
    for(int i=0; i < nThreads; ++i)
        frameTV[i] = new ActiveSetTrunc2D(w, h);

}

TruncFrameProjector::~TruncFrameProjector(void)
{
    delete[] iter;

    for(int i=0; i < nThreads; ++i)
        delete frameTV[i];

    delete[] frameTV;
}

// size(y) = size(s) = size(u) = w x h; size(W1) = h x w-1; size(W2) = w * h-1
void TruncFrameProjector::parallelProject(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, double epsilon)
{
    size_t i, j, l;

    size_t * iter = this->iter;

    #pragma omp parallel shared(y, s, W1, W2, u, iter, epsilon) private(i, j) default(none)
    {
        #pragma omp for
        for(l=0; l < nThreads; ++l)
        {
            double **y_l = y+(l*w);
            double **s_l = s+(l*w);
            const double *const*u_l = u+(l*w);
            const double *const*W1_l = W1+(l*h);
            const double *const*W2_l = W2+(l*w);

            iter[l] = project(y_l, s_l, W1_l, W2_l, u_l, epsilon, l);
            for(int i=0; i < w; ++i)
                for(int j=0; j < h; ++j)
                    s_l[i][j] += u_l[i][j];
        }
    }
}

size_t TruncFrameProjector::project(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, double epsilon, size_t threadID)
{
    frameTV[threadID]->initialize(y, s, W1, W2, u, epsilon);
    size_t iter = frameTV[threadID]->optimize();
    return iter;
}
