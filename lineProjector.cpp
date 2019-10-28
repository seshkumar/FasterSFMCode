#include <iostream>
#include "lineProjector.h"

LineProjector::LineProjector(size_t nThreads, const size_t *n)
              :nThreads(nThreads)
{
    this->n = new size_t[nThreads]();
    std::copy(n, n+nThreads, this->n);

    iter   = new size_t[nThreads];
    lineTV = new ActiveSetChain_base*[nThreads];
    for(int i=0; i < nThreads; ++i)
        lineTV[i] = new ActiveSetChain_base(n[i]);

}

LineProjector::~LineProjector(void)
{
    delete[] n;
    delete[] iter;

    for(int i=0; i < nThreads; ++i)
        delete lineTV[i];

    delete[] lineTV;
}

void LineProjector::parallelProject(double **y, double **s, const double *const*W, const double *const*u)
{
    size_t i, j;

    size_t * iter = this->iter;

    #pragma omp parallel shared(y, s, W, u, iter) private(i, j) default(none)
    {
        #pragma omp for
        for(i=0; i < nThreads; ++i)
        {
            iter[i] = project(y[i], s[i], W[i], u[i], i);
            for(int j=0; j < n[i]; ++j)
                s[i][j] += u[i][j];
        }
    }
}

size_t LineProjector::project(double *y, double *s, const double *W, const double *u, size_t threadID)
{
    lineTV[threadID]->initialize(y, s, W, u);
    size_t iter = lineTV[threadID]->optimize();
    return iter;
}
