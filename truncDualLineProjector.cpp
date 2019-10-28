#include <iostream>
#include "truncDualLineProjector.h"

TruncDualLineProjector::TruncDualLineProjector(size_t nThreads, const size_t *n)
                       :TruncLineProjector(nThreads, n)
{
    t = new double*[nThreads]();
    for(int i=0; i<nThreads; ++i)
        t[i] = new double[n[i]]();

}

TruncDualLineProjector::~TruncDualLineProjector(void)
{
    for(int i=0; i<nThreads; ++i)
        delete[] t[i];
    delete[] t;
}

void TruncDualLineProjector::parallelProject(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u, double epsilon)
{
    double **t   = this->t;
    size_t  *n   = this->n;
    size_t *iter = this->iter;

    size_t i, j;
    #pragma omp parallel shared(y, s, tp, W, u, t, n, epsilon, iter) private(i, j) default(none)
    {
        #pragma omp for 
        for(i=0; i < nThreads; ++i)
        {
            for(j=0; j < n[i]; ++j)
                t[i][j] = u[i][j] + tp[i][j];
            iter[i] = project(y[i], s[i], W[i], t[i], epsilon, i);
            for(j=0; j < n[i]; ++j)
                s[i][j] += tp[i][j];
        }
    }
}
