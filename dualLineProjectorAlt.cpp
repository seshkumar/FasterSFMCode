#include <iostream>
#include "dualLineProjectorAlt.h"

DualLineProjectorAlt::DualLineProjectorAlt(size_t nThreads, const size_t *n)
                  :LineProjector(nThreads, n)
{
    t = new double*[nThreads]();
    for(int i=0; i<nThreads; ++i)
        t[i] = new double[n[i]]();

}

DualLineProjectorAlt::~DualLineProjectorAlt(void)
{
    for(int i=0; i<nThreads; ++i)
        delete[] t[i];
    delete[] t;
}

void DualLineProjectorAlt::parallelProject(double **y, double **s, const double *const*tp, const double *const*W, const double *const*u)
{
    double **t   = this->t;
    size_t  *n   = this->n;
    size_t *iter = this->iter;

    size_t i, j;
    #pragma omp parallel shared(y, s, tp, W, u, t, n, iter) private(i, j) default(none)
    {
        #pragma omp for 
        for(i=0; i < nThreads; ++i)
        {
            for(j=0; j < n[i]; ++j)
                t[i][j] = u[i][j] - tp[i][j];
            iter[i] = project(y[i], s[i], W[i], t[i], i);
            for(j=0; j < n[i]; ++j)
                s[i][j] += tp[i][j];
        }
    }
}

void DualLineProjectorAlt::parallelReflect(double **y, double **s, const double * const *tp, const double * const *u, const double * const *W)
{
    double **t = this->t;
    size_t  *n = this->n;
    size_t i, j;
    #pragma omp parallel  shared(x, s, y, u, W, t, n) private(i, j) default(none)
    {
        // Projecting y onto S_B(W+u) is equivalent to projecting (y+u) onto S_B(W)
        #pragma omp for
        for(i=0; i < nThreads; ++i)
        {
            for(j=0; j < n[i]; ++j)
                t[i][j] = u[i][j] - tp[i][j];
            iter[i] = project(y[i], s[i], W[i], t[i], i);

            // displacement from the original vector (t -x - u = y - x)
            for(int j=0; j < n[i]; ++j)
            {
                s[i][j] += tp[i][j];
                //reflection step
                y[i][j] = 2*s[i][j] - tp[i][j];
            }
        }
    }
}
