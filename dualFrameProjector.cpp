#include <iostream>
#include "dualFrameProjector.h"

DualFrameProjector::DualFrameProjector(const size_t nThreads, const size_t w, const size_t h)
                   :FrameProjector(nThreads, w, h)
{
    t    = new double*[nThreads*w]();
    t[0] = new double [nThreads*w*h]();

    for(int i=1; i < nThreads*w; ++i)
        t[i] = t[0] + i*h;
}

DualFrameProjector::~DualFrameProjector(void)
{
    delete[] t[0];
    delete[] t;
}

void DualFrameProjector::parallelProject(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u)
{
    size_t l;
    size_t *iter = this->iter;

    double **t = this->t;
    double **y_l;
    double **s_l;
    double **t_l;

    #pragma omp parallel shared(y, s, tp, W1, W2, t, u, iter) private(l, y_l, s_l, t_l) default(none)
    {
        #pragma omp for
        for(l=0; l < nThreads; ++l)
        {
            y_l  = y  + l*w;
            s_l  = s  + l*w;
            t_l  = t  + l*w;

            const double *const* tp_l = tp + l*w;
            const double *const* u_l  = u  + l*w;
            const double *const* W1_l = W1 + l*h;
            const double *const* W2_l = W2 + l*w;

            for(int i=0; i < w; ++i)
                for(int j=0; j < h; ++j)
                    t_l[i][j] = u_l[i][j] + tp_l[i][j];

            iter[l] = project(y_l, s_l, W1_l, W2_l, t_l, l);

            for(int i=0; i < w; ++i)
                for(int j=0; j < h; ++j)
                    s_l[i][j] += tp_l[i][j];
        }
    }

}

void DualFrameProjector::parallelProjectS(double **y, double **s, const double *const*tp, const double *const*W1, const double *const*W2, const double *const*u)
{
    size_t l;
    size_t * iter = this->iter;

    double **t = this->t;

    #pragma omp parallel shared(y, s, tp, W1, W2, t, u, iter) private(l) default(none)
    {
        #pragma omp for
        for(l=0; l < nThreads; ++l)
        {
            double **y_l = y+(l*w);
            double **s_l = s+(l*w);
            double **t_l = t+(l*w);
            const double *const*u_l = u+(l*w);
            const double *const*tp_l = tp+(l*w);
            const double *const*W1_l = W1+(l*h);
            const double *const*W2_l = W2+(l*w);

            for(int i=0; i < w; ++i)
                for(int j=0; j < h; ++j)
                    t_l[i][j] = u_l[i][j] + tp_l[i][j];

            iter[l] = project(y_l, s_l, W1_l, W2_l, t_l, l);

            for(int i=0; i < w; ++i)
                for(int j=0; j < h; ++j)
                    s_l[i][j] += t_l[i][j];
        }
    }
}
