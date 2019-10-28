#ifndef _ACTIVESET_TRUNC_2D_H_
#define _ACTIVESET_TRUNC_2D_H_

#include <cmath>
#include <vector>
#include "activeSet2D_base.h"
#include "graph.h"

typedef Graph<double, double, double> GraphType;

inline double cval(size_t **x, const double *const*u, const double *const* W1, const double *const* W2, size_t w, size_t h)
{
    double val = 0.0f;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            val -= u[i][j]*x[i][j];
            if(i < w-1)
                val+= W1[j][i]*abs(int(x[i][j] - x[i+1][j]));
            if(j < h-1)
                val+= W2[i][j]*abs(int(x[i][j] - x[i][j+1]));
        }
            
    return val;
};

class ActiveSetTrunc2D:public ActiveSet2D_base
{
    public:
	ActiveSetTrunc2D(const size_t w, const size_t h);
        void initialize(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u, double epsilon);
        int optimize();
        ~ActiveSetTrunc2D();

   private:
        void   truncateWeights(const double *const* W1, const double *const* W2, const double *const* u);
        double getPrimalDual(size_t **p, double **d, const double *const* W1, const double *const* W2, const double *const* u, double lambda);
        double truncatedTVNorm(void);
        double epsilon;
        double primal1, dual1;
        double primal2, dual2;
    
        GraphType *g;

        size_t **A1;
        size_t **A2;
        size_t **V;

        double **u_m;
        double **W1_m;
        double **W2_m;

        double **u_n;
        double **W1_n;
        double **W2_n;

        double **s1;
        double **s2;
};

#endif
