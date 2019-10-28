#ifndef _ACTIVESET_TRUNC_CHAIN_H_
#define _ACTIVESET_TRUNC_CHAIN_H_

#include <cmath>
#include <vector>
#include "activeSetChain_base.h"
#include "graph.h"

typedef Graph<double, double, double> GraphType;

inline double cval(size_t *x, const double *u, const double *W, size_t length)
{
    double val = 0.0f;
    for(int j=0; j < length; ++j)
    {
        val -= u[j]*x[j];
        if(j < length-1)
            val+= W[j]*abs(int(x[j] - x[j+1]));
    }
            
    return val;
};

class ActiveSetTruncChain:public ActiveSetChain_base
{
    public:
	ActiveSetTruncChain(const size_t length);
        void initialize(double *y, double *s, const double *W, const double *u, double epsilon);
        int optimize();
        ~ActiveSetTruncChain();

   private:
        void   truncateWeights(const double *W, const double *u);
        double getPrimalDual(size_t *p, double *d, const double *W, const double *u, double lambda);
        double epsilon;
        double primal1, dual1;
        double primal2, dual2;
    
        GraphType *g;

        size_t *A1;
        size_t *A2;
        size_t *V;

        double *u_m;
        double *W_m;

        double *u_n;
        double *W_n;

        double *s1;
        double *s2;
};

#endif
