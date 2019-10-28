#ifndef _DR_OPTIMIZER_4N2D_H_
#define _DR_OPTIMIZER_4N2D_H_

#include <iostream>
#include <cstddef>
#include <vector>
#include "dualLineProjector.h"
#include "graph4N2D.h"
#include "optInfo.h"

class node
{
    public:
        node(double v, size_t x, size_t y): s(v), i(x), j(y) {}
        double s;
        size_t i;
        size_t j;
};

inline bool compare(const node &n1, const node &n2)
{ 
    return (n1.s > n2.s);
}

class DROptimizer4N2D
{
    public:
        // size(W1) = (h, w-1) size(W2) = (w, h-1) size(u) = (w, h)
        DROptimizer4N2D(Graph4N2D &g);
        DROptimizer4N2D(double **W1, double **W2, double **u, double **y, const size_t *dims);
        ~DROptimizer4N2D();
        void optimize(const size_t maxIterations, std::vector<OptInfo> & stats);
        void optimize(double ** start, const size_t maxIterations, std::vector<OptInfo> & stats);

    private:
        void initialize(double **W1, double **W2, double **u, double **y, const size_t *dims);
        void parallelReflectLambda(double ** tb, double ** ta);

        // size(x) = (h, w) + (w, h) 
        // size(u) = (w, h)
        double duplicate(double  ** x, double ** u);

        void   computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> & stats);
        double greedyTV(void);

        std::vector<node> nodes;

        DualLineProjector *projL;
        size_t w;
        size_t h;
        size_t nThreads;
        
        size_t * n;
        const double **W;
        double **u2D;
        double **s;
        double **ta;
        double **tb;

        double **u;
        double **y;
        size_t maxIter;

        bool **C;
};

#endif
