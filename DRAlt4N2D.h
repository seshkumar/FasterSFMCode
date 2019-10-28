#ifndef _DR_ALT_4N2D_H_
#define _DR_ALT_4N2D_H_

#include <iostream>
#include <cstddef>
#include <vector>
#include "dualLineProjector.h"
#include "graph4N2D.h"
#include "optInfo_mex.h"
//#include "optInfo.h"

class node
{
    public:
        node(double v, size_t x, size_t y): s(v), i(x), j(y){}
        double s;
        size_t i;
        size_t j;
};

inline bool compare(const node &n1, const node &n2)
{ 
    return (n1.s > n2.s);
}

class DRAlt4N2D
{
    public:
        // size(W1) = (h, w-1) size(W2) = (w, h-1) size(u) = (w, h)
        DRAlt4N2D(Graph4N2D &g);
        DRAlt4N2D(double **y, double **s, double **W1, double **W2, double **u, size_t *dims);
        ~DRAlt4N2D();
        void optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);

    private:
        void initialize(double **y, double **s, double **W1, double **W2, double **u, size_t *dims);
        void parallelAccumulate(size_t dir);
        
        // size(x) = (w, d*h) + (h, w*d) + (d, w*h) 
        // size(u) = (w, h, d)
        double duplicate(double  ** x, double ** u);

        void computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> & stats);
        double greedyTV(void);
        double threshold;
        double mean;

        std::vector<node> nodes;

        DualLineProjector *projL1;
        DualLineProjector *projL2;
        size_t w;
        size_t h;
        size_t nThreads;
        
        size_t * n;
        const double **W;
        double **u2D;
        double **z;
        double **ta;
        double **tb;
        
        bool **optimal;
        double **u;
        double **s;
        double **y;
        size_t maxIter;

        bool **C;
};

#endif
