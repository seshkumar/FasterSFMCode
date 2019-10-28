#ifndef _CYCLIC_BCD_TRUNC_TV_4N2D_H_
#define _CYCLIC_BCD_TRUNC_TV_4N2D_H_

#include <iostream>
#include <cstddef>
#include <vector>
#include "truncDualLineProjector.h"
#include "graph4N2D.h"
#include "optInfo.h"

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

class CyclicBCDTruncTV4N2D
{
    public:
        // size(W1) = (h, w-1) size(W2) = (w, h-1) size(u) = (w, h)
        CyclicBCDTruncTV4N2D(Graph4N2D &g, double epsilon);
        CyclicBCDTruncTV4N2D(double **y, double **s, double **W1, double **W2, double **u, const size_t *dims, double epsilon);
        ~CyclicBCDTruncTV4N2D();
        void optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var_sqrt(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_acc(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var_acc(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var_sqrt_acc(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        //void optimize(double ** start, const size_t maxIterations, std::vector<OptInfo> & stats);

    private:
        void initialize(double **y, double **s, double **W1, double **W2, double **u, const size_t *dims, double epsilon);
        void parallelAccumulate(size_t dir, double beta = 0.0f);
        
        // size(x) = (w, d*h) + (h, w*d) + (d, w*h) 
        // size(u) = (w, h, d)
        double duplicate(double  ** x, double ** u);

        void computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> & stats);
        double greedyTV(void);
        double threshold;
        double mean;

        std::vector<node> nodes;

        TruncDualLineProjector *projL1;
        TruncDualLineProjector *projL2;
        size_t w;
        size_t h;
        size_t nThreads;
        
        size_t * n;
        const double **W;
        double **u2D;
        double **z;
        double **ta;
        double **tb;

        double epsilon;
        
        bool **optimal;
        double **u;
        double **s;
        double **y;
        size_t maxIter;

        bool **C;
};

#endif
