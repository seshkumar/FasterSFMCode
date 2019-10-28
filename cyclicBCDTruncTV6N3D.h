#ifndef _CYCLIC_BCD_TRUNC_TV_6N3D_H_
#define _CYCLIC_BCD_TRUNC_TV_6N3D_H_

#include <iostream>
#include <cstddef>
#include <vector>
#include "truncDualLineProjector.h"
#include "graph6N3D.h"
//#include "optInfo_mex.h"
#include "optInfo.h"

class node
{
    public:
        node(double v, size_t x, size_t y, size_t z): s(v), i(x), j(y), k(y) {}
        double s;
        size_t i;
        size_t j;
        size_t k;
};

inline bool compare(const node &n1, const node &n2)
{ 
    return (n1.s > n2.s);
}

class CyclicBCDTruncTV6N3D
{
    public:
        // size(W1) = (d*h, w-1) size(W2) = (w*d, h-1) size(W3) = (w*h, d-1) size(u) = (w, h, d)
        CyclicBCDTruncTV6N3D(Graph6N3D &g, double epsilon);
        CyclicBCDTruncTV6N3D(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims, double epsilon);
        ~CyclicBCDTruncTV6N3D();
        void optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);
        void optimize_var_sqrt(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname);

    private:
        void initialize(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims, double epsilon);
        void parallelAccumulate(double ** tb, double ** ta, size_t dir);
        
        // size(x) = (w, d*h) + (h, w*d) + (d, w*h) 
        // size(u) = (w, h, d)
        double duplicate(double  ** x, double *** u);

        void computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> & stats);
        double greedyTV(void);
        double threshold;
        double mean;

        std::vector<node> nodes;

        TruncDualLineProjector *projL1;
        TruncDualLineProjector *projL2;
        TruncDualLineProjector *projL3;
        size_t w;
        size_t d;
        size_t h;
        size_t nThreads;
       
        double epsilon;

        size_t * n;
        const double **W;
        double **u2D;
        double **s;
        double **z;
        double **ta;
        double **tb;
        
        bool ***optimal;
        double ***u;
        double ***y;
        size_t maxIter;

        bool ***C;
};

#endif
