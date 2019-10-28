#ifndef _ACTIVESET_2D_BASE_H_
#define _ACTIVESET_2D_BASE_H_

#include <vector>
#include "graph4N2D.h"
#include "graph.h"
#include "optInfo.h"

typedef Graph<double, double, double> GraphType;

class ActiveSet2D_base
{
    public:
        // size(W1) = (h, w-1) size(W2) = (w, h-1) size(u) = (w, h)
        // dims[0] = w; dims[1] = h;
        ActiveSet2D_base(const size_t w, const size_t h);
        ~ActiveSet2D_base();

        void initialize(double **y, double **s, const double *const*W1, const double *const*W2, const double *const*u);
        size_t optimize(void);

        double **s;
        double **y;
        size_t **A;

        int iter;

        std::vector<double> a;
        std::vector<double> c;
        std::vector<double> d;
    
    protected:
        void calculateMarginals(void);
        void calculateMarginalsGreedy(void);
        void cleanup(void);
        double cutValue(void);
        double cutValue(size_t p, size_t &npi);

        void wpav(void);
        void decouple(void);
        void checkOptimality(void);
        void splitSets(void);
        double frobeniusNorm(void);
        template <class T> double cutValue(T **x);
        template <class T> double cutValue(T **x, const double *const*u, const double *const*W1, const double *const*W2);
        void initializeParitionList();
   
        size_t w;
        size_t h;

        const double *const*u;
        const double *const*W1;
        const double *const*W2;

        // used for checking optimality
        double **uD;
        double **W1D;
        double **W2D;

        size_t **C;
        size_t **D;
        size_t **S;

        size_t np;

        std::vector< std::vector<int> > pList;

        double cutVal;

        GraphType *gr;

        std::vector<double> v;
        std::vector<double> p;

        double *m;
        double *wp;
};

#endif
