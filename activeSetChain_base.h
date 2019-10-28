#ifndef _ACTIVESET_CHAIN_BASE_H_
#define _ACTIVESET_CHAIN_BASE_H_

#include <vector>
#include "graph.h"

typedef Graph<double, double, double> GraphType;

class ActiveSetChain_base
{
    public:
        ActiveSetChain_base(const size_t length);
        ~ActiveSetChain_base();

        void initialize(double *y, double *s, const double *W, const double *u);
        size_t optimize(void);

        double *s;
        double *y;
        size_t *A;

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
        template <class T> double cutValue(T *x);
        template <class T> double cutValue(T *x, const double *u, const double *W);
        void initializeParitionList();
   
        size_t length;

        const double *W;
        const double *u;

        // used for checking optimality
        double *uD;
        double *WD;

        size_t *C;
        size_t *D;
        size_t *S;

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
