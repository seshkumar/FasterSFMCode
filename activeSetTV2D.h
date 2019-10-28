#ifndef _ACTIVESET_TV_2D_H_
#define _ACTIVESET_TV_2D_H_

#include <vector>
#include "graph4N2D.h"
#include "graph.h"

typedef Graph<double, double, double> GraphType;

class ActiveSetTV2D
{
    public:
        ActiveSetTV2D(Graph4N2D &g);
        ActiveSetTV2D(double **W1, double **W2, double **u, double **y, const size_t *dims);
        void initialize(double **W1, double **W2, double **u, double **y, const size_t *dims);
        ~ActiveSetTV2D();

        void calculateMarginals(void);
        void calculateMarginalsGreedy(void);
        double cutValue(void);
        double cutValue(size_t p, size_t &npi);

        void wpav(void);
        void decouple(void);
        void checkOptimality(void);
        void splitSets(void);
        void optimize(void);
        double frobeniusNorm(void);
        template <class T> double cutValue(T **x);

        std::vector<double> a;
        std::vector<double> c;
        std::vector<double> d;
    
    private:
        size_t w;
        size_t h;

        size_t np;

        double **W1;
        double **W2;
        double **u;
        double **y;

        // used for checking optimality
        double **uD;
        double **W1D;
        double **W2D;

        double **s;
        double **x;

        size_t **A;
        size_t **C;
        size_t **D;
        size_t **S;

        std::vector< std::vector<int> > pList;

        double cutVal;

        GraphType *gr;

        std::vector<double> v;
        std::vector<double> p;

        double *m;
        double *wp;

        int iter;
};

#endif
