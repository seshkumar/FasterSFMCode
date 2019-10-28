#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "activeSetTruncChain.h"
#include "contraction.h"
#include "restriction.h"

ActiveSetTruncChain::ActiveSetTruncChain(const size_t length)
                    :ActiveSetChain_base(length)
{
    g = new GraphType(length, 3*length);

    g->add_node(length);

    for(size_t j=0; j < length; ++j)
    {
        g->add_tweights(j, 0.0f, 0.0f);
        if(j < length-1)
            g->add_edge(j, j+1, 0.0f, 0.0f);
    } 

    A1 = new size_t[length];
    A2 = new size_t[length];
    V  = new size_t[length];

    s1 = new double[length];
    s2 = new double[length];

    u_m = new double[length];
    W_m = new double[length-1];

    u_n = new double[length];
    W_n = new double[length-1];

}

void ActiveSetTruncChain::initialize(double *y, double *s, const double *W, const double *u, double eps)
{
    for(size_t i=0; i < length; ++i)
        V[i] = 1;

    std::memset(A1, 0, length*sizeof(size_t));
    std::memset(A2, 0, length*sizeof(size_t));

    std::memset(s1, 0, length*sizeof(double));
    std::memset(s2, 0, length*sizeof(double));

    std::memset(u_m, 0, length*sizeof(double));
    std::memset(W_m, 0, (length-1)*sizeof(double));

    std::memset(u_n, 0, length*sizeof(double));
    std::memset(W_n, 0, (length-1)*sizeof(double));

    std::memset(A, 0, length*sizeof(size_t));
    std::memset(C, 0, length*sizeof(size_t));
    std::memset(D, 0, length*sizeof(size_t));
    std::memset(S, 0, length*sizeof(size_t));

    std::memset(uD, 0, length*sizeof(double));
    std::memset(WD, 0, (length-1)*sizeof(double));

    this->epsilon = eps;

    this->y  = y;
    this->s  = s;

    np = 1;

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    initializeParitionList();
    truncateWeights(W, u);
    splitSets();

    this->W  = W_n;
    this->u  = u_n;

    calculateMarginalsGreedy();
}

void ActiveSetTruncChain::truncateWeights(const double *W, const double *u)
{    
    std::memset(A1, 0, length*sizeof(size_t));
    std::memset(A2, 0, length*sizeof(size_t));
    std::memset(C , 0, length*sizeof(size_t));
    std::memset(D , 0, length*sizeof(size_t));

    // A1 is equivalent to A_+
    // A2 is equivalent to A_-
    // A_+ \subset A_-
    // A1 \subset A2

    dual1 = getPrimalDual(A1, s1, W, u, epsilon);
    dual2 = getPrimalDual(A2, s2, W, u, -epsilon);

    primal1 = cval(A1, u, W, length);
    primal2 = cval(A2, u, W, length);

    assert(fabs(primal1 - dual1) < 1e-1);
    assert(fabs(primal2 - dual2) < 1e-1);

    //std::cout << "cut1 : Primal = " << primal1 << "; Dual =  " << dual1 << "; Gap = " << primal1 - dual1 << std::endl;
    //std::cout << "cut2 : Primal = " << primal2 << "; Dual =  " << dual2 << "; Gap = " << primal2 - dual2 << std::endl;

    for(size_t j=0; j < length; ++j)
    {
        size_t k = A[j];

       if(A2[j] - A1[j] == 1)
        {
            C[j] = k;
            D[j] = 1;
            c[k] += (j+1);
        }
        else
        {
            C[j] = np+2;
        }

        a[k] += (j+1);
    }

    // F: 2^{A2 - A1} \to R
    restriction(u_m, W_m, u, W, A2,  V, length);
    contraction(u_n, W_n, u_m, W_m, A1, A2, length);
}



double ActiveSetTruncChain::getPrimalDual(size_t *primal, double *dual, const double *W, const double *u, double lambda)
{
    double dualCost = 0.0f;
    GraphType::arc_id a = g->get_first_arc();
    for(size_t j=0; j < length; ++j)
    {
        g->set_trcap(j, - u[j] + lambda);
        if( j < length - 1)
        {
            g->set_rcap(a, W[j]);
            a = g->get_next_arc(a);

            g->set_rcap(a, W[j]);
            a = g->get_next_arc(a);
        }
    }

    g->maxflow();

    for(size_t j=0; j < length; ++j)
    {
        dual[j] = g->get_trcap(j) - lambda;
        if(g->what_segment(j))
        {
            primal[j] = 1;
            dualCost += dual[j];
        }
    }
    return dualCost;
}
int ActiveSetTruncChain::optimize(void)
{
    iter = 1;
    while(true)
    {
        wpav();
        decouple();
        checkOptimality();
        //printf("Iter = %d; Partitions= %ld; TV = %e; epsilon= %e\n", iter, np, cutValue()+frobeniusNorm(), cutVal);
        if(cutVal > -1e-6)
            break;
        splitSets();
        calculateMarginalsGreedy();
        iter++;
    }

    for(size_t j=0; j < length; ++j)
    {
        if(A1[j] == 1)
        {
            y[j] = 1;
            s[j] = s1[j];
        }
        else if(A2[j] == 0) // V - A2 = 1
        {
            y[j] = -1;
            s[j] = s2[j];
        }
        else if(A2[j] - A1[j] == 1) // A2 - A1
        {
                //assert(abs(y[j]) <= epsilon);
                y[j] = y[j]/epsilon;
        }
    }

    cleanup();
    return iter+2;
}

ActiveSetTruncChain::~ActiveSetTruncChain()
{
    delete[] A1;
    delete[] A2;
    delete[] V;

    delete[] s1;
    delete[] s2;

    delete[] u_m;
    delete[] W_m;

    delete[] u_n;
    delete[] W_n;

    g->reset();
    delete g;

}
