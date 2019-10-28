#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "activeSetTrunc2D.h"
#include "contraction.h"
#include "restriction.h"

ActiveSetTrunc2D::ActiveSetTrunc2D(const size_t w, const size_t h)
                 :ActiveSet2D_base(w, h)
{

    g = new GraphType(w*h, (w*h + (w-1)*h + w*(h-1))*2);

    g->add_node(w*h);

    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            g->add_tweights(i*h + j, 0, 0);
            if(i < w-1)
                g->add_edge(i*h + j, (i+1)*h + j, 0, 0);
            if(j < h-1)
                g->add_edge(i*h + j, i*h + j + 1, 0, 0);
        }

    A1 = new size_t*[w];
    A2 = new size_t*[w];
    V  = new size_t*[w];

    s1 = new double*[w];
    s2 = new double*[w];

    u_m = new double*[w];
    W1_m = new double*[h];
    W2_m = new double*[w];

    u_n = new double*[w];
    W1_n = new double*[h];
    W2_n = new double*[w];

    A1[0] = new size_t[w*h];
    A2[0] = new size_t[w*h];
    V[0] = new size_t[w*h];

    s1[0] = new double[w*h];
    s2[0] = new double[w*h];

    u_m[0] = new double[w*h];
    u_n[0] = new double[w*h];

    W1_m[0] = new double[h*(w-1)];
    W1_n[0] = new double[h*(w-1)];

    W2_m[0] = new double[w*(h-1)];
    W2_n[0] = new double[w*(h-1)];

    for(size_t i=1; i < w; ++i)
    {
        A1[i] = A1[0] + i*h;
        A2[i] = A2[0] + i*h;
        V [i] = V [0] + i*h;

        s1[i] = s1[0] + i*h;
        s2[i] = s2[0] + i*h;

        u_m[i] = u_m[0] + i*h;
        u_n[i] = u_n[0] + i*h;

        W2_m[i] = W2_m[0] + i*(h-1);
        W2_n[i] = W2_n[0] + i*(h-1);
    }

    for(size_t i=1; i < h; ++i)
    {
        W1_m[i] = W1_m[0] + i*(w-1);
        W1_n[i] = W1_n[0] + i*(w-1);
    }
}

void ActiveSetTrunc2D::initialize(double **y, double **s, const double *const*W1, const double *const*W2, const double *const* u, double eps)
{
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            V[i][j] = 1;

    std::memset(A1[0], 0, w*h*sizeof(size_t));
    std::memset(A2[0], 0, w*h*sizeof(size_t));

    std::memset(s1[0], 0, w*h*sizeof(double));
    std::memset(s2[0], 0, w*h*sizeof(double));

    std::memset(u_m[0], 0, w*h*sizeof(double));
    std::memset(u_n[0], 0, w*h*sizeof(double));

    std::memset(W1_m[0], 0, h*(w-1)*sizeof(double));
    std::memset(W1_n[0], 0, h*(w-1)*sizeof(double));

    std::memset(W2_m[0], 0, w*(h-1)*sizeof(double));
    std::memset(W2_n[0], 0, w*(h-1)*sizeof(double));

    std::memset(A[0], 0, w*h*sizeof(size_t));
    std::memset(C[0], 0, w*h*sizeof(size_t));
    std::memset(D[0], 0, w*h*sizeof(size_t));
    std::memset(S[0], 0, w*h*sizeof(size_t));

    std::memset(uD[0], 0, w*h*sizeof(double));
    std::memset(W1D[0],0, h*(w-1)*sizeof(double));
    std::memset(W2D[0], 0,w*(h-1)*sizeof(double));

    this->epsilon = eps;

    this->y  = y;
    this->s  = s;

    np = 1;

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    initializeParitionList();
    truncateWeights(W1, W2, u);
    splitSets();

    this->W1 = W1_n;
    this->W2 = W2_n;
    this->u  = u_n;

    calculateMarginalsGreedy();
}

void ActiveSetTrunc2D::truncateWeights(const double *const* W1, const double *const* W2, const double *const* u)
{    
    // A1 is equivalent to A_+
    // A2 is equivalent to A_-
    // A_+ \subset A_-
    // A1 \subset A2

    dual1 = getPrimalDual(A1, s1, W1, W2, u, epsilon);
    dual2 = getPrimalDual(A2, s2, W1, W2, u, -epsilon);

    primal1 = cval(A1, u, W1, W2, w, h);
    primal2 = cval(A2, u, W1, W2, w, h);

    assert(fabs(primal1 - dual1) < 1e-1);
    assert(fabs(primal2 - dual2) < 1e-1);

    //std::cout << "cut1 : Primal = " << primal1 << "; Dual =  " << dual1 << "; Gap = " << primal1 - dual1 << std::endl;
    //std::cout << "cut2 : Primal = " << primal2 << "; Dual =  " << dual2 << "; Gap = " << primal2 - dual2 << std::endl;

    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            size_t k = A[i][j];

           if(A2[i][j] - A1[i][j] == 1)
            {
                C[i][j] = k;
                D[i][j] = 1;
                c[k] += (i*h + j + 1);
            }
            else
            {
                C[i][j] = np+2;
            }

            a[k] += (i*h + j + 1);
        }

    // F: 2^{A2 - A1} \to R
    restriction2D(u_m, W1_m, W2_m,   u,   W1,   W2, A2, w, h);
    contraction2D(u_n, W1_n, W2_n, u_m, W1_m, W2_m, A1, w, h);
}



double ActiveSetTrunc2D::getPrimalDual(size_t **primal, double **dual, const double *const* W1, const double *const* W2, const double *const*u, double lambda)
{
    double dualCost = 0.0f;
    GraphType::arc_id a = g->get_first_arc();
    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            g->set_trcap(i*h + j, - u[i][j] + lambda);
            if( i < w - 1)
            {
                g->set_rcap(a, W1[j][i]);
                a = g->get_next_arc(a);

                g->set_rcap(a, W1[j][i]);
                if(i*h + j < (w*h)-1)
                    a = gr->get_next_arc(a);
            }
            if(j < h-1)
            {
                g->set_rcap(a, W2[i][j]);
                a = g->get_next_arc(a);

                g->set_rcap(a, W2[i][j]);
                if(i*h + j < (w*h)-1)
                    a = g->get_next_arc(a);
            }
        }

    g->maxflow();

    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            dual[i][j] = g->get_trcap(i*h + j) - lambda;
            if(g->what_segment(i*h + j))
            {
                primal[i][j] = 1;
                dualCost += dual[i][j];
            }
        }
    return dualCost;
}

double ActiveSetTrunc2D::truncatedTVNorm(void)
{
    double result = 0.0;
    for(int i = 0; i < w; ++i)
    {
        for(int j = 0; j < h; ++j)
        {
            double value = y[i][j];
            if(value < epsilon)
                result += value * value;
        }
    }
    return 0.5*result/epsilon;
}

int ActiveSetTrunc2D::optimize(void)
{
    iter = 1;
    while(true)
    {
        wpav();
        decouple();
        checkOptimality();
        //printf("Iter = %d; Partitions= %ld; TV = %e; epsilon= %e\n", iter, np, cutValue()+truncatedTVNorm(), cutVal);
        if(cutVal > -1e-6)
            break;
        splitSets();
        calculateMarginalsGreedy();
        iter++;
    }

    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            if(A1[i][j] == 1)
            {
                y[i][j] = 1;
                s[i][j] = s1[i][j];
            }
            else if(A2[i][j] == 0) // V - A2 = 1
            {
                y[i][j] = -1;
                s[i][j] = s2[i][j];
            }
            else if(A2[i][j] - A1[i][j] == 1) // A2 - A1
            {
                    //XXX: Pontential bug
                    //assert(abs(y[i][j]) <= epsilon);
                    //if(abs(y[i][j]) > epsilon)
                    //{
                    //    std::cerr << "i = " << i << "; j = " << j << std::endl;
                    //    std::cerr << "A1[i][j] = " << A1[i][j] << std::endl;
                    //    std::cerr << "A2[i][j] = " << A2[i][j] << std::endl;
                    //    if(i < w-1)
                    //    {
                    //        std::cerr << "A1[i+1][j] = " << A1[i+1][j] << std::endl;
                    //        std::cerr << "A2[i+1][j] = " << A2[i+1][j] << std::endl;
                    //    }
                    //    if(j < h-1)
                    //    {
                    //        std::cerr << "A1[i][j+1] = " << A1[i][j+1] << std::endl;
                    //        std::cerr << "A2[i][j+1] = " << A2[i][j+1] << std::endl;
                    //    }
                    //    std::cerr << "s1[i][j] = " << s1[i][j] << std::endl;
                    //    std::cerr << "s2[i][j] = " << s2[i][j] << std::endl;
                    //    std::cerr << "y[i][j] = "  << y[i][j] << std::endl;
                    //    std::cerr << "epsilon = "  << epsilon << std::endl;
                    //    std::cerr << std::endl;

                    //}
                    if(y[i][j] > epsilon)
                        y[i][j] = epsilon;
                    else if(y[i][j] < -epsilon)
                        y[i][j] = -epsilon;

                    y[i][j] = y[i][j]/epsilon;
            }
        }

    cleanup();
    return iter+2;
}

ActiveSetTrunc2D::~ActiveSetTrunc2D()
{
    //cleanup();
    delete[] A1[0];
    delete[] A2[0];
    delete[] V [0];

    delete[] A1;
    delete[] A2;
    delete[] V;

    delete[] s1[0];
    delete[] s2[0];

    delete[] s1;
    delete[] s2;

    delete[] u_m[0];
    delete[] u_n[0];

    delete[] u_m;
    delete[] u_n;

    delete[] W1_m[0];
    delete[] W2_m[0];

    delete[] W1_n[0];
    delete[] W2_n[0];

    delete[] W1_m;
    delete[] W2_m;

    delete[] W1_n;
    delete[] W2_n;

    g->reset();
    delete g;

}
