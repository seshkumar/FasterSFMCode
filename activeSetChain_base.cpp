#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "activeSetChain_base.h"

//h = length   = dims[1]
ActiveSetChain_base::ActiveSetChain_base(const size_t length)
{
    // Start Memory Allocation
    this->length   = length;

    uD = new double[length];
    WD = new double[length-1]; 

    A  = new size_t[length];
    C  = new size_t[length];
    D  = new size_t[length];
    S  = new size_t[length];

    gr = new GraphType(length, 3*length);

    gr->add_node(length);

    for(size_t j=0; j < length; ++j)
    {
        gr->add_tweights(j, 0, 0);
        if(j < length-1)
            gr->add_edge(j, j + 1, 0, 0);
    }
}

void ActiveSetChain_base::initialize(double *y, double *s, const double *W, const double *u)
{
    // Initializing the variables
    this->W  = W;
    this->u  = u;
    this->y  = y;
    this->s  = s;

    std::memset(A, 0, length*sizeof(size_t));
    std::memset(C, 0, length*sizeof(size_t));
    std::memset(D, 0, length*sizeof(size_t));
    std::memset(S, 0, length*sizeof(size_t));

    std::memset(uD, 0, length*sizeof(double));
    std::memset(WD, 0, (length-1)*sizeof(double));

    np = 1;

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    initializeParitionList();
    calculateMarginalsGreedy();
}

void ActiveSetChain_base::initializeParitionList(void)
{
    //PList(Partition List) initialization related to initialization of A
    std::vector<int> list;
    for(int i=0; i < np; ++i)
        pList.push_back(list);

    for(int j=0; j < length; ++j)
    {
        int k = A[j];
        pList[k].push_back(j);
    }

    a.clear();
    c.clear();
    d.clear();
    for(int i=0; i < np; ++i)
    {
        a.push_back(0.0f);
        c.push_back(0.0f);
        d.push_back(0.0f);
    }

    m   = new double[np];
    wp  = new double[np];
}


ActiveSetChain_base::~ActiveSetChain_base()
{
    delete[] A;
    delete[] C;
    delete[] D;
    delete[] S;

    delete[] uD;
    delete[] WD;

    gr->reset();
    delete gr;

    a.clear();
    c.clear();
    d.clear();

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    v.clear();
    p.clear();
}

void ActiveSetChain_base::calculateMarginals(void)
{
    delete[] m;
    delete[] wp;

    m  = new double[np];
    wp = new double[np];

    size_t n0 = 0;
    double fb0 = cutValue(0, n0);

    m[0]  = -fb0/n0;
    wp[0] = n0/2.0f;

    #pragma omp parallel    
    { 
        #pragma omp for
        for(int i=1; i < np; ++i)
        {
            size_t ni, nj;

            double fbi = cutValue(i, ni);    
            double fbj = cutValue(i-1, nj);

            m[i] = ((fbj - fbi)/(ni-nj)); 
            wp[i] = ((ni-nj)/2.0f);
        }
    }
}

// Function updates m and wp
void ActiveSetChain_base::calculateMarginalsGreedy(void)
{
    delete[] m;
    delete[] wp;
    
    m  = new double[np];
    wp = new double[np];

    std::memset(S, 0, length*sizeof(size_t));

    if(pList.size() != np)
    {
        std::cerr << "PList size neq np :"  << pList.size() << " " << np << std::endl;
        
    }

    for(int l=0; l < pList.size(); ++l)
    {
        size_t ni = pList[l].size();
        double gain = 0.0f;
        for(int k=0; k < pList[l].size(); ++k)
        {
            int v = pList[l][k];
            int j = v;
    
            gain -= u[v];
            if( j > 0)
            {
                size_t j1 = j - 1;
                gain += (1.0f-2*S[j1])*W[j1];
            }
    
            if(j < length-1)
            {
                size_t j1 = j + 1;
                gain += (1.0f-2*S[j1])*W[j];
            }
            S[j] = 1;
        } 
        m[l]  = -gain/ni;
        wp[l] = ni/2.0f;
    }
}

double ActiveSetChain_base::cutValue(void)
{
    return cutValue(this->y);
}

double ActiveSetChain_base::cutValue(size_t p, size_t &ni)
{
    double val = 0.0f;
    ni = 0;
    for(int j=0; j < length; ++j)
    {
        int xj = (A[j] <= p)? 1 : 0;
        ni += xj;
        val-= xj*u[j];

        if(j < length-1)
        {
            int xj1 = (A[j+1] <= p)? 1 : 0;
            val += W[j]*(abs(xj - xj1));
        }
    }
    return val;
}

template <class T> double ActiveSetChain_base::cutValue(T *x)
{
    return cutValue(x, this->u, this->W);            
}

template <class T> double ActiveSetChain_base::cutValue(T *x, const double *u, const double *W)
{
    double val = 0.0f;
    for(int j=0; j < length; ++j)
    {
        val -= u[j]*x[j];
        if(j < length-1)
            val+= W[j]*fabs(x[j] - x[j+1]);
    }
            
    return val;
}

size_t ActiveSetChain_base::optimize(void)
{
    iter = 1;
    while(true)
    {
        wpav();
        decouple();
        checkOptimality();
        //printf("Iter = %d; Partitions= %ld; TV = %e; epsilon= %e\n", iter, np, cutValue()+frobeniusNorm(), cutVal);
        //printf("Iter = %d; Partitions= %ld; epsilon= %e\n", iter, np, cutVal);
        if(cutVal > -1e-6)
            break;
        splitSets();
        calculateMarginalsGreedy();
        iter++;
    }
    cleanup();
    return iter;
}

void ActiveSetChain_base::cleanup(void)
{
    delete[] m;
    delete[] wp;
}

// Run wpava and update partitions, weights and dual variables
void ActiveSetChain_base::wpav(void)
{
    std::vector<size_t> index (np, 0);
    std::vector<double> weight (np, 0.0f);
    
    v.clear();
    p.clear();
    for(int i=0; i < np; ++i)
    {
        v.push_back(0.0f);
        p.push_back(0);
    }

    int ci = 0;
    index[ci] = 0;
    weight[ci] = wp[0];
    v[ci] = -m[0];

    for(int j=1; j < np; ++j)
    {
        ci++;
        
        index[ci] = j;
        weight[ci] = wp[j];
        v[ci]      = -m[j];
       
        while(ci >=1 && v[ci-1] >= v[ci])
        {
            double nw = weight[ci-1] + weight[ci];
            v[ci-1] = ((weight[ci-1]*v[ci-1] + weight[ci]*v[ci]))/nw;
            weight[ci-1] = nw;
            ci--;
        }
    } 

    int n_t = np-1;
    np = 0;
    while(n_t >= 0)
    {
        p[ci] = ci;
        for(int j = index[ci]; j <= n_t; ++j)
        {
            p[j] = ci;
            np   = (p[j] > np)? p[j]: np;
            v[j] = v[ci];
        }
        n_t = index[ci] - 1;
        ci--;
    }

    np++;


    index.clear();
    weight.clear();

    a.clear();
    c.clear();
    d.clear();
    for(int i=0; i < np; ++i)
    {
        a.push_back(0.0f);
        c.push_back(0.0f);
        d.push_back(0.0f);
    }


    //std::cout << "In wpav: A = " ;
    for(size_t j=0; j <length; ++j)
    {
        size_t k = A[j]; 
        y[j] = -v[k];
        s[j] = v[k];
        A[j] = p[k];
        size_t k_new = p[k];
        //std::cout << A[j] << " " ;
        
        a[k_new] += (j+1);
    }
    //std::cout << std::endl;
}

double ActiveSetChain_base::frobeniusNorm(void)
{
    double result = 0.0;
    for(int j = 0; j < length; ++j)
    {
        double value = y[j];
        result += value * value;
    }
    return 0.5*result;
}

void ActiveSetChain_base::decouple(void)
{
    //std::cout << "Decoupling started" << std::endl;
    std::memcpy(uD, u, length*sizeof(double));
    for(size_t j=0; j < length-1; ++j)
    {
        size_t aj  = A[j];

        size_t aj1 = A[j+1]; 
        if(aj > aj1)
        {
            WD[j]   = 0;
            uD[j]  += W[j];
            uD[j+1]-= W[j];
        }
        else if(aj < aj1)
        {
            WD[j]    = 0;
            uD[j]   -= W[j];
            uD[j+1] += W[j];
        }
        else
            WD[j] = W[j];

    }
}

void ActiveSetChain_base::checkOptimality(void)
{

     GraphType::arc_id a = gr->get_first_arc();
     for(size_t j=0; j < length; ++j)
     {
         gr->set_trcap(j, -(s[j] + uD[j]));

         if(j < length-1)
         {
             gr->set_rcap(a, WD[j]);
             a = gr->get_next_arc(a);

             gr->set_rcap(a, WD[j]);
             a = gr->get_next_arc(a);
         }
     }
     gr->maxflow();

    std::memset(D, 0, length*sizeof(size_t));
    //std::cout << "In checkOptimality : D = " ;
    for(size_t j=0; j < length; ++j)
    {
        if(gr->what_segment(j))
        {
            size_t k = A[j];
            C[j] = k;
            D[j] = 1;
            c[k] += (j+1);
        }
        else
        {
            C[j] = np + 2; 
        }
        //std::cout << D[j] << " ";
    }
    //std::cout << std::endl;

    cutVal = 0.0f;
    //std::cout << "In checkOptimality : C = " ;
    for(int j=0; j < length; ++j)
    {
        //std::cout << C[j] << " ";
        cutVal -= (s[j]+uD[j])*(D[j]);
        if(j < length-1)
            cutVal+= WD[j]*abs(int(D[j] - D[j+1]));
    }
    //std::cout << std::endl;

    //for(int i=0; i < np ; ++ i)
    //    std::cout << "c[" << i << "] = " << c[i] << std::endl;

}

void ActiveSetChain_base::splitSets(void)
{
    //for(int i=0; i < np ; ++ i)
    //    std::cout << "a[" << i << "] = " << a[i] << " ";
    //std::cout << std::endl;
    //for(int i=0; i < np ; ++ i)
    //    std::cout << "c[" << i << "] = " << c[i] << " ";
    //std::cout << std::endl;
    size_t np_new = 0;
    for(size_t i=0; i < np ; ++i)
    {
        d[i] = np_new;
        if(c[i] != 0 && a[i] != c[i])
            np_new += 2;
        else
            np_new ++;
    } 
    //std::cout << "In splitSets: np_new = "  << np_new << std::endl;

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    std::vector<int> list;
    list.clear();
    for(int i=0; i < np_new; ++i)
        pList.push_back(list);
   
    for(size_t j=0; j<length; ++j)
    {
        int k = A[j];
        if( c[k] != 0 && a[k] != c[k])
        {
            if(C[j] == A[j])
            {
                A[j] = d[k];
                pList[d[k]].push_back(j);
            }
            else
            {
                A[j] = d[k]+1;
                pList[d[k]+1].push_back(j);
            }
        }
        else
        {
            A[j] = d[k];
            pList[d[k]].push_back(j);
        }
    }
    np = np_new;
}
