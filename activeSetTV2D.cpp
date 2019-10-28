#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "activeSetTV2D.h"

ActiveSetTV2D::ActiveSetTV2D(Graph4N2D &g)
{
    initialize(g.W1, g.W2, g.u, g.y, g.dims);
}

ActiveSetTV2D::ActiveSetTV2D(double **W1, double **W2, double **u, double **y, const size_t *dims)
{
    initialize(W1, W2, u, y, dims);
}

void ActiveSetTV2D::initialize(double **W1, double **W2, double **u, double **y, const size_t *dims)
{
    w = dims[0];
    h = dims[1];


    this->W1 = W1;
    this->W2 = W2;
    this->u  = u;
    this->y  = y;

    x   = new double*[w];
    s   = new double*[w];
    uD  = new double*[w];
    W1D = new double*[h];
    W2D = new double*[w];
   
    A = new size_t*[w];
    C = new size_t*[w];
    D = new size_t*[w];
    S = new size_t*[w];

    x  [0] = new double[w*h];
    s  [0] = new double[w*h];
    uD [0] = new double[w*h];
    W1D[0] = new double[h*(w-1)];
    W2D[0] = new double[w*(h-1)]; 

    A[0]  = new size_t[w*h];
    C[0]  = new size_t[w*h];
    D[0]  = new size_t[w*h];
    S[0]  = new size_t[w*h];

    for (size_t i=1 ; i < w ; ++i)
    {
        x [i] = x [0] + i*h;
        s [i] = s [0] + i*h;
        uD[i] = uD[0] + i*h;
        W2D[i]= W2D[0]+ i*(h-1);

        A[i] = A[0] + i*h;
        C[i] = C[0] + i*h;
        D[i] = D[0] + i*h;
        S[i] = S[0] + i*h;
    }

    for(size_t i=1 ; i < h; ++i)
        W1D[i] = W1D[0] + i*(w-1);

    np = 1;
    std::memset(A[0] , 0, w*h*sizeof(size_t));
    std::memset(C[0] , 0, w*h*sizeof(size_t));
    std::memset(D[0] , 0, w*h*sizeof(size_t));
    std::memset(S[0] , 0, w*h*sizeof(size_t));

    //XXX: PList initialization related to initialization of A
    std::vector<int> list;
    for(int i=0; i < np; ++i)
        pList.push_back(list);

    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            int k = A[i][j];
            pList[k].push_back(i*h + j);
        }

    std::memset(x [0] , 0, w*h*sizeof(double));
    std::memset(s [0] , 0, w*h*sizeof(double));
    std::memset(uD[0] , 0, w*h*sizeof(double));
    std::memset(W1D[0], 0, h*(w-1)*sizeof(double));
    std::memset(W2D[0], 0, w*(h-1)*sizeof(double));

    gr = new GraphType(w*h, (w*h + (w-1)*h + w*(h-1))*2);

    gr->add_node(w*h);

    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            gr->add_tweights(i*h + j, -(s[i][j] + uD[i][j]), 0);
            if(i < w-1)
                gr->add_edge(i*h + j, (i+1)*h + j, W1D[j][i], W1D[j][i]);
            if(j < h-1)
                gr->add_edge(i*h + j, i*h + j + 1, W2D[i][j], W2D[i][j]);
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

    //calculateMarginals();
    calculateMarginalsGreedy();
}

ActiveSetTV2D::~ActiveSetTV2D()
{
    delete[]  x[0];
    delete[]  s[0];
    delete[] uD[0];
    delete[] W1D[0];
    delete[] W2D[0];

    delete[] A[0];
    delete[] C[0];
    delete[] D[0];
    delete[] S[0];

    delete[] A;
    delete[] C;
    delete[] D;
    delete[] S;

    delete[] x;
    delete[] s;
    delete[] uD;
    delete[] W1D;
    delete[] W2D;

    gr->reset();
    delete gr;

    a.clear();
    c.clear();
    d.clear();

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    delete[] m;
    delete[] wp;

    v.clear();
    p.clear();
}

void ActiveSetTV2D::calculateMarginals(void)
{
    delete[] m;
    delete[] wp;

    m  = new double[np];
    wp = new double[np];

    size_t n0 = 0;
    double fb0 = cutValue(0, n0);

    m[0]  = -fb0/n0;
    wp[0] = n0/2.0f;

    //#pragma omp parallel    
    { 
        //#pragma omp for
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

void ActiveSetTV2D::calculateMarginalsGreedy(void)
{
    delete[] m;
    delete[] wp;
    
    m  = new double[np];
    wp = new double[np];

    std::memset(S[0] , 0, w*h*sizeof(size_t));

    assert(pList.size() == np);
    //for(int i=0; i < w; ++i)
    //{ 
    //    for(int j=0; j < h; ++j)
    //        std::cout << A[i][j] << " " ;
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;

    //for(int l=0; l < pList.size(); ++l)
    //{
    //    std::cout << "Component - " << l << " " << pList[l].size() << "--";
    //    for(int k=0; k < pList[l].size(); ++k)
    //        std::cout << pList[l][k] << " ";
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;
    

    for(int l=0; l < pList.size(); ++l)
    {
        size_t ni = pList[l].size();
        double gain = 0.0f;
        for(int k=0; k < pList[l].size(); ++k)
        {
            int v = pList[l][k];
            int i = v/h;
            int j = v%h;
    
            gain -= u[i][j];
            if(i > 0)
            {
                size_t i1 = i - 1;
                gain += (1.0f-2*S[i1][j])*W1[j][i1];
            }
    
            if( i < w-1)
            {
                size_t i1 = i + 1;
                gain += (1.0f-2*S[i1][j])*W1[j][i];
            }
    
            if( j > 0)
            {
                size_t j1 = j - 1;
                gain += (1.0f-2*S[i][j1])*W2[i][j1];
            }
    
            if(j < h-1)
            {
                size_t j1 = j + 1;
                gain += (1.0f-2*S[i][j1])*W2[i][j];
            }
            S[i][j] = 1;
        } 
        m[l]  = -gain/ni;
        wp[l] = ni/2.0f;
    }
}

double ActiveSetTV2D::cutValue(void)
{
    return cutValue(this->x);
}

double ActiveSetTV2D::cutValue(size_t p, size_t &ni)
{
    double val = 0.0f;
    ni = 0;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            int xij = (A[i][j] <= p)? 1 : 0;
            ni += xij;
            val -= xij*u[i][j];

            if(i < w-1)
            {
                int xi1j = (A[i+1][j] <= p)? 1 : 0;
                val += W1[j][i]*(abs(xij - xi1j));
            }

            if(j < h-1)
            {
                int xij1 = (A[i][j+1] <= p)? 1 : 0;
                val += W2[i][j]*(abs(xij - xij1));
            }
        }
    return val;
}

template <class T> double ActiveSetTV2D::cutValue(T **x)
{
    double val = 0.0f;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            val -= u[i][j]*x[i][j];
            if(i < w-1)
                val+= W1[j][i]*fabs(x[i][j] - x[i+1][j]);
            if(j < h-1)
                val+= W2[i][j]*fabs(x[i][j] - x[i][j+1]);
        }
            
    return val;
}


void ActiveSetTV2D::optimize(void)
{
    iter = 1;
    while(true)
    {
        wpav();
        decouple();
        checkOptimality();
        printf("Iter = %d; Partitions= %ld; TV = %.10e; gap = %e; epsilon= %e\n", iter, np, cutValue()+frobeniusNorm(), cutValue() + 2*frobeniusNorm(), cutVal);
        //printf("Iter = %d; Partitions= %ld; epsilon= %e\n", iter, np, cutVal);
        if(cutVal > -1e-6)
            break;
        //if(iter > 16)
        //    break;
        splitSets();
        calculateMarginalsGreedy();
        iter++;
    }

    //for(int i=0; i < w; ++i)
    //{
    //    for(int j=0; j < h; ++j)
    //        printf("%e ", x[i][j]);
    //    printf("\n");
    //}
    //printf("\n");
}

// Run wpava and update partitions, weights and dual variables
void ActiveSetTV2D::wpav(void)
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

    for(size_t i=0; i < w; ++i)
    {
        for(size_t j=0; j <h; ++j)
        {
            size_t k = A[i][j]; 
            x[i][j] = -v[k];
            s[i][j] = v[k];
            A[i][j] = p[k];
            //std::cout << "Partition : " << k  << " size : " << a[k].size() << " capacity :" << a[k].capacity() << std::endl;
            //a[k].push_back(j*w + i);
            a[k] += (i*h + j);
        }
        //std::cout << std::endl;
    }
    //std::cout << std::endl;

    //std::cout << "Pava and updating partitions done" << std::endl;
}

double ActiveSetTV2D::frobeniusNorm(void)
{
    double result = 0.0;
    for(int i = 0; i < w; ++i)
    {
        for(int j = 0; j < h; ++j)
        {
            double value = x[i][j];
            result += value * value;
        }
    }
    return 0.5*result;
}

void ActiveSetTV2D::decouple(void)
{
    //std::cout << "Decoupling started" << std::endl;
    std::memcpy(uD[0], u[0], w*h*sizeof(double));
    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            size_t aij  = A[i][j];

            if(i < w-1)
            {
                size_t ai1j = A[i+1][j]; 
                if(aij > ai1j)
                {
                    W1D[j][i]  = 0;
                    uD[i][j]   += W1[j][i];
                    uD[i+1][j] -= W1[j][i]; 
                }
                else if(aij < ai1j)
                {
                    W1D[j][i]  = 0;
                    uD[i][j]   -= W1[j][i];
                    uD[i+1][j] += W1[j][i];
                }
                else
                    W1D[j][i] = W1[j][i];
            }
            
            if(j < h-1)
            {
                size_t aij1 = A[i][j+1]; 
                if(aij > aij1)
                {
                    W2D[i][j]  = 0;
                    uD[i][j]  += W2[i][j];
                    uD[i][j+1]-= W2[i][j];
                }
                else if(aij < aij1)
                {
                    W2D[i][j]   = 0;
                    uD[i][j]   -= W2[i][j];
                    uD[i][j+1] += W2[i][j];
                }
                else
                    W2D[i][j] = W2[i][j];
            }

        }
    //for(int i=0; i < w; ++i)
    //{
    //    for(int j=0; j < h; ++j)
    //        std::cout << uD[i][j] << " " ;
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;

    //for(int i=0; i < w-1; ++i)
    //{
    //    for(int j=0; j < h; ++j)
    //        std::cout << W1D[j][i] << " ";
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;

    //for(int i=0; i < w; ++i)
    //{
    //    for(int j=0; j < h-1; ++j)
    //        std::cout << W2D[i][j] << " ";
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;

    //std::cout << "Decoupling done" << std::endl;

}

void ActiveSetTV2D::checkOptimality(void)
{
    //TODO: Find alternative way other to reset edge weights
//    gr->reset();
//    gr->add_node(w*h);
//
//    for(size_t i=0; i < w; ++i)
//        for(size_t j=0; j < h; ++j)
//        {
//            gr->add_tweights(i*h + j, -(s[i][j] + uD[i][j]), 0);
//            if(i < w-1)
//                gr->add_edge(i*h + j, (i+1)*h + j, W1D[j][i], W1D[j][i]);
//            if(j < h-1)
//                gr->add_edge(i*h + j, i*h + j + 1, W2D[i][j], W2D[i][j]);
//        }

    GraphType::arc_id a = gr->get_first_arc();
    for(size_t i=0; i < w; ++i)
        for(size_t j=0; j < h; ++j)
        {
            gr->set_trcap(i*h + j, -(s[i][j] + uD[i][j]));

            if(i < w-1)
            {
                gr->set_rcap(a, W1D[j][i]);

                a = gr->get_next_arc(a);

                gr->set_rcap(a, W1D[j][i]);

                if(i*h + j < (w*h)-1)
                    a = gr->get_next_arc(a);
            }
            if(j < h-1)
            {
                gr->set_rcap(a, W2D[i][j]);

                a = gr->get_next_arc(a);

                gr->set_rcap(a, W2D[i][j]);
    
                if(i*h + j < (w*h)-1)
                    a = gr->get_next_arc(a);
            }
        }

    gr->maxflow();

    std::memset(D[0] , 0, w*h*sizeof(size_t));
    for(size_t i=0; i < w; ++i)
    {
        for(size_t j=0; j < h; ++j)
        {
            if(gr->what_segment(i*h + j))
            {
                //std::cout << gr->what_segment(j*w + i) << " " ;
                size_t k = A[i][j];
                C[i][j] = k;
                D[i][j] = 1;
                //c[k].push_back(j*w + i);
                c[k] += (i*h + j);
            }
            else
            {
            //    std::cout << "0" << " " ;
                C[i][j] = np+1; 
            }
        }
        //std::cout << std::endl;
    }
    //std::cout << std::endl;

    cutVal = 0.0f;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            cutVal -= (s[i][j]+uD[i][j])*(D[i][j]);
            if(i < w-1)
                cutVal+= W1D[j][i]*abs(D[i][j] - D[i+1][j]);
            if(j < h-1)
                cutVal+= W2D[i][j]*abs(D[i][j] - D[i][j+1]);
        }

}
void ActiveSetTV2D::splitSets(void)
{
    size_t np_new = 0;
    for(size_t i=0; i < np ; ++i)
    {
        d[i] = np_new;
        if(c[i] != 0 && a[i] != c[i])
            np_new += 2;
        else
            np_new ++;
        
        //std::set_difference(a[i].begin(), a[i].end(), c[i].begin(), c[i].end(), std::inserter(d[i], d[i].begin()));
        
        //std::cout << "**** " << i  << " ****" << std::endl;
        //for(size_t j=0; j< a[i].size(); ++j)
        //    std::cout << a[i][j] << " " ;
        //std::cout << std::endl;

        //for(size_t j=0; j< c[i].size(); ++j)
        //    std::cout << c[i][j] << " " ;
        //std::cout << std::endl;

        //for(size_t j=0; j< d[i].size(); ++j)
        //    std::cout << d[i][j] << " " ;
        //std::cout << std::endl;

        //for(size_t j=0; j < a[i].size(); ++j)
        //{
        //    int k = a[i][j];
        //    A[k%w][k/w] = np_new;
        //}
        //if(!c[i].empty() && !d[i].empty())
        //{
        //    for(size_t j=0; j < d[i].size(); ++j)
        //    {
        //        int k = d[i][j];
        //        A[k%w][k/w] = np_new + 1;
        //    }
        //    np_new += 2;
        //}
        //else
        //    np_new++;
    } 

    for(int i=0; i < pList.size(); ++i)
        pList[i].clear();
    pList.clear();

    std::vector<int> list;
    list.clear();
    for(int i=0; i < np_new; ++i)
        pList.push_back(list);
   
    //for(size_t i=0; i < np; ++i)
    //    std::cout << a[i] << " " << c[i] << " " << d[i] <<  std::endl;
    //std::cout << std::endl; 
    for(size_t i=0; i< w; ++i)
        for(size_t j=0; j<h; ++j)
        {
            int k = A[i][j];
            if( c[k] != 0 && a[k] != c[k])
            {
                if(C[i][j] == A[i][j])
                {
                    A[i][j] = d[k];
                    pList[d[k]].push_back(i*h + j);
                }
                else
                {
                    A[i][j] = d[k]+1;
                    pList[d[k]+1].push_back(i*h + j);
                }
            }
            else
            {
                A[i][j] = d[k];
                pList[d[k]].push_back(i*h + j);
            }
        }
    np = np_new;

    //std::cout << std::endl;
    //for(size_t i=0; i < w; ++i)
    //{
    //    for(size_t j=0; j < h; ++j)
    //        std::cout << A[i][j] << " ";
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //std::cout << "New Partitions = " << np << std::endl;
}
