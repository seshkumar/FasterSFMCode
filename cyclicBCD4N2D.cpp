#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <limits>
//#include <algorithm>
#include <parallel/algorithm>
#include "cyclicBCD4N2D.h"

CyclicBCD4N2D::CyclicBCD4N2D(Graph4N2D &g)
{
    initialize(g.y, g.s, g.W1, g.W2, g.u, g.dims);
}

CyclicBCD4N2D::CyclicBCD4N2D(double **y, double **s, double **W1, double **W2, double **u, size_t *dims)
{
    initialize(y, s, W1, W2, u, dims);
}

void CyclicBCD4N2D::initialize(double **y, double **s, double **W1, double **W2, double **u, size_t *dims)
{
    w = dims[0];
    h = dims[1];
    nThreads = h + w;

    for(size_t i=0; i < w*h; ++i)
    {
        node n(0.0, i%w, i/w);
        nodes.push_back(n);
    }

    C = new bool*[h];
    for(size_t j=0; j < h; ++j)
        C[j] = new bool[w]();

    optimal = new bool*[w]();
    for(size_t i=0; i < w; ++i)
        optimal[i] = new bool[h]();

    n = new size_t[nThreads];
    size_t offset = 0;
    for(int i=offset; i < offset+h; ++i)
        n[i] = w;
    offset += h;
    for(int i=offset; i < offset+w; ++i)
        n[i] = h;

    this->u = u;
    this->y = y;
    this->s = s;

    projL1 = new DualLineProjector(h, n  );
    projL2 = new DualLineProjector(w, n+h);

    W = new const double*[nThreads]();

    std::copy(W1, W1+h, W   );
    std::copy(W2, W2+w, W+h );

    u2D  = new double*[nThreads]();
    z    = new double*[nThreads]();
    ta   = new double*[nThreads]();
    tb   = new double*[nThreads]();
    for(int i=0; i < nThreads; ++i)
    {
        u2D[i]  = new double[n[i]]();
        z[i]    = new double[n[i]]();
        ta[i]   = new double[n[i]]();
        tb[i]   = new double[n[i]]();
    }

    float m = duplicate(u2D, u);

    
    // divide the unary vector by 3 and initializing tb with mean
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
            {
                u2D[  j][i] /= 2;
                u2D[h+i][j] /= 2;

                tb[  j][i] = -m;
                tb[h+i][j] = -m;
            }
    }
}

CyclicBCD4N2D::~CyclicBCD4N2D()
{
    delete[] n;
    delete projL1;
    delete projL2;
    delete[] W;
    for(int i=0; i < nThreads; ++i)
    {
        delete[] u2D[i];
        delete[] z[i];  
        delete[] ta[i];
        delete[] tb[i];
    }
    for(size_t j=0; j < h; ++j)
        delete[] C[j];
    delete[] C;

    for(size_t i=0; i < w; ++i)
        delete[] optimal[i];
    delete[] optimal;

    delete[] u2D;
    delete[] z;
    delete[] ta;
    delete[] tb;
    
    nodes.clear();

    u = NULL;
}

void CyclicBCD4N2D::parallelAccumulate(size_t dir)
{
    size_t w = this->w;
    size_t h = this->h;
 
    double **ta1 = this->ta;
    double **ta2 = ta1+h;
    double **tb1 = this->tb;
    double **tb2 = tb1+h;
  
    size_t i, j;

    switch(dir)
    {
        case 1:
            //#pragma omp parallel shared(tb1, ta2, i, j, w, h) default(none)
            {
                //#pragma omp for
                for(i = 0; i < w; ++i)
                    for(j = 0; j < h; ++j)
                        tb1[j][i] = -ta2[i][j];
            }
            break;
        case 2:
            //#pragma omp parallel shared(tb2, ta1, i, j, w, h) default(none)
            {
                //#pragma omp for
                for(j = 0; j < h; ++j)
                    for(i = 0; i < w; ++i)
                        tb2[i][j] = -ta1[j][i];
            }
            break;
    }
}

void CyclicBCD4N2D::computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> &stats)
{
    double **x1 = this->ta;
    double **x2 = x1+h;

    double **y = this->y;
    double **s = this->s;
    double **u = this->u;

    bool **C = this->C;

    double dValueSfm = 0.0;
    double dValue    = 0.0;
    double pValue    = 0.0;

    double val = 0.0;

    int w = this->w;
    int h = this->h;
    
    std::vector<node>::iterator it = nodes.begin();

    double yij;
    double sij;

    int i, j, idx;
//    #pragma omp parallel shared(x1, x2, s1, s2, y, u, C, it, w, h, d) private(i, j, val, idx) reduction(+: pValue, dValueSfm) reduction(-:dValue) default(none)
    {
//        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
            {
                C[j][i] = false;

                val = x1[j][i] + x2[i][j];

                y[i][j] = yij = -val;
                s[i][j] = sij = val;

                if(sij < 0)
                    dValueSfm += sij;

                dValue -= (sij*sij)/2;

                pValue += (yij*yij/2 - yij*u[i][j]);


                int idx = j*w + i;
                (it+idx)->s = yij;
                (it+idx)->i = i;
                (it+idx)->j = j;
            }
    }

    //std::cout << std::endl  << "y = " << std::endl;
    //for(size_t i=0; i < w; ++i)
    //{
    //    for(size_t j=0; j < h; ++j)
    //        std::cout << y[i][j] << " ";
    //    std::cout << std::endl;
    //}

    //std::cout << std::endl  << "s = " << std::endl;
    //for(size_t i=0; i < w; ++i)
    //{
    //    for(size_t j=0; j < h; ++j)
    //        std::cout << s[i][j] << " ";
    //    std::cout << std::endl;
    //}



    __gnu_parallel::sort(nodes.begin(), nodes.end(), compare);
    //std::sort(nodes.begin(), nodes.end(), compare);

    double pValueSfm = greedyTV();

    double diff = 0.0;
    const double **W1 = this->W;
    const double **W2 = W1+h;
    double p1, p2;
    //#pragma omp parallel shared(s1, s2, s3, y, W1, W2, W3, w, h, d) private(i, j, diff, p1, p2, p3) reduction(+: pValue, jDistNum, jDistDen) default(none)
    {
        //#pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
            {
                if(i < w-1)
                    pValue += fabs(y[i][j] - y[i+1][j])*W1[j][i];
                if(j < h-1)
                    pValue += fabs(y[i][j] - y[i][j+1])*W2[i][j];
            }
    }

    stats.push_back(OptInfo(iter, pValue, dValue, pValueSfm, dValueSfm));
}

double CyclicBCD4N2D::greedyTV(void)
{
    const double **W1 = W;
    const double **W2 = W1+h;

    double **u  = this->u;

    double fmin = std::numeric_limits<double>::max();
    double f    = 0.0;
    for(std::vector<node>::iterator it= nodes.begin(); it!=nodes.end(); ++it)
    {
        size_t x = it->i;
        size_t y = it->j;
        double gain = -u[x][y]; 
        if(x > 0)
        {
            size_t x1 = x-1;
            gain += (1-2*C[y][x1])*W1[y][x1];
        }
        if( x < w-1)
        {
            size_t x1 = x+1;
            gain += (1-2*C[y][x1])*W1[y][x];
        }
        if( y > 0)
        {
            size_t y1 = y-1;
            gain += (1-2*C[y1][x])*W2[x][y1];
        }
        if(y < h-1)
        {
            size_t y1 = y+1;
            gain += (1-2*C[y1][x])*W2[x][y];
        }
        C[y][x] = true;

        f   += gain;
        //fmin = (f < fmin)? f: fmin;
        if(f < fmin)
        {
            fmin = f;
            threshold = it->s;
        }
    }
    return fmin;
}

void CyclicBCD4N2D::optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname)
{   
    size_t w = this->w;
    size_t h = this->h;
 
    double **ta1 = ta;
    double **ta2 = ta1+h;

    double **z1 = this->z;
    double **z2 = z1+h;

    double **tb1 = tb;
    double **tb2 = tb1+h;

    const double **W1 = this->W;
    const double **W2 = W1+h;

    double **u2D1 = this->u2D;
    double **u2D2 = u2D1+h;

    double **w1 = z;
    double **w2 = w1+h;

    maxIter = maxIterations;

    size_t nIter1;
    size_t nIter2;

    for(int iter =0; iter < maxIter; ++iter)
    {

        nIter1 = 0;
        nIter2 = 0;

        projL1->parallelProject(z1, ta1, tb1, W1, u2D1);
        parallelAccumulate(2);

        for(size_t i=0; i < h; ++i)
            nIter1 = (nIter1 > projL1->iter[i])? nIter2: projL1->iter[i];
            
        projL2->parallelProject(z2, ta2, tb2, W2, u2D2);
        parallelAccumulate(1);

        for(size_t i=0; i < w; ++i)
            nIter2 = (nIter2 > projL2->iter[i])? nIter2: projL2->iter[i];

        computeGapsUpdateSolution(iter, stats);
        stats.back().nSFM = nIter1 + nIter2;
        double primal_tv = stats.back().primalS;
        double dual_tv = stats.back().dualS;
        double primal_sfm = stats.back().primalD;
        double dual_sfm   = stats.back().dualD;
        std::cerr << std::setprecision(10) << "Iter :" << iter  <<  " " << primal_sfm  << " " << dual_sfm << " " << primal_sfm - dual_sfm << " " << stats.back().nSFM << std::endl;
        if(abs(primal_sfm - dual_sfm) < 0.1)
            break;
    }
}

// size(x) =  (h, w) + (w, h) size(u) = (w, h)
double CyclicBCD4N2D::duplicate(double **x, double ** u)
{
    double **x1 = x;
    double **x2 = x1+h;

    double mean = 0.0;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
        {
            double val = u[i][j];

            x1[j][i] = val;
            x2[i][j] = val;

            mean += val;
        }
    mean /= (w*h);
    return mean;
}
