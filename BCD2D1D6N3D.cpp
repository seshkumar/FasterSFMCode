#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <limits>
//#include <parallel/algorithm>
#include <algorithm>
#include "BCD2D1D6N3D.h"

BCD2D1D6N3D::BCD2D1D6N3D(Graph6N3D &g)
{
    initialize(g.y, g.W1, g.W2, g.W3, g.u, g.dims);
}

BCD2D1D6N3D::BCD2D1D6N3D(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims)
{
    initialize(y, W1, W2, W3, u, dims);
}

void BCD2D1D6N3D::initialize(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims)
{
    w = dims[0];
    h = dims[1];
    d = dims[2];

    nThreadsF = d;
    nThreadsL = w*h;

    for(size_t i=0; i < w*h*d; ++i)
    {
        node n(0.0, i%w, (i/w)%h, i/(w*h));
        nodes.push_back(n);
    }

    C = new bool**[d];
    for(size_t k=0; k < d; ++k)
    {
        C[k] = new bool*[h];
        for(size_t j=0; j < h; ++j)
            C[k][j] = new bool[w]();
    }

    optimal = new bool**[w]();
    for(size_t i=0; i < w; ++i)
    {
        optimal[i] = new bool*[h]();
        for(size_t j=0; j < h; ++j)
            optimal[i][j] = new bool[d]();
    }

    n = new size_t[nThreadsL];
    for(int i=0; i < nThreadsL; ++i)
        n[i] = d;

    this->u = u;
    this->y = y;

    projF = new DualFrameProjector(nThreadsF, w, h);
    projL = new DualLineProjector (nThreadsL, n);

    W = new const double*[d*h + w*d + w*h];

    std::copy(W1, W1+d*h, W            );
    std::copy(W2, W2+w*d, W + d*h      );
    std::copy(W3, W3+w*h, W + d*h + w*d);

    lines = d*w + h*w;
    m = new size_t[lines];
    for(size_t i=0   ; i < d*w  ; ++i)
        m[i] = h;
    for(size_t i=d*w ; i < lines; ++i)
        m[i] = d;

    u2D  = new double*[lines]();
    s    = new double*[lines]();
    z    = new double*[lines]();
    ta   = new double*[lines]();
    tb   = new double*[lines]();

    u2D[0] = new double[2*w*d*h];
    s  [0] = new double[2*w*d*h];
    z  [0] = new double[2*w*d*h];
    ta [0] = new double[2*w*d*h];
    tb [0] = new double[2*w*d*h];

    for(int i=1; i < lines; ++i)
    {
        u2D[i]  = u2D[i-1] + m[i];
        s[i]    = s[i-1]   + m[i];
        z[i]    = z[i-1]   + m[i];
        ta[i]   = ta[i-1]  + m[i];
        tb[i]   = tb[i-1]  + m[i];
    }

    // duplicate returns mean of u
    mean = duplicate(u2D, u);

    // divide the unary vector by 3 and initializing tb with mean
    double **u1 = this->u2D;
    double **u2 = u1 + w*d;

    double **tb1 = this->tb;
    double **tb2 = tb1 + w*d;
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
                for(int k=0; k < d; ++k)
                {
                    u1[k*w + i][j] /= 2;
                    u2[j*w + i][k] /= 2;

                    tb1[k*w + i][j] = -mean;
                    tb2[j*w + i][k] = -mean;
                }
    }
}

BCD2D1D6N3D::~BCD2D1D6N3D()
{
    delete[] n;
    delete[] m;

    delete projF;
    delete projL;

    delete[] W;
    for(size_t k=0; k < d; ++k)
    {
        for(size_t j=0; j < h; ++j)
            delete[] C[k][j];
        delete[] C[k];
    }
    delete[] C;

    for(size_t i=0; i < w; ++i)
    {
        for(size_t j=0; j < h; ++j)
            delete[] optimal[i][j];
        delete[] optimal[i];
    }
    delete[] optimal;

    delete[] u2D[0];
    delete[] s[0];
    delete[] z[0];
    delete[] ta[0];
    delete[] tb[0];
 

    delete[] u2D;
    delete[] s;
    delete[] z;
    delete[] ta;
    delete[] tb;
    
    nodes.clear();

    u = NULL;
}

void BCD2D1D6N3D::parallelAccumulate(size_t dir, double beta)
{
    size_t w = this->w;
    size_t h = this->h;
    size_t d = this->d;
 
    double **ta1 = ta;
    double **ta2 = ta1+(w*d);
    double **tb1 = tb;
    double **tb2 = tb1+(w*d);
  
    size_t i, j, k;

    switch(dir)
    {
        case 1:
            for(i = 0; i < w; ++i)
                for(j = 0; j < h; ++j)
                    for(k = 0; k < d; ++k)
                        tb1[k*w + i][j] = -ta2[j*w + i][k];
            break;
        case 2:
            for(j = 0; j < h; ++j)
                for(i = 0; i < w; ++i)
                    for(k = 0; k < d; ++k)
                        tb2[j*w + i][k] = -ta1[k*w + i][j];
            break;
	 case 3:
            for(i = 0; i < w; ++i)
                for(j = 0; j < h; ++j)
                    for(k = 0; k < d; ++k)
                        tb1[k*w + i][j] = -(ta2[j*w + i][k] + beta*(ta2[j*w + i][k] + tb1[k*w + i][j]));
	    break;


    }
}

void BCD2D1D6N3D::computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> &stats)
{
    double **x1 = this->ta;
    double **x2 = x1+(w*d);

    double ***y = this->y;
    double ***u = this->u;

    bool ***C = this->C;

    double dValueSfm = 0.0f;
    double dValue    = 0.0f;
    double pValue    = 0.0f;

    double val = 0.0;

    int w = this->w;
    int h = this->h;
    int d = this->d;
    
    std::vector<node>::iterator it = nodes.begin();

    double yijk;
    double sijk;

    int i, j, k, idx;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            for(int k=0; k < d; ++k)
            {
                C[k][j][i] = false;

                val = (x1[k*w + i][j] + x2[j*w + i][k]);

                y[i][j][k] = yijk = -val;
                sijk = val;

                if(sijk < 0)
                    dValueSfm += sijk;

                dValue -= (sijk*sijk)/2;

                pValue += (yijk*yijk/2 - yijk*u[i][j][k]);

                // Converting to Primal
                y[i][j][k]     = -val;

                int idx = k*w*h + j*w + i;
                (it+idx)->s = yijk;
                (it+idx)->i = i;
                (it+idx)->j = j;
                (it+idx)->k = k;
            }

    //__gnu_parallel::sort(nodes.begin(), nodes.end(), compare);
    std::sort(nodes.begin(), nodes.end(), compare);

    double pValueSfm = greedyTV();

    const double **W1 = this->W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            for(int k=0; k < d; ++k)
            {
                if(i < w-1)
                    pValue += labs(y[i][j][k] - y[i+1][j][k])*W1[k*h + j][i];
                if(j < h-1)
                    pValue += labs(y[i][j][k] - y[i][j+1][k])*W2[k*w + i][j];
                if(k < d-1)
                    pValue += labs(y[i][j][k] - y[i][j][k+1])*W3[j*w + i][k];
            }


    stats.push_back(OptInfo(iter, pValue, dValue, pValueSfm, dValueSfm));
}

double BCD2D1D6N3D::greedyTV(void)
{
    const double **W1 = W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);

    double ***u  = this->u;

    double fmin = std::numeric_limits<double>::max();
    double f    = 0.0;
    for(std::vector<node>::iterator it= nodes.begin(); it!=nodes.end(); ++it)
    {
        size_t x = it->i;
        size_t y = it->j;
        size_t z = it->k;
        double gain = -u[x][y][z]; 
        if(x > 0)
        {
            size_t x1 = x-1;
            gain += (1-2*C[z][y][x1])*W1[z*h + y][x1];
        }
        if( x < w-1)
        {
            size_t x1 = x+1;
            gain += (1-2*C[z][y][x1])*W1[z*h + y][x];
        }
        if( y > 0)
        {
            size_t y1 = y-1;
            gain += (1-2*C[z][y1][x])*W2[z*w + x][y1];
        }
        if(y < h-1)
        {
            size_t y1 = y+1;
            gain += (1-2*C[z][y1][x])*W2[z*w + x][y];
        }
        if( z > 0)
        {
            size_t z1 = z-1;
            gain += (1-2*C[z1][y][x])*W3[y*w + x][z1];
        }
        if( z < d-1)
        {
            size_t z1 = z+1;
            gain += (1-2*C[z1][y][x])*W3[y*w + x][z];
        }
        C[z][y][x] = true;

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

void BCD2D1D6N3D::optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname)
{   
    size_t w = this->w;
    size_t h = this->h;
    size_t d = this->d;

    double **z1 = this->z;
    double **z2 = z1+(w*d);
 
    double **ta1 = this->ta;
    double **ta2 = ta1+(w*d);

    double **tb1 = this->tb;
    double **tb2 = tb1+(w*d);

    const double **W1 = this->W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);

    double **u2D1 = this->u2D;
    double **u2D2 = u2D1+(w*d);

    maxIter = maxIterations;
    
    size_t nIterF, nIterL;
    for(int iter =0; iter < maxIter; ++iter)
    {
        nIterF = 0;
        nIterL = 0;

        projF->parallelProject(z1, ta1, tb1, W1, W2, u2D1);
        parallelAccumulate(2);

        for(size_t i=0; i < nThreadsF; ++i)
            nIterF += projF->iter[i];
                                                      
        projL->parallelProject(z2, ta2, tb2, W3, u2D2);
        parallelAccumulate(1);
                                                      
        for(size_t i=0; i < nThreadsL; ++i)
            nIterL += projL->iter[i];

        computeGapsUpdateSolution(iter, stats);
        
        stats.back().nSFM = nIterF/nThreadsF;

        double primal_sfm = stats.back().primalD;
        double dual_sfm   = stats.back().dualD;
        std::cerr << std::setprecision(10) << "Iter :" << iter  <<  " " << primal_sfm  << " " << dual_sfm << " " << primal_sfm - dual_sfm << " " << stats.back().nSFM << std::endl;
        if(primal_sfm - dual_sfm < 1e-2)
            break;

    }
}

void BCD2D1D6N3D::optimize_acc(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname)
{   
    size_t w = this->w;
    size_t h = this->h;
    size_t d = this->d;

    double **z1 = this->z;
    double **z2 = z1+(w*d);
 
    double **ta1 = this->ta;
    double **ta2 = ta1+(w*d);

    double **tb1 = this->tb;
    double **tb2 = tb1+(w*d);

    const double **W1 = this->W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);

    double **u2D1 = this->u2D;
    double **u2D2 = u2D1+(w*d);

    maxIter = maxIterations;
    
    double beta;

    size_t nIterF, nIterL;
    for(int iter =0; iter < maxIter; ++iter)
    {
        nIterF = 0;
        nIterL = 0;

	beta = (iter - 1.0)/(iter + 2.0);

        projF->parallelProject(z1, ta1, tb1, W1, W2, u2D1);
        parallelAccumulate(2);

        for(size_t i=0; i < nThreadsF; ++i)
            nIterF = (nIterF > projF->iter[i])? nIterF: projF->iter[i];
                                                      
        projL->parallelProject(z2, ta2, tb2, W3, u2D2);
        parallelAccumulate(3, beta);
                                                      
        for(size_t i=0; i < nThreadsL; ++i)
            nIterL = (nIterL > projL->iter[i])? nIterL: projL->iter[i];

        computeGapsUpdateSolution(iter, stats);
        
        stats.back().nSFM = nIterF;

        double primal_sfm = stats.back().primalD;
        double dual_sfm   = stats.back().dualD;
        std::cerr << std::setprecision(10) << "Iter :" << iter  <<  " " << primal_sfm  << " " << dual_sfm << " " << primal_sfm - dual_sfm << " " << stats.back().nSFM << std::endl;
        if(primal_sfm - dual_sfm < 1e-2)
            break;

    }
}

//void BCD2D1D6N3D::optimize(double ** start, const size_t maxIterations, std::vector<OptInfo> & stats)
//{
//    // Initialize the tb with the warm start
//    double **s1 = start;
//    double **s2 = s1+(d*h);
//    double **s3 = s2+(w*d);
//
//    double **tb1 = tb;
//    double **tb2 = tb1+(d*h);
//    double **tb3 = tb2+(w*d);
//
//    for(int i=0; i < d*h; ++i)
//        std::copy(s1[i], s1[i]+w, tb1[i]);
//    for(int i=0; i < w*d; ++i)
//        std::copy(s2[i], s2[i]+h, tb2[i]);
//    for(int i=0; i < w*h; ++i)
//        std::copy(s3[i], s3[i]+h, tb3[i]);
//
//    optimize(maxIterations, stats, NULL);
//}

//TODO : check indexing consistency with mex
// size(x) = (w*d, h) + (w*h, d) size(u) = (w, h, d)
double BCD2D1D6N3D::duplicate(double **x, double *** u)
{
    double **x1 = x;
    double **x2 = x1+(w*d);

    double mean = 0.0;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            for(int k=0; k < d; ++k)
            {
                double val = u[i][j][k];
                x1[k*w + i][j] = val;
                x2[j*w + i][k] = val;

                mean += val;
            }
    mean /= (w*d*h);
    return mean;
}
