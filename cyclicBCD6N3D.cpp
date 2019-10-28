#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <parallel/algorithm>
//#include <algorithm>
#include "cyclicBCD6N3D.h"

CyclicBCD6N3D::CyclicBCD6N3D(Graph6N3D &g)
{
    initialize(g.y, g.W1, g.W2, g.W3, g.u, g.dims);
}

CyclicBCD6N3D::CyclicBCD6N3D(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims)
{
    initialize(y, W1, W2, W3, u, dims);
}

void CyclicBCD6N3D::initialize(double ***y, double **W1, double **W2, double **W3, double ***u, const size_t *dims)
{
    w = dims[0];
    h = dims[1];
    d = dims[2];
    nThreads = d*h + w*d + w*h;

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

    n = new size_t[nThreads];
    size_t offset = 0;
    for(int i=offset; i < offset+d*h; ++i)
        n[i] = w;
    offset += d*h;
    for(int i=offset; i < offset + w*d; ++i)
        n[i] = h;
    offset += w*d;
    for(int i=offset; i < offset + w*h; ++i)
        n[i] = d;

    this->u = u;
    this->y = y;

    projL1 = new DualLineProjector(d*h, n                );
    projL2 = new DualLineProjector(w*d, n + (d*h)        );
    projL3 = new DualLineProjector(w*h, n + (d*h) + (w*d));

    W = new const double*[nThreads];

    std::copy(W1, W1+d*h, W            );
    std::copy(W2, W2+w*d, W + d*h      );
    std::copy(W3, W3+w*h, W + d*h + w*d);

    u2D  = new double*[nThreads];
    s    = new double*[nThreads];
    z    = new double*[nThreads];
    ta   = new double*[nThreads];
    tb   = new double*[nThreads];
    for(int i=0; i < nThreads; ++i)
    {
        u2D[i]  = new double[n[i]];
        s[i]    = new double[n[i]];
        z[i]    = new double[n[i]];
        ta[i]   = new double[n[i]];
        tb[i]   = new double[n[i]];
    }

    // duplicate returns mean of u
    mean = duplicate(u2D, u);

    // divide the unary vector by 3 and initializing tb with mean
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
                for(int k=0; k < d; ++k)
                {
                    u2D[            k*h + j][i] /= 3;
                    u2D[      d*h + k*w + i][j] /= 3;
                    u2D[w*d + d*h + j*w + i][k] /= 3;

                    tb[            k*h + j][i] = mean/3;
                    tb[      d*h + k*w + i][j] = mean/3;
                    tb[w*d + d*h + j*w + i][k] = mean/3;
                }
    }
}

CyclicBCD6N3D::~CyclicBCD6N3D()
{
    delete[] n;
    delete projL1;
    delete projL2;
    delete projL3;
    delete[] W;
    for(int i=0; i < nThreads; ++i)
    {
        delete[] u2D[i];
        delete[] s[i];  
        delete[] z[i];  
        delete[] ta[i];
        delete[] tb[i];
    }
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


    delete[] u2D;
    delete[] s;
    delete[] z;
    delete[] ta;
    delete[] tb;
    
    nodes.clear();

    u = NULL;
}

void CyclicBCD6N3D::parallelAccumulate(double ** tb, double ** ta, size_t dir)
{
    size_t w = this->w;
    size_t h = this->h;
    size_t d = this->d;
 
    double **ta1 = ta;
    double **ta2 = ta1+(d*h);
    double **ta3 = ta2+(w*d);
    double **tb1 = tb;
    double **tb2 = tb1+(d*h);
    double **tb3 = tb2+(w*d);
  
    size_t i, j, k;

    switch(dir)
    {
        case 1:
            #pragma omp parallel shared(tb1, ta2, ta3, i, j, k, w, h, d) default(none)
            {
                #pragma omp for
                for(i = 0; i < w; ++i)
                    for(j = 0; j < h; ++j)
                        for(k = 0; k < d; ++k)
                            tb1[k*h + j][i] = -(ta2[k*w + i][j] + ta3[j*w + i][k]);
            }
            break;
        case 2:
            #pragma omp parallel shared(tb2, ta1, ta3, i, j, k, w, h, d) default(none)
            {
                #pragma omp for
                for(j = 0; j < h; ++j)
                    for(i = 0; i < w; ++i)
                        for(k = 0; k < d; ++k)
                            tb2[k*w + i][j] = -(ta1[k*h + j][i] + ta3[j*w + i][k]);
            }
            break;
        case 3:
            #pragma omp parallel shared(tb3, ta1, ta2, i, j, k, w, h, d) default(none)
            {
                #pragma omp for
                for(k = 0; k < d; ++k)
                    for(i = 0; i < w; ++i)
                        for(j = 0; j < h; ++j)
                            tb3[j*w + i][k] = -(ta1[k*h + j][i] + ta2[k*w + i][j]);
            }
            break;
    }
}

void CyclicBCD6N3D::computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> &stats)
{
    double **x1 = tb;
    double **x2 = x1+(d*h);
    double **x3 = x2+(w*d);

    double **s1 = s;
    double **s2 = s1+(d*h);
    double **s3 = s2+(w*d);

    double ***y = this->y;
    double ***u = this->u;

    bool ***C = this->C;

    double dValueSfm = 0.0;
    double dValue    = 0.0;
    double pValue    = 0.0;

    double val = 0.0;

    int w = this->w;
    int h = this->h;
    int d = this->d;
    
    std::vector<node>::iterator it = nodes.begin();

    int i, j, k, idx;
    #pragma omp parallel shared(x1, x2, x3, s1, s2, s3, y, u, C, it, w, h, d) private(i, j, k, val, idx) reduction(+: pValue, dValueSfm) reduction(-:dValue) default(none)
    {
        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
                for(int k=0; k < d; ++k)
                {
                    C[k][j][i] = false;

                    val = (x1[k*h + j][i] + x2[k*w + i][j] + x3[j*w + i][k]);

                    if(val < 0)
                        dValueSfm += val;

                    dValue -= (val*val)/2;

                    pValue += (val*u[i][j][k] + (val*val)/2);

                    // Converting to Primal
                    y[i][j][k]     = -val;
                    s1[k*h + j][i] = val;
                    s2[k*w + i][j] = val;
                    s3[j*w + i][k] = val;

                    int idx = k*w*h + j*w + i;
                    (it+idx)->s = -val;
                    (it+idx)->i = i;
                    (it+idx)->j = j;
                    (it+idx)->k = k;
                }
    }

    __gnu_parallel::sort(nodes.begin(), nodes.end(), compare);
    //std::sort(nodes.begin(), nodes.end(), compare);

    double pValueSfm = greedyTV();

    int jDistNum = 0;
    int jDistDen = 0;

    double diff = 0.0;
    const double **W1 = this->W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);
    double p1, p2, p3;
    #pragma omp parallel shared(s1, s2, s3, y, W1, W2, W3, w, h, d) private(i, j, diff, p1, p2, p3) reduction(+: pValue, jDistNum, jDistDen) default(none)
    {
        #pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
                for(int k=0; k < d; ++k)
                {
                    jDistNum += (y[i][j][k] >= threshold && optimal[i][j][k])?1:0;
                    jDistDen += (y[i][j][k] >= threshold || optimal[i][j][k])?1:0;

                    p1 = 0.0;
                    p2 = 0.0;
                    p3 = 0.0;
                    if(i < w-1)
                    {
                        diff = abs(s1[k*h + j][i] - s1[k*h +j][i+1]);
                        p1  += (W1[k*h + j][i] * diff);
                    }
                    if(j < h-1)
                    {
                        diff = abs(s2[k*w +i][j] - s2[k*w + i][j+1]);
                        p2  += (W2[k*w +i][j] * diff);
                    }
                    if(k < d-1)
                    {
                        diff = abs(s3[j*w +i][k] - s3[j*w + i][k+1]);
                        p3  += (W3[j*w +i][k] * diff);
                    }
                    pValue += (p1 + p2 + p3);
                }
    }


    stats.push_back(OptInfo(iter, pValue, dValue, pValueSfm, dValueSfm, jDistNum, jDistDen));
}

double CyclicBCD6N3D::greedyTV(void)
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

void CyclicBCD6N3D::optimize(const size_t maxIterations, std::vector<OptInfo> & stats, char * fname)
{   
    size_t w = this->w;
    size_t h = this->h;
    size_t d = this->d;

    double **z1 = z;
    double **z2 = z1+(d*h);
    double **z3 = z2+(w*d);
 
    double **ta1 = ta;
    double **ta2 = ta1+(d*h);
    double **ta3 = ta2+(w*d);

    double **tb1 = tb;
    double **tb2 = tb1+(d*h);
    double **tb3 = tb2+(w*d);

    const double **W1 = this->W;
    const double **W2 = W1+(d*h);
    const double **W3 = W2+(w*d);

    double **u2D1 = this->u2D;
    double **u2D2 = u2D1+(d*h);
    double **u2D3 = u2D2+(w*d);

    maxIter = maxIterations;
    
    std::ifstream f(fname);
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            for(int k=0; k < d; ++k)
                f >> optimal[i][j][k];
    f.close();

    size_t nIter1, nIter2, nIter3;
    for(int iter =0; iter < maxIter; ++iter)
    {
        nIter1 = 0;
        nIter2 = 0;
        nIter3 = 0;

        parallelAccumulate(ta, tb, 1);
        projL1->parallelProject(z1, tb1, ta1, W1, u2D1);

        for(size_t i=0; i < d*h; ++i)
            nIter1 += projL1->iter[i];
                                                      
        parallelAccumulate(ta, tb, 2);                
        projL2->parallelProject(z2, tb2, ta2, W2, u2D2);
                                                      
        for(size_t i=0; i < w*d; ++i)
            nIter2 += projL2->iter[i];

        parallelAccumulate(ta, tb, 3);                
        projL3->parallelProject(z3, tb3, ta3, W3, u2D3);
        
        for(size_t i=0; i < w*h; ++i)
            nIter3 += projL3->iter[i];

        computeGapsUpdateSolution(iter, stats);
        
        stats.back().nSFM = (nIter1 + nIter2 + nIter3)/nThreads;

        double primal_sfm = stats.back().primalD;
        double dual_sfm   = stats.back().dualD;
        //std::cerr << std::setprecision(10) << "Iter :" << iter  <<  " " << primal_sfm  << " " << dual_sfm << " " << primal_sfm - dual_sfm << " " << stats.back().nSFM << std::endl;
        if(abs(primal_sfm - dual_sfm) < 1e-3)
            break;
    }
}

void CyclicBCD6N3D::optimize(double ** start, const size_t maxIterations, std::vector<OptInfo> & stats)
{
    // Initialize the tb with the warm start
    double **s1 = start;
    double **s2 = s1+(d*h);
    double **s3 = s2+(w*d);

    double **tb1 = tb;
    double **tb2 = tb1+(d*h);
    double **tb3 = tb2+(w*d);

    for(int i=0; i < d*h; ++i)
        std::copy(s1[i], s1[i]+w, tb1[i]);
    for(int i=0; i < w*d; ++i)
        std::copy(s2[i], s2[i]+h, tb2[i]);
    for(int i=0; i < w*h; ++i)
        std::copy(s3[i], s3[i]+h, tb3[i]);

    optimize(maxIterations, stats, NULL);
}

//TODO : check indexing consistency with mex
// size(x) = (w, d*h) + (h, w*d) + (d, w*h) size(u) = (w, h, d)
double CyclicBCD6N3D::duplicate(double **x, double *** u)
{
    double **x1 = x;
    double **x2 = x1+(d*h);
    double **x3 = x2+(w*d);

    double mean = 0.0;
    for(int i=0; i < w; ++i)
        for(int j=0; j < h; ++j)
            for(int k=0; k < d; ++k)
            {
                double val = u[i][j][k];
                x1[k*h + j][i] = val;
                x2[k*w + i][j] = val;
                x3[j*w + i][k] = val;

                mean += val;
            }
    mean /= (w*d*h);
    return mean;
}
