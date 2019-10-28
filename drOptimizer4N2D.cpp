#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <parallel/algorithm>
#include "drOptimizer4N2D.h"

DROptimizer4N2D::DROptimizer4N2D(Graph4N2D &g)
{
    initialize(g.W1, g.W2, g.u, g.y, g.dims);
}

DROptimizer4N2D::DROptimizer4N2D(double **W1, double **W2, double **u, double **y, const size_t *dims)
{
    initialize(W1, W2, u, y, dims);
}

void DROptimizer4N2D::initialize(double **W1, double **W2, double **u, double **y, const size_t *dims)
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
    for(size_t i=0; i < h; ++i)
      C[i] = new bool[w]; 

    n = new size_t[nThreads]();
    size_t offset = 0;
    for(int i=offset; i < offset+h; ++i)
        n[i] = w;
    offset += h;
    for(int i=offset; i < offset+w; ++i)
        n[i] = h;

    this->u = u;
    this->y = y;

    projL = new DualLineProjector(nThreads, n);

    W = new const double*[nThreads]();

    std::copy(W1, W1+h, W   );
    std::copy(W2, W2+w, W+h );

    u2D  = new double*[nThreads]();
    s    = new double*[nThreads]();
    ta   = new double*[nThreads]();
    tb   = new double*[nThreads]();
    for(int i=0; i < nThreads; ++i)
    {
        u2D[i]  = new double[n[i]]();
        s[i]    = new double[n[i]]();
        ta[i]   = new double[n[i]]();
        tb[i]   = new double[n[i]]();
    }

    // convert3Dto2D returns mean of u
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

                tb[  j][i] = m/2;
                tb[h+i][j] = m/2;
            }
    }
}

DROptimizer4N2D::~DROptimizer4N2D()
{
    delete[] n;
    delete projL;
    delete[] W;
    for(int i=0; i < nThreads; ++i)
    {
        delete[] u2D[i];
        delete[] s[i];  
        delete[] ta[i];
        delete[] tb[i];
    }
    for(int i=0; i < h; ++i)
        delete[] C[i];

    delete[] C;
    delete[] u2D;
    delete[] s;
    delete[] ta;
    delete[] tb;

    nodes.clear();
    
    u = NULL;
}

void DROptimizer4N2D::parallelReflectLambda(double ** tb, double ** ta)
{
    size_t w = this->w;
    size_t h = this->h;

    double **ta1 = ta;
    double **ta2 = ta1+h;
    double **tb1 = tb;
    double **tb2 = tb1+h;

    size_t i, j;
    #pragma omp parallel
    {
        #pragma omp for
        for(i=0; i < w; ++i)
            for(j=0; j < h; ++j)
            {
                double meant =  (ta1[j][i] + ta2[i][j])/2;
                tb1[j][i] = (tb1[j][i] + ta1[j][i])*0.5 - meant;
                tb2[i][j] = (tb2[i][j] + ta2[i][j])*0.5 - meant;
            }
    }
}

void DROptimizer4N2D::computeGapsUpdateSolution(const size_t iter, std::vector<OptInfo> & stats)
{
    double **x1 = s;
    double **x2 = x1+h; 

    double **y = this->y;
    double **u = this->u;

    bool **C = this->C;

    const double **W = this->W;

    double dValueSfm = 0.0;
    double dValue    = 0.0;
    double pValue    = 0.0;

    double val       = 0.0;
    
    int w = this->w;
    int h = this->h;

    std::vector<node>::iterator it = nodes.begin();

    int i, j, idx;
    //#pragma omp parallel shared(x1, x2, y, u, C, it, w, h) private(i, j, val, idx) reduction(+: pValue, dValueSfm) reduction(-:dValue) default(none)
    {
        //#pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
            {
                //Initialize cut for greedyTV
                C[j][i] = false;

                val = (x1[j][i] + x2[i][j]);

                if(val < 0)
                    dValueSfm += val;

                dValue -= (val*val)/2;
                
                pValue += val*u[i][j] * (val*val)/2; 

                // Converting to Primal
                y[i][j]  = -val;
                x1[j][i] = val;
                x2[i][j] = val;

                int idx = j*w + i;
                (it+idx)->s = -val;
                (it+idx)->i = i;
                (it+idx)->j = j;            
            }
    }


    double diff = 0.0;
    const double **W1 = this->W;
    const double **W2 = W1+h;
    double p1, p2;

    //#pragma omp parallel shared(x1, x2, W1, W2, w, h) private(i, j, diff, p1, p2) reduction(+: pValue) default(none)
    {
        //#pragma omp for
        for(int i=0; i < w; ++i)
            for(int j=0; j < h; ++j)
            {
                p1 = 0.0;
                p2 = 0.0;
                if(i < w-1)
                {
                    diff = abs(x1[j][i] - x1[j][i+1]);
                    p1  += (W1[j][i] * diff);
                }
                if(j < h-1)
                {
                    diff = abs(x2[i][j] - x2[i][j+1]);
                    p2  += (W2[i][j] * diff);
                }
                pValue += (p1 + p2);
            }
    }

    __gnu_parallel::sort(nodes.begin(), nodes.end(), compare);

    double pValueSfm = greedyTV();

    stats.push_back(OptInfo(iter, pValue, dValue, pValueSfm, dValueSfm));
}

double DROptimizer4N2D::greedyTV(void)
{
    const double **W1 = W;
    const double **W2 = W+h;
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
        fmin = (f < fmin)? f: fmin;
    }
    return fmin;
}

void DROptimizer4N2D::optimize(const size_t maxIterations, std::vector<OptInfo> & stats)
{
    maxIter = maxIterations;

    //reflection on base polytopes 
    projL->parallelReflect(ta, s, tb, u2D, W);
    for(int iter =0; iter < maxIter; ++iter)
    {
        computeGapsUpdateSolution(iter, stats);
        double primal_sfm = stats.back().primalD;
        double dual_sfm   = stats.back().dualD;
        std::cout << std::setprecision(10) << "Iter :" << iter  <<  " " << primal_sfm  << " " << dual_sfm << std::endl;
        if(primal_sfm - dual_sfm < 0.00001)
            break;

        //reflection on lambda polytope
        parallelReflectLambda(tb, ta);

        //reflection on base polytopes 
        projL->parallelReflect(ta, s, tb, u2D, W);
    }
}

void DROptimizer4N2D::optimize(double ** start, const size_t maxIterations, std::vector<OptInfo> & stats)
{
    // Initialize the tb with the warm start
    duplicate(tb, start);
    optimize(maxIterations, stats);
}

// size(x) = (h, w) + (w, h) 
// size(u) = (w, h)
double DROptimizer4N2D::duplicate(double **x, double **u)
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
