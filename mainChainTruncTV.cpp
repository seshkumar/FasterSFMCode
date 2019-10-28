#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/time.h>
#include "optInfo.h"
#include "graph4N2D.h"
#include "activeSetTruncChain.h"

inline double cval(double **x, double **u, double **W, size_t nThreads, size_t length)
{
    double val = 0.0f;
    for(int i=0; i < nThreads; ++i)
        for(int j=0; j < length; ++j)
        {
            val -= u[i][j]*x[i][j];
            if(j < length-1)
                val+= W[i][j]*fabs(x[i][j] - x[i][j+1]);
        }
            
    return val;
};

inline double frobeniusNorm(double **y, size_t nThreads, size_t length)
{
    double result = 0.0;
    for(int i = 0; i < nThreads; ++i)
    {
        for(int j = 0; j < length; ++j)
        {
            double value = y[i][j];
            result += value * value;
        }
    }
    return 0.5*result;
};
inline double getDual(double ** s, size_t nThreads, size_t length, double epsilon)
{
    double dual = 0.0f;
    for(int i = 0; i < nThreads; ++i)
    {
        for(int j = 0; j < length; ++j)
        {
            double val = fabs(s[i][j]);
            if(val <= epsilon)
                dual -= (val*val)/(2*epsilon);
            else
                dual -= (val - epsilon/2);
        }
    }
    return dual; 
}

inline double cval(double *x, double *u, double *W, size_t length)
{
    double val = 0.0f;
    for(int j=0; j < length; ++j)
    {
        val -= u[j]*x[j];
        if(j < length-1)
            val+= W[j]*fabs(x[j] - x[j+1]);
    }
            
    return val;
};

inline double frobeniusNorm(double *y, size_t length)
{
    double result = 0.0;
    for(int j = 0; j < length; ++j)
    {
        double value = y[j];
        result += value * value;
    }
    return 0.5*result;
};

inline double getDual(double *s, size_t length, double epsilon)
{
    double dual = 0.0f;
    for(int j = 0; j < length; ++j)
    {
        double val = fabs(s[j]);
        if(val <= epsilon)
            dual -= (val*val)/(2*epsilon);
        else
            dual -= (val - epsilon/2);
    }
    return dual; 
}






int main(int argc, char ** argv)
{
    timeval t0, t1, t2, t3;

    gettimeofday(&t0,0);

    Graph4N2D g;
    g.readDimacs(argv[1]);

    gettimeofday(&t1,0);

    double epsilon = std::stod(argv[2]);

    size_t nThreads = g.dims[0];
    size_t length   = g.dims[1];

    ActiveSetTruncChain **d = new ActiveSetTruncChain*[g.dims[0]];

    #pragma omp parallel    
    { 
        #pragma omp for
        for(int i=0; i < nThreads; ++i)
            d[i] = new ActiveSetTruncChain(g.dims[1]);
    }

    int maxiter=100;
    if(argc >=3)
    {
        maxiter = atoi(argv[2]);
    } 
    std::vector<OptInfo> stats;

    gettimeofday(&t2,0);

    #pragma omp parallel    
    { 
        #pragma omp for
       for(int i=0; i < nThreads; ++i)
    //int i = 2;
    //std::cout << "u = ";
    //for(int j=0; j < length; ++j)
    //    std::cout << g.u[i][j] << " ";
    //std::cout << std::endl;
    //std::cout << "W = ";
    //for(int j=0; j < length-1; ++j)
    //    std::cout << g.W2[i][j] << " ";
    //std::cout << std::endl;

        {
            //std::cout << "Thread = " << i << std::endl;
            d[i]->initialize(g.y[i], g.s[i], g.W2[i], g.u[i], epsilon);
            d[i]->optimize();
        }
    }
    gettimeofday(&t3,0);

    for(int i=0; i < nThreads; ++i)
        delete d[i];
    delete[] d;
    print(stats);
    stats.clear();    
    for(int i=0; i < nThreads; ++i)
    {
        double primal = cval(g.y[i], g.u[i], g.W2[i], length) + frobeniusNorm(g.y[i], length)*epsilon;
        double dual   = getDual(g.s[i], length, epsilon);
        assert(fabs(primal - dual) < 1e-6);

        //std::cout << "w[" << i << "] = " ; 
        //for(int j=0; j < length; ++j)
        //    std::cout << g.y[i][j] << " ";
        //std::cout << std::endl;
        //std::cout << "s[" << i << "] = " ; 
        //for(int j=0; j < length; ++j)
        //    std::cout << g.s[i][j] << " ";
        //std::cout << std::endl;

    }

    //double primal = cval(g.y, g.u, g.W2, g.dims[0], g.dims[1]) + frobeniusNorm(g.y, g.dims[0], g.dims[1])*epsilon;
    //double dual   = getDual(g.s, g.dims[0], g.dims[1], epsilon);
    ////assert(fabs(primal - dual) < 1e-6);

    //std::cout << "Primal = " << primal << "; Dual = " << dual << "; Gap = " << primal - dual << std::endl;

    std::cout << "Reading time : "      << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << std::endl;
    std::cout << "Initializing time : " << (t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec << std::endl;
    std::cout << "Optimizing time : "   << (t3.tv_sec-t2.tv_sec)*1000000 + t3.tv_usec-t2.tv_usec << std::endl;
}

