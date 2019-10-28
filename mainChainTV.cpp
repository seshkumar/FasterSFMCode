#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include "optInfo.h"
#include "graph4N2D.h"
#include "activeSetChain_base.h"
int main(int argc, char ** argv)
{
    timeval t0, t1, t2, t3;

    gettimeofday(&t0,0);

    Graph4N2D g;
    g.readDimacs(argv[1]);

    std::cout << "XXXXX   Filename = " << argv[1] << std::endl;

    gettimeofday(&t1,0);

    size_t nThreads = g.dims[0];
    ActiveSetChain_base **d = new ActiveSetChain_base*[g.dims[0]];
    #pragma omp parallel    
    { 
        #pragma omp for
        for(int i=0; i < nThreads; ++i)
            d[i] = new ActiveSetChain_base(g.dims[1]);
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
        {
            d[i]->initialize(g.y[i], g.s[i], g.W2[i], g.u[i]);
            d[i]->optimize();
        }
    }
    gettimeofday(&t3,0);

    for(int i=0; i < nThreads; ++i)
        delete d[i];
    delete[] d;
    print(stats);
    stats.clear();    

    std::cout << "Reading time : "      << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << std::endl;
    std::cout << "Initializing time : " << (t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec << std::endl;
    std::cout << "Optimizing time : "   << (t3.tv_sec-t2.tv_sec)*1000000 + t3.tv_usec-t2.tv_usec << std::endl;
}

