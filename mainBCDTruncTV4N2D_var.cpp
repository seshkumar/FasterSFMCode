#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/time.h>
#include "cyclicBCDTruncTV4N2D.h"
#include "optInfo.h"

int main(int argc, char ** argv)
{
    timeval t0, t1, t2, t3;

    gettimeofday(&t0,0);

    Graph4N2D g;
    g.readDimacs(argv[1]);

    gettimeofday(&t1,0);

    double epsilon = std::stod(argv[2]);

    CyclicBCDTruncTV4N2D d(g, epsilon);
    int maxiter=5000;
    std::vector<OptInfo> stats;

    gettimeofday(&t2,0);

    d.optimize_var(maxiter, stats, argv[2]);

    gettimeofday(&t3,0);

    printAll(stats);
    stats.clear();

    //std::cout << "Reading time : "      << (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec << std::endl;
    //std::cout << "Initializing time : " << (t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec << std::endl;
    //std::cout << "Optimizing time : "   << (t3.tv_sec-t2.tv_sec)*1000000 + t3.tv_usec-t2.tv_usec << std::endl;
}

