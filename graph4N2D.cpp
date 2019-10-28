#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cstddef>
#include <cstring>
#include "graph4N2D.h"
#include "flowNetwork.h"
#include "readDimacsMax2D.h"

void Graph4N2D::initialize(size_t w, size_t h)
{
    dims = new size_t[2];
    dims[0] = w;
    dims[1] = h;

    this->w = w;
    this->h = h;

    u = new double*[w];
    y = new double*[w];
    s = new double*[w];

    W1 = new double*[h];
    W2 = new double*[w];

    u[0] = new double[w*h];
    y[0] = new double[w*h];
    s[0] = new double[w*h];

    W1[0] = new double[h*(w-1)];
    W2[0] = new double[w*(h-1)];

    P    = new size_t*[w];
    P[0] = new size_t[w*h];

    np = 1;

    for (size_t i=1; i < w; ++i)
    {
        u[i] = u[0] + i*h;
        y[i] = y[0] + i*h;
        s[i] = s[0] + i*h;
        P[i] = P[0] + i*h;
        W2[i] = W2[0] + i*(h-1);
    } 
    
    for (size_t i=1; i < h; ++i)
        W1[i] = W1[0] + i*(w-1);

    std::memset(u[0], 0, w*h*sizeof(double));
    std::memset(y[0], 0, w*h*sizeof(double));
    std::memset(s[0], 0, w*h*sizeof(double));
    std::memset(W1[0], 0, h*(w-1)*sizeof(double));
    std::memset(W2[0], 0, w*(h-1)*sizeof(double));
    std::memset(P[0], 0, w*h*sizeof(size_t));
}

void Graph4N2D::readWeightsFromFile(char * fname)
{
    std::ifstream inputStream(fname);
    if(!inputStream.good())
    {
        std::cerr << "Unable to read file : " << fname << std::endl;
        return;
    }

    double wf;
    double hf;

    inputStream >> wf;
    inputStream >> hf;

    std::cout << wf << " " << hf << std::endl;

    initialize(static_cast<int>(wf), static_cast<int>(hf));

    size_t w = static_cast<size_t>(wf);
    size_t h = static_cast<size_t>(hf);

    for(size_t i=0; i < w*h; ++i)
            inputStream >> u[i%w][i/w];

    for(size_t i=0; i < h*(w-1); ++i)
            inputStream >> W1[i/(w-1)][i%(w-1)];

    for(size_t i=0; i < w*(h-1); ++i)
            inputStream >> W2[i/(h-1)][i%(h-1)];

    inputStream.close();
}

void Graph4N2D::indextoSubscript(int &x, int &y, int v)
{
    x = v%w;
    y = v/w;
}

void Graph4N2D::print (double **W, size_t w, size_t h)
{
    for(size_t j=0; j < h; ++j)
    {
        for(size_t i=0; i < w; ++i)
        {
            std::cout<< std::setw(4) << W[j][i] << " " ;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Graph4N2D::readDimacs(char * fname)
{
    flowNetwork f;
    readDimacsMax2D(fname, f);
    initialize(f.width, f.height);

    int nd = 0;
    for(size_t i=0; i < f.edges; ++i)
    {
        int x1, y1;
        int x2, y2;
        // Assuming 0 and 1 and source and sink
        indextoSubscript(x1, y1, f.u[i]-2);
        indextoSubscript(x2, y2, f.v[i]-2);

        if(!f.c[i])
            continue;

        if(f.u[i] == 0)
        {
            u[x2][y2] = f.c[i];
        }
        else if(f.v[i] == 1)
        {
            u[x1][y1] = -f.c[i];
        }
        else if(f.v[i] == f.u[i] + 1)
        {
            assert(y1 == y2);
            W1[y1][x1] = f.c[i];          
        }
        else if(f.v[i] == f.u[i] + w)
        {
            assert(x1 == x2);
            W2[x1][y1] = f.c[i];          
        }
        else if(f.u[i] == f.v[i] + 1 )
        {
            assert(y1 == y2);
            W1[y2][x2] = f.c[i];          
        }
        else if(f.u[i] == f.v[i] + w)
        {
            assert(x1 == x2);
            W2[x2][y2] = f.c[i];          
        }
        else
            std::cerr << "Something wrong in edge reading" << std::endl;
    }  
}

Graph4N2D::~Graph4N2D(void)
{
    delete[] u[0];
    delete[] y[0];
    delete[] s[0];
    delete[] P[0];

    delete[] W1[0];
    delete[] W2[0];

    delete[] u;
    delete[] y;
    delete[] s;
    delete[] P;
    delete[] W1;
    delete[] W2;

    delete[] dims;
}
