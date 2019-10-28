#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstddef>
#include <cmath>
#include "graph6N3D.h"
#include "flowNetwork.h"
#include "readDimacsMax3D.h"

void Graph6N3D::initialize(size_t w, size_t h, size_t d)
{
    dims = new size_t[3];
    dims[0] = w;
    dims[1] = h;
    dims[2] = d;

    this->w = w;
    this->h = h;
    this->d = d; 

    u = new double**[w];
    y = new double**[w];
    for (size_t i=0; i < w; ++i)
    {
        u[i] = new double*[h];
        y[i] = new double*[h];
        for(size_t j=0; j < h; ++j)
        {
            u[i][j] = new double[d];
            y[i][j] = new double[d];
                for(size_t k=0; k < d;++k)
                {
                    u[i][j][k] = 0.0;
                    y[i][j][k] = 0.0;
                }
        }
    } 

    W1 = new double*[h*d]();
    W2 = new double*[w*d]();
    W3 = new double*[w*h]();

    W1[0] = new double[(w-1)*h    *d    ]();
    W2[0] = new double[w    *(h-1)*d    ]();
    W3[0] = new double[w    *h    *(d-1)]();
    for (size_t i=1; i < h*d; ++i)
        W1[i] = W1[i-1] + (w-1);
    for(size_t i=1; i < w*d; ++i)
        W2[i] = W2[i-1] + (h-1);
    for(size_t i=1; i < w*h; ++i)
        W3[i] = W3[i-1] + (d-1);

}

void Graph6N3D::indextoSubscript(int &x, int &y, int &z, int v)
{
    x = v%w;
    y = (v/w)%h;
    z = v/(w*h);
}

void Graph6N3D::print (double **W, size_t w, size_t h)
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

void Graph6N3D::print (double ***u, size_t w, size_t h, size_t d)
{
    for(size_t k=0; k < d; ++k)
    {
        for(size_t j=0; j < h; ++j)
        {
            for(size_t i=0; i < w; ++i)
            {
                std::cout << std::setw(4) << u[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void Graph6N3D::readDimacs(char * fname, double gamma)
{
    flowNetwork f;
    readDimacsMax3D(fname, f);
    initialize(f.width, f.height, f.depth);

    for(size_t i=0; i < f.edges; ++i)
    {
        int x1, y1, z1;
        int x2, y2, z2;
        // Assuming 0 and 1 and source and sink
        indextoSubscript(x1, y1, z1, f.u[i]-2);
        indextoSubscript(x2, y2, z2, f.v[i]-2);

        if(!f.c[i])
            continue;

        if(f.u[i] == 0)
        {
            u[x2][y2][z2] = f.c[i];
        }
        else if(f.v[i] == 1)
        {
            u[x1][y1][z1] = -f.c[i];
        }
	else
	{
		// for calculating epsilon 
		y[x1][y1][z1] += f.c[i];
		y[x2][y2][z2] += f.c[i];
		// ends here
        	if(f.v[i] == f.u[i] + 1)
        	{
        	    assert(y1 == y2 && z1 == z2);
        	    W1[z1*h + y1][x1] = round(gamma*f.c[i]);          
        	}
        	else if(f.v[i] == f.u[i] + w)
        	{
        	    assert(x1 == x2 && z1 == z2);
        	    W2[z1*w + x1][y1] = round(gamma*f.c[i]);          
        	}
        	else if(f.v[i] == f.u[i] + w*h)
        	{
        	    assert(x1 == x2 && y1 == y2);
        	    W3[y1*w + x1][z1] = round(gamma*f.c[i]);
        	}
        	else if(f.u[i] == f.v[i] + 1)
        	{
        	    assert(y1 == y2 && z1 == z2);
        	    W1[z2*h + y2][x2] = round(gamma*f.c[i]);          
        	}
        	else if(f.u[i] == f.v[i] + w)
        	{
        	    assert(x1 == x2 && z1 == z2);
        	    W2[z2*w + x2][y2] = round(gamma*f.c[i]);          
        	}
        	else if(f.u[i] == f.v[i] + w*h)
        	{
        	    assert(x1 == x2 && y1 == y2);
        	    W3[y2*w + x2][z2] = round(gamma*f.c[i]);
        	}
        	else
        	    std::cerr << "Something wrong in edge reading" << std::endl;
	}
    }  
}

double Graph6N3D::calculateEps(void)
{
    double eps;
    for(size_t k=0; k < d; ++k)
        for(size_t j=0; j < h; ++j)
            for(size_t i=0; i < w; ++i)
	    {
		eps += y[i][j][k]*y[i][j][k];
		y[i][j][k] = 0;
	    }
    eps = 2*sqrt(eps/(w*h*d));
    return eps;
}

Graph6N3D::~Graph6N3D(void)
{
    delete[] W1[0];
    delete[] W2[0];
    delete[] W3[0];
    delete[] W1;
    delete[] W2;
    delete[] W3;

    for (size_t i=0; i < w; ++i)
    {
        for (size_t j=0; j < h; ++j)
        {
            delete[] u[i][j];
            delete[] y[i][j];
        }
        delete[] u[i];
        delete[] y[i];
    }
    delete[] u;
    delete[] y;
    delete[] dims;
}
