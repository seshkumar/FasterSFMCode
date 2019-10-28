#ifndef _GRAPH_6N3D_H_
#define _GRAPH_6N3D_H_

#include <cstddef>

class Graph6N3D
{
    public:
        double  **W1;
        double  **W2;
        double  **W3;
        double  ***u;
        double  ***y;
        size_t *dims;

        Graph6N3D(void){}
        void initialize(size_t w, size_t h, size_t d);
        ~Graph6N3D(void);
        void readDimacs(char * fname, double gamma=1.0);
        double calculateEps(void);

    private:
        size_t w;
        size_t h;
        size_t d;

        void indextoSubscript(int &x, int &y, int &z, int v);
        void print (double  **W, size_t w, size_t h);
        void print (double ***u, size_t w, size_t h, size_t d);
};
#endif
