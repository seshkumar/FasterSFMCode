#ifndef _GRAPH_4N2D_H_
#define _GRAPH_4N2D_H_

#include <cstddef>

class Graph4N2D
{
    public:
        double  **W1;
        double  **W2;
        double  **u;
        double  **y;
        double  **s;
        size_t *dims;

        size_t **P;
        size_t np;

        Graph4N2D(void){}
        void initialize(size_t w, size_t h);
        ~Graph4N2D(void);
        void readDimacs(char * fname);
        void readWeightsFromFile(char * fname);

    private:
        size_t w;
        size_t h;

        void indextoSubscript(int &x, int &y, int v);
        void print (double  **W, size_t w, size_t h);
};
#endif
