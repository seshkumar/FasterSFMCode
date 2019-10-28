#include <cstring>
#include "flowNetwork.h"

flowNetwork::flowNetwork(void)
{
    width  = 0;
    height = 0;
    depth  = 0;

    nodes = 0;
    edges = 0;

    terminal = new int[2];
}

flowNetwork::~flowNetwork(void)
{
    delete[] terminal;
    
    delete[] u;
    delete[] v;
    delete[] c;
}

void flowNetwork::initializeEdgeArrays(void)
{
    u = new int[edges];
    v = new int[edges];
    c = new int[edges];
}

