#ifndef FLOW_NETWORK_H_
#define FLOW_NETWORK_H_

class flowNetwork
{
    public:
        flowNetwork(void);
        ~flowNetwork(void);

        void initializeEdgeArrays(void);

        int width ;
        int height; 
        int depth ;

        int *u;
        int *v;
        int *c;

        int nodes;
        int edges;

        int *terminal;
};

#endif
