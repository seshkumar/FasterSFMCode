#include <fstream>
#include <iostream>
#include <cstring>

#include "readDimacsMax3D.h"

using namespace std;

void readDimacsMax3D( char * filename, flowNetwork & f)
{
    ifstream inputStream(filename);
    if(!inputStream.good())
    {
        cerr << "Unable to read file :" << filename << endl;
        return;
    }
    
    bool error = false;

    int eIndex = 0;


    while(!inputStream.eof())
    {
        char c;
        inputStream >> c;
        switch(c)
        {
            case 'c':
            {
                char str[100];
                inputStream >> str;
                if(!strcmp(str, "regulargrid"))
                    inputStream >> f.width >> f.height >> f.depth;
                else if (!strcmp(str, "p") | !strcmp(str, "n") | !strcmp(str, "a"))
                    inputStream.unget();
                else
                    while(!inputStream.eof() && inputStream.get() != '\n');
                break;
            }
            case 'p':
            {
                char str[100];
                inputStream >> str >> f.nodes >> f.edges;
                if(strcmp(str, "max"))
                    error = true;
                else
                    f.initializeEdgeArrays();
                break;
            }
            case 'n':
            {
                int n;
                char type;
                inputStream >> n >> type;
                switch (type)
                {
                    case 's': //source
                        f.terminal[0] = n;
                        break;
                    case 't': //sink
                        f.terminal[1] = n;
                        break;
                }
                break;
            }
            case 'a':
            {
                if (f.edges == 0 | f.nodes == 0)
                    error = true;
                else if(eIndex < f.edges)
                {
                    inputStream >> f.u[eIndex] >> f.v[eIndex] >> f.c[eIndex];
                    // To be consistent with array notation in C/C++
                    f.u[eIndex] -= 1;
                    f.v[eIndex] -= 1;
                    ++eIndex;
                }
                break;
            }
            default:
                error = true;
        } 
        
        if(error)
        {
            cerr << "ERROR: input file " << filename << "not correct" << endl;
            return;
        }
    }
    
    if(!f.terminal[0] || !f.terminal[1])
    {
        cerr << "ERROR: source node not defined" << filename << endl;
    }
}

