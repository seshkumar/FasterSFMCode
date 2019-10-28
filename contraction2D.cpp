#include <cstring>
#include "contraction.h"

void contraction2D(double **uC, double **W1C, double **W2C,const double *const*u, const double *const*W1, const double *const*W2, const size_t *const*S, const size_t w, const size_t h)
{
    std::memcpy(uC[0], u[0], w*h*sizeof(double));

    std::memset(W1C[0], 0, (w-1)*h*sizeof(double));
    std::memset(W2C[0], 0, w*(h-1)*sizeof(double));

    for(int i=0; i < w-1; ++i)
        for(int j=0; j < h-1; ++j)
        {
            if(S[i][j])
            {
                uC[i][j] = 0;
                //if(S[i][j+1])  //SV
                if(!S[i][j+1]) //SScV
                    uC[i][j+1] += W2[i][j]; 

                //if(S[i+1][j])  //SH
                if(!S[i+1][j]) //SScH
                    uC[i+1][j]  += W1[j][i];
            }
            else
            {
                if(S[i][j+1]) //ScSV
                    uC[i][j] += W2[i][j]; 
                else          //ScV
                    W2C[i][j] = W2[i][j];

                if(S[i+1][j]) //ScSH
                    uC[i][j]  += W1[j][i];
                else          //ScH
                    W1C[j][i] = W1[j][i];

            }
        }
}
