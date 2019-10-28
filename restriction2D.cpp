#include <cstring>
#include "restriction.h"

void restriction2D(double **uR, double **W1R, double **W2R,const double *const*u, const double *const*W1, const double *const*W2, const size_t *const*S, const size_t w, const size_t h)
{
    std::memcpy(uR[0], u[0], w*h*sizeof(double));

    std::memset(W1R[0], 0, (w-1)*h*sizeof(double));
    std::memset(W2R[0], 0, w*(h-1)*sizeof(double));

    for(int i=0; i < w-1; ++i)
        for(int j=0; j < h-1; ++j)
        {
            if(S[i][j])
            {
                if(S[i][j+1]) //SV
                    W2R[i][j] = W2[i][j];
                else          //SScV
                    uR[i][j] -= W2[i][j]; 

                if(S[i+1][j])  //SH
                    W1R[j][i] = W1[j][i];
                else           //SScH
                    uR[i][j]  -= W1[j][i];
            }
            else
            {
                uR[i][j] = 0;
                if(S[i][j+1]) //ScSV
                    uR[i][j+1] -= W2[i][j]; 
                //else          //ScV
                //{
                //}

                if(S[i+1][j]) //ScSH
                    uR[i+1][j]  -= W1[j][i];
                //else          //ScH
                //{
                //}

            }
        }
}
