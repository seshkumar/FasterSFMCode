#include "contraction.h"

void contraction(double * uC, double * wC, const double * u, const double * w, const size_t * S, const size_t * V, const size_t n)
{
    for (int i = 0; i < n; ++i)
        uC[i] = S[i] ? 0 : u[i];
    
    for (int i = 0; i < n-1; ++i)
    {
        double w_i = w[i];
        wC[i] = (S[i] && S[i+1]) ? 0 : w_i;

        int ds_i = S[i+1] - S[i];
        if(S[i+1] && !S[i])      // left case
        {
            uC[i] += w_i;
            wC[i] =  0;
        }
        else if(!S[i+1] && S[i]) // right case
        {
            uC[i+1] += w_i;
            wC[i]   = 0;
        }

    }
}
