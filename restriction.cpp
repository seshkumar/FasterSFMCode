#include "restriction.h"

void restriction(double * uR, double * wR, const double * u, const double * w, const size_t * S, const size_t * V, const size_t n)
{
    size_t * Sc = new size_t[n];
    for (int i = 0; i < n ; ++i)
    {
        Sc[i] = V[i] - S[i];
        uR[i] = Sc[i] ? 0 : u[i];
    }
    for (int i = 0; i < n-1; ++i)
    {
        double w_i = w[i];
        wR[i] = (Sc[i] && Sc[i+1]) ? 0 : w_i;
        
        if(Sc[i+1] && !Sc[i])      // left case
        {
            uR[i] -= w_i;
            wR[i] = 0;
        }
        else if(!Sc[i+1] && Sc[i]) // right case
        {
            uR[i+1] -= w_i;
            wR[i] = 0;
        }
    }
    delete[] Sc;
}
