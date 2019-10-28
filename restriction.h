#ifndef __RESTRICTION_H_
#define __RESTRICTION_H_

#include <cstddef>

void restriction2D(double **uR, double **W1R, double **W2R,const double *const*u, const double *const*W1, const double *const*W2, const size_t *const*S, const size_t w, const size_t h);
void restriction(double * uR, double * wR, const double * u, const double * w, const size_t * S, const size_t * V, const size_t n);

#endif
