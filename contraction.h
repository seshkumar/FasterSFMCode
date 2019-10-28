#ifndef __CONTRACTION_H_
#define __CONTRACTION_H_

#include <cstddef>

void contraction2D(double **uC, double **W1C, double **W2C,const double *const*u, const double *const*W1, const double *const*W2, const size_t *const*S, const size_t w, const size_t h);
void contraction(double * uC, double * wC, const double * u, const double * w, const size_t * S, const size_t * V, const size_t n);

#endif
