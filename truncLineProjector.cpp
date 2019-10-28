#include <iostream>
#include "truncLineProjector.h"

TruncLineProjector::TruncLineProjector(size_t nThreads, const size_t *n)
              :nThreads(nThreads)
{
    this->n = new size_t[nThreads]();
    std::copy(n, n+nThreads, this->n);

    iter   = new size_t[nThreads];
    lineTV = new ActiveSetTruncChain*[nThreads];
    for(int i=0; i < nThreads; ++i)
        lineTV[i] = new ActiveSetTruncChain(n[i]);

}

TruncLineProjector::~TruncLineProjector(void)
{
    delete[] n;
    delete[] iter;

    for(int i=0; i < nThreads; ++i)
        delete lineTV[i];

    delete[] lineTV;
}

void TruncLineProjector::parallelProject(double **y, double **s, const double *const*W, const double *const*u, double epsilon)
{
    size_t i, j;

    size_t * iter = this->iter;

    #pragma omp parallel shared(y, s, W, u, iter, epsilon) private(i, j) default(none)
    {
        #pragma omp for
        for(i=0; i < nThreads; ++i)
        {
            iter[i] = project(y[i], s[i], W[i], u[i], epsilon, i);
            for(int j=0; j < n[i]; ++j)
                s[i][j] += u[i][j];
        }
    }
}

size_t TruncLineProjector::project(double *y, double *s, const double *W, const double *u, double epsilon, size_t threadID)
{
    lineTV[threadID]->initialize(y, s, W, u, epsilon);
    size_t iter = lineTV[threadID]->optimize();
    return iter;
}

//// x is the labels
//// y is the input
//// lambda are the weights
//void TruncLineProjector::project(double *x, const double *y, const double *lambda, const size_t n) 
//{
//    /* Minorant and minorant slopes */
//    double mn(0), mx(0);
//    /* Relative height of majorant and minorant slopes at the current points w.r.t. the tube center */
//    double mnHeight(0), mxHeight(0);
//    /* Last break points of minorant and majorant */
//    int mnBreak(0), mxBreak(0);
//    /* Last break point of taut string */
//    int lastBreak(0);
//    /* Auxiliary variables */
//    int i, j;
//        
//    /* Starting point */
//    mnHeight = mxHeight = 0;
//    mn = -lambda[0] + y[0];
//    mx = lambda[0] + y[0];
//    lastBreak = -1;
//    mnBreak = mxBreak = 0;
//        
//    /* Proceed along string */
//    i = 0;
//    while ( i < n ) {
//        /* Loop over all points except the last one, that needs special care */
//        while ( i < n-1 ) {
//            #ifdef DEBUG
//                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//            #endif
//            
//            /* Update height of minorant slope w.r.t. tube center */
//            /* This takes into account both the slope of the minorant and the change in the tube center */
//            mnHeight += mn - y[i];
//        
//            /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
//            /* Majorant is r + lambda (except for last point), which is computed on the fly */   
//            if ( lambda[i] < mnHeight ) {
//                #ifdef DEBUG
//                    fprintf(DEBUG_FILE,"CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,lambda[i],mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//                #endif
//                /* Break segment at last minorant breaking point */
//                i = mnBreak + 1;
//                /* Build valid segment up to this point using the minorant slope */
//                for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
//                    x[j] = mn;
//                /* Start new segment after the break */
//                lastBreak = mnBreak;
//                /* Build first point of new segment, which can be done in closed form */
//                mn = y[i] + lambda[i-1] - lambda[i];  // When computing the slopes we need to account for the differences in lambdas
//                mx = y[i] + lambda[i-1] + lambda[i]; // Here too
//                mxHeight = lambda[i];
//                mnHeight = -lambda[i];
//                mnBreak = mxBreak = i;
//                i++;
//                continue;
//            }
//            
//            /* Update height of minorant slope w.r.t. tube center */
//            /* This takes into account both the slope of the minorant and the change in the tube center */
//            mxHeight += mx - y[i];
//            
//            /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
//            /* Minorant is r - lambda (except for last point), which is computed on the fly */
//            if ( -lambda[i] > mxHeight ) {
//                #ifdef DEBUG
//                    fprintf(DEBUG_FILE,"FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,-lambda[i],mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//                #endif
//                /* If violated, break segment at last majorant breaking point */
//                i = mxBreak + 1;
//                /* Build valid segment up to this point using the majorant slope */
//                for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
//                    x[j] = mx;
//                /* Start new segment after the break*/
//                lastBreak = mxBreak;
//                /* Build first point of new segment, which can be done in closed form */
//                mx = y[i] - lambda[i-1] + lambda[i]; // When computing the slopes we need to account for the differences in lambdas
//                mn = y[i] - lambda[i-1] - lambda[i]; // Here too
//                mxHeight = lambda[i];
//                mnHeight = -lambda[i];
//                mnBreak = mxBreak = i;
//                i++;
//                continue;
//            }
//            
//            /* No violations at this point */
//
//            /* Check if proyected majorant height is above ceiling */
//            if ( mxHeight >= lambda[i] ) {
//                /* Update majorant slope */
//                mx += ( lambda[i] - mxHeight ) / ( i - lastBreak );
//                /* Get correct majorant height (we are touching it!) */
//                mxHeight = lambda[i];
//                /* This is a possible majorant breaking point */
//                mxBreak = i;
//            }
//            
//            /* Check if proyected minorant height is under actual minorant */
//            if ( mnHeight <= -lambda[i] ) {
//                /* Update minorant slope */
//                mn += ( -lambda[i] - mnHeight ) / ( i - lastBreak );
//                /* Compute correct minorant height (we are touching it!) */
//                mnHeight = -lambda[i];
//                /* This is a possible minorant breaking point */
//                mnBreak = i;
//            }
//
//            #ifdef DEBUG
//                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//            #endif
//            
//            /* At this point: no violations, so keep up building current segment */
//            i++;
//        }
//        
//        /* Special case i == n-1 (last point) */
//        /* We try to validate the last segment, and if we can, we are finished */
//        /* The code is essentially the same as the one for the general case, 
//           the only different being that here the tube ceiling and floor are both 0 */
//        
//        /* Update height of minorant slope w.r.t. tube center */
//        /* This takes into account both the slope of the minorant and the change in the tube center */
//        mnHeight += mn - y[i];
//    
//        /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
//        /* Majorant is 0 at this point */   
//        if ( 0 < mnHeight ) {
//            #ifdef DEBUG
//                fprintf(DEBUG_FILE,"ENDING CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//            #endif
//            /* Break segment at last minorant breaking point */
//            i = mnBreak + 1;
//            /* Build valid segment up to this point using the minorant slope */
//            for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
//                x[j] = mn;
//            /* Start new segment after the break */
//            lastBreak = mnBreak;
//            /* Go back to main loop, starting a new segment */
//            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
//            mn = y[i] + lambda[i-1] - (i==n-1 ? 0 : lambda[i]);  // When computing the slopes we need to account for the differences in lambdas
//            mx = y[i] + lambda[i-1] + (i==n-1 ? 0 : lambda[i]); // Here too
//            mxHeight = mnHeight = -lambda[i-1];
//            mnBreak = mxBreak = i;
//            continue;
//        }
//            
//        /* Update height of minorant slope w.r.t. tube center */
//        /* This takes into account both the slope of the minorant and the change in the tube center */
//        mxHeight += mx - y[i];
//        
//        /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
//        /* Minorant is 0 at this point */
//        if ( 0 > mxHeight ) {
//            #ifdef DEBUG
//                fprintf(DEBUG_FILE,"ENDING FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//            #endif
//            /* If violated, break segment at last majorant breaking point */
//            i = mxBreak + 1;
//            /* Build valid segment up to this point using the majorant slope */
//            for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
//                x[j] = mx;
//            /* Start new segment after the break*/
//            lastBreak = mxBreak;
//            /* Go back to main loop, starting a new segment */
//            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
//            mx = y[i] - lambda[i-1] + (i==n-1 ? 0 : lambda[i]); // When computing the slopes we need to account for the differences in lambdas
//            mn = y[i] - lambda[i-1] - (i==n-1 ? 0 : lambda[i]); // Here too
//            mxHeight = mnHeight = lambda[i-1];
//            mnBreak = mxBreak = i;
//            continue;
//        }
//        
//        /* No violations at this point */
//        
//        /* Check if proyected minorant height is under actual minorant */
//        if ( mnHeight <= 0 ) {
//            /* Update minorant slope */
//            mn += ( - mnHeight ) / ( i - lastBreak );
//        }
//
//        #ifdef DEBUG
//            fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
//        #endif
//        
//        /* At this point: we are finished validating last segment! */
//        i++;
//    }
//    
//    /* Build last valid segment */
//    for ( i = lastBreak+1 ; i < n ; i++ )
//        x[i] = mn;
//}
