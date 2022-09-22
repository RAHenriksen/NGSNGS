#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::reverse
#include <iostream>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>
#include <random>
#include <map>
#include <math.h>
#include <pthread.h>
#include "RandSampling.h"
#include "mrand.h"

#define LENS 4096
#define MAXBINS 100

<<<<<<< HEAD
int BinarySearch_fraglength(double* SearchArray,int low, int high, double key){
    int ans = 0; 
    while (low <= high) {
        int mid = low + (high - low + 1) / 2;
        double midVal = SearchArray[mid];
 
        if (midVal < key) {
            ans = mid;
            low = mid + 1;
        }
        else if (midVal > key) {

            high = mid - 1;
        }
        else if (midVal == key) {
 
            high = mid - 1;
        }
    }
 
    return ans+1;
}

=======
>>>>>>> development
//! Allocate workspace for random-number sampling.
ransampl_ws* ransampl_alloc( int n )
{
    ransampl_ws *ws;
    if ( !(ws = (ransampl_ws*) malloc( sizeof(ransampl_ws) )) ||
         !(ws->alias = (int*) malloc( n*sizeof(int) )) ||
         !(ws->prob = (double*) malloc( n*sizeof(double) )) ) {
        fprintf( stderr, "ransampl: workspace allocation failed\n" );
        exit(ENOMEM); // ENOMEM
    }
    ws->n = n;
    return ws;
}

//! Initialize workspace by precompute alias tables from given probabilities.
void ransampl_set( ransampl_ws *ws, const double* p )
{
    const int n = ws->n;
    int i, a, g;

    // Local workspace:
    double *P;
    int *S, *L;
    if ( !(P = (double*) malloc( n*sizeof(double) ) ) ||
         !(S = (int*) malloc( n*sizeof(int) ) ) ||
         !(L = (int*) malloc( n*sizeof(int) ) ) ) {
        fprintf( stderr, "ransampl: temporary allocation failed\n" );
        exit(ENOMEM);
    }

    // Normalise given probabilities:
    double sum=0;
    for ( i=0; i<n; ++i ) {
        if( p[i]<0 ) {
            fprintf( stderr, "ransampl: invalid probability p[%i]<0\n", i );
            exit(EINVAL);
        }
        sum += p[i];
    }
    if ( !sum ) {
        fprintf( stderr, "ransampl: no nonzero probability\n" );
        exit(EINVAL);
    }
    for ( i=0; i<n; ++i )
        P[i] = p[i] * n / sum;

    // Set separate index lists for small and large probabilities:
    int nS = 0, nL = 0;
    for ( i=n-1; i>=0; --i ) {
        // at variance from Schwarz, we revert the index order
        if ( P[i]<1 )
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    // Work through index lists
    while ( nS && nL ) {
        a = S[--nS]; // Schwarz's l
        g = L[--nL]; // Schwarz's g
        ws->prob[a] = P[a];
        ws->alias[a] = g;
        P[g] = P[g] + P[a] - 1;
        if ( P[g] < 1 )
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while ( nL )
        ws->prob[ L[--nL] ] = 1;

    while ( nS )
        // can only happen through numeric instability
        ws->prob[ S[--nS] ] = 1;

    // Cleanup:
    free( P );
    free( S );
    free( L );
}

int ransampl_draw2( ransampl_ws *ws,double r1, double r2)
{   
    const int i = (int) (ws->n * r1);
    return r2 < ws->prob[i] ? i : ws->alias[i];
}

//! Free the random-number sampling workspace.
void ransampl_free( ransampl_ws *ws )
{
    free( ws->alias );
    free( ws->prob );
    free( ws );
}
