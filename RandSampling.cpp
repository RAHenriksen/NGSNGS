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

//! Allocate workspace for random-number sampling.
ransampl_ws* ransampl_alloc( int n ){
    /*
    ransampl_alloc - Allocates memory for ransampl_ws` structure which is used for efficient sampling from discrete probability distributions using the alias method

    @param n The number of elements in the distribution.
    @return A pointer to the allocated `ransampl_ws` structure. On failure, the function exits the program with an ENOMEM error code.
    */

    ransampl_ws *ws;
    // Attempt to allocate memory for the structure and its components.
    if ( !(ws = (ransampl_ws*) malloc( sizeof(ransampl_ws) )) ||
         !(ws->alias = (int*) malloc( n*sizeof(int) )) ||
         !(ws->prob = (double*) malloc( n*sizeof(double) )) ) {
            
        fprintf( stderr, "ransampl: workspace allocation failed\n" );
        exit(ENOMEM);
    }

    // Set the number of elements in the workspace.
    ws->n = n;
    return ws;
}

void ransampl_set( ransampl_ws *ws, const double* p ){
    /*
    ransampl_set - Initializes the random sampling by precomputing alias tables and probability arrays.

    @param ws A pointer to the `ransampl_ws` structure, previously allocated with `ransampl_alloc`.
    @param p  An array of probabilities for each outcome. The length of the array must match the size (`n`) specified in the allocated memory.
    
    */

    //! Initialize workspace by precompute alias tables from given probabilities.
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

int ransampl_draw2( ransampl_ws *ws,double r1, double r2){
    /*
    ransampl_draw2 - Draws a sample from the probability distribution using two uniform random numbers.

    @param ws A pointer to the `ransampl_ws` structure, previously allocated with `ransampl_alloc`.
    @param r1 A uniform random number in the range [0, 1) used to select the initial index.
    @param r2 A uniform random number in the range [0, 1) used to decide between the index and its alias.
    
    @return An integer representing the sampled index from the probability distribution.
    */

    int i = (int) (ws->n * r1);
    return r2 < ws->prob[i] ? i : ws->alias[i];
}

//! Free the random-number sampling workspace.
void ransampl_free( ransampl_ws *ws ){
    /*
    ransampl_free - frees the allocated memory for the random-number sampling structure

    @param ws A pointer to the `ransampl_ws` structure, previously allocated with `ransampl_alloc`.
    */

    free( ws->alias );
    free( ws->prob );
    free( ws );
}
