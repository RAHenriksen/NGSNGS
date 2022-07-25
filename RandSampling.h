#ifndef RANDSAMPLING_H
#define RANDSAMPLING_H

#include "NGSNGS_func.h"
#include "mrand.h"

int BinarySearch_fraglengthold(double* SearchArray,int low, int high, double key);

typedef struct{
    int n;
    int* alias;
    double* prob;
} ransampl_ws;

ransampl_ws* ransampl_alloc(int n);

void ransampl_set( ransampl_ws *ws, const double *p );

int ransampl_draw2( ransampl_ws *ws,double r1, double r2); //added below function to make it threadsafe tsk 23dec 2021

void ransampl_free( ransampl_ws *ws );

#endif /* NGSNGSFUNC_H */