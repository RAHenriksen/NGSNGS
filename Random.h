#ifndef RANDOM_H
#define RANDOM_H

int __drand48_iterate (unsigned short int xsubi[3], struct drand48_data *buffer);

int __erand48_r (unsigned short int xsubi[3], struct drand48_data *buffer,double *result);

int __srand48_r (long int seedval, struct drand48_data *buffer);

int drand48_r (struct drand48_data *buffer, double *result);

int rand_r (unsigned int *seed);
#endif