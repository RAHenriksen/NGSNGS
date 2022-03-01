#ifndef NGSNGS_RANDOM_H
#define NGSNGS_RANDOM_H

/*
 https://code.woboq.org/userspace/glibc/stdlib/stdlib.h.html#drand48_data
 https://code.woboq.org/userspace/glibc/stdlib/rand_r.c.html#rand_r
 https://code.woboq.org/userspace/glibc/stdlib/drand48_r.c.html#drand48_r
 https://code.woboq.org/userspace/glibc/stdlib/erand48_r.c.html#47
 https://code.woboq.org/userspace/glibc/sysdeps/ieee754/ieee754.h.html#ieee754_double
 */
struct drand48_data
  {
    unsigned short int __x[3];        /* Current state.  */
    unsigned short int __old_x[3]; /* Old state.  */
    unsigned short int __c;        /* Additive const. in congruential formula.  */
    unsigned short int __init;        /* Flag for initializing.  */
    __extension__ unsigned long long int __a;        /* Factor in congruential
                                                   formula.  */
  };

int drand48_r (struct drand48_data *buffer, double *result);

int srand48_r (long int seedval, struct drand48_data *buffer);
#endif
