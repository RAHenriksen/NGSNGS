#ifndef MRAND_H
#define MRAND_H

#include <random>
#include <cstdlib>

typedef struct{
  int type;
  #if defined(__linux__) || defined(__unix__)
    struct drand48_data buf0;
  #elif defined(__APPLE__) || defined(__MACH__) 
    std::random_device rd;
    std::default_random_engine eng;
    std::uniform_real_distribution<float> distr;
  #endif
}mrand_t;

//type=0 is drand48 style
//type=1 is osx compatible
mrand_t *mrand_alloc(int type_a, long int seedval);
double mrand_pop(mrand_t *mr);
#endif