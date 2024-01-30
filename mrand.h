#ifndef MRAND_H
#define MRAND_H

#include <random>
#include <iostream>
#include <cstdlib>
#include <climits>

/*
  struct Ran {
  Ullong u,v,w;
  Ran(Ullong j) : v(4101842887655102017LL), w(1) {
  u = j ^ v; int64();
  v = u; int64();
  w = v; int64();
  }
  inline Ullong int64() {
  u = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
  w = 4294957665U*(w & 0xffffffff) + (w >> 32);
  Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
  return (x + v) ^ w;
  }
  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
  };
*/
typedef struct{
  int type;
#if defined(__linux__) || defined(__unix__)
  struct drand48_data buf0;
#endif
  std::random_device rd;
  std::default_random_engine eng;
  std::uniform_real_distribution<double> distr;
  std::uniform_int_distribution<> distrInt;
  unsigned int rand_r_seed;
  unsigned short xsubi[3];
  //below are states for nr. After comparing it should be removed.
  unsigned long long nr_uvw[3];
  double nr_inv_rec;
  unsigned long long nr_int64(){
    unsigned long long u,v,w;
    u= nr_uvw[0];v=nr_uvw[1];w=nr_uvw[2];
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    nr_uvw[0] = u ;nr_uvw[1]=v;nr_uvw[2]=w;
    return (x + v) ^ w;
  }
}mrand_t;

//type=0 is drand48 style
//type=1 is osx compatible
//type=2 is rand_r
//type=3 is erand, the most recent member of the family.
//type=4 is ran from nr
mrand_t *mrand_alloc(int type_a, long int seedval);
double mrand_pop(mrand_t *mr);
long mrand_pop_long(mrand_t *mr);
int Random_geometric_k(const double p,mrand_t *mr);
#endif
