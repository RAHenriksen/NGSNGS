#include <random>
#include <cstdlib>

typedef struct{
  int type;
  struct drand48_data buf0;
  std::random_device rd;
  std::default_random_engine eng;
  std::uniform_real_distribution<float> distr;
}mrand_t;

//type=0 is drand48 style
//type=1 is osx compatible
mrand_t *mrand_alloc(int type_a, long int seedval);
double mrand_pop(mrand_t *mr);
