#include <stdlib.h>

typedef struct{
  int type;
  struct drand48_data buf0;
}mrand_t;

//type=0 is drand48 style
mrand_t *mrand_alloc(int type_a, long int seedval);

