#ifndef SIM_FRAGMENT_H
#define SIM_FRAGMENT_H
typedef struct{
  int type;
  int fixlength;
  
  
}sim_fragment;

sim_fragment *sim_fragment_alloc(int type,double par1, double par2, char *fname);

int getFragmentLength(sim_fragment *sf);

#endif
