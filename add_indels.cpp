#include "add_indels.h"
char *bas = "ACGTN";
char *bas2 = "QWERY";

int add_indel(mrand_t *mr,char *frag,int readlength,double *pars){
  int beg = 0;
  int end = strlen(frag);//remember this is the length, not the index of the last element.
  int ops[2] ={0,0};
  while(beg<end-1){
    //fprintf(stderr,"beg[%d]: %s len:%lu\n",beg,frag,strlen(frag));
    //do dont anything if we are in the internal part of the fragmenth that is not getting sequenced
    if(end>2*readlength && beg > readlength && ((end-beg)>readlength)){
      beg++;
      continue;
    }
    if(mrand_pop(mr)<pars[0]){
      ops[0]++;
      int len = Random_geometric_k(pars[2],mr);
      //fprintf(stderr,"Having insertion[%d]: %d end-1: %d\n",beg,len,end-1);
      for(int i=end-1;i>=beg;i--){
        //fprintf(stderr,"%d -> %d\n",i,i+len);
        frag[i+len] = frag[i];
      }
      for(int i=0;i<len;i++){
        frag[beg+i] = bas[mrand_pop_long(mr) %5]; //bas[mrand_pop_long(mr) %5];
        //fprintf(stderr,"Setting insertion: %d\n",beg+i);
      }
      beg += len;
      end += len;
      continue;
    }
    if(mrand_pop(mr)<pars[1]){
      ops[1]++;
      int len = Random_geometric_k(pars[3],mr);
      //fprintf(stderr,"Having deletion[%d]: %d\n",beg,len);
      if(len>1){
        //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
        //2) Deletion will cover the rest of the read. 
        if(len+beg>end){
          //case 2)
          frag[beg+1] = '\0';
          end = beg+1;
          //should make debug statement to validate correctness of the thinking machines.
        }else{
          //case 1) we dont need to check for running of the data
          for(int i=0;i<len;i++)
            frag[beg+i] = frag[beg+i+len];
          end -= len;
          frag[end] = '\0';
        }
      }
      beg++;
      continue;
    }
    //fprintf(stderr,"Nothing happens[%d]\n",beg);
    beg++;
  }
  return ops[0]+ops[1];
}


#ifdef __WITH_MAIN__
int main(int argc,char **argv){

  int readlength = 100;
  mrand_t *mr = mrand_alloc(2,777);
  char fragment[1024];
  double pars[4] = {0.05,0.1,0.1,0.2};
  while(1){
    memset(fragment,'\0',1024);
    int flen = mrand_pop_long(mr) % 32;
    for(int i=0;i<flen;i++)
      fragment[i] = bas[mrand_pop_long(mr) %5];
    fprintf(stderr,"1) %s len: %lu\n",fragment,strlen(fragment));
    add_indel(mr,fragment,100,pars);
    fprintf(stderr,"2) %s len: %lu\n",fragment,strlen(fragment));
    break;
  }
  return 0;
}
#endif
