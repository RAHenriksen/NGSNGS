#include "add_indels.h"
#include "fasta_sampler.h"
const char *bas = "ACGTN";
const char *bas2 = "QWERY";

void add_indel(mrand_t *mr,char *frag,int readlength,double *pars,char *INDEL_INFO,int* ops){
  int beg = 0;
  int end = strlen(frag);
  //int ops[2] ={0,0};
  int fragbefore = strlen(frag);
  
  char DelOps[512] = "\0";
  char InsOps[512] = "\0";
  int offsetDel = 0;//sprintf(DelOps, "");
  int offsetIns = 0;//sprintf(InsOps, "");

  int i=0;
  int cumlen = 0;
  int offset = 0;
  while(beg<end-1){
    if(end>2*readlength && beg > readlength && ((end-beg)>readlength)){
      beg++;
      continue;
    }
    if(mrand_pop(mr)<pars[0]){
      ops[0]++;
      int len = Random_geometric_k(pars[2],mr)+1;
      offsetIns+=sprintf(InsOps+offsetIns,"[%d,%d:%d],",beg-cumlen,beg,len);
      cumlen = cumlen+len;
      for(int i=end-1;i>=beg;i--){
        frag[i+len] = frag[i];
      }
      for(int i=0;i<len;i++){
        frag[beg+i] =  bas[mrand_pop_long(mr) %5];
      }
      beg += len;
      end += len;
      continue;
    }
    if(mrand_pop(mr)<pars[1]){
      ops[1]++;
      int len = Random_geometric_k(pars[3],mr)+1;
      offsetDel+=sprintf(DelOps+offsetDel,"[%d,%d:%d],",beg+offset,beg-cumlen+offset,len);
      cumlen = cumlen+len;
        //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
        //2) Deletion will cover the rest of the read. 
        if(len+beg>end){
          //case 2)
          frag[beg+1] = '\0';
          end = beg+1;
        }else{
          //case 1) we dont need to check for running of the data
          memmove(frag + beg, frag+(beg+len), strlen(frag) - len);
          end -= len;
          frag[end] = '\0';
        }
      i++;      
      beg++;
      offset++;
      continue;
    }
    beg++;
  }
  if (ops[0] == 0){
    snprintf(INDEL_INFO,1024,"%d\tNA\t%d\t%s\t%d\t%ld",ops[0],ops[1],DelOps,fragbefore,strlen(frag));
  }
  else if (ops[1] == 0){
    snprintf(INDEL_INFO,1024,"%d\t%s\t%d\tNA\t%d\t%ld",ops[0],InsOps,ops[1],fragbefore,strlen(frag));
  }
  else
  {
    snprintf(INDEL_INFO,1024,"%d\t%s\t%d\t%s\t%d\t%ld",ops[0],InsOps,ops[1],DelOps,fragbefore,strlen(frag));
  }
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
  char INDEL_INFO[512];
  int readlength = 100;
  mrand_t *mr = mrand_alloc(2,777);
  char fragment[1024];
  double pars[4] = {0.05,0.0,0.9,0.0};
  while(1){
    memset(fragment,'\0',1024);
    int flen = 100; //mrand_pop_long(mr) % 32;
    for(int i=0;i<flen;i++)
      fragment[i] = bas[mrand_pop_long(mr) %5];
    fprintf(stderr,"1) %s len: %lu\n",fragment,strlen(fragment));
    add_indel(mr,fragment,100,pars,INDEL_INFO);
    fprintf(stderr,"2) %s len: %lu\n",fragment,strlen(fragment));
    break;
  }
  return 0;
}
#endif
