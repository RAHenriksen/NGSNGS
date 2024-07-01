#include "add_indels.h"
#include "fasta_sampler.h"
#include <htslib/kstring.h>

const char *bas = "ACGTN";
const char *bas2 = "QWERY";

void add_indel(mrand_t *mr,char *frag,int readlength,double *pars,char *INDEL_INFO,int* ops){
  int beg = 0;
  int end = strlen(frag);
  //int ops[2] ={0,0};
  int fragbefore = strlen(frag);
  
  char DelOps[512] = "\0";
  char InsOps[512] = "\0";
  int offsetDel = 0;
  int offsetIns = 0;

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
      offsetIns+=snprintf(InsOps+offsetIns,512,"[%d,%d:%d],",beg-cumlen,beg,len);
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
      offsetDel+=snprintf(DelOps+offsetDel,512,"[%d,%d:%d],",beg+offset,beg-cumlen+offset,len);
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

void add_indel_amplicon_fa(mrand_t *mr,kstring_t* seq,double *pars,int* ops){
  int beg = 0;
  int end = seq->l;
  //int ops[2] ={0,0};
  int fragbefore = seq->l;

  int i=0;

  while(beg<end-1){
    if(end>2*fragbefore && beg > fragbefore && ((end-beg)>fragbefore)){
      beg++;
      continue;
    }

    //insertions
    if(mrand_pop(mr)<pars[0]){
      ops[0]++;
      int len = Random_geometric_k(pars[2],mr)+1;

      // Ensure there is enough capacity for the insertion
      if (seq->l + len > seq->m) {
        size_t new_capacity = seq->l + len;
        char *new_seq = (char *)realloc(seq->s, new_capacity + 1);
        if (!new_seq) {
          fprintf(stderr, "Memory reallocation failed\n");
          free(seq->s);
          return;
        }
        seq->s = new_seq;
        seq->m = new_capacity;
      }
      
      // Shift the existing characters to the right
      memmove(seq->s + beg + len, seq->s + beg, seq->l - beg + 1);

      // Insert new characters
      for (int i = 0; i < len; i++) {
        seq->s[beg + i] = bas[mrand_pop_long(mr) %5];
      }

      //memset(seq->s + beg, 'Q', len);
      //memset(qual->s + beg, 'B', len);

      // Update the lengths
      seq->l += len;
      beg += len;
      end += len;
      continue;
    }

    // Handle deletion
    if(mrand_pop(mr)<pars[1]){
      ops[1]++;
      int len = Random_geometric_k(pars[3],mr)+1;
      //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
      //2) Deletion will cover the rest of the read. 
      if(len+beg>=end){
        //case 2)
        seq->s[beg] = '\0';
        seq->l = beg;
        end = beg+1;
      }
      else{
        //case 1) we dont need to check for running of the data
        memmove(seq->s + beg, seq->s + beg + len, seq->l - beg - len + 1);
        end -= len;
        seq->l -= len;
      }
      continue;
    }
    beg++;
  }

  // Null-terminate the sequences
  seq->s[seq->l] = '\0';
}

void add_indel_amplicon_fqbam(mrand_t *mr,kstring_t* seq,kstring_t* qual,double *pars,int* ops,int ErrProbTypeOffset){
  int beg = 0;
  int end = seq->l;
  //int ops[2] ={0,0};
  int fragbefore = seq->l;

  int i=0;
  int qscoresval;

  while(beg<end-1){
    if(end>2*fragbefore && beg > fragbefore && ((end-beg)>fragbefore)){
      beg++;
      continue;
    }

    //insertions
    if(mrand_pop(mr)<pars[0]){
      ops[0]++;
      int len = Random_geometric_k(pars[2],mr)+1;

      // Ensure there is enough capacity for the insertion
      if (seq->l + len > seq->m) {
        size_t new_capacity = seq->l + len;
        char *new_seq = (char *)realloc(seq->s, new_capacity + 1);
        char *new_qual = (char *)realloc(qual->s, new_capacity + 1);
        if (!new_seq || !new_qual) {
          fprintf(stderr, "Memory reallocation failed\n");
          free(seq->s);
          free(qual->s);
          return;
        }
        seq->s = new_seq;
        qual->s = new_qual;
        seq->m = qual->m = new_capacity;
      }
      
      // Shift the existing characters to the right
      memmove(seq->s + beg + len, seq->s + beg, seq->l - beg + 1);
      memmove(qual->s + beg + len, qual->s + beg, qual->l - beg + 1);

      // Insert new characters
      for (int i = 0; i < len; i++) {
        seq->s[beg + i] = bas[mrand_pop_long(mr) %5];
        int quality_val = (int)(mrand_pop(mr) * 39 + 1);
        qual->s[beg + i] = (char)(quality_val + ErrProbTypeOffset);
      }

      //memset(seq->s + beg, 'Q', len);
      //memset(qual->s + beg, 'B', len);

      // Update the lengths
      seq->l += len;
      qual->l += len;
      beg += len;
      end += len;
      continue;
    }

    // Handle deletion
    if(mrand_pop(mr)<pars[1]){
      ops[1]++;
      int len = Random_geometric_k(pars[3],mr)+1;
      //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
      //2) Deletion will cover the rest of the read. 
      if(len+beg>=end){
        //case 2)
        seq->s[beg] = '\0';
        qual->s[beg] = '\0';
        seq->l = beg;
        qual->l = beg;
        end = beg+1;
      }
      else{
        //case 1) we dont need to check for running of the data
        memmove(seq->s + beg, seq->s + beg + len, seq->l - beg - len + 1);
        memmove(qual->s + beg, qual->s + beg + len, qual->l - beg - len + 1);
        end -= len;
        seq->l -= len;
        qual->l -= len;
      }
      continue;
    }
    beg++;
  }

  // Null-terminate the sequences
  seq->s[seq->l] = '\0';
  qual->s[qual->l] = '\0';
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
  char INDEL_INFO[1024];
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


      /*
      if (seq_intermediate.l + len > qual_intermediate.l) {
          // Reallocate memory to increase the capacity
          size_t new_capacity = seq_intermediate.l + len;
          char *new_seq = (char *)realloc(seq_intermediate.s, new_capacity + 1);  // +1 for the null terminator
          char *new_qual = (char *)realloc(qual_intermediate.s, new_capacity + 1);  // +1 for the null terminator
          seq_intermediate.s = new_seq;
          seq_intermediate.l = new_capacity;
          seq_intermediate.m = new_capacity;
          seq->l = new_capacity;
          seq->m = new_capacity;

          qual_intermediate.s = new_qual;
          qual_intermediate.l = new_capacity;
          qual_intermediate.m = new_capacity;
          qual->l = new_capacity;
          qual->m = new_capacity;

      }*/