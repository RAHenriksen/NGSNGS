#include "add_indels.h"
#include "fasta_sampler.h"
const char *bas = "ACGTN";
const char *bas2 = "QWERY";

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
      int len = Random_geometric_k(pars[2],mr)+1;
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
      int len = Random_geometric_k(pars[3],mr)+1;
      //fprintf(stderr,"Having deletion[%d]: %d\n",beg,len);
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
      beg++;
      continue;
    }
    //fprintf(stderr,"Nothing happens[%d]\n",beg);
    beg++;
  }
  return ops[0]+ops[1];
}


/*while(beg<end-1){
    //fprintf(stderr,"beg[%d]: %s len:%lu\n",beg,frag,strlen(frag));
    //do dont anything if we are in the internal part of the fragmenth that is not getting sequenced
    if(end>2*readlength && beg > readlength && ((end-beg)>readlength)){
      beg++;
      continue;
    }
    if(mrand_pop(mr)<pars[0]){
      ops[0]++;
      int len = Random_geometric_k(pars[2],mr)+1;
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
*/

void add_indel_fs(fasta_sampler *fs,mrand_t *mr,double *pars){
  for(int j=0;j<fs->nref;j++){
    //fprintf(stderr,"ref %s \t chr number %d \n",fs->seqs_names[j],j);
    int beg = 0;
    int end = fs->seqs_l[j];
    //fs->seqs[j]
    //fprintf(stderr,"BEG %d \t END %d \n",beg,end);
    int len_value = -1;
    while(beg<end-1){
      //fprintf(stderr,"beg[%d]: %s len:%lu\n",beg,frag,strlen(frag));
      //do dont anything if we are in the internal part of the fragmenth that is not getting sequenced
      if(end>2*150 && beg > 150 && ((end-beg)>150)){
        beg++;
        continue;
      }
      if(mrand_pop(mr)<pars[0]){
        int len = Random_geometric_k(pars[2],mr)+1;
        //fprintf(stderr,"Having insertion[%d]: %d end-1: %d\n",beg,len,end-1);
        for(int i=end-1;i>=beg;i--){
          //fprintf(stderr,"%d -> %d\n",i,i+len);
          fs->seqs[j][i+len] = fs->seqs[j][i];
        }
        for(int i=0;i<len;i++){
          fs->seqs[j][beg+i] = bas[mrand_pop_long(mr) %5]; //bas[mrand_pop_long(mr) %5];
          //fprintf(stderr,"Setting insertion: %d\n",beg+i);
        }
        beg += len;
        end += len;
        continue;
      }
      if(mrand_pop(mr)<pars[1]){
        int len = Random_geometric_k(pars[3],mr)+1;
        //fprintf(stderr,"Having deletion[%d]: %d\n",beg,len);
          //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
          //2) Deletion will cover the rest of the read. 
          if(len+beg>end){
            //case 2)
            fs->seqs[j][beg+1] = '\0';
            end = beg+1;
            //should make debug statement to validate correctness of the thinking machines.
          }else{
            //case 1) we dont need to check for running of the data
            for(int i=0;i<len;i++)
              fs->seqs[j][beg+i] = fs->seqs[j][beg+i+len];
            end -= len;
            fs->seqs[j][end] = '\0';
          }      
        beg++;
        continue;
      }
      //fprintf(stderr,"Nothing happens[%d]\n",beg);
      beg++;
    }
  }
}

/*void add_indel_fs(fasta_sampler *fs,mrand_t *mr,double *pars){
  fprintf(stderr,"ISNIDE ADD_INDELS \n");
  fprintf(stderr,"indels Param 1 %f 2 %f 3 %f 4 %f\n",pars[0],pars[1],pars[2],pars[3]);
  fprintf(stderr,"test %f\n",mrand_pop(mr));
  fprintf(stderr,"test2 %d\n",abs(mrand_pop_long(mr)%5));
  for(int j=0;j<fs->nref;j++){
    int beg = 0;
    int end = fs->seqs_l[j];//remember this is the length, not the index of the last element.
    while(beg<end-1){
      //fprintf(stderr,"beg[%d]: %lu len:%lu\n",beg,fs->seqs_l[j],strlen(fs->seqs[j]));
      double dtemp = mrand_pop(mr);
      fprintf(stderr,"dtemp %f\n",dtemp);
      if(dtemp<pars[0]){
        int len = Random_geometric_k(pars[2],mr)+1;
        //fprintf(stderr,"Having insertion[%d]: %d end-1: %d\n",beg,len,end-1);
        for(int i=end-1;i>=beg;i--){
          //fprintf(stderr,"%d -> %d\n",i,i+len);
          fs->seqs[j][i+len] = fs->seqs[j][i];
        }
        fprintf(stderr,"BEFORE BAS \n");
        for(int i=0;i<len;i++){
          fprintf(stderr,"test %d\n",i);
          fs->seqs[j][beg] = bas[3];//bas[mrand_pop_long(mr) %5]; //bas[mrand_pop_long(mr) %5];
          fprintf(stderr,"Setting insertion: %d\n",beg+i);
        }
        fprintf(stderr,"AFTER BAS \n");
        beg += len;
        end += len;
        continue;
      }
      if(mrand_pop(mr)<pars[1]){
        int len = Random_geometric_k(pars[3],mr)+1;
        fprintf(stderr,"Having deletion[%d]: %d\n",beg,len);
          //two cases: 1) is that deletion is in the middle of the fragment and we should shift all data to the left
          //2) Deletion will cover the rest of the read. 
          if(len+beg>end){
            //case 2)
            fs->seqs[j][beg+1] = '\0';
            end = beg+1;
            //should make debug statement to validate correctness of the thinking machines.
          }
          else{
            //case 1) we dont need to check for running of the data
            for(int i=0;i<len;i++)
              fs->seqs[j][beg+i] = fs->seqs[j][beg+i+len];
            end -= len;
            fs->seqs[j][end] = '\0';
          }      
        beg++;
        continue;
      }
      //fprintf(stderr,"Nothing happens[%d]\n",beg);
      beg++;
    }
  }
}*/

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
