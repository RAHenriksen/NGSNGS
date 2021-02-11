#include <cstdio>
#include <htslib/kstring.h>

void minskriver(FILE *fp, char *str){
  fprintf(fp,"%s",str);
}


int main(){
  printf("hej verden\n");
  fprintf(stdout,"hej verden2\n");
  fprintf(stderr,"Detter er inadsfasdf debgug beskegagsd\n");
  minskriver(stdout,"besked1\n");
  minskriver(stderr,"besked2\n");

  FILE *fp = fopen("minfil.txt","wb");
  minskriver(fp,"besked3\n");
  minskriver(fp,"besked1\n");
  fclose(fp);
  return 0;
}
