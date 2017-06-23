#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <inttypes.h>
#include <unistd.h>


int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: kliass k11file\n");
    return 1;
  }
  for (int i = 2; i < argc; i++){
    char* rg = argv[i];
    D_(0, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (!strcmp(rg, "-v")) DEBUG = 1;
      if (!strcmp(rg, "-vv")) DEBUG = 2;
    }
  }
  kmerHolder* kh = readIn(argv[1]);

  //seqCollection* sc = allTraces(&kh);
  kcLL* s1 = nextTrace(&kh);
  while (s1){
    char* yo = getTraceSeq(&kh, &s1);
    printf("%s\n", yo);
    free(yo);
    resetKcLL(&s1);
    s1 = nextTrace(&kh);
  }
  resetKcLL(&s1);
  destroyKh(&kh);
  //summarize(kh->ms); exit(0);
  /*kcLL* ft = nextTrace(&kh);
  uint32_t i = 0;
  while (ft){
    i++;
    printf(">trace%lu\n", (LUI) i);
    char* seq = getTraceSeq(&kh, &ft);
    printf("%s\n", seq);
    if (ft) resetKcLL(&ft);
    ft = nextTrace(&kh);
    free(seq);
  }*/
}
