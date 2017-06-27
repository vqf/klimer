#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

int main(int argc, char** argv){
  char* k11File = NULL;
  for (int i = 1; i < argc; i++){
    char* rg = argv[i];
    D_(1, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (!strcmp(rg, "-v")) DEBUG = 1;
      if (!strcmp(rg, "-vv")) DEBUG = 2;
    }
    else{
      k11File = rg;
    }
  }
  if (!k11File){
    fprintf(stderr, "Use: %s k11_index_file [-v]\n", argv[0]);
    return 1;
  }
  printf("Reading %s...\n", k11File);
  kmerHolder* kh = readIn(k11File);
  printf("Done\n");
  printf("Please enter one sequence at a time\n");
  char* seq = (char*) calloc(0xFFFF, sizeof(char));
  while (scanf("%s", seq) > 0){
    seqCollection* result = searchSeq(&kh, seq);
    printSeqCollection(&kh, &result);
    printf("\n");
  }
  printf("Cleaning...\n");
  free(seq);
  destroyKh(&kh);
  return 0;
}
