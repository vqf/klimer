#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

int main(int argc, char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: %s k11_index_file sequence [-v]\n", argv[0]);
    return 1;
  }
  char* k11File = NULL;
  char* seq     = NULL;
  for (int i = 1; i < argc; i++){
    char* rg = argv[i];
    D_(0, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (!strcmp(rg, "-v")) DEBUG = 1;
      if (!strcmp(rg, "-vv")) DEBUG = 2;
    }
    else{
      if (k11File){
        seq = rg;
      }
      else{
        k11File = rg;
      }
    }
  }

  kmerHolder* kh = readIn(k11File);
  D_(0, "Read %s\n", k11File);
  seqCollection* result = searchSeq(&kh, seq);
  printSeqCollection(&kh, &result);
  destroyKh(&kh);
  return 0;
}
