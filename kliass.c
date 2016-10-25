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
  kmerHolder* kh = readIn(argv[1]);
  kcLL* ft = nextTrace(&kh);
  printf("%s\n", getTraceSeq(&kh, &ft));
  ft = nextTrace(&kh);
  printf("%s\n", getTraceSeq(&kh, &ft));
}
