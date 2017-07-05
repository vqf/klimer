#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/readSeqFiles.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"

#define KMERLENGTH 11


int main(int argc,char** argv){
  char* k11File = NULL;
  uint8_t userLength = 0;
  for (int i = 1; i < argc; i++){
    char* rg = argv[i];
    D_(1, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (EQ(rg, "-v")) DEBUG = 1;
      if (EQ(rg, "-vv")) DEBUG = 2;
      if (EQ(rg, "-k")){
        i++;
        userLength = (uint8_t) atoi(argv[i]);
      }
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
  uint8_t klength = kh->kmerSize + 1;
  if (userLength > klength) klength = userLength;
  if (klength <= kh->kmerSize){
    printf("This version only generates histograms for k-mer lengths larger than %"
           PRIu8"\n", kh->kmerSize);
    exit(0);
  }
  printf("Generating histogram with k-mer length %" PRIu8 "\n", klength);
  printHistogram(&kh, klength);
  destroyKh(&kh);
}
