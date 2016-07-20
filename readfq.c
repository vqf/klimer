#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../myseqfunctions/kmer.h"
#include "../myseqfunctions/readFastq.h"

#define KMERSIZE 11


int main(int argc,char** argv){
  if (argc < 3){
    fprintf(stderr, "Use: readfq file1 file2\n");
    return 1;
  }
  char* fName1 = argv[1];
  char* fName2 = argv[2];
  //fprintf(stderr, "%s\n", fName);
  fastqReader *fq1 = initFastqReader(fName1);
  fastqReader *fq2 = initFastqReader(fName2);
  if (!fq1 || !fq2){
    fprintf(stderr, "Could not create fqReaders from %s and %s\n", fName1, fName2);
    return 1;
  }
  kmerHolder *km = initKmer(KMERSIZE);
  while (getNextFqRead(fq1) && getNextFqRead(fq2)){
    char s1 = getNextFqBase(fq1);
    char s2 = getNextFqBase(fq2);
    while (s1){
      updateKmer(&km, &s1, addRelationship);
      s1 = getNextFqBase(fq1);
    }
    resetTrace(&km);
  }
  summarize(km->ms);
  printf("Sizeof kmerConnector: %lui\n", sizeof(memstruct));
  return 0;
}

