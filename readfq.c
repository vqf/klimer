#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/readFastq.h"
#include <inttypes.h>
#include <unistd.h>

#define KMERSIZE 11

fastqReader* fq1;

/*void signal_handler(int signal){
  if (signal == SIGABRT){
    uint8_t zoff = doTell(fq1);
    uint64_t sz = fSize(fq1);
    printf("Read about %d%% of %" PRIu64 "\n", zoff, sz);
  }
}*/

int main(int argc,char** argv){
  if (argc < 5){
    fprintf(stderr, "Use: readfq file1 file2 infile outfile\n");
    return 1;
  }
  char* fName1 = argv[1];
  char* fName2 = argv[2];
  char* infile = argv[3];
  char* outfile = argv[4];
  //fprintf(stderr, "%s\n", fName1); exit(0);
  fq1 = initFastqReader(fName1);
  fastqReader *fq2 = initFastqReader(fName2);
  if (!fq1 || !fq2){
    fprintf(stderr, "Could not create fqReaders from %s and %s\n", fName1, fName2);
    return 1;
  }
  kmerHolder *km = NULL;
  if (infile && (access(infile, R_OK) != -1)){
    km = readIn(infile);
  }
  else{
    km = initKmer(KMERSIZE, 4);
  }
  while (getNextFqRead(&fq1) && getNextFqRead(&fq2)){
    char s1 = getNextFqBase(&fq1);
    char s2 = getNextFqBase(&fq2);
    while (s1){
      updateKmer(&km, &s1, addRelationship);
      s1 = getNextFqBase(&fq1);
    }
    resetTrace(&km);
    while (s2){
      updateKmer(&km, &s2, addRelationship);
      s2 = getNextFqBase(&fq2);
    }
    resetTrace(&km);
  }
  writeOut(&km, outfile);
  destroyKh(&km);
  destroyFastqReader(&fq1);
  destroyFastqReader(&fq2);
  return 0;
}

