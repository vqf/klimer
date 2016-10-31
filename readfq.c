
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
  if (argc < 4){
    fprintf(stderr, "Use: readfq fastq_file infile outfile\n");
    return 1;
  }
  char* fName1 = argv[1];
  //char* fName2 = argv[2];
  char* infile = argv[2];
  char* outfile = argv[3];
  //fprintf(stderr, "%s\n", fName1); exit(0);
  fq1 = initFastqReader(fName1);
  //fastqReader *fq2 = initFastqReader(fName2);
  if (!fq1){
    fprintf(stderr, "Could not create fqReaders from %s\n", fName1);
    return 1;
  }
  kmerHolder *km = NULL;
  if (infile && (access(infile, R_OK) != -1)){
    km = readIn(infile);
  }
  else{
    km = initKmer(KMERSIZE, 4);
  }
  while (getNextFqRead(&fq1)){
    char s1 = getNextFqBase(&fq1);
    //char s2 = getNextFqBase(&fq2);
    while (s1){
      updateKmer(&km, &s1, addRelationship);
      s1 = getNextFqBase(&fq1);
    }
    resetTrace(&km);
    /*while (s2){x
      updateKmer(&km, &s2, addRelationship);
      s2 = getNextFqBase(&fq2);
    }
    resetTrace(&km);*/
  }
  writeOut(&km, outfile);
  summarize(km->ms);
  destroyKh(&km);
  destroyFastqReader(&fq1);
  //destroyFastqReader(&fq2);
  return 0;
}

