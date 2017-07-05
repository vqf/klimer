#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/readSeqFiles.h"
#include "myseqfunctions/kmerRead.h"

//#define KMERLENGTH 15

uint32_t tellUserEvery = 10000;
uint32_t canonizeEvery = 0;

void print_usage(){
  printf("buildMap [opts] infile [outfile]");
}

int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: %s fasta_file [outfile] [-v] [-k kmerLength]\n", argv[0]);
    return 1;
  }
  uint8_t kmerLength = 11;
  char* infile  = argv[1];
  char* outfile = (char*) calloc(80, sizeof(char));

  sprintf(outfile, "%s.k11", argv[1]);
  for (int i = 2; i < argc; i++){
    char* rg = argv[i];
    D_(0, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (EQ(rg, "-v")) DEBUG = 1;
      if (EQ(rg, "-vv")) DEBUG = 2;
      if (EQ(rg, "-k")){
        i++;
        kmerLength = (uint8_t) atoi(argv[i]);
      }
      if (EQ(rg, "-c")){
        i++;
        canonizeEvery = (uint32_t) atoi(argv[i]);
      }
    }
    else{
      outfile = rg;
    }
  }
  time_t start = time(NULL);
  seqReader* sf = newSeqReader(infile);
  kmerHolder* kh = initKmer(kmerLength, 4);
  uint32_t counter = 0;
  uint32_t counter2 = 0;
  uint32_t nseq = 0;
  while (getNextRead(&sf)){
    D_(1, "Reading %s\n", chromName(&sf));
    char fst = getNextBase(&sf);
    while (fst){
      updateKmer(&kh, &fst, addRelationship);
      fst = getNextBase(&sf);
    }
    resetTrace(&kh);
    counter++;
    counter2++;
    nseq++;
    if (counter >= tellUserEvery){
      counter = 0;
      time_t now = time(NULL);
      D_(0, "%d\t%ld\n", nseq, now-start);
    }
    if (canonizeEvery > 0 && counter2 >= canonizeEvery){
      counter2 = 0;
      time_t now = time(NULL);
      _canonize(&kh);
      D_(0, "%d\t%ld\n", nseq, now-start);
    }
  }
  writeOut(&kh, outfile);
  E_(2, summarize(kh->ms);)
  destroySeqReader(&sf);
  free(outfile);
  destroyKh(&kh);
  time_t now = time(NULL);
  D_(0, "Done in %ld seconds\n", now-start);
  TRACE;
  return 0;
}
