#include <stdlib.h>
#include <time.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/readSeqFiles.h"
#include "myseqfunctions/kmerRead.h"

#define KMERLENGTH 11
#define TELLUSEREVERY 10000

int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: %s fasta_file [outfile]\n", argv[0]);
    return 1;
  }
  char* infile  = argv[1];
  char* outfile = (char*) calloc(80, sizeof(char));
  sprintf(outfile, "%s.k11", argv[1]);
  if (argc >= 3){
    outfile = argv[2];
  }
  time_t start = time(NULL);
  seqReader* sf = newSeqReader(infile);
  kmerHolder* kh = initKmer(KMERLENGTH, 4);
  uint32_t counter = 0;
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
    nseq++;
    if (counter >= TELLUSEREVERY){
      counter = 0;
      time_t now = time(NULL);
      //_canonize(&kh);
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
