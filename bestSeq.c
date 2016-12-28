#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include "myseqfunctions/readFastq.h"
#include <inttypes.h>
#include <unistd.h>


int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: bestSeq fastq_file\n");
    return 1;
  }
  char* fName1 = argv[1];
  kmerHolder *km = initKmer(11, 4);
  fastqReader* fq = initFastqReader(fName1);
  if (!fq){
    fprintf(stderr, "Could not create fqReaders from %s\n", fName1);
    return 1;
  }
  while (getNextFqRead(&fq)){
    char s1 = getNextFqBase(&fq);
    while (s1){
      updateKmer(&km, &s1, addRelationship);
      s1 = getNextFqBase(&fq);
    }
    resetTrace(&km);
  }
  seqCollection* sc = allTraces(&km);
  seqCollection* best = bestSeqCollection(&sc);
  printSeq(&km, &best);
  destroyKh(&km);
  destroyFastqReader(&fq);
}
