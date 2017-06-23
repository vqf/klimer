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
  if (argc < 2){
    fprintf(stderr, "Use: kmerHist fastqFile [fromk11]\n");
    return 1;
  }
  char* infile  = argv[1];
  if (argc == 2){
    char* outfile = (char*) calloc(80, sizeof(char));
    sprintf(outfile, "%s.k11", argv[1]);
    seqReader* sf = newSeqReader(infile);
    kmerHolder* kh = initKmer(KMERLENGTH, 4);
    while (getNextRead(&sf)){
      D_(1, "Reading %s\n", chromName(&sf));
      char fst = getNextBase(&sf);
      while (fst){
        updateKmer(&kh, &fst, addCount);
        fst = getNextBase(&sf);
      }
      resetTrace(&kh);
    }
    writeOut(&kh, outfile);
    printHistogram(&kh);
    free(outfile);
    destroySeqReader(&sf);
    destroyKh(&kh);
  }
  else{
    kmerHolder* kh = readIn(infile);
    memstruct* ms = kh->ms;
    for (uint32_t i = 0; i < ms->nPos; i++){
      kmerConnector* kc = getConnector(&ms, i, 0);
      if (kc){
        char* seq = calloc(kh->kmerSize + 1, sizeof(char));
        pos2seq(&ms, i, seq);
        printf("%s\t%lu\n", seq, (LUI) kc->n);
        free(seq);
      }
    }
    //printHistogram(&kh);
    destroyKh(&kh);
  }
}
