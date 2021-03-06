#include <assert.h>
#include <stdarg.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


void addRawSeq(kmerHolder** khp, int num_args, char** arrp){
  kmerHolder* km = *khp;
  for(int i = 0; i < num_args; i++) {
    char* dtext = arrp[i];
    uint16_t l = strlen(dtext);
    for (int j = 0; j < l; j++){
      updateKmer(&km, &dtext[j], addRelationship);
    }
    resetTrace(&km);
  }
}

#define nseqs 6

int main(int argc, char** argv){
  char* fileout = "/users/vqf/desktop/tmp.txt";
  for (int i = 1; i < argc; i++){
    char* rg = argv[i];
    //D_(0, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (EQ(rg, "-v")) DEBUG = 1;
      if (EQ(rg, "-vv")) DEBUG = 2;
    }
    else{
      fileout = rg;
    }
  }
  kmerHolder *km = initKmer(11, 4);
  //uint32_t p = seq2pos(km->ms, "AAATGCCGGAT");
  char* seqs[nseqs] = {
    "AGCATTACTTCAAGCGGCGCATCACCAGATGTAAACTGTTTTT",
    "AGTGTCACCAGCAATACTTCAAGCGGCGCATCACCAGAAAGGGCAGGCGCTTCTCCCTGACCTTCTGGACGTTCAGCATTACGAATATATTC",
    "GTTCTTCGGTGGTCATTTTATTGCGTAATCTTTTCATATTATTTGCTATCGCGAATCCCTAAGCGTTGTAACAATCCGATAATCCCTTCCCGCAGAAGCAAGGATGTAAACTGTTTTTGTTCAACTGGCGTCATCTCTTTCGCCAGTGTCACCAGCAATACTTCAAGCGGCGCATCACCAGAAAGGGCAGGCGCTTCTCCCTGACCTTCTGGACGTTCAGCATTACGAATATATTCACGAACCTGTTCATT",
    "GATGTAAACTGTTTTTGTTCAACTGGCGTCATCTCTTTCGCCAGTGTCACCAGCATTACTTCAAGCGGCGCATCACCAGAAAGGGCAGGCGCTTCTCCCTGACCTTCTGGACGTTCAGCATTACGAATATATTCACGAACCTGTTCATTGACGTGAACCAGTCGGGCTTTGCCACCCTGGACGCCAGGTTTTGGTGACGTTGTCCAGCCTTCCTTGCGTACCCATTTATTAATGGTCTGGCGGCTATAGC",
    "GCGTTGTAACAATCCGATAATCCCTTCCCGCAGAAGCAAGGATGTAAACTGTTTTTGTTCAACTGGCGTCATCTCTTTCGCCAGTGTCACCAGCAATACTTCAAGCGGCGCATCACCAGAAAGGGCAGGCGCTTCTCCCTGACCTTCTGGACGTTCAGCATTACGAATATATTCACGAACCTGTTCATTGACGTGAACCAGTCGGGCTTTGCCACCCTGGACGCCAGGTTTTGGTGACGTTGTCCAGCCTT",
    "GTTCTTCGGTGGTCATTTTATTGCGTAATCTTTTCATATTATTTGCTATCGCGAATCCCTAAGCGTTGTAACAATCCGATAATCCCTTCCCGCAGAAGCAAGGATGTAAACTGTTTTTGTTCAACTGGCGTCATCTCTTTCGCCAGTGTCACCAGCAATACTTCAAGCGGCGCATCACCAGAAAGGGCAGGCGCTTCTCCCTGACCTTCTGGACGTTCAGCATTACGAATATATTCACGAACCTGTTCATT"
  };
  addRawSeq(&km, nseqs, seqs);
  writeOut(&km, fileout);
  char* t = (char*) calloc(12, sizeof(char));
  pos2seq(&km->ms, (uint32_t) (7*km->ms->nPos / 8), t);
  summarize(km->ms);
  printf("%s\n", t);
  free(t);
  //printf("Written\n");
  //kmerHolder* khr = readIn("/users/vqf/desktop/tmp.txt");
  //printf("Read\n");
  //summarize(khr->ms);
  //drawMs(km->ms, "/Users/vqf/Desktop/delme.txt");
  /*uint32_t init = seq2pos(km->ms, "GAAAAAAAAAA");
  kcLL* t = followTrace(&km, init, 1);
  char* yo = getTraceSeq(&km, &t);
  printf("%s\n", yo);
  resetKcLL(&t);
  free(yo);*/
  /*kcLL* ft = nextTrace(&km);
  if (ft){
    char* seq = getTraceSeq(&km, &ft);
    printf("%s\n", seq);
    free(seq);
    resetKcLL(&ft);
  }*/
  destroyKh(&km);
  return 0;
}
