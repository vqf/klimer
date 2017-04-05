#include <assert.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>



int main(int argc,char** argv){
  srand((unsigned) time(NULL));
  char* fileout = argv[1];
  if (argc < 2){
    fileout = "/users/vqf/desktop/tmp.txt";
  }
  kmerHolder *km = initKmer(11, 4);
  char* dtext = "GTAAAAAAAAAAAAAAACGTACGTACGAGCT";
  uint8_t l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  dtext = "GACGAGCTAACGCGTACGTAGAAAAAAAAAAAAAAAACGTACGTACGAGCTGTACGTACGTTCGTACGAGCGCGTACGTAGAAAAAAAAAAAA";
  l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  dtext = "AAAAAAAAAAAAAAGTACGTAGAAAAAAAAAAAAACCCGTACGTAAAAAA";

  l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  dtext = "CCGTACGTACGTTCGTACGAGCTACGTTCGTAAAACGTACGTACGCGTACGTACGTTCGAAAAAAAAAAAAAGACGAGCT";
  l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  dtext = "GTACGAGCGCGTACGTAGAAAAAAAAAAAAAAAAAAGTACGTAGAAAAA";
  l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  summarize(km->ms);
  writeOut(&km, fileout);
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
