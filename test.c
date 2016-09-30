#include <assert.h>
#include "myseqfunctions/kmer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>



int main(int argc,char** argv){
  srand((unsigned) time(NULL));
  kmerHolder *km = initKmer(11);
  char* dtext = "GTAAAAAAAAAAAAAAACGTACGTACGAGCT";
  uint8_t l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);

  dtext = "AAAAAAAAAAAACGTACGTACGAGCTTAGG";
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
  dtext = "GAAAAAAAAAAAAACCCGTACGTACGTTCGTACGAGCTACGTTCGTAAAA";
  l = strlen(dtext);
  for (int i = 0; i < l; i++){
    updateKmer(&km, &dtext[i], addRelationship);
  }
  resetTrace(&km);
  summarize(km->ms);
  //drawMs(km->ms, "/Users/vqf/Desktop/delme.txt");
  destroyKh(&km);
}
