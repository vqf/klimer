#include <stdio.h>
#include <stdlib.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/readFastas.h"
//#include "myseqfunctions/readFastq.h"

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
  FILE* fp = fopen(infile, "r");
  fastaReader* fr = newFastaReader(fp, 0);
  kmerHolder* kh = initKmer(11, 4);
  char fst = getNextBase(fr);
  updateKmer(&kh, &fst, addRelationship);
  while (fst){
    //printf("%c", fst);
    fst = getNextBase(fr);
    updateKmer(&kh, &fst, addRelationship);
    if (fr->newchr){
      printf(">%s\n", fr->cname);
      resetKmer(&kh);
    }
  }
  printf("\n");
  writeOut(&kh, outfile);
  destroyFastaReader(&fr);
  free(outfile);
  return 0;
}
