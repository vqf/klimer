#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/mydmacros.h"
#include "myseqfunctions/kmerIO.h"

#define OUT "merge_filt500.k11"
#define FILTER 500

int main(int argc,char** argv){
  if (argc < 2){
    fprintf(stderr, "Use: kmerHist k11_file1, [k11_file2, ...]\n");
    return 1;
  }
  kmerHolder* kh = readIn(argv[1]);
  for (int i = 2; i < argc; i++){
    D_(0, "Reading %s\n", argv[i]);
    kmerHolder* tkh = readIn(argv[i]);
    //mergeKh(&kh, &tkh);
    filterMergeKh(&kh, &tkh, FILTER);
    destroyKh(&tkh);
  }
  writeOut(&kh, OUT);
  destroyKh(&kh);
}
