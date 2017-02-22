#include <stdio.h>
#include <stdlib.h>
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include "myseqfunctions/mydmacros.h"




int main(int argc,char** argv){
  if (argc < 3){
    D_(0, "Please provide at least two k11 files\n");
    exit(1);
  }
  uint8_t n = (uint8_t) (argc - 1);
  kmerHolder* khv[n];
  uint8_t i = 0;
  for (i = 1; i <= n; i++){
    khv[i] = readIn(argv[i]);
  }
  for (i = 1; i <= n; i++){
    destroyKh(&khv[i]);
  }
  return 0;
}
