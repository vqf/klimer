#include <assert.h>
#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


int main(int argc,char** argv){
  kmerHolder* kh = readIn("/users/vqf/desktop/tmp.txt");
  summarize(kh->ms);
  destroyKh(&kh);
  return 0;
}
