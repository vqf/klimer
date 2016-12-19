#include <assert.h>
//#include "myseqfunctions/kmer.h"
//#include "myseqfunctions/kmerIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "myseqfunctions/traceIdList.h"


int main(int argc,char** argv){
  tIdList* a = newTIdList(7);
  addPosInTrace(&a, 8);
  insertInTIdList(&a, 3);
  addPosInTrace(&a, 9);
  addPosInTrace(&a, 19);
  printTIdList(a);
  return 0;
}
