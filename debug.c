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

  tIdList* l = newTIdList(5);
  insertInTIdList(&l, 7);
  insertInTIdList(&l, 1);
  delTIdFromList(&l, 2);
  printTIdList(l); printf("\n");
  tIdList* m = _getTrace(&l, 1, 0);
  delTIdFromList(&l, 1);
  printTIdList(l);
  return 0;
}
