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
  tIdList* n = _getTrace(&l, 5, 0);
  SET(n, FIRST_IN_TRACE);
  n = _getTrace(&l, 1, 0);
  SET(n, FIRST_IN_TRACE);
  n = _getTrace(&l, 7, 0);
  SET(n, FIRST_IN_TRACE);
  printTIdList(n);
  traceVessel* v = traceFirst(&l);
  delTraceFromVessel(&v, 7);
  delTraceFromVessel(&v, 1);
  delTraceFromVessel(&v, 5);
  printf("\n");
  printTraceLL((traceLL*) v);
  return 0;
}
