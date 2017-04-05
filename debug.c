#include <assert.h>
//#include "myseqfunctions/kmer.h"
//#include "myseqfunctions/kmerIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "myseqfunctions/traceIdList.h"
#include "myseqfunctions/readFastas.h"



int main(int argc,char** argv){

  tIdList* n = newTIdList(1741);
  addPosInTrace(&n, 238);
  tIdList* second = addTrace(&n, 3276);
  addPosInTrace(&second, 1);
  SET(second, FIRST_IN_TRACE);
  SET(second->posInTrace, FIRST_IN_TRACE);
  traceVessel* tv = traceFirst(&n);
  printTraceLL((traceLL*) tv);

  printTIdList(n);
  return 0;
}
