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
  tSet* t = newTSet(3, 5);
  unshiftTSet(&t, 7, 8);
  /*tSet* nav = t;
  unshiftTSet(&t, 9, 8);
  tSet* nav2 = t;
  unshiftTSet(&t, 10, 8);*/
  exciseTSet(&t, &t->next);
  exciseTSet(&t, &t);
  printTSet(t);
  destroyTSet(&t);


  tIdList* l = newTIdList(5);
  insertInTIdList(&l, 7);
  insertInTIdList(&l, 1);
  printTIdList(l); printf("\n");
  tIdList* n = _getTrace(&l, 1, 0);
  printTIdList(n);
  return 0;
}
