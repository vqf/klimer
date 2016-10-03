#include <assert.h>
#include "myseqfunctions/kmer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


int main(int argc,char** argv){
  tIdList* one = newTIdList(0);
  tIdList* two = newTIdList(6);
  insertInTIdList(&one, 2, destroyCircular);
  insertInTIdList(&one, 6, destroyCircular);
  insertInTIdList(&two, 2, destroyCircular);
  insertInTIdList(&two, 0, destroyCircular);
  tIdList* tmp = _getTrace(&one, 2);
  setAsFirst(&tmp, tmp);
  printTIdList(one);
  printTIdList(two);
  printf("isInc: %d\n", isIncluded(&one, &two));
}
