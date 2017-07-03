#include "myseqfunctions/kmer.h"
#include "myseqfunctions/kmerIO.h"
#include "myseqfunctions/kmerRead.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>


#define REPS 100000

int main(int argc, char** argv){
  char* k11File = NULL;
  bool MEASURETIME = false;
  for (int i = 1; i < argc; i++){
    char* rg = argv[i];
    D_(1, "%s\n", rg);
    if (rg[0] == 0x2d){
      if (EQ(rg, "-v")) DEBUG = 1;
      if (EQ(rg, "-vv")) DEBUG = 2;
      if (EQ(rg, "-t")) MEASURETIME = true;
    }
    else{
      k11File = rg;
    }
  }
  if (!k11File){
    fprintf(stderr, "Use: %s k11_index_file [-v]\n", argv[0]);
    return 1;
  }
  printf("Reading %s...\n", k11File);
  kmerHolder* kh = readIn(k11File);
  printf("Done\n");
  printf("Please enter one sequence at a time\n");
  char* seq = (char*) calloc(0xFFFF, sizeof(char));
  while (scanf("%s", seq) > 0){
    uint32_t times = 1;
    time_t ts = 0;
    time_t te = 0;
    bool done = false;
    if (MEASURETIME){
      times = REPS;
      ts = time(NULL);
    }
    for (uint32_t i = 0; i < times; i++){
      seqCollection* result = searchSeq(&kh, seq);
      if (!done){
        printSeqCollection(&kh, &result);
        printf("\n");
      }
      destroySeqCollection(&result);
      done = true;
    }
    if (MEASURETIME){
      te = time(NULL);
      printf("Done %" PRIu32 " times in %"PRIu32"\n", times, (uint32_t)((te - ts)));
    }
  }
  printf("Cleaning...\n");
  free(seq);
  destroyKh(&kh);
  return 0;
}
