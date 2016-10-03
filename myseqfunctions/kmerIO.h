#ifndef KMERIO_H_INCLUDED
#define KMERIO_H_INCLUDED

#include "kmer.h"

typedef struct wKc{
  uint32_t pos;
  uint32_t dest;
  LISTTYPE tId;
} wKc;


uint8_t writeOut(memstruct* ms, char* fname){
  FILE* fout = fopen(fname, "w");
  uint32_t cpos = 0;
  if (fout){
    uint32_t emtpy = NOKMER;
    // NOKMER on initial memory chunk (tId is zero)
    for (int i = 0; i < ms->nPos; i++){
      size_t chSize = fwrite(&emtpy, sizeof(uint32_t), 1, fout);
      cpos += chSize;
    }
    fclose(fout);
  }
  else{
    printf("Could not open %s\n", fname);
    return 0x01;
  }
  return 0x00;
}


#endif // KMERIO_H_INCLUDED
