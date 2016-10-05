#ifndef KMERIO_H_INCLUDED
#define KMERIO_H_INCLUDED

#define _LARGEFILE64_SOURCE

#include "kmer.h"



uint8_t writeOut(memstruct** msp, char* fname){
  memstruct* ms = *msp;
  FILE* fout = fopen(fname, "w");
  if (fout){
    uint32_t i = 0;
    for (i = 0; i < ms->nPos; i++){
      if (ms->kmerArray[i] && ms->kmerArray[i]->n > 0){
        kmerConnector* kc = ms->kmerArray[i];
        tIdList* ids = kc->idflags;
        while (kc){
          fwrite(&kc->dest, sizeof(uint32_t), 1, fout);
          uint32_t lid = tIdListLength(&kc->idflags);
          fwrite(&lid, sizeof(uint32_t), 1, fout);
          while (ids){
            fwrite(&ids->trace.n, sizeof(LISTTYPE), 1, fout);
            fwrite(&ids->trace.flag, sizeof(uint8_t), 1, fout);
            fwrite(&ids->trace.nReads, sizeof(uint32_t), 1, fout);
            kcLL* loops = (kcLL*) ids->trace.circular;
            uint32_t ll = lengthKcll(&loops);
            fwrite(&ll, sizeof(uint32_t), 1, fout);
            while (loops){
              fwrite(&loops->kc->dest, sizeof(uint32_t), 1, fout);
              loops = loops->next;
            }
            ids = ids->next;
          }
          kc = kc->next;
        }
      }
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
