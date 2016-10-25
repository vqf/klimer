#ifndef KMERIO_H_INCLUDED
#define KMERIO_H_INCLUDED

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif /* _LARGEFILE64_SOURCE */
#define MAGIC 0xDABAFA00
#include "kmer.h"


//#include "debugFileIO.h"


typedef struct kcInfo{
  uint32_t orig;
  uint32_t dest;
} kcInfo;

uint32_t* rebuildCirc(uint32_t l, uint32_t* arr){
  uint32_t* result = (uint32_t*) calloc(l + 1, sizeof(uint32_t));
  result[0] = l;
  uint32_t i = 0;
  for (i = 0; i < l; i++){
    result[i+1] = arr[i];
    D_(2, "Adding uint32_t %lu\n", (LUI) arr[i]);
  }
  free(arr);
  arr = NULL;
  D_(2, "Rebuilt array");
  return result;
}

void rebuildFromCirc(kmerHolder** khp){
  memstruct* ms = (*khp)->ms;
  D_(2, "Rebuilding circ traces\n");
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    if (kc && kc->n > 0){
      while (kc){
        tIdList* ptr = kc->idflags;
        uint32_t orig = kc->dest;
        while (ptr){
          if (ptr->trace.circular){
            kcLL* circ = NULL;
            uint32_t* ns = (uint32_t*) ptr->trace.circular;
            uint32_t lns = ns[0];
            uint32_t j = 0;
            for (j = 0; j < lns; j++){
              uint32_t dest = ns[j + 1];
              kmerConnector* kc = getConnector(&ms, orig, dest);
              D_(2, "Adding kc from %lu to %lu\n", (LUI) orig, (LUI) dest);
              kcpush(&circ, &kc);
            }
            free(ptr->trace.circular); ptr->trace.circular = NULL;
            ptr->trace.circular = circ;
          }
          ptr = ptr->next;
        }
        kc = kc->next;
      }
    }
  }
}

uint8_t writeOut(kmerHolder** khp, char* fname){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  FILE* fout = fopen(fname, "w");
  if (fout){
    D_(2, "Writing header\n")
    uint32_t magic = MAGIC;
    fwrite(&magic, sizeof(uint32_t), 1, fout);
    fwrite(&kh->kmerSize, sizeof(uint8_t), 1, fout);
    fwrite(&kh->nBases, sizeof(uint8_t), 1, fout);
    uint32_t i = 0;
    for (i = 0; i < ms->nPos; i++){
      if (ms->kmerArray[i] && ms->kmerArray[i]->n > 0){
        kmerConnector* kc = ms->kmerArray[i];
        D_(2, "Writing kcs from %lu\n", (LUI) i);
        while (kc){
          fwrite(&i, sizeof(uint32_t), 1, fout);
          D_(2, "...to %lu\n", (LUI) kc->dest);
          fwrite(&kc->dest, sizeof(uint32_t), 1, fout);
          uint32_t lid = tIdListLength(&kc->idflags);
          D_(2, "Has %lu traces\n", (LUI) lid);
          fwrite(&lid, sizeof(uint32_t), 1, fout);
          tIdList* ids = kc->idflags;
          while (ids){
            D_(2, "-New trace\n");
            fwrite(&ids->trace.n, sizeof(LISTTYPE), 1, fout);
            fwrite(&ids->trace.flag, sizeof(uint8_t), 1, fout);
            fwrite(&ids->trace.nReads, sizeof(uint32_t), 1, fout);
            kcLL* loops = (kcLL*) ids->trace.circular;
            uint32_t ll = lengthKcll(&loops);
            D_(2, "Has %lu loops\n", (LUI) ll);
            fwrite(&ll, sizeof(uint32_t), 1, fout);
            while (loops){
              fwrite(&loops->kc->dest, sizeof(uint32_t), 1, fout);
              loops = loops->next;
            }
            D_(2, "-Next trace\n");
            ids = ids->next;
          }
          D_(2, "NextKc\n");
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

kmerHolder* readIn(char* file){
  FILE* fin = fopen(file, "r");
  if (fin){
    uint32_t magic;
    uint8_t  kmerLength;
    uint8_t  nBases;
    D_(2, "Reading Header\n");
    fread(&magic, sizeof(uint32_t), 1, fin);
    fread(&kmerLength, sizeof(uint8_t), 1, fin);
    fread(&nBases, sizeof(uint8_t), 1, fin);
    if (magic != MAGIC){
      fprintf(stderr, "Your magic is incorrect, the file might be corrupted\n");
    }
    kmerHolder* result = initKmer(kmerLength, nBases);
    uint32_t orig;
    uint32_t dest;
    uint32_t iSize = 0;  // idList size
    uint32_t lSize = 0;  // loop size
    size_t success = fread(&orig, sizeof(uint32_t), 1, fin);
    success = fread(&dest, sizeof(uint32_t), 1, fin);
    while (success){
      D_(2, "Reading kc from %lu to %lu\n", (LUI) orig, (LUI) dest);
      uint32_t totalReads = 0;
      kmerConnector* thisKc = getConnector(&result->ms, orig, dest);
      tIdList* idPtr = thisKc->idflags;
      uint32_t i = 0;
      if (fread(&iSize, sizeof(uint32_t), 1, fin)){
        D_(2, "Has %lu traces\n", (LUI) iSize);
        for (i = 0; i < iSize; i++){
          LISTTYPE n;
          uint8_t flag;
          uint32_t nReads;
          D_(2, "-New trace\n");
          fread(&n, sizeof(LISTTYPE), 1, fin);
          fread(&flag, sizeof(uint8_t), 1, fin);
          fread(&nReads, sizeof(uint32_t), 1, fin);
          tIdList* thisId = newTIdList(n);
          thisId->trace.flag = flag;
          thisId->trace.nReads = nReads;
          totalReads += nReads;
          uint32_t* circ = NULL;
          if (fread(&lSize, sizeof(uint32_t), 1, fin)){
            D_(2, "Has %lu loops\n", (LUI) lSize);
            if (lSize){
              circ = (uint32_t*) calloc(lSize, sizeof(uint32_t));
              fread(circ, sizeof(uint32_t), lSize, fin);
              thisId->trace.circular = rebuildCirc(lSize, circ);
            }
          }
          if (idPtr){
            idPtr->next = thisId;
            idPtr = idPtr->next;
          }
          else{
            thisKc->idflags = thisId;
            idPtr = thisId;
          }
        }
      }
      thisKc->n = totalReads;
      success = fread(&orig, sizeof(uint32_t), 1, fin);
      success = fread(&dest, sizeof(uint32_t), 1, fin);
    }
    fclose(fin);
    rebuildFromCirc(&result);
    cleanTraceStatus(&result);
    result->ms->status->current = 0;
    result->ms->status->cId = 0;
    return result;
  }
  else{
    fprintf(stderr, "Could not open %s\n", file);
    return NULL;
  }
}


#endif // KMERIO_H_INCLUDED
