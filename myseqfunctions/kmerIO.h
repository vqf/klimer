#ifndef KMERIO_H_INCLUDED
#define KMERIO_H_INCLUDED

#define _LARGEFILE64_SOURCE
#define MAGIC 0xDABAFA00
#include "kmer.h"

//#include "debugFileIO.h"


typedef struct fHeader{
  uint32_t magic;
  uint8_t  kmerLength;
  uint8_t  nBases;
} fHeader;

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
    //printf("%lu\n", result[i+1]);
  }
  free(arr);
  arr = NULL;
  return result;
}

void rebuildFromCirc(kmerHolder** khp){
  memstruct* ms = (*khp)->ms;
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    if (kc && kc->n > 0){
      while (kc){
        tIdList* ptr = kc->idflags;
        while (ptr){
          if (ptr->trace.circular){
            kcLL* circ = NULL;
            uint32_t* ns = (uint32_t*) ptr->trace.circular;
            uint32_t lns = ns[0];
            uint32_t j = 0;
            for (j = 0; j < lns; j++){
              uint32_t dest = ns[j + 1];
              kmerConnector* kc = getConnector(&ms, i, dest);
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
    uint32_t magic = MAGIC;
    fwrite(&magic, sizeof(uint32_t), 1, fout);
    fwrite(&kh->kmerSize, sizeof(uint8_t), 1, fout);
    fwrite(&kh->nBases, sizeof(uint8_t), 1, fout);
    uint32_t i = 0;
    for (i = 0; i < ms->nPos; i++){
      if (ms->kmerArray[i] && ms->kmerArray[i]->n > 0){
        kmerConnector* kc = ms->kmerArray[i];
        while (kc){
          fwrite(&i, sizeof(uint32_t), 1, fout);
          fwrite(&kc->dest, sizeof(uint32_t), 1, fout);
          uint32_t lid = tIdListLength(&kc->idflags);
          fwrite(&lid, sizeof(uint32_t), 1, fout);
          tIdList* ids = kc->idflags;
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

kmerHolder* readIn(char* file){
  FILE* fin = fopen(file, "r");
  if (fin){
    fHeader fh;
    fread(&fh, sizeof(fHeader), 1, fin);
    if (fh.magic != MAGIC){
      fprintf(stderr, "Your magic is incorrect, the file might be corrupted\n");
    }
    kmerHolder* result = initKmer(fh.kmerLength, fh.nBases);
    kcInfo coords;
    uint32_t iSize = 0;  // idList size
    uint32_t lSize = 0;  // loop size
    size_t success = fread(&coords, sizeof(kcInfo), 1, fin);
    while (success){
      uint32_t totalReads = 0;
      kmerConnector* thisKc = getConnector(&result->ms, coords.orig, coords.dest);
      tIdList* idPtr = thisKc->idflags;
      uint32_t i = 0;
      if (fread(&iSize, sizeof(uint32_t), 1, fin)){
        for (i = 0; i < iSize; i++){
          LISTTYPE n;
          uint8_t flag;
          uint32_t nReads;
          fread(&n, sizeof(LISTTYPE), 1, fin);
          fread(&flag, sizeof(uint8_t), 1, fin);
          fread(&nReads, sizeof(uint32_t), 1, fin);
          tIdList* thisId = newTIdList(n);
          thisId->trace.flag = flag;
          thisId->trace.nReads = nReads;
          totalReads += nReads;
          uint32_t* circ = NULL;
          if (fread(&lSize, sizeof(uint32_t), 1, fin)){
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
      success = fread(&coords, sizeof(kcInfo), 1, fin);
    }
    fclose(fin);
    rebuildFromCirc(&result);
    return result;
  }
  else{
    fprintf(stderr, "Could not open %s\n", file);
    return NULL;
  }
}


#endif // KMERIO_H_INCLUDED
