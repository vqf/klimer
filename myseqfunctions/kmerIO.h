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


uint8_t writeOut(kmerHolder** khp, char* fname){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  _canonize(khp);
  FILE* fout = fopen(fname, "wb");
  char* headerInfo = "";
  if (fout){
    D_(2, "Writing header\n")
    uint32_t magic = MAGIC;
    fwrite(&magic, sizeof(uint32_t), 1, fout);
    //
    uint32_t hL = (uint32_t) (strlen(headerInfo));
    fwrite(&hL, sizeof(uint32_t), 1, fout);
    fwrite(headerInfo, sizeof(char), hL, fout);
    //
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
          fwrite(&kc->n, sizeof(uint16_t), 1, fout);
          uint32_t lid = tIdListLength(&kc->idflags);
          D_(2, "Has %lu traces\n", (LUI) lid);
          fwrite(&lid, sizeof(uint32_t), 1, fout);
          tIdList* ids = kc->idflags;
          while (ids){
            D_(2, "-New trace\n");
            fwrite(&ids->trace.n, sizeof(LISTTYPE), 1, fout);
            fwrite(&ids->trace.flag, sizeof(uint8_t), 1, fout);
            fwrite(&ids->trace.nReads, sizeof(uint32_t), 1, fout);
            uint32_t ll = tIdListLength(&ids->posInTrace);
            D_(2, "Has %lu Pos\n", (LUI) ll);
            fwrite(&ll, sizeof(uint32_t), 1, fout);
            tIdList* loops = ids->posInTrace;
            while (loops){
              fwrite(&loops->trace.n, sizeof(LISTTYPE), 1, fout);
              fwrite(&loops->trace.flag, sizeof(uint8_t), 1, fout);
              fwrite(&loops->trace.nReads, sizeof(uint32_t), 1, fout);
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
  FILE* fin = fopen(file, "rb");
  if (fin){
    uint32_t magic;
    uint32_t hL;
    char* headerInfo = NULL;
    uint8_t  kmerLength;
    uint8_t  nBases;
    D_(1, "Reading Header\n");
    fread(&magic, sizeof(uint32_t), 1, fin);
    //
    fread(&hL, sizeof(uint32_t), 1, fin);
    fread(headerInfo, sizeof(char), hL, fin);
    //
    fread(&kmerLength, sizeof(uint8_t), 1, fin);
    fread(&nBases, sizeof(uint8_t), 1, fin);
    if (magic != MAGIC){
      fprintf(stderr, "Your magic is incorrect, the file might be corrupted\n");
    }
    D_(1, "KmerLength: %" PRIu8 ", NBases: %" PRIu8 "\n", kmerLength, nBases);
    kmerHolder* result = initKmer(kmerLength, nBases);
    uint32_t orig = 0;
    uint32_t dest = 0;
    uint32_t iSize = 0;  // idList size
    uint32_t lSize = 0;  // loop size
    size_t success = fread(&orig, sizeof(uint32_t), 1, fin);
    success = fread(&dest, sizeof(uint32_t), 1, fin);
    while (success){
      D_(2, "Reading kc from %lu to %lu\n", (LUI) orig, (LUI) dest);
      uint32_t totalReads = 0;
      kmerConnector* thisKc = getConnector(&result->ms, orig, dest);
      fread(&thisKc->n, sizeof(uint16_t), 1, fin);
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
          tIdList* posPtr = thisId->posInTrace;
          if (fread(&lSize, sizeof(uint32_t), 1, fin)){
            D_(2, "Has %lu posInTrace\n", (LUI) lSize);
            uint32_t j = 0;
            for (j = 0; j < lSize; j++){
              LISTTYPE pn;
              uint8_t pflag;
              uint32_t preads;
              fread(&pn, sizeof(LISTTYPE), 1, fin);
              fread(&pflag, sizeof(uint8_t), 1, fin);
              fread(&preads, sizeof(uint32_t), 1, fin);
              tIdList* thisPos = newTIdList(pn);
              thisPos->trace.flag   = pflag;
              thisPos->trace.nReads = preads;
              if (posPtr){
                posPtr->next = thisPos;
                posPtr = posPtr->next;
              }
              else{
                thisId->posInTrace = thisPos;
                posPtr = thisId->posInTrace;
              }
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
      //thisKc->n = totalReads;
      success = fread(&orig, sizeof(uint32_t), 1, fin);
      success = fread(&dest, sizeof(uint32_t), 1, fin);
    }
    fclose(fin);
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
