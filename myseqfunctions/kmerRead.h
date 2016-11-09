#ifndef KMERREAD_H_INCLUDED
#define KMERREAD_H_INCLUDED

#include "kmer.h"
#include "traceIdList.h"
#include "stdio.h"

typedef struct kcLL2D{
  kcLL* kcll;
  struct kcLL2D* next;
} kcLL2D;

kcLL2D* _newkcLL2D(kcLL** kclp){
  kcLL* kcl = *kclp;
  kcLL2D* result = (kcLL2D*) calloc(1, sizeof(kcLL2D));
  result->kcll = kcl;
  result->next = NULL;
  return result;
}

void pushKcLL(kcLL2D** kclp, kcLL** toaddp){
  kcLL2D* kcl = *kclp;
  if (kcl){
    while (kcl->next){
      kcl = kcl->next;
    }
    kcl->next = _newkcLL2D(toaddp);
  }
  else{
    kcl = _newkcLL2D(toaddp);
  }
  *kclp = kcl;
}

void destroyKcLL2D(kcLL2D** todelp){
  kcLL2D* todel = *todelp;
  while (todel){
    kcLL2D* tmp = todel;
    todel = todel->next;
    free(tmp);
  }
  todelp = NULL;
}


char getLastBase(kmerHolder** khp, uint32_t pos){
  kmerHolder* kh = *khp;
  uint8_t nb = kh->nBases;
  uint8_t md = (uint8_t) (pos % nb);
  char result = kh->ms->bi->valBase[md];
  return result;
}

kcLL* followTrace(kmerHolder** khp, uint32_t pos, LISTTYPE tid){
  kcLL* result = NULL;
  memstruct* ms = (*khp)->ms;
  kmerConnector* kc = getKcWithTId(&ms, pos, tid);
  tIdList* thisTrace = _getTrace(&kc->idflags, tid);
  while (kc && !IS(thisTrace, LAST_IN_TRACE)){

    E_(1,
      char* seq2 = calloc(12, sizeof(char));
      pos2seq(&ms, kc->dest, seq2);
      printKmerConnector(kc, seq2);
      printTIdList(kc->idflags);
      free(seq2);
      //getc(stdin)
    );

    kcpush(&result, &kc);
    kc = nextKc(&ms, &kc, tid);
  }
  result->pos = pos;
  return result;
}

kcLL* nextTrace(kmerHolder** khp){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  kcLL* result = NULL;
  uint32_t i = 0;//ms->status->current;
  while (i < ms->nPos){
    kmerConnector* kc = ms->kmerArray[i];
    while (kc && kc->n > 0){
      tIdList* l = kc->idflags;
      while (l){
        if (!IS(l, IN_USE) && IS(l, FIRST_IN_TRACE)){
          SET(l, IN_USE);
          D_(2, "Found starting trace at %lu\n", (LUI) i);
          ms->status->current = i;
          result = followTrace(khp, i, l->trace.n);
          return result;
        }
        l = l->next;
      }
      kc = kc->next;
    }
    i++;
  }
  return result;
}

char* getTraceSeq(kmerHolder** khp, kcLL** kcp){
  kmerHolder* kh = *khp;
  kcLL* kcll = *kcp;
  kcLL* tmp = kcll;
  uint32_t l = 0;
  SCALAR(tmp, l);
  uint32_t cl = (uint32_t) (kh->kmerSize + l + 1);
  char* result = (char*) calloc((size_t) cl, sizeof(char));
  pos2seq(&kh->ms, kcll->pos, result);
  uint32_t j = (uint32_t) (kh->kmerSize);
  while (kcll){
    result[j] = getLastBase(khp, kcll->kc->dest);
    j++;
    kcll = kcll->next;
  }
  return result;
}


#endif // KMERREAD_H_INCLUDED
