#ifndef KMERREAD_H_INCLUDED
#define KMERREAD_H_INCLUDED

#include "kmer.h"
#include "traceIdList.h"



kmerConnector* getKcWithTId(kmerHolder** khp, uint32_t pos, LISTTYPE tid){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  kmerConnector* kc = ms->kmerArray[pos];
  while (kc){
    if (_getTrace(&kc->idflags, tid)){
      return kc;
    }
    kc = kc->next;
  }
  if (!kc){
    D_(0, "Could not find trace\n");
  }
  return NULL;
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
  kmerConnector* kc = getKcWithTId(khp, pos, tid);
  tIdList* thisTrace = _getTrace(&kc->idflags, tid);
  while (kc && !IS(thisTrace, LAST_IN_TRACE)){
    kcpush(&result, &kc);
    tIdList* tidl = _getTrace(&kc->idflags, tid);
    bool goOn = true;
    if (goOn && tidl->trace.circular){
      kcLL* c = (kcLL*) tidl->trace.circular;
      while (c && ISKC(c, IN_USE)){
        c = c->next;
      }
      if (c){
        SETKC(c, IN_USE);
        kc = c->kc;
      }
      else{
        goOn = false;
        kc = getKcWithTId(khp, kc->dest, tid);
      }
    }
    else{
      kc = getKcWithTId(khp, kc->dest, tid);
    }
    thisTrace = _getTrace(&kc->idflags, tid);
    D_(2, "Seq %lu\n", (LUI) kc->dest);
  }
  result->pos = pos;
  return result;
}

char* getTraceSeq(kmerHolder** khp, kcLL** kcp){
  kmerHolder* kh = *khp;
  kcLL* kcll = *kcp;
  kcLL* tmp = kcll;
  uint32_t l = 0;
  SCALAR(tmp, l);
  uint32_t cl = (uint32_t) kh->kmerSize + l;
  char* result = calloc(cl, sizeof(char));
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
