#ifndef KMERREAD_H_INCLUDED
#define KMERREAD_H_INCLUDED

#include "kmer.h"
#include "traceIdList.h"
#include "stdio.h"


typedef struct traceInfo{
  LISTTYPE n;
  kmerConnector* first;
  bool extreme;
  tIdList* downstreamCompat;
  struct traceInfo* next;
} traceInfo;

typedef struct traces{
  LISTTYPE currentUID;
  traceInfo* traceList;
} traces;

typedef struct seqCollection{
  char* traceName;
  uint32_t nReads;
  char* seq;
  struct seqCollection* next;
} seqCollection;

seqCollection* newSeqCollection(){
  seqCollection* result = calloc(1, sizeof(seqCollection));
  result->traceName = NULL;
  result->seq = NULL;
  result->next = NULL;
  return result;
}

void pushSeq(seqCollection** scp, char** tNamep, uint32_t nReads, char** seqp){
  seqCollection* sc = *scp;
  char* tName = *tNamep;
  char* seq   = *seqp;
  while (sc && sc->next){
    sc = sc->next;
  }
  seqCollection* nxt = newSeqCollection();
  nxt->traceName = tName;
  nxt->seq = seq;
  nxt->nReads = nReads;
  sc->next = nxt;
}

void destroySeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  while(sc){
    seqCollection* nxt = sc->next;
    free(sc->traceName);
    free(sc->seq);
    free(sc);
    sc = nxt;
  }
  scp = NULL;
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
    if (kc) thisTrace = _getTrace(&kc->idflags, tid);
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

seqCollection* allSeqs(kmerHolder** khp){
}

traceInfo* newTraceInfo(LISTTYPE n, kmerConnector** first){
  traceInfo* result = calloc(1, sizeof(traceInfo));
  result->n = n;
  result->first = *first;
  result->extreme = true;
  result->downstreamCompat = NULL;
  result->next = NULL;
  return result;
}

void destroyTraceInfo(traceInfo** tip){
  traceInfo* ti = *tip;
  traceInfo* p = ti;
  while (p){
    p = ti->next;
    destroyTIdList(&ti->downstreamCompat, destroyCircular);
    free(ti);
    ti = p;
  }
  free(tip);
}

traces* startTraces(){
  traces* result = calloc(1, sizeof(traces));
  result->currentUID = 0;
  result->traceList = NULL;
  return result;
}

void pushTraceInfo(traces** tp, LISTTYPE n, kmerConnector** first){
  traces* t = *tp;
  traceInfo* ti = newTraceInfo(n, first);
  traceInfo* tmp = t->traceList;
  while(tmp->next){
    tmp = tmp->next;
  }
  tmp->next = ti;
}

void destroyTraces(traces** tp){
  traces* t = *tp;
  destroyTraceInfo(&t->traceList);
  free(t);
  free(tp);
}

/*
void consolidateTraces(kmerHolder** khp){
  kmerHolder* kh = *khp;
  memstruct* ms = kh->ms;
  traces* tr = startTraces();
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    kmerConnector* kc = ms->kmerArray[i];
    tIdList* il = kc->idflags;
    while (il){
      if (IS(il, FIRST_IN_TRACE)){
        kmerConnector* nxt = nextKc(&ms, &kc, il->trace.n);
      }
      il = il->next;
    }
  }
}
*/


#endif // KMERREAD_H_INCLUDED
