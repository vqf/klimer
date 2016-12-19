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
  kcLL* trace;
  struct seqCollection* next;
} seqCollection;

seqCollection* newSeqCollection(){
  seqCollection* result = calloc(1, sizeof(seqCollection));
  result->trace = NULL;
  result->next = NULL;
  return result;
}

seqCollection* pushSeq(seqCollection** scp, kcLL** tracep){
  // returns a pointer to the last element
  seqCollection* sc = *scp;
  while (sc && sc->next){
    sc = sc->next;
  }
  seqCollection* nxt = newSeqCollection();
  nxt->trace = kcCopy(tracep);
  sc->next = nxt;
  return nxt;
}

void clearTraceUse(kcLL** kclp){
  kcLL* kcl = *kclp;
  while (kcl){
    kmerConnector* kc = kcl->kc;
    if (kc && ISKC(kc, IN_USE)){
      UNSETKC(kc, IN_USE);
    }
    kcl = kcl->next;
  }
}

void destroySeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  while(sc){
    seqCollection* nxt = sc->next;
    free(sc->trace);
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
  if (kc) kcpush(&result, &kc);
  if (result) result->pos = pos;
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

seqCollection* copySeqCollection(seqCollection** scp){
  seqCollection* result = newSeqCollection();
  seqCollection* sc = *scp;
  while (sc){
    kcLL* trace = kcCopy(&sc->trace);
    pushSeq(&result, &trace);
    sc = sc->next;
  }
  return result;
}

seqCollection* forkSeqCollections(seqCollection** sc1p, seqCollection** sc2p){
  seqCollection* result = newSeqCollection();
  seqCollection* sc1 = *sc1p;
  seqCollection* sc2 = *sc2p;
  if (!sc2){
    result = copySeqCollection(sc1p);
    return result;
  }
  if (sc1 && sc1->trace){
    while(sc2){
      if (sc2->trace){
        D_(1, "To trace %lu\n", (LUI) sc2->trace->kc->dest);
        kcLL* start = kcCopy(&sc1->trace);
        kcLL* goon  = kcCopy(&sc2->trace);
        start->last->next = goon;
        pushSeq(&result, &start);
      }
      sc2 = sc2->next;
    }
  }
  else{
    D_(1, "Copy op\n");
    result = copySeqCollection(&sc2);
  }
  return result;
}

seqCollection* recursivelyFollowTrace(kmerHolder** khp, kmerConnector** kcp, LISTTYPE tid, uint32_t startPos){
  seqCollection* result = newSeqCollection();
  memstruct* ms = (*khp)->ms;
  kmerConnector* kc = *kcp;
  kcpush(&result->trace, &kc);
  kc = nextKc(&ms, &kc, tid);
  if (!kc) return result;
  result->trace->pos = startPos;
  tIdList* thisTrace = _getTrace(&kc->idflags, tid);
  tIdList* extra = NULL;
  bool finishIt = IS(thisTrace, LAST_IN_TRACE);
  while (kc && thisTrace && !finishIt){
    kcpush(&result->trace, &kc);
    kc = nextKc(&ms, &kc, tid);
    if (kc){
      thisTrace = _getTrace(&kc->idflags, tid);
      //
      tIdList* tmpextra = traceFirst(&kc->idflags, destroyCircular);
      intersectTIdLists(&extra, kc->idflags, destroyCircular);
      if (IS(thisTrace, LAST_IN_TRACE)){ // Shall we go on?
        D_(1, "Last in trace\n");
        while (extra){
          D_(1, "Following\n");
          if (!IS(extra, IN_USE)){
            seqCollection* addme = recursivelyFollowTrace(khp, &kc, extra->trace.n, 0);
            if (addme){
              seqCollection* old = result;
              D_(1, "Forking\n");
              result = forkSeqCollections(&result, &addme);
              destroySeqCollection(&old);
              destroySeqCollection(&addme);
            }
          }
          extra = extra->next;
          //D_(1, "Passing to trace %lu\n", (LUI) extra->trace.n);
        }
        finishIt = true;
      }
      if (tmpextra){
        mergeTIdLists(&extra, tmpextra, destroyCircular);
      }
    }
  }
  seqCollection* pointer = result;
  while (pointer && pointer->trace){
    clearTraceUse(&pointer->trace);
    pointer = pointer->next;
  }
  return result;
}

seqCollection* allTraces(kmerHolder** khp){
  seqCollection* result = newSeqCollection();
  seqCollection* pointer = result;
  memstruct* ms = (*khp)->ms;
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    D_(2, "Reading %lu\n", (LUI) i);
    kmerConnector* tc = ms->kmerArray[i];
    while (tc){
      tIdList* ti = tc->idflags;
      while (ti){
        if (IS(ti, FIRST_IN_TRACE)){
          D_(1, "Found starting trace at %lu\n", (LUI) i);
          SET(ti, IN_USE);
          seqCollection* newseqs = recursivelyFollowTrace(khp, &tc, ti->trace.n, i);
          pointer->next = newseqs;
          pointer = newseqs;
        }
        ti = ti->next;
      }
      tc = tc->next;
    }
  }
  return result;
}

uint16_t getTotalNumberOfReads(kcLL** kclp){
  uint16_t result = 0;
  kcLL* tmp = *kclp;
  while (tmp){
    result += tmp->kc->n;
    tmp = tmp->next;
  }
  return result;
}

seqCollection* bestSeqCollection(seqCollection** scp){
  seqCollection* sc = *scp;
  seqCollection* result = sc;
  uint16_t top = getTotalNumberOfReads(&sc->trace);
  while (sc && sc->next){
    sc = sc->next;
    uint16_t n = getTotalNumberOfReads(&sc->trace);
    if (n > top){
      top = n;
      result = sc;
    }
  }
  return result;
}

void printSeq(kmerHolder** khp, seqCollection** scp){
  seqCollection* sc = *scp;
  if (sc->trace){
    char* seq = getTraceSeq(khp, &sc->trace);
    printf("%s\n", seq);
    free(seq);
  }
}

void printSeqCollection(kmerHolder** khp, seqCollection** scp){
  seqCollection* sc = *scp;
  uint32_t i = 0;
  while (sc){
    i++;
    printf(">trace%lu\n", (LUI) i);
    printSeq(khp, scp);
    sc = sc->next;
  }
}

#endif // KMERREAD_H_INCLUDED
