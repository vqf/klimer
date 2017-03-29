#ifndef KMERREAD_H_INCLUDED
#define KMERREAD_H_INCLUDED

#include "kmer.h"
#include "traceIdList.h"
#include "stdio.h"




char getLastBase(kmerHolder** khp, uint32_t pos){
  kmerHolder* kh = *khp;
  uint8_t nb = kh->nBases;
  uint8_t md = (uint8_t) (pos % nb);
  char result = kh->ms->bi->valBase[md];
  return result;
}

kcLL* followTrace(kmerHolder** khp, uint32_t pos, LISTTYPE tid, LISTTYPE posInTrace){
  kcLL* result = NULL;
  memstruct* ms = (*khp)->ms;
  kmerConnector* kc = getKcWithTId(&ms, pos, tid, posInTrace);
  if (!kc){
    D_(0, "Error at kc %lu, id %lu, pos %lu\n", (LUI) pos, (LUI) tid, (LUI) posInTrace);
    SEQ(&ms, ms->kmerArray[pos]);
    X_;
  }
  tIdList* thisTrace = _getTrace(&kc->idflags, tid, posInTrace);
  while (kc && !isLastInTrace(&thisTrace, posInTrace)){
    E_(1,
      char* seq2 = calloc(12, sizeof(char));
      pos2seq(&ms, kc->dest, seq2);
      printKmerConnector(kc, seq2);
      printTIdList(kc->idflags);
      free(seq2);
      //getc(stdin)
    );
    kcpush(&result, &kc, tid, posInTrace);
    kc = nextKc(&ms, &kc, tid, posInTrace);
    if (!kc){
      D_(0, "Trace %lu interrupted\n", (LUI) tid);
    }
    posInTrace++;
    if (kc) thisTrace = _getTrace(&kc->idflags, tid, posInTrace);
  }
  if (kc) kcpush(&result, &kc, tid, posInTrace);
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
        if (IS(l, FIRST_IN_TRACE)){
          tIdList* start = isInTIdList(&l->posInTrace, 1);
          if (!IS(start, IN_USE)){
            SET(start, IN_USE);
            kcpush(&result, &kc, l->trace.n, 0);
            result->posInTrace = i;
            D_(2, "Found starting trace at %lu\n", (LUI) i);
            ms->status->current = i;
            kcLL* f = followTrace(khp, i, l->trace.n, 1);
            result->next = f;
            return result;
          }
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
  pos2seq(&kh->ms, kcll->posInTrace, result);
  uint32_t j = (uint32_t) ((kh->kmerSize) - 1);
  while (kcll){
    result[j] = getLastBase(khp, kcll->kc->dest);
    j++;
    kcll = kcll->next;
  }
  return result;
}


seqCollection* recursivelyFollowTrace(kmerHolder** khp, kmerConnector** kcp, LISTTYPE tid, LISTTYPE posInTrace){
  seqCollection* result = newSeqCollection();
  memstruct* ms = (*khp)->ms;
  kmerConnector* kc = *kcp;
  kcpush(&result->trace, &kc, tid, posInTrace);
  kc = nextKc(&ms, &kc, tid, posInTrace);
  posInTrace++;
  if (!kc) return result;
  tIdList* thisTrace = _getTrace(&kc->idflags, tid, posInTrace);
  while (kc && thisTrace){
    kcpush(&result->trace, &kc, tid, posInTrace);

    kc = nextKc(&ms, &kc, tid, posInTrace);
    posInTrace++;
    if (kc) thisTrace = _getTrace(&kc->idflags, tid, posInTrace);
  }
  return result;
}

seqCollection* allTraces(kmerHolder** khp){
  seqCollection* result = newSeqCollection();
  seqCollection* pointer = result;
  memstruct* ms = (*khp)->ms;
  uint32_t i = 0;
  for (i = 0; i < ms->nPos; i++){
    D_(3, "Reading %lu\n", (LUI) i);
    kmerConnector* tc = ms->kmerArray[i];
    while (tc){
      tIdList* ti = tc->idflags;
      while (ti){
        if (IS(ti, FIRST_IN_TRACE)){
          D_(1, "Found starting trace at %lu\n", (LUI) i);
          seqCollection* newseqs = recursivelyFollowTrace(khp, &tc, ti->trace.n, 1);
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
